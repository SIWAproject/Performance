import pandas as pd
from re import findall
from boto3 import client
from dataclasses import dataclass
from importlib import reload
import s3utils

reload(s3utils)
from s3utils import *


@dataclass
class Performance:
    filename_input: str
    sheet_name: str
    bucket_name: str
    experiment: str
    folder: str
    parameters: dict
    overwrite: bool


    def __post_init__(self, overwrite=False):
        s3_client = client(
            "s3",
            aws_access_key_id=AWS_ACCES_KEY_ID,
            aws_secret_access_key=AWS_SECRET_ACCES_KEY,
        )
        folders(self.experiment, s3_client, self.bucket_name, self.folder)
        self.overwrite = overwrite
        self.raw = pd.read_excel(self.filename_input, sheet_name=self.sheet_name)
        # subir RAW a S3 - SIEMPRE
        raw_to_s3(
            pd.read_excel(self.filename_input, sheet_name=None),
            s3_client,
            self.bucket_name,
            self.experiment,
        )
        main = self.extract(["Pen", "Trt", "Sex"]).rename(columns={"Trts": "Treatment"})
        # buscar pen y house del animal file
        animal_file = getAnimalFile(s3_client).query("Project==@self.experiment")[
            ["AnimalID", "House", "Pen"]
        ]
        main = pd.merge(main, animal_file, on=("Pen"), how="left").convert_dtypes()
        for variable in self.parameters.keys():
            print(f"Processing {variable}")
            df = self.ETL(**self.parameters[variable])
            if df.shape[0] == 0:
                continue
            if variable.startswith(("adg", "fcr", "adfi")) and variable.endswith(
                "phase"
            ):
                var_name = variable.split("_")[0].upper()
                df_extra = self.extract([f"StoF{var_name}"]).rename(
                    columns={
                        f"StoF{var_name}": self.parameters[variable]["variable"][0]
                    }
                )
                df_extra["phase"] = "StoF"
                df = pd.concat([df, df_extra])
            df = pd.merge(main, df, on="main_id").drop(
                columns=["main_id", "Sex", "Trt"]
            )
            df["House"] = 1
            df_to_S3(
                df,
                s3_client,
                self.bucket_name,
                self.experiment,
                self.folder,
                variable,
                overwrite=self.overwrite,
            )
        df = pd.merge(
            main, self.extract(["MortCount", "Mortality"]), on="main_id"
        ).drop(columns=["main_id", "Sex", "Trt"])
        df["House"] = 1
        df_to_S3(
            df,
            s3_client,
            self.bucket_name,
            self.experiment,
            self.folder,
            "mortality",
            overwrite=self.overwrite,
        )

    def extract(self, columns: list) -> pd.DataFrame:
        """Se extraen unas columnas en particular de df raw

        Args:
            columns (list): Nombres de las columnas que desean ser extraidas

        Returns:
            pd.DataFrame: dataframe con las columnas que se deseaban mas main_id
        """
        return (
            self.raw.loc[:, columns].reset_index().rename(columns={"index": "main_id"})
        )

    def transform(self, row: list, columns: list) -> pd.DataFrame:
        """Se transforma una fila de raw data en un DataFrame

        Args:
            row (list): fila que se desea transformar
            columns (list): nombres de las columnas del dataframe creado

        Returns:
            pd.DataFrame:
        """
        df = pd.DataFrame(row[1:], columns=columns)
        df["main_id"] = row[0]
        df[self.name] = self.periods
        return df

    def process(self, variable: str, cols: list) -> pd.DataFrame:
        """Concatenacion de cada una de las filas convertidas a DataFrames

        Args:
            variable (str): nombre de la columna de la data transformada
            cols (list): nombres de las columnas que van a ser procesadas

        Returns:
            pd.DataFrame:
        """
        df = pd.concat(
            [self.transform(row, variable) for row in self.extract(cols).to_numpy()]
        ).reset_index(drop=True)
        return df.convert_dtypes()

    def ETL(self, variable: list, cols_lambda, name: str):
        cols = list(filter(cols_lambda, self.raw.columns))
        self.periods = list(map(lambda x: int(findall("[0-9]+", x)[0]), cols))
        self.name = name
        df = self.process(
            variable,
            cols,
        )
        return df


@dataclass
class Treatments:
    filename_input: str
    sheet_name: str
    bucket_name: str
    experiment: str
    folder: str
    filename_output: str
    overwrite: bool

    def __post_init__(self, overwrite=False):
        self.raw = pd.read_excel(self.filename_input, sheet_name=self.sheet_name)
        self.overwrite = overwrite

        treatments = pd.concat(
            [self.transform(phase) for phase in self.raw.Phase.unique()]
        )
        s3_client = client(
            "s3",
            aws_access_key_id=AWS_ACCES_KEY_ID,
            aws_secret_access_key=AWS_SECRET_ACCES_KEY,
        )
        folders(self.experiment, s3_client, self.bucket_name, self.folder)
        df_to_S3(
            treatments,
            s3_client,
            self.bucket_name,
            self.experiment,
            self.folder,
            self.filename_output,
            self.overwrite,
        )

    def transform(self, phase: int) -> pd.DataFrame:
        """Se buscan las filas de una fase particular del experimento y se
        transforman los datos con melt(deshacer una tabla pivoteada)


        Args:
            phase (int): int de la phase que se desea transformar

        Returns:
            pd.DataFrame:
        """

        df = pd.melt(
            self.raw.query("Phase==@phase"),
            id_vars=self.raw.columns[0],
            value_vars=self.raw.columns[1:-1],
            var_name="Treatment",
            value_name="grams",
        )
        df["Phase"] = phase
        return df


@dataclass
class Results:
    filename_input: str
    sheet_name: str
    bucket_name: str
    experiment: str
    folder: str
    filename_output: str
    overwrite: bool

    def __post_init__(self, overwrite=False):
        self.setVariables()
        self.ETL()
        self.overwrite = overwrite

    def setVariables(self) -> None:
        """Se definen otros atributos que no son seteados por el constructor de la clase"""
        self.raw = pd.read_excel(self.filename_input, sheet_name=self.sheet_name)
        self.treatments = pd.Series(
            list(map(lambda x: int(x[-1]), self.raw.columns[1:-1]))
        ).unique()

    def ETL(self) -> None:
        s3_client = client(
            "s3",
            aws_access_key_id=AWS_ACCES_KEY_ID,
            aws_secret_access_key=AWS_SECRET_ACCES_KEY,
        )
        folders(self.experiment, s3_client, self.bucket_name, self.folder)

        df = pd.concat([self.process(phase) for phase in self.raw.Phase.unique()])
        df_to_S3(
            df,
            s3_client,
            self.bucket_name,
            self.experiment,
            self.folder,
            self.filename_output,
            self.overwrite,
        )

    def transform(self, treatment: int) -> pd.DataFrame:
        """Se seleccionan las columnas de un tratamiento y se deshace
        la tabla pivoteada usando pandas.melt()

        Args:
            treatment (int): entero que representa un tratamiento

        Returns:
            pd.DataFrame:
        """
        data = self.df_phase.loc[
            :,
            [self.df_phase.columns[0]]
            + list(
                filter(
                    lambda col: col.endswith(str(treatment)),
                    self.df_phase.columns[1:-1],
                )
            ),
        ]
        data.columns = list(map(lambda col: col.split()[0].title(), data.columns))
        melt = pd.melt(
            data,
            id_vars=data.columns[0],
            value_vars=data.columns[1:],
            var_name="Category",
            value_name="Value",
        )
        melt["Treatment"] = treatment
        return melt

    def process(self, phase: int) -> pd.DataFrame:
        """Concatenacion de DataFrames de cada uno de los tratamientos

        Args:
            phase (int): fase del experimento que se desea procesar

        Returns:
            pd.DataFrame:
        """
        self.df_phase = self.raw.query("Phase==@phase")
        df = pd.concat([self.transform(treatment) for treatment in self.treatments])
        df["Phase"] = phase
        return df
