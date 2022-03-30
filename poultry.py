import numpy as np
import boto3
import pandas as pd
from dataclasses import dataclass
from importlib import reload
import s3utils

reload(s3utils)
from s3utils import *


@dataclass
class Treatment:
    filename: str
    sheet_name: str
    bucket_name: str
    experiment: str
    folder_aws: str
    filename_output: str
    overwrite: bool

    def __post_init__(self, overwrite=False) -> None:
        self.setVariables()
        self.transform()
        self.setTreatments()
        self.setDays()
        self.loadToS3()
        self.overwrite = overwrite

    def setVariables(self) -> None:
        """Se carga raw data, se definen los nombres de las columnas del archivo de salida
        y se identifican los indices de cada fase
        """
        self.raw = pd.read_excel(self.filename, sheet_name=self.sheet_name, header=None)
        self.name_cols = self.raw.iloc[1, 2:].unique()
        self.idx_phases = np.where(self.raw.iloc[0].isna() == False)[0][2:]

    def transform(self) -> None:
        """Se transforma la seccion de cada fase en un DataFrame y luego se concatenan"""
        steps = np.concatenate(
            (np.diff(self.idx_phases), [self.raw.shape[1] - self.idx_phases[-1]])
        )
        self.df = pd.concat(
            [self.extract(i, step) for i, step in zip(self.idx_phases, steps)]
        ).reset_index(drop=True)

    def extract(self, i: int, step: int) -> pd.DataFrame:
        """Se extra una seccion de una fase en especifico

        Args:
            i (int): indice de la primera columna de la seccion
            step (int): i+step retorna el indice de la ultima columna de la seccion

        Returns:
            pd.DataFrame: region extraida de raw data
        """
        df = self.raw.iloc[2:, i : i + step]
        df.columns = self.name_cols
        return df

    def setTreatments(self) -> None:
        """Se setea la columna 'Treatment' en 'df'"""
        self.treatments = self.raw.iloc[2:, 0].fillna(0).values.astype(int)
        self.phases = list(self.raw.iloc[0, self.idx_phases])
        self.df["Treatment"] = list(self.treatments) * len(self.phases)

    def setDays(self) -> None:
        """Se setea la columna 'days' en 'df'"""
        days = []
        for day in map(lambda stage: stage.split()[-1], self.phases):
            days += [day] * len(self.treatments)
        self.df["days"] = days
        self.df.days = self.df.days.astype(str).str.rstrip("d")

    def loadToS3(self) -> None:
        """Se carga 'df' a S3"""
        s3_client = client(
            "s3",
            aws_access_key_id=AWS_ACCES_KEY_ID,
            aws_secret_access_key=AWS_SECRET_ACCES_KEY,
        )
        folders(self.experiment, s3_client, self.bucket_name, self.folder_aws)

        df_to_S3(
            self.df,
            s3_client,
            self.bucket_name,
            self.experiment,
            self.folder_aws,
            self.filename_output,
            self.overwrite,
        )
        phase2days = pd.DataFrame(
            list(map(lambda x: x.split(), self.phases)), columns=["phase", "days"]
        )
        df_to_S3(
            phase2days,
            s3_client,
            self.bucket_name,
            self.experiment,
            self.folder_aws,
            "phases2days",
            self.overwrite,
        )


@dataclass
class DataToStat:
    filename: str
    sheet_name: str
    bucket_name: str
    experiment: str
    folder_aws: str
    parameters: dict
    phase2days: dict
    overwrite: bool

    def __post_init__(self, overwrite=False) -> None:
        self.raw = pd.read_excel(self.filename, sheet_name=self.sheet_name)
        self.overwrite = overwrite
        s3_client = client(
            "s3",
            aws_access_key_id=AWS_ACCES_KEY_ID,
            aws_secret_access_key=AWS_SECRET_ACCES_KEY,
        )
        folders(self.experiment, s3_client, self.bucket_name, self.folder_aws)
        raw_to_s3(
            pd.read_excel(self.filename, sheet_name=None),
            s3_client,
            self.bucket_name,
            self.experiment,
        )
        treatments = self.extract(["House", "Pen", "Trts"]).rename(
            columns={"Trts": "Treatment"}
        )
        if "Jaulas" in treatments.House.unique():
            treatments.House = treatments.House.str.replace("Jaulas", "5").astype(int)
        animal_file = getAnimalFile(s3_client).query("Project==@self.experiment")[
            ["AnimalID", "House", "Pen"]
        ]
        main = pd.merge(
            treatments, animal_file, on=("House", "Pen"), how="left"
        ).convert_dtypes()
        for variable in self.parameters.keys():
            print(f"Processing {variable}")
            df = self.ETL(self.phase2days, **self.parameters[variable])
            if df.shape[0] == 0:
                continue
            df = pd.merge(main, df, on="main_id").drop(columns=["main_id"])
            df_to_S3(
                df,
                s3_client,
                self.bucket_name,
                self.experiment,
                self.folder_aws,
                variable,
                self.overwrite,
            )

    def extract(self, columns: list) -> pd.DataFrame:
        """Se extraen unas columnas en especifico de raw data y se agrega la columna 'main_id'

        Args:
            columns (list): columnas que se quieren extraer

        Returns:
            pd.DataFrame:
        """
        return (
            self.raw.loc[:, columns].reset_index().rename(columns={"index": "main_id"})
        )

    def transform(
        self,
        row: list,
    ) -> pd.DataFrame:
        """Se transforma una fila de raw data en un DataFrame
        donde se adicionan las columnas 'days y 'main'

        Args:
            row (list): indice de la fila que se desea transformar

        Returns:
            pd.DataFrame:
        """
        df = pd.DataFrame(
            [row[i : i + self.step] for i in range(1, len(row), self.step)],
            columns=self.name_cols,
        )
        df["main_id"] = row[0]
        df["days"] = self.days
        return df

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """Se concatena la lista de DataFrames, donde cada uno de estos es la fila
        que fue transformada.

        Args:
            data (pd.DataFrame):

        Returns:
            pd.DataFrame:
        """
        df = pd.concat([self.transform(row) for row in data.to_numpy()]).reset_index(
            drop=True
        )
        return df.convert_dtypes()

    def ETL(
        self, phase2days: dict, name_cols: list, startswith: str, p2d: bool
    ) -> pd.DataFrame:
        """

        Args:
            phase2days (dict): diccionario que permite convertir una fase en el rango de dias
            name_cols (list): nombre de las columnas del DataFrame retornado por esta funcion
            startswith (str): nombre de la variable de performance que se desea transformar
            p2d (bool): booleano que indica si se debe convertir de phase a dias(por ejemplo en count)

        Returns:
            pd.DataFrame:
        """
        cols = list(filter(lambda col: col.startswith(startswith), self.raw.columns))
        self.name_cols = name_cols
        self.step = len(name_cols)
        self.days = [cols[stage].split()[1] for stage in range(0, len(cols), self.step)]
        df = self.process(self.extract(cols))
        if p2d:
            df.days.replace(phase2days, inplace=True)
        df.days = df.days.astype(str).str.rstrip("d")
        return df
