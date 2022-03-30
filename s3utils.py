from io import StringIO, BytesIO
from pandas import read_csv, ExcelWriter
from boto3 import client, resource


AWS_ACCES_KEY_ID = "AKIATNNTOXHYYITGFJP5"
AWS_SECRET_ACCES_KEY = "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"


def check_if_key_exists(key, bucket):
    s3_resource = resource(
        "s3",
        aws_access_key_id=AWS_ACCES_KEY_ID,
        aws_secret_access_key=AWS_SECRET_ACCES_KEY,
    )
    my_bucket = s3_resource.Bucket(bucket)
    # key = os.path.relpath(complete_path, bucket) En windows no lee igual
    if list(my_bucket.objects.filter(Prefix=key)):
        return True
    else:
        return False


def folders(experiment, s3_client, bucket_name, folder):
    experiments_folders = list(
        map(
            lambda folder: folder["Prefix"][:-1],
            s3_client.list_objects_v2(Bucket=bucket_name, Prefix="", Delimiter="/").get(
                "CommonPrefixes"
            ),
        )
    )
    if experiment not in experiments_folders:
        s3_client.put_object(Bucket=bucket_name, Key=(experiment + "/"))
        if folder != "":
            s3_client.put_object(Bucket=bucket_name, Key=(experiment + f"/{folder}/"))
    else:
        folders_in_experiment = s3_client.list_objects_v2(
            Bucket=bucket_name, Prefix=experiment + "/", Delimiter="/"
        ).get("CommonPrefixes")
        if folders_in_experiment == None or experiment + f"/{folder}" not in list(
            map(lambda folder: folder["Prefix"][:-1], folders_in_experiment)
        ):
            if folder != "":
                s3_client.put_object(
                    Bucket=bucket_name, Key=(experiment + f"/{folder}/")
                )


def df_to_S3(df, s3_client, bucket_name, experiment, folder, file_name, overwrite):
    if folder != "":
        path = f"/{folder}/{file_name}.csv"
    else:
        path = f"/{file_name}.csv"
    csv_buf = StringIO()
    df.to_csv(csv_buf, index=False)
    csv_buf.seek(0)
    if not check_if_key_exists((experiment + path), bucket_name):
        print(f"Uploading  to s3://{bucket_name}{(experiment + path)}")
        s3_client.put_object(
            Bucket=bucket_name, Body=csv_buf.getvalue(), Key=(experiment + path)
        )
    elif check_if_key_exists((experiment + path), bucket_name) and overwrite:
        print(f"Overwriting  to s3://{bucket_name}{(experiment + path)}")
        s3_client.put_object(
            Bucket=bucket_name, Body=csv_buf.getvalue(), Key=(experiment + path)
        )
    elif check_if_key_exists((experiment + path), bucket_name) and not overwrite:
        print(f"File already exists: {(experiment + path)} and overwrite FALSE")


def raw_to_s3(data, s3_client, bucket_name, experiment):
    """Upload raw file to S3 with the name = EXPERIMENT/raw.xlsx"""
    with BytesIO() as output:
        with ExcelWriter(output, engine="xlsxwriter") as writer:
            for name, df in data.items():
                df.to_excel(writer, sheet_name=name, index=False)
        data = output.getvalue()
    print(f"Uploading  to s3://{bucket_name}{experiment}/raw.xlsx")
    s3_client.put_object(Bucket=bucket_name, Body=data, Key=f"{experiment}/raw.xlsx")


def getAnimalFile(s3_client):
    nombre = s3_client.list_objects_v2(Bucket="siwaothers", Prefix="all_animalfile")[
        "Contents"
    ][0]["Key"]
    df_string = (
        s3_client.get_object(Bucket="siwaothers", Key=nombre)["Body"]
        .read()
        .decode("utf-8")
    )
    return read_csv(StringIO(df_string), sep="\t")
