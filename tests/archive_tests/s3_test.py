import sys
import boto3
from botocore.client import Config
from s3path import PureS3Path, register_configuration_parameter, S3Path
import os

sys.path.append(".")

def s3_init():

    minio_bucket_path = PureS3Path('/')
    s3 = boto3.resource(
                's3',
                endpoint_url=os.environ.get("MINIO_URL", 'http://localhost:9000'),
                aws_access_key_id=os.environ.get("MINIO_ACCESS_KEY"),
                aws_secret_access_key=os.environ.get("MINIO_SECRET_KEY"),
                config=Config(signature_version='s3v4'),
                region_name='us-east-1')

    register_configuration_parameter(minio_bucket_path, resource=s3)
    return s3

if __name__ == "__main__":

    s3 = s3_init()
        
    from corems.mass_spectrum.input.massList import ReadMassList
    from corems.mass_spectrum.calc.Calibration import MzDomainCalibration

    filepath = "1/ESI_NEG_ESFA.ascii"
    s3path = S3Path('/fticr-data/' + filepath)

    mass_list_reader = ReadMassList(s3path)
    #mass_list_reader = ReadMassList(s3path, header_lines=7, isCentroid=False, isThermoProfile=True)
    mass_spectrum = mass_list_reader.get_mass_spectrum(-1)
    
    s3path = S3Path('/fticr-data/') / "1/SRFA.ref"

    mass_spectrum.filter_by_noise_threshold()

    MzDomainCalibration(mass_spectrum, s3path).run()

