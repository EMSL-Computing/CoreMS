import sys
import boto3
from botocore.client import Config
from s3path import PureS3Path, register_configuration_parameter, S3Path
import os
from minio import Minio

sys.path.append(".")

def s3_init():

    minio_bucket_path = PureS3Path('/')
    s3 = boto3.resource(
                's3',
                endpoint_url=os.environ.get("MINIO_URL", 'http://localhost:64977'),
                aws_access_key_id=os.environ.get("MINIO_ACCESS_KEY"),
                aws_secret_access_key=os.environ.get("MINIO_SECRET_KEY"),
                config=Config(signature_version='s3v4'),
                region_name='us-east-1')

    register_configuration_parameter(minio_bucket_path, resource=s3)
    return s3

def minio_init():

    minio = Minio(
            os.environ.get("MINIO_URL", 'localhost:64977').replace('http://', ''),
            access_key=os.environ.get("MINIO_ACCESS_KEY"),
            secret_key=os.environ.get("MINIO_SECRET_KEY"),
            secure=False
                )

    return minio

def check_create_buckets(minio, buckets_list):

    buckets = ['fticr-data', 'gcms-data']
    for bucket in buckets:
        if not minio.bucket_exists(bucket):
            minio.make_bucket(bucket)

s3 = s3_init()

if __name__ == "__main__":
    
    from corems.mass_spectrum.input.massList import ReadMassList

    filepath = "1/NEG_ESI_SRFA_CoreMS.xlsx"
    s3path = S3Path('/fticr-data/' + filepath)

    mass_list_reader = ReadMassList(s3path)
    #mass_list_reader = ReadMassList(s3path, header_lines=7, isCentroid=False, isThermoProfile=True)
    mass_spectrum = mass_list_reader.get_mass_spectrum(-1)
    print(mass_spectrum)



    '''
    from corems.mass_spectra.input.andiNetCDF import ReadAndiNetCDF
    filepath = "1/28-Combined_AA_S.cdf"
    s3path = S3Path('/gcms-data/' + filepath)
    reader_gcms = ReadAndiNetCDF(s3path)
    '''

    '''
    from corems.transient.input.brukerSolarix import ReadBrukerSolarix
    from matplotlib import pyplot
    filepath = "1/20190205_WK_SRFA_opt_000001.d/"
    s3path = S3Path('/fticr-data/' + filepath)

    with ReadBrukerSolarix(s3path) as bruker_transient:
        mass_spectrum_obj = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=True)
        mass_spectrum_obj.plot_profile_and_noise_threshold()
        #pyplot.show()
    ''' 


    ''' 
    from corems.mass_spectra.input.boosterHDF5 import ReadHDF_BoosterMassSpectra
    from matplotlib import pyplot
    filepath = "1/ESFA_100k_9767-13548_chB.A_re_pc_CoAddAll_mFT.h5/"
    s3path = S3Path('/fticr-data/' + filepath)

    mass_list_reader = ReadHDF_BoosterMassSpectra(s3path)
    ''' 

    ''' 
    from corems.mass_spectra.input.brukerSolarix import ReadBruker_SolarixTransientMassSpectra

    filepath = "1/NEG_ESI_SRFA_Auto.d/"
    s3path = S3Path('/fticr-data/' + filepath)

    read_lcms = ReadBruker_SolarixTransientMassSpectra(s3path)
    read_lcms.start()
    read_lcms.join()
    ''' 

    from corems.encapsulation.factory.parameters import MSParameters
    from corems.mass_spectra.input import rawFileReader
    import matplotlib.pyplot as plt

    filepath = "1/SRFA_NEG_ESI_ORB.raw/"
    s3path = S3Path('/fticr-data/' + filepath)

        # change parameters here
    MSParameters.mass_spectrum.threshold_method = 'relative_abundance'
    MSParameters.mass_spectrum.relative_abundance_threshold = 1

    # creates the parser obj
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(s3path)

    # sums all the mass spectra
    mass_spectrum = parser.get_average_mass_spectrum_in_scan_range()

    # sums scans in selected range
    mass_spectrum = parser.get_average_mass_spectrum_in_scan_range(first_scan=1, last_scan=5)

    scans_list = [1,4,6,9]
    # sums scans in selected range
    mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans_list)

    mass_spectrum.plot_mz_domain_profile()
    mass_spectrum.plot_profile_and_noise_threshold()

    #print("polarity", mass_spectrum.polarity)
    plt.savefig("test.png")