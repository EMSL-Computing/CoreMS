__author__ = "Yuri E. Corilo"
__date__ = "Nov 11, 2019"


from io import StringIO, BytesIO
from pathlib import Path
from copy import deepcopy

from pandas import read_csv, read_pickle, read_excel
import chardet
from pandas.core.frame import DataFrame
from s3path import S3Path

from corems.encapsulation.factory.processingSetting import DataInputSetting
from corems.encapsulation.factory.parameters import default_parameters
from corems.encapsulation.constant import Labels
from corems.encapsulation.input.parameter_from_json import load_and_set_parameters_class, load_and_set_parameters_ms

class MassListBaseClass:
    '''
    # The MassListBaseClass object reads mass list data types and returns the mass spectrum obj

    # Parameters
    ----------
    ## file_location : Path or S3Path
        file_location is full data path
    ## * delimiter: str
        Delimeter to read text based files ("," "\t", " ", "  ", etc)
    ## **data_type : str
        The keyword argument data_type is used to determine what type of file to read:
        pandas(pickle), txt(.csv, txt, asci, etc) or HDF5
    ## **isCentroid : bool
        The keyword argument isCentroid is used to determine the mass spectrum data structure:
        if set to False will assume profile mode and will attempt to peak pick
    ## **header_lines : int
        Line count to skip: should include the columns labels line
     ## **isThermoProfile : bool
        The keyword argument is_thermo_profile is used to change the number of expected columns to only m/z and intensity
        S/N and RP will be calculated based on the data. i.e signal to noise might not be accurate because the lack of noise data 
    # Attributes
    ----------
    ## _expected_columns : set
       The file has to have a least the values inside the set.

    For label translation see:
                corems.encapsulation.settings.input.InputSetting 
        or add your labels to the SettingsCoreMS.json file and parse the settings

    '''

    def __init__(self, file_location, isCentroid=True, analyzer='Unknown', instrument_label='Unknown',
                 sample_name=None, header_lines=0, isThermoProfile=False):

        self.file_location = Path(file_location) if isinstance(file_location, str) else file_location
        
        if not self.file_location.exists():
            raise FileExistsError("File does not exist: %s" % file_location)

        # (newline="\n")

        self.header_lines = header_lines

        if isThermoProfile:
            
            self._expected_columns = {Labels.mz, Labels.abundance}

        else:

            self._expected_columns = {Labels.mz, Labels.abundance, Labels.s2n, Labels.rp}


        self._delimiter = None

        self.isCentroid = isCentroid

        self._data_type = None

        self.analyzer = analyzer

        self.instrument_label = instrument_label

        self.sample_name = sample_name

        self._parameters = deepcopy(DataInputSetting())

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, instance_DataInputSetting):
        self._parameters = instance_DataInputSetting

    def set_parameter_from_json(self, parameters_path):
        self._parameters = load_and_set_parameters_class(
            'DataInput', self.parameters, parameters_path=parameters_path)

    @property
    def data_type(self):
        return self._data_type

    @data_type.setter
    def data_type(self, data_type):
        self._data_type = data_type

    @property
    def delimiter(self):
        return self._delimiter

    @delimiter.setter
    def delimiter(self, delimiter):
        self._delimiter = delimiter

    def encoding_detector(self, file_location):
        with file_location.open('rb') as rawdata:
            result = chardet.detect(rawdata.read(10000))
        return result['encoding']

    def set_data_type(self):

        if self.file_location.suffix == '.csv':

            self.data_type = 'txt'
            self.delimiter = ','

        elif self.file_location.suffix == '.txt':

            self.data_type = 'txt'
            self.delimiter = '\t'

        elif self.file_location.suffix == '.xlsx':

            self.data_type = 'excel'

        elif self.file_location.suffix == '.ascii':

            self.data_type = 'txt'
            self.delimiter = '  '

        elif self.file_location.suffix == '.pkl':

            self.data_type = 'dataframe'

        elif self.file_location.suffix == '.pks':

            self.data_type = 'pks'
            self.delimiter = '          '
            self.header_lines = 9
            
        else:
            raise TypeError(
                "Data type could not be automatically recognized for %s; please set data type and delimiter manually." % self.file_location.name)

    def get_dataframe(self):

        if not self.data_type or not self.delimiter:

            self.set_data_type()

        if isinstance(self.file_location, S3Path):
            # data = self.file_location.open('rb').read()
            data = BytesIO(self.file_location.open('rb').read())
        
        else:
            data = self.file_location

        if self.data_type == 'txt':
            
            dataframe = read_csv(data,  skiprows= self.header_lines, delimiter=self.delimiter,
                                 encoding=self.encoding_detector(self.file_location), engine='python')

        elif self.data_type == 'pks':
            
            names=["m/z", "I", "Scaled Peak Height", "Resolving Power", "Frequency", 'S/N']
            
            clean_data = []
            
            with self.file_location.open() as maglabfile:
                for i in  maglabfile.readlines()[8:-1]:
                    
                    clean_data.append(i.split())
            
            dataframe = DataFrame(clean_data, columns=names)
            
            #dataframe = read_csv(data,  skiprows= self.header_lines, delimiter=self.delimiter,
            #                     encoding=self.encoding_detector(self.file_location), engine='python',
            #                     header=None)
            #dataframe.columns = names
            #print(dataframe)
            #dataframe = read_csv(data,  skiprows= self.header_lines, delimiter=self.delimiter,
            #                     encoding=self.encoding_detector(self.file_location), engine='python')

        elif self.data_type == 'dataframe':

            dataframe = read_pickle(data)

        elif self.data_type == 'excel':

            dataframe = read_excel(data)

        else:

            raise TypeError('Data type %s is not supported' % self.data_type)

        return dataframe

    def load_settings(self, mass_spec_obj, output_parameters):

        import json
        import warnings

        settings_file_path = self.file_location.with_suffix('.json')

        if settings_file_path.exists():

            self._parameters = load_and_set_parameters_class(
                'DataInput', self._parameters, parameters_path=settings_file_path)

            load_and_set_parameters_ms(
                mass_spec_obj, parameters_path=settings_file_path)

        else:

            warnings.warn(
                "auto settings loading is enabled but could not locate the file:  %s. Please load the settings manually" % settings_file_path)

        # TODO this will load the setting from SettingCoreMS.json
        # coreMSHFD5 overrides this function to import the attrs stored in the h5 file
        #loaded_settings = {}
        #loaded_settings['MoleculaSearch'] = self.get_scan_group_attr_data(scan_index,  time_index, 'MoleculaSearchSetting')
        #loaded_settings['MassSpecPeak'] = self.get_scan_group_attr_data(scan_index,  time_index, 'MassSpecPeakSetting')

        #loaded_settings['MassSpectrum'] = self.get_scan_group_attr_data(scan_index, time_index, 'MassSpectrumSetting')
        #loaded_settings['Transient'] = self.get_scan_group_attr_data(scan_index, time_index, 'TransientSetting')

    def get_output_parameters(self, polarity, scan_index=0):

        # TODO pull attrs from json settings file in load_settings function MassSpecAttrs group and analyzer, instrument_label and sample_name
        from copy import deepcopy

        output_parameters = default_parameters(self.file_location)

        if self.isCentroid:
            output_parameters['label'] = Labels.corems_centroid
        else:
            output_parameters['label'] = Labels.bruker_profile

        output_parameters['analyzer'] = self.analyzer

        output_parameters['instrument_label'] = self.instrument_label

        output_parameters['sample_name'] = self.sample_name

        output_parameters["Aterm"] = None

        output_parameters["Bterm"] = None

        output_parameters["Cterm"] = None

        output_parameters["polarity"] = polarity

        '''scan_number and rt will be need to lc ms'''

        output_parameters["mobility_scan"] = 0

        output_parameters["mobility_rt"] = 0

        output_parameters["scan_number"] = scan_index

        output_parameters["rt"] = 0

        return output_parameters

    def clean_data_frame(self, dataframe):

        # print(dataframe.head())

        for column_name in dataframe.columns:

            expected_column_name = self.parameters.header_translate.get(
                column_name)
            if expected_column_name not in self._expected_columns:

                #dataframe = dataframe.drop(column_name, axis=1)
                del dataframe[column_name]
        # return dataframe

    def check_columns(self, header_labels):

        # print(header_labels)
        #print( self.parameters.header_translate.keys())
        #inverted_name_dict = {value: key for key, value in self.parameters.header_translate.items()}

        # print(inverted_name_dict)

        found_label = set()

        for label in header_labels:

            if not label in self._expected_columns:

                user_column_name = self.parameters.header_translate.get(label)

                if user_column_name in self._expected_columns:

                    found_label.add(user_column_name)

            else:
                found_label.add(label)

        not_found = self._expected_columns - found_label

        if len(not_found) > 0:

            raise Exception(
                "Please make sure to include the columns %s" % ', '.join(not_found))
