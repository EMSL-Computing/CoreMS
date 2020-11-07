__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

import sys
sys.path.append('.')
from pathlib import Path
from threading import Thread
from corems.mass_spectrum.input.massList import ReadCoremsMasslist
from corems.mass_spectra.factory.LC_Class import LCMSBase

class ReadCoremsMassSpectraText(ReadCoremsMasslist, Thread):
    
    def __init__(self, file_location, analyzer='Unknown', instrument_label='Unknown'):
        
        """
         # Parameters
		----------
        file_location: text,  pathlib.Path(), or s3path.S3Path 
            Path object from pathlib containing the file location
        """
        
        if  isinstance(file_location, str):
			# if obj is a string it defaults to create a Path obj, pass the S3Path if needed
            file_location = Path(file_location)

        if not file_location.exists():
            raise FileNotFoundError("%s not found" % file_location)
        
        if not file_location.suffix == '.corems':
            
            raise TypeError("%s is not a valid CoreMS file" % file_location)
        
        Thread.__init__(self)
        
        ReadCoremsMasslist.__init__(self, file_location)
        
        self.lcms = LCMSBase(self.file_location, analyzer=analyzer,instrument_label=instrument_label)

    def get_scans_filename(self):
        
        all_other = self.file_location.glob('*_scan*[!.json]')

        scans_filepath = [(file_path_obj.stem.split('scan')[1], file_path_obj)  for file_path_obj in all_other]
        
        scans_filepath.sort(key=lambda m: int(m[0]))

        return scans_filepath
        
    def set_filepath_datatype_and_delimiter(self, file_path_obj):

        self.file_location = file_path_obj

        if file_path_obj.suffix == '.pkl':

            self.data_type == 'dataframe'
            
        else:

            if file_path_obj.suffix == '.csv':
                self.data_type == 'txt'    
                self.delimiter = ','    
            
            elif file_path_obj.suffix == '.xlsx':
                self.data_type == 'excel'    
                self.delimiter = ',' 
            
            elif file_path_obj.suffix == '.txt':
                self.data_type == 'txt'    
                self.delimiter = '\t'
                print('WARNING using tab as delimiter')     
            else:
                raise NotImplementedError('%s data not yet supported ' % file_path_obj.suffix)
     
    def import_mass_spectra(self):
        
        list_rt, list_tic, list_scan = list(), list(), list()
        
        for scan_number, file_path_obj in self.get_scans_filename():
            
            self.set_filepath_datatype_and_delimiter(file_path_obj)
            
            mass_spec = self.get_mass_spectrum(int(scan_number))

            list_scan.append(int(scan_number))

            list_rt.append(mass_spec.rt)

            list_tic.append(mass_spec.tic)
            
            self.lcms.add_mass_spectrum(mass_spec)

        self.lcms.retention_time = list_rt
        self.lcms.tic_list = list_tic
        self.lcms.scans_number = list_scan
     
    def run(self):
        '''creates the lcms obj'''

        self.import_mass_spectra()
            
    def get_lcms_obj(self):
        
        if self.lcms:
            
            return self.lcms
        
        else:
            
            raise Exception("returning a empty lcms class")

if __name__ == "__main__":
    
    file_location = Path.cwd() / "NEG_ESI_SRFA_CoreMS.corems"
    
    all_other = file_location.glob('*_scan*[!.json]')

    all_json = file_location.glob('*_scan*.json')

    names_scans = [filename.stem.split('scan')  for filename in all_other]
    
    print(names_scans[0]) 
    #names_scans.sort(key=lambda m: int(m[1]))
    
    