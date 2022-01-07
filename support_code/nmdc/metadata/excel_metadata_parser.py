from dataclasses import asdict, dataclass, field
import hashlib
import json
import os
from pathlib import Path
from typing import Dict, List, Tuple
import glob
from datetime import datetime

from openpyxl import load_workbook

import jsonschema
from jsonschema import validate

@dataclass
class BioSample:
    
    id: str
    env_broad_scale: float
    env_local_scale: float
    env_medium: float
    env_package: str
    geo_loc_name: str
    lat_lon: str
    collection_date: str
    name: str

@dataclass
class OmicsProcessing:
    
    id: str
    has_input: List[str]
    has_output: List[str]
    omics_type: str
    part_of: List[str] = field(default_factory=list)
    processing_institution: str =  "Environmental Molecular Science Laboratory"
    type: str = "nmdc:OmicsProcessing"

@dataclass
class DataObject:
    
    id: str
    name: str
    file_size_bytes: int
    md5_checksum: str
    url: str
    was_generated_by: str
    data_object_type: str = "Direct Infusion FT ICR-MS Raw Data"
    description: str = "Raw 21T Direct Infusion Data"
    type: str = "nmdc:DataObject"
    

@dataclass
class SampleDataMap:
    name: str
    id: str
    dataset_names: List[str]
    dataset_ids: List[int]

class Metadata_Mapping():

    def __init__(self, file_path, results_path):

        self.wb = load_workbook(filename=file_path)
        self.results_dir = results_path
        self.biosample_dataset_map, self.dataset_biosample_map, self.datasetname_id_map = self.get_biosample_dataset_map()
        self.biosamples = self.create_biosamples()

    def get_biosample_dataset_map(self) -> Tuple[Dict[str,SampleDataMap], Dict[str,str], Dict[str,str]] :
        
        first_sheet = 'ID Map'

        id_map = self.wb[first_sheet]

        ids = id_map['A']
        sample_name = id_map['B']
        methanol_dataset = id_map['N']
        methanol_dataset_id = id_map['O']
        cloroform_dataset = id_map['P']
        cloroform_dataset_id = id_map['Q']
        
        sample_name_id_map = {}
        
        datasetname_sampe_map = {}

        datasetname_id_map = {}

        for x in range(1, len(ids)):

            if ids[x].value:
                dataset_names = methanol_dataset[x].value.replace('\n', '').split(";") + cloroform_dataset[x].value.replace('\n', '').split(";")
                dataset_names = list(filter(None, dataset_names))
                dataset_ids = methanol_dataset_id[x].value.replace('\n', '').split(";") + cloroform_dataset_id[x].value.replace('\n', '').split(";")
                dataset_ids = list(filter(None, dataset_ids))
                id = ids[x].value.replace('EMSL:', 'emsl:')
                
                sample_data_map = SampleDataMap(
                                                name=sample_name[x].value,
                                                id=id,
                                                dataset_names=dataset_names,
                                                dataset_ids=dataset_ids)

                sample_name_id_map[sample_name[x].value] = sample_data_map
                
                #print(sample_data_map)
                
                for i, dataset_name in enumerate(dataset_names):
                    
                    datasetname_sampe_map[dataset_name] = sample_name[x].value

                    datasetname_id_map[dataset_name] = dataset_ids[i]
                #print(ids[x].value.replace('EMSL:', 'emsl:'), sample_name[x].value)

        return sample_name_id_map, datasetname_sampe_map, datasetname_id_map

    def create_biosamples(self) -> Dict[str, BioSample]:

        first_sheet = 'Metadata'

        sample_metadata = self.wb[first_sheet]

        sample_names = sample_metadata['B']
        env_packages = sample_metadata['C']
        geo_loc_names = sample_metadata['F']
        lat_lons = sample_metadata['G']
        collection_dates = sample_metadata['H']
        env_broad_scales = sample_metadata['M']
        env_local_scales = sample_metadata['N']
        env_mediums = sample_metadata['O']

        biosample_objs = {}
        
        for x in range(1, len(sample_names)):
            
            name = sample_names[x].value

            if name in self.biosample_dataset_map.keys(): 
                
                env_broad_scale = {"has_raw_value": env_broad_scales[x].value}
                env_local_scale = {"has_raw_value": env_local_scales[x].value}
                env_medium = {"has_raw_value": env_mediums[x].value}
                env_package = {"has_raw_value": env_packages[x].value}
                geo_loc_name = {"has_raw_value": geo_loc_names[x].value}
                lat_lon = {"has_raw_value": lat_lons[x].value.replace(',', ''), 
                           "latitude": float(lat_lons[x].value.split(',')[0]),
                           "longitude": float(lat_lons[x].value.split(',')[1].strip()) }
                collection_date = {"has_raw_value": str(collection_dates[x].value)}

                biosample = BioSample(id=self.biosample_dataset_map.get(name).id,
                                        env_broad_scale= env_broad_scale,
                                        env_local_scale= env_local_scale,
                                        env_medium= env_medium,
                                        env_package= env_package,
                                        geo_loc_name= geo_loc_name,
                                        lat_lon= lat_lon,
                                        collection_date= collection_date,
                                        name= name)
                
                biosample_objs[name] = biosample
        
        return biosample_objs

    def get_processed_dataset_name(self):
        '''list data products processed, assumes directory cointains csv data table'''
        for filename in os.listdir(self.results_dir):

            dataset_name = (filename.replace('.csv', ''))

            yield  dataset_name

    def dump_biosample_set(self):

        biosamples_processed = set()
        
        for dataset_name in self.get_processed_dataset_name():
            
            if dataset_name in self.dataset_biosample_map.keys():

                biosamples_processed.add(self.dataset_biosample_map.get(dataset_name))

        
        biosample_set = {"biosample_set": []}
        
        for biosample_name in biosamples_processed:
        
            biosample_obj = self.biosamples.get(biosample_name)
            biosample_set.get("biosample_set").append(asdict(biosample_obj))
        
        valid_json, message = validate_json(biosample_set)

        if valid_json:
            nom_biosample_set_path = Path("spruce_nom_biosample_set.json")

            with nom_biosample_set_path.open('w') as nom_json_out:
                nom_json_out.write(json.dumps(biosample_set, indent=1))  

    def dump_omics_processing_set(self):

        omics_processing_set = {"omics_processing_set": []}
        
        for dataset_name in self.get_processed_dataset_name():
            
            dataset_id = self.datasetname_id_map.get(dataset_name)
            biosample_name = self.dataset_biosample_map.get(dataset_name)
            biosample_obj = self.biosamples.get(biosample_name)
            
            if biosample_obj:

                omics_processing_obj = OmicsProcessing(
                    omics_type={'has_raw_value':  "Organic Matter"},
                    has_input = [biosample_obj.id],
                    has_output = ['emsl:output_{}'.format(dataset_id)],
                    id = 'emsl:{}'.format(dataset_id),
                    part_of = ["gold:Gs0110138"]
                )
                
                omics_processing_set.get("omics_processing_set").append(asdict(omics_processing_obj))
        
        valid_json, message = validate_json(omics_processing_set)
        
        if valid_json:
        
            nom_omics_processing_set_path = Path("spruce_nom_omics_processing_set.json")

            with nom_omics_processing_set_path.open('w') as nom_json_out:
                nom_json_out.write(json.dumps(omics_processing_set, indent=1))  
    
    def dump_data_objs_raw_data(self, rawdata_path_dir):

        data_obj_set = {"data_object_set": []}
        
        for dataset_name in self.get_processed_dataset_name():
            
            dataset_id = self.datasetname_id_map.get(dataset_name)

            rawdata_path = rawdata_path_dir / dataset_name
            
            rawdata_path = rawdata_path.with_suffix('.raw')
            
            print(rawdata_path)

            dataobj = DataObject(
                id = 'emsl:output_{}'.format(dataset_id),
                name = rawdata_path.name,
                file_size_bytes = rawdata_path.stat().st_size,
                md5_checksum = hashlib.md5(rawdata_path.open('rb').read()).hexdigest(),
                url = "{}/{}/{}/{}".format("https://nmdcdemo.emsl.pnnl.gov", "nom", "raw", rawdata_path.name),
                was_generated_by = 'emsl:{}'.format(dataset_id),
                
            )

            data_obj_set.get('data_object_set').append(asdict(dataobj))

        valid_json, message = validate_json(data_obj_set)
        
        if valid_json:
        
            nom_rawdata_objset_path = Path("spruce_nom_raw_data_object_set.json")


            with nom_rawdata_objset_path.open('w') as nom_json_out:
                nom_json_out.write(json.dumps(data_obj_set, indent=1))      
        
    def dump_analysis_activity_set(self, analysis_activity_path):

        data_obj_set = {"nom_analysis_activity_set": []}

        for dataset_name in self.get_processed_dataset_name():
            
            analysis_obj_path = analysis_activity_path / dataset_name
            
            analysis_obj_path = analysis_obj_path.with_suffix('.json')

            with analysis_obj_path.open('r') as tmp_json:
                analysis_obj =  json.load(tmp_json)
                analysis_obj['has_calibration'] = 'false'
                data_obj_set.get('nom_analysis_activity_set').append(analysis_obj)

        valid_json, message = validate_json(data_obj_set)
        
        if valid_json:

            nom_rawdata_objset_path = Path("spruce_nom_analysis_activity_set.json")
            with nom_rawdata_objset_path.open('w') as nom_json_out:
                nom_json_out.write(json.dumps(data_obj_set, indent=1))    

    def validate_data_products_set(self, filepath):
        
        with open(filepath, 'r') as file:
            data = json.load(file)
            validate_json(data)

    def create_registration_objs(self, file_paths):


        for name in glob.glob(str(file_paths) + '/*.json'):
            path = Path(name)
            url = "{}/{}/{}/{}/{}".format("https://nmdcdemo.emsl.pnnl.gov", "nom", "registration", 'spruce', path.name),
            
            registration_obj = {
                "description": "Spruce NOM ",
                "name": path.name,
                 "access_methods": [ 
                    {"access_url": {'url': url}
                    }
                 ],
                 'checksums' : [
                     {
                    "checksum": hashlib.sha256(path.open('rb').read()).hexdigest(),
                    "type": "sha256"
                    }
                 ],
                 'size': path.stat().st_size,
                 "created_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
            
            
            print(registration_obj)
            
            nom_rawdata_objset_path = Path(str(path.stem)+"_registration.json")
            
            with nom_rawdata_objset_path.open('w') as nom_json_out:
                nom_json_out.write(json.dumps(registration_obj, indent=1))    

def get_schema():
    import requests
    """This function loads the given schema available"""
    url = 'https://raw.githubusercontent.com/microbiomedata/nmdc-schema/main/jsonschema/nmdc.schema.json'
    resp = requests.get(url)
    schema = json.loads(resp.text)
    
    #schema = json.loads(data)
    
    return schema

def validate_json(json_data):
    """REF: https://json-schema.org/ """
    # Describe what kind of json you expect.
    execute_api_schema = get_schema()

    try:
        validate(instance=json_data, schema=execute_api_schema)
    except jsonschema.exceptions.ValidationError as err:
        print(err)
        error_file = Path('error.txt')
        with open('error.txt', 'w') as file:
            file.write(str(err))
        err = "Given JSON data is InValid"
        return False, err

    message = "Given JSON data is Valid"
    return True, message

if __name__ == '__main__':
    
    metadata_sheet_path = Path("/Users/eber373/OneDrive - PNNL/Documents/Data/FT_ICR_MS/Spruce_Data/SPRUCE_Schadt-Wilson_NMDC_SampleMetadata_Subset.xlsx")
    
    results_path = Path("/Users/eber373/OneDrive - PNNL/Documents/Data/FT_ICR_MS/Spruce_Data/results")

    rawdata_path_dir = Path("/Users/eber373/OneDrive - PNNL/Documents/Data/FT_ICR_MS/Spruce_Data/raw")

    analysis_activity_path = Path("/Users/eber373/OneDrive - PNNL/Documents/Data/FT_ICR_MS/Spruce_Data/metadata")

    data_products_path = Path("spruce_ftms_nom_data_products_set.json")

    registration_path = Path("/Users/eber373/OneDrive - PNNL/Documents/Data/FT_ICR_MS/Spruce_Data/registration")

    metadata_mapping = Metadata_Mapping(metadata_sheet_path, results_path)
    
    #metadata_mapping.dump_biosample_set()
    #metadata_mapping.dump_omics_processing_set()    
    #metadata_mapping.dump_data_objs_raw_data(rawdata_path_dir)
    #metadata_mapping.dump_analysis_activity_set(analysis_activity_path)
    metadata_mapping.validate_data_products_set(data_products_path)
    #metadata_mapping.create_registration_objs(registration_path)