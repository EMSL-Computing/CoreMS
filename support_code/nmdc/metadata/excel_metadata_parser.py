from dataclasses import asdict, dataclass, field
import hashlib
import json
import os
from pathlib import Path
from typing import Dict, List, Tuple
import shutil

from openpyxl import load_workbook

@dataclass
class BioSample:
    
    uid: str
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
    
    uid: str
    has_input: List[str]
    has_output: List[str]
    part_of: List[str] = field(default_factory=list)
    omics_type: str = "Organic Matter"
    processing_institution: str =  "Environmental Molecular Science Laboratory"
    type: str = "nmdc:OmicsProcessing"

@dataclass
class DataObject:
    
    uid: str
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
    uid: str
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
                uid = ids[x].value.replace('EMSL:', 'emsl:')
                
                sample_data_map = SampleDataMap(
                                                name=sample_name[x].value,
                                                uid=uid,
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
        geo_loc_name = sample_metadata['F']
        lat_lon = sample_metadata['G']
        collection_date = sample_metadata['H']
        env_broad_scale = sample_metadata['M']
        env_local_scale = sample_metadata['N']
        env_medium = sample_metadata['O']

        biosample_objs = {}
        
        for x in range(1, len(sample_names)):
            
            name = sample_names[x].value

            if name in self.biosample_dataset_map.keys(): 
                
                biosample = BioSample(uid=self.biosample_dataset_map.get(name).uid,
                                        env_broad_scale= env_broad_scale[x].value,
                                        env_local_scale= env_local_scale[x].value,
                                        env_medium= env_medium[x].value,
                                        env_package= env_packages[x].value,
                                        geo_loc_name= geo_loc_name[x].value,
                                        lat_lon= lat_lon[x].value,
                                        collection_date= str(collection_date[x].value),
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
        
        
        nom_biosample_set_path = Path("nom_biosample_set.json")

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
                    has_input = [biosample_obj.uid],
                    has_output = ['emsl:output_{}'.format(dataset_id)],
                    uid = ['emsl:{}'.format(dataset_id)],
                    part_of = ["gold:Gs0110138"]
                )
                
                omics_processing_set.get("omics_processing_set").append(asdict(omics_processing_obj))
        
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
                uid = 'emsl:output_{}'.format(dataset_id),
                name = rawdata_path.name,
                file_size_bytes = rawdata_path.stat().st_size,
                md5_checksum = hashlib.md5(rawdata_path.open('rb').read()).hexdigest(),
                url = "{}/{}/{}/{}".format("https://nmdcdemo.emsl.pnnl.gov", "nom", "raw", rawdata_path.name),
                was_generated_by = 'emsl:{}'.format(dataset_id),
                
            )

            data_obj_set.get('data_object_set').append(asdict(dataobj))

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
                data_obj_set.get('nom_analysis_activity_set').append(analysis_obj)

        
        nom_rawdata_objset_path = Path("spruce_nom_analysis_activity_set.json")
        with nom_rawdata_objset_path.open('w') as nom_json_out:
            nom_json_out.write(json.dumps(data_obj_set, indent=1))    

if __name__ == '__main__':
    
    metadata_sheet_path = Path("/Users/eber373/OneDrive - PNNL/Documents/Data/FT_ICR_MS/Spruce_Data/SPRUCE_Schadt-Wilson_NMDC_SampleMetadata_Subset.xlsx")
    
    results_path = Path("/Users/eber373/OneDrive - PNNL/Documents/Data/FT_ICR_MS/Spruce_Data/results")

    rawdata_path_dir = Path("/Users/eber373/OneDrive - PNNL/Documents/Data/FT_ICR_MS/Spruce_Data/raw")

    analysis_activity_path = Path("/Users/eber373/OneDrive - PNNL/Documents/Data/FT_ICR_MS/Spruce_Data/metadata")

    metadata_mapping = Metadata_Mapping(metadata_sheet_path, results_path)
    
    #metadata_mapping.dump_biosample_set()
    #metadata_mapping.dump_omics_processing_set()    
    #metadata_mapping.dump_data_objs_raw_data(rawdata_path_dir)
    metadata_mapping.dump_analysis_activity_set(analysis_activity_path)