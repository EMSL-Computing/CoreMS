from dataclasses import dataclass, field, asdict
from datetime import datetime
import hashlib
from json import dumps
from pathlib import Path
from typing import List
import yaml

from linkml_runtime.dumpers import json_dumper
import nmdc_schema.nmdc as nmdc
import oauthlib
import requests_oauthlib

@dataclass
class NomAnalysisActivity:
    codebase_url:str = "https://github.com/microbiomedata/enviroMS"
    cluster_name:str = "EMSL-RZR"
    nom_21T_instrument_name: str = "21T Agilent"
    nom_12T_instrument_name: str = "12T_FTICR_B"
    
@dataclass
class OmicsProcessing:
    nom_omics_processing_type:str = "Organic Matter Characterization"
    nom_omics_processing_description:str = "High resolution MS spectra only"
    nom_21T_instrument_name: str = "21T Agilent"
    nom_12T_instrument_name: str = "12T_FTICR_B"

@dataclass
class DataObject:
    nom_raw_data_object_type:str = "Direct Infusion FT ICR-MS Raw Data"
    nom_raw_data_object_description:str = "Raw 21T Direct Infusion Data"
    nom_dp_data_object_type:str = "FT ICR-MS Analysis Results"
    nom_dp_data_object_description:str = "EnviroMS FT ICR-MS natural organic matter workflow molecular formula assignment output details",
    
@dataclass
class NMDC_Types: 
    
    BioSample:str = "nmdc:BioSample"
    OmicsProcessing:str = "nmdc:OmicsProcessing"
    NomAnalysisActivity:str = "nmdc:NomAnalysisActivity"
    DataObject:str = "nmdc:DataObject"

@dataclass
class NMDC_Mint:
    
    schema_class: dict = field(default_factory= lambda: {
        'schema': None,
    })
    how_many:int = 1

    @property
    def __dict__(self):
        return asdict(self)

    @property
    def json(self):
        return dumps(self.__dict__)

def mint_nmdc_id(type:NMDC_Types, how_many:int = 1) -> List[str]: 
    
    config = yaml.safe_load(open('./config.yaml','r'))
    client = oauthlib.oauth2.BackendApplicationClient(client_id=config['client_id'])
    oauth = requests_oauthlib.OAuth2Session(client=client)
    
    token = oauth.fetch_token(token_url='https://api.microbiomedata.org/token',
                              client_id=config['client_id'], 
                              client_secret=config['client_secret'])

    nmdc_mint_url = "https://api.microbiomedata.org/pids/mint"
    
    payload = NMDC_Mint(type, how_many)
    
    #response = s.post(nmdc_mint_url, data=payload.json, )
    #list_ids = response.json()
    print(payload.json)
    response = oauth.post(nmdc_mint_url, data=payload.json)
    list_ids = response.json()
    print(list_ids)
    return list_ids

def get_data_object(file_path:Path, base_url:str, was_generated_by:str,
                data_object_type:str, description:str) -> nmdc.DataObject:
    
    nmdc_id = mint_nmdc_id({'id': NMDC_Types.DataObject})[0]

    data_dict = {
                'id': nmdc_id,
                "name": file_path.name,
                "file_size_bytes": file_path.stat().st_size,
                "md5_checksum": hashlib.md5(file_path.open('rb').read()).hexdigest(),
                "url": base_url + str(file_path.name),
                "was_generated_by": was_generated_by,
                "data_object_type": data_object_type,
                "description": description,
                "type": "nmdc:DataObject"
                } 

    data_object = nmdc.DataObject(**data_dict)

    return data_object

def get_omics_processing(file_path:Path, instrument_name:str, sample_id:str,
                raw_data_id:str, omics_type:str,  
                description:str, project_id:str) -> nmdc.OmicsProcessing:
    
    nmdc_id = mint_nmdc_id({'id': NMDC_Types.OmicsProcessing})[0]
    
    data_dict = {
                'id': nmdc_id,
                "name": file_path.stem,
                "instrument_name": instrument_name,
                "has_input": [sample_id],
                "has_output": [raw_data_id],
                "omics_type": {"has_raw_value": omics_type},
                "part_of": project_id,
                "processing_institution": "EMSL",
                "description": description,
                "type": "nmdc:OmicsProcessing"
                } 

    omicsProcessing = nmdc.OmicsProcessing(**data_dict)

    return omicsProcessing

def get_nom_analysis_activity(cluster_name:str, code_repository_url:str, 
                          raw_data_id:str, data_product_id:str, 
                          has_calibration:bool,  omics_processing_id:str, 
                          instrument_name:str) -> nmdc.NomAnalysisActivity:
    
    nmdc_id = mint_nmdc_id({'id': NMDC_Types.NomAnalysisActivity})[0]
    
    data_dict = {
                'id': nmdc_id,
                "execution_resource": cluster_name,
                "git_url": code_repository_url,
                "has_input": [raw_data_id],
                "has_output": [data_product_id],
                "has_calibration": has_calibration,
                "ended_at_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "started_at_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "used": instrument_name,
                "was_informed_by": omics_processing_id,
                "type": "nmdc:NomAnalysisActivity"
                } 

    nomAnalysisActivity = nmdc.NomAnalysisActivity(**data_dict)

    return nomAnalysisActivity

def start_nmdc_database() -> nmdc.Database:
    return nmdc.Database()

def create_nmdc_metadata(raw_data_path:Path, data_product_path:Path, base_url:str,
                         nmdc_study_id:str, nom_metadata_db:nmdc.Database,
                         biosample_id=None):

    if not biosample_id:
        biosample_id = mint_nmdc_id({'id': NMDC_Types.BioSample})[0]
        ''' needs to finish the logic for creating biosamples, this will fail because it is missing some required fields'''
        bioSample =  nmdc.BioSample(id=biosample_id)
        biosample_id = bioSample.id

    omicsProcessing = get_omics_processing(raw_data_path,
                                           OmicsProcessing.nom_12T_instrument_name,
                                           biosample_id, None, 
                                           OmicsProcessing.nom_omics_processing_type,
                                           OmicsProcessing.nom_omics_processing_description,
                                           nmdc_study_id
                                           )
    
    rawDataObject = get_data_object(raw_data_path, base_url + 'nom/grow/raw/', 
                                    was_generated_by=omicsProcessing.id, 
                                    data_object_type =DataObject.nom_raw_data_object_type,
                                    description =DataObject.nom_raw_data_object_description)
    
    nomAnalysisActivity = get_nom_analysis_activity(NomAnalysisActivity.cluster_name,
                                                NomAnalysisActivity.codebase_url,
                                                rawDataObject.id, None, False, 
                                                omicsProcessing.id,
                                                NomAnalysisActivity.nom_12T_instrument_name)

    dataProductDataObject = get_data_object(data_product_path, base_url + 'nom/grow/results/', 
                                    was_generated_by=nomAnalysisActivity.id, 
                                    data_object_type =DataObject.nom_dp_data_object_type,
                                    description =DataObject.nom_dp_data_object_description)
    
    #circular dependencies : great! 
    nomAnalysisActivity.has_input = rawDataObject.id
    nomAnalysisActivity.has_output = dataProductDataObject.id
    omicsProcessing.has_output = rawDataObject.id

    nom_metadata_db.data_object_set.append(rawDataObject)
    nom_metadata_db.nom_analysis_activity_set.append(nomAnalysisActivity)
    nom_metadata_db.omics_processing_set.append(omicsProcessing)
    nom_metadata_db.data_object_set.append(dataProductDataObject)

def dump_nmdc_database(ndmc_database:nmdc.Database, output_filepath:str):
    
    json_dumper.dump(ndmc_database, output_filepath)
    