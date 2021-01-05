''' NMDC Schema

 "ChemicalEntity": {
         "additionalProperties": false,
         "description": "An atom or molecule that can be represented with a chemical formula. Include lipids, glycans, natural products, drugs. There may be different terms for distinct acid-base forms, protonation states",
         "properties": {
            "chemical_formula": {
               "description": "A generic grouping for miolecular formulae and empirican formulae",
               "type": "string"
            },
            "has_raw_value": {
               "description": "The value that was specified for an annotation in raw form, i.e. a string. E.g. \"2 cm\" or \"2-4 cm\"",
               "type": "string"
            },
            "inchi": {
               "type": "string"
            },
            "inchi_key": {
               "type": "string"
            },
            "smiles": {
               "description": "A string encoding of a molecular graph, no chiral or isotopic information. There are usually a large number of valid SMILES which represent a given structure. For example, CCO, OCC and C(O)C all specify the structure of ethanol.",
               "items": {
                  "type": "string"
               },
               "type": "array"
            },
            "term": {
               "$ref": "#/definitions/OntologyClass",
               "description": "pointer to an ontology class"
            },
            "was_generated_by": {
               "type": "string"
            }
         },
         "required": [
            "inchi"
         ],
         "title": "ChemicalEntity",
         "type": "object"
      },

"MetaboliteQuantification": {
         "additionalProperties": false,
         "description": "This is used to link a metabolomics analysis workflow to a specific metabolite",
         "properties": {
            "highest_similarity_score": {
               "description": "TODO: Yuri to fill in",
               "type": "number"
            },
            "metabolite_quantified": {
               "description": "the specific metabolite identifier",
               "type": "string"
            }
         },
         "required": [],
         "title": "MetaboliteQuantification",
         "type": "object"
      },

"MetabolomicsAnalysisActivity": {
         "additionalProperties": false,
         "description": "",
         "properties": {
            "activity_id": {
               "type": "string"
            },
            "ended_at_time": {
               "type": "string"
            },
            "execution_resource": {
               "description": "Example: NERSC-Cori",
               "type": "string"
            },
            "git_url": {
               "description": "Example: https://github.com/microbiomedata/mg_annotation/releases/tag/0.1",
               "type": "string"
            },
            "has_input": {
               "description": "An input to a process.",
               "items": {
                  "type": "string"
               },
               "type": "array"
            },
            "has_metabolite_quantifications": {
               "items": {
                  "$ref": "#/definitions/MetaboliteQuantification"
               },
               "type": "array"
            },
            "has_output": {
               "description": "An output biosample to a processing step",
               "items": {
                  "type": "string"
               },
               "type": "array"
            },
            "started_at_time": {
               "type": "string"
            },
            "used": {
               "description": "The instrument used to collect the data used in the analysis",
               "type": "string"
            },
            "was_associated_with": {
               "$ref": "#/definitions/Agent"
            },
            "was_informed_by": {
               "type": "string"
            }
         },
         "required": [
            "activity_id"
         ],
         "title": "MetabolomicsAnalysisActivity",
         "type": "object"
      },
'''
from datetime import datetime
from pathlib import Path
import hashlib
import uuid
import json

from .dms_mapping import get_mapping


def save_nmdc_metadata(in_file_path, gcms_obj, out_file_path, dataset_mapping):

   in_file_path = Path(in_file_path)
   out_file_path = Path(out_file_path)

   now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
   data_id = get_dataid(in_file_path, dataset_mapping)
   output_id = get_md5(out_file_path)
   was_informed_by = get_was_informed_by(in_file_path, dataset_mapping)
   activity_id = "{}:{}".format("nmdc", uuid.uuid4().hex)

   metabolomics_analysis_activity = {
        
        "activity_id": activity_id,
        "ended_at_time" : now,
        "execution_resource" : "EMSL-RZR",
        "git_url" : "https://github.com/microbiomedata/metaMS",
        "has_input": data_id,
        "has_output": output_id,
        "started_at_time": now,
        "used": "Agilent_GC_MS",
        "was_informed_by": was_informed_by,
        "has_metabolite_quantifications": get_metabolites_objs(gcms_obj)
    }
   
   with open(out_file_path.with_suffix('.json'), 'w') as metadata_output:
      metadata_output.write(json.dumps(metabolomics_analysis_activity, indent=1))    
   
   add_metabolomics_data_product(output_id, out_file_path, activity_id)

def add_metabolomics_data_product(output_id, out_file_path, activity_id):
      
      data_obj = [{

         "id": output_id,
         "name": out_file_path.name,
         "description": "MetaMS GC-MS metabolomics output detail CSV file",
         "file_size_bytes": out_file_path.stat().st_size,
         "md5_checksum": hashlib.md5(out_file_path.open('rb').read()).hexdigest(),
         "url": "https://data.corems.emsl.pnnl.gov",
         "was_generated_by": activity_id
      }]

      metabolomics_data_path = Path("gcms_metabolomics_data_products.json")
      
      if metabolomics_data_path.exists():
      
         with  metabolomics_data_path.open('r') as metabolomics_data_products:
            products = json.load(metabolomics_data_products)
            products.extend(data_obj)

      else:
            products = data_obj
      
      with metabolomics_data_path.open('w') as metabolomics_data_products:
         metabolomics_data_products.write(json.dumps(products, indent=1))  
   
def get_dataid(in_file_path, dataset_mapping):
    
    data_id = dataset_mapping.get(in_file_path.stem).get("data_id")
    return "{}:{}_{}".format("emsl", "output", data_id)

def get_was_informed_by(in_file_path, dataset_mapping):
    
    data_id = dataset_mapping.get(in_file_path.stem).get("data_id")
    return "{}:{}".format("emsl", data_id)

def get_md5(out_file_path):
    bytes_io = out_file_path.open('rb').read()

    return "{}:{}".format("nmdc", hashlib.md5(bytes_io).hexdigest())

def get_metabolites_objs(gcms_obj):
    
    all_metabolites = []
    
    for metabolite in gcms_obj.metabolites_data:
        
        metabolite_quantification = {}
        metabolite_quantification["highest_similarity_score"] = metabolite.get("highest_similarity_score")
        metabolite_quantification["metabolite_quantified"] = metabolite.get("inchi")
        metabolite_quantification["inchi"] = metabolite.get("inchi")
        metabolite_quantification["inchi_key"] = metabolite.get("inchi_key")
        metabolite_quantification["chebi"] = metabolite.get("chebi")
        metabolite_quantification["kegg"] = metabolite.get("kegg")
        all_metabolites.append(metabolite_quantification)
    
    return all_metabolites

def create_nmdc_metadata(in_file_path, gcms_obj, out_file_path):

   dms_file = 'db/GC-MS Metabolomics Experiments to Process Final.xlsx'
   dataset_mapping = get_mapping(dms_file)
   #print(dataset_mapping.get('GCMS_Blank_08_GC-01_20150622'))
   save_nmdc_metadata(in_file_path, gcms_obj, out_file_path, dataset_mapping)

if __name__ == "__main__":
   pass 
   dms_file = 'db/GC-MS Metabolomics Experiments to Process Final.xlsx'
   get_mapping(dms_file)
   #print(dataset_mapping.get('GCMS_Blank_08_GC-01_20150622'))
   create_nmdc_metadata(in_file_path, gcms_obj, out_file_path, dataset_mapping)