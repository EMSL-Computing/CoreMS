
import hashlib
import json
import os
from pathlib import Path

from minio import Minio

import sys
sys.path.append("./")

from doc.examples.nmdc.NMDC_Metadata import DMS_Mapping


def fix_data_registration_nom():
    
    #registration_path = Path("db/gcms_metabolomics_data_products.json")
    registration_path = Path("results/ftms_nom_data_products.json")

    metadata_path = Path("results/ftms_nom_metadata_products.json")

    metabolomics_data_path = Path(registration_path)
    data_dir = Path("db/results")

    if metabolomics_data_path.exists():

        url_list = []
        with metabolomics_data_path.open('r') as metabolomics_data_products:
            products = json.load(metabolomics_data_products)
            for product in products:

                #name = product.get("name")
                #out_file_path = data_dir / name
                #md5_checksum = hashlib.md5(out_file_path.open('rb').read()).hexdigest()
                #product["md5_checksum"] = md5_checksum
                #product["id"] = "{}:{}".format("nmdc", md5_checksum)
                #product["file_size_bytes"] = out_file_path.stat().st_size
                #product["has_input"] = [product["has_input"]]
                #product["has_output"] = [product["has_output"]]
                
                url_list.append(product.get("url").replace('.csv', ".json").replace('results', 'metadata'))

        #with registration_path.open('w') as metabolomics_data_products:
        #    metabolomics_data_products.write(json.dumps(products, indent=1))

        with metadata_path.open('w') as metabolomics_data_products:
            metabolomics_data_products.write(json.dumps(url_list, indent=1))

def fix_data_registration():
    
    
    registration_path = Path("results/ftms_nom_data_products.json")
    #metadata_path = Path("db/gcms_metabolomics_metadata_products.json")

    metabolomics_data_path = Path(registration_path)
    data_dir = Path("results/")

    if metabolomics_data_path.exists():

        url_list = []
        with metabolomics_data_path.open('r') as metabolomics_data_products:
            products = json.load(metabolomics_data_products)
            for product in products:

                #name = product.get("name")
                #out_file_path = data_dir / name
                #md5_checksum = hashlib.md5(out_file_path.open('rb').read()).hexdigest()
                #product["md5_checksum"] = md5_checksum
                #product["id"] = "{}:{}".format("nmdc", md5_checksum)
                product["type"] = "nmdc:DataObject"
                #product["file_size_bytes"] = out_file_path.stat().st_size
                #url_list.append(product.get("url").replace('.csv', ".json").replace('results', 'metadata'))

        with registration_path.open('w') as metabolomics_data_products:
            metabolomics_data_products.write(json.dumps(products, indent=1))

        #with metadata_path.open('w') as metabolomics_data_products:
        #    metabolomics_data_products.write(json.dumps(url_list, indent=1))

def fix_metadata_registration():
    
    #from copy import deepcopy
    cvs_file_dir = Path("results/data")
    metadata_file_dir = Path("results/metadata")
    newmetada_file_dir = Path("results/fix_metadata")
    
    dms_file_path = Path("db/NOM Data to Process.xlsx")
    
    dms_dataset_mapping = DMS_Mapping(dms_file_path).get_mapping()

    directory = cvs_file_dir
    for filename in os.listdir(directory):

        csv_path = Path(cvs_file_dir/filename)
        #md5_checksum = hashlib.md5(csv_path.open('rb').read()).hexdigest()
        
        metadata_path = metadata_file_dir / Path(csv_path.name).with_suffix('.json')
        
        with metadata_path.open('r') as metabolomics_data_products:
            
            metadata_dict = json.load(metabolomics_data_products)
            
            metadata_dict["has_output"] = [metadata_dict["has_output"]]

            metadata_dict["has_input"] = [metadata_dict["has_input"]]

            print(metadata_path.stem)    
            mapping = dms_dataset_mapping.get(metadata_path.stem)

            metadata_dict["used"] = mapping.get("instrument_name")

            #id_dict = {"id": deepcopy(metadata_dict["activity_id"]) }

            #del metadata_dict["activity_id"]

            #id_dict.update(metadata_dict)

        new_metadata = newmetada_file_dir / metadata_path.name
        
        with new_metadata.open('w') as metabolomics_data_products:

            metabolomics_data_products.write(json.dumps(metadata_dict, indent=1))

def upload_data():
    
    minio = Minio(
            "admin.nmdcdemo.emsl.pnl.gov",
            access_key=os.environ.get("dMUv0sYh3K"),
            secret_key=os.environ.get("xkHUzeZ8KlXaDSXI5bfoheWqqFddsLYjDE8784yB"),
            secure=True
                )

    filepath = Path("results/ftms_nom_data_products.json")
    s3_path = "{}/{}".format("data", filepath.name)
    
    minio.fput_object('nom', s3_path, str(filepath))
    
if __name__ == "__main__":
    
    #fix_metadata_registration()
    #fix_data_registration_nom()
    #upload_data()
    fix_data_registration()
    