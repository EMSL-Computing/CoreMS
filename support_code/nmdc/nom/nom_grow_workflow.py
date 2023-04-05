from dataclasses import dataclass
import json
import os
from pathlib import Path
import shutil
from zipfile import ZipFile

from openpyxl import load_workbook

from nom_workflow import run_nmdc_workflow
import nmdc_metadata_gen

@dataclass
class EMSL_Metadata:
    data_path: str
    dms_dataset_id: str
    myemsl_link: str
    sample_name: str
    nmdc_biosample_id: str
    gold_name: str
    nmdc_study: str
    env_broad_scale: str
    env_local_scale: str
    env_medium: str
    img_id: str
    gold_biosample_id: str

def parse_metadata(metadata_file_path:Path) -> EMSL_Metadata:
        
        wb = load_workbook(filename=metadata_file_path)

        first_sheet = wb.sheetnames[0]

        full_list_worksheet = wb[first_sheet]

        data_name = full_list_worksheet['A']
        dms_dataset_id = full_list_worksheet['B']
        myemsl_link = full_list_worksheet['C']
        sample_name = full_list_worksheet['D']
        nmdc_biosample_id = full_list_worksheet['E']
        gold_name = full_list_worksheet['F']
        nmdc_study = full_list_worksheet['G']
        env_broad_scale = full_list_worksheet['H']
        env_local_scale = full_list_worksheet['I']
        env_medium = full_list_worksheet['J']
        img_id = full_list_worksheet['K']
        gold_biosample_id = full_list_worksheet['L']
                                
        for x in range(1, len(full_list_worksheet['A'])):

            metadata = EMSL_Metadata(data_path = Path(data_name[x].value), 
                                 dms_dataset_id = dms_dataset_id[x].value,
                                 myemsl_link = myemsl_link[x].value,
                                 sample_name = sample_name[x].value,
                                 nmdc_biosample_id = nmdc_biosample_id[x].value,
                                 gold_name = gold_name[x].value,
                                 nmdc_study = nmdc_study[x].value,
                                 env_broad_scale = env_broad_scale[x].value,
                                 env_local_scale = env_local_scale[x].value,
                                 env_medium = env_medium[x].value,
                                 img_id = img_id[x].value,
                                 gold_biosample_id = gold_biosample_id[x].value
                                )
            
            yield metadata

def run_nom_nmdc_data_processing():
    
    file_ext = '.d' 
    data_dir = Path("/Users/eber373/Library/CloudStorage/OneDrive-PNNL/Desktop/data/nmdc_data/GROW /Freshwater/")
    metadata_file_path = Path("/Users/eber373/Library/CloudStorage/OneDrive-PNNL/Desktop/data/nmdc_data/GROW /Freshwater/surface_water_metadata.xlsx")

    raw_dir_zip = data_dir / Path("raw_zip/")
    raw_dir_zip.mkdir(parents=True, exist_ok=True)

    results_dir = data_dir / Path("results/")
    results_dir.mkdir(parents=True, exist_ok=True)
    
    failed_files = results_dir / "nom_failed_files.json"
    
    registration_dir = data_dir / 'registration'
    registration_dir.mkdir(parents=True, exist_ok=True)
    registration_file = registration_dir / 'surface_water.json'
    
    field_strength = 12
    
    ref_calibration_path = Path("db/Hawkes_neg.ref")
    failed_list = []
    
    nmdc_database = nmdc_metadata_gen.start_nmdc_database()

    for each_data in parse_metadata(metadata_file_path):

        raw_file_path = data_dir / each_data.data_path / each_data.data_path.with_suffix(file_ext)    

        print(raw_file_path)
        
        issue, ms = run_nmdc_workflow((raw_file_path, ref_calibration_path, field_strength))
        
        try:
            if ms:

                if raw_file_path.suffix =='.d':
                    #raw_file_to_upload_path = zip_d_folder(raw_file_path)    
                    raw_file_to_upload_path = Path(raw_dir_zip / raw_file_path.stem)
                    shutil.make_archive(raw_file_to_upload_path , 'zip', raw_file_path)
                else:
                    raw_file_to_upload_path = raw_file_path

                result_file_name = Path(raw_file_path.name)
                output_file_path = results_dir / result_file_name.with_suffix('.csv')
                ms.to_csv(output_file_path, write_metadata=False)
                
                nmdc_metadata_gen.create_nmdc_metadata(raw_file_to_upload_path.with_suffix('.zip'), 
                                                        output_file_path,
                                                        "https://nmdcdemo.emsl.pnnl.gov/", 
                                                        each_data.nmdc_study,
                                                        nmdc_database,
                                                        each_data.nmdc_biosample_id)

            else:
                
                print(issue)
                failed_list.append(raw_file_path)

        except Exception as inst:

            print(type(inst))    # the exception instance
            print(inst.args)     # arguments stored in .args
            print(inst)
            failed_list.append(str(raw_file_path))
    
    nmdc_metadata_gen.dump_nmdc_database(nmdc_database, registration_file)    

    with failed_files.open('w+') as json_file:

        json_file.write(json.dumps(failed_list, indent=1))

if __name__ == "__main__":

    # run_multiprocess()
    # cpu_percents = monitor(target=run_multiprocess)
    # print(cpu_percents)
    run_nom_nmdc_data_processing()
    # file_location = get_filename()
    # if file_location:
    #    run_assignment(file_location)

