from api_info_retriever import ApiInfoRetriever
from pathlib import Path
from dataclasses import asdict
from datetime import datetime
from metadata_parser import BiosampleIncludedMetadata, MetadataParser, NmdcTypes
from nom_workflow import run_nmdc_workflow
from linkml_runtime.dumpers import json_dumper
from tqdm import tqdm

import logging
import shutil
import nmdc_schema.nmdc as nmdc
import numpy as np
import pandas as pd
import hashlib
import json
import yaml
import oauthlib
import requests_oauthlib

# TODO: Update script to for Sample Processing - has_input for MassSpectrometry will have to be changed to be a processed sample id - not biosample id
# TODO: og_url in ApiInfoGetter in metadata_gen_supplement.py to be regular url once Berkeley is integrated

class MetadataGenerator:
   
    """
    A class to generate metadata for natural organic matter (NOM) analysis.

    Parameters
    ----------
    metadata_file : str
        The path to the metadata file.
    data_dir : str
        The directory containing the data files.
    ref_calibration_path : str
        The path to the reference calibration file.
    raw_data_object_type : str
        The type of raw data object.
    processed_data_object_type : str
        The type of processed data object.
    database_dump_json_path : str
        The path to the database dump JSON file.
    execution_resource : str
        The execution resource used in the nom analysis.
    field_strength : int
        The field strength for the analysis.
    workflow_version : str
        The version of the workflow.
    minting_client_config_path : str
        The path to the minting client configuration file.
    """

    def __init__(self, metadata_file: str, data_dir: str, ref_calibration_path: str,
                 raw_data_object_type: str, processed_data_object_type: str,
                 database_dump_json_path: str, execution_resource: str,
                 field_strength: int, workflow_version: str,
                 minting_client_config_path: str):
        self.metadata_file = metadata_file
        self.data_dir = data_dir
        self.ref_calibration_path = ref_calibration_path
        self.raw_data_object_type = raw_data_object_type
        self.processed_data_object_type = processed_data_object_type
        self.database_dump_json_path = database_dump_json_path
        self.execution_resource = execution_resource
        self.field_strength = field_strength
        self.workflow_version = workflow_version
        self.minting_client_config_path = minting_client_config_path
        self.analyte_category = "nom"
        self.processing_institution = "EMSL"
        self.raw_data_category = "instrument_data"
        self.base_url = "https://nmdcdemo.emsl.pnnl.gov/"
        self.workflow_analysis_name = "NOM Analysis"
        self.workflow_description = ("Natural Organic Matter analysis of raw mass "
                                      "spectrometry data.")
        self.workflow_git_url = "https://github.com/microbiomedata/enviroMS"
        self.processed_data_category = "processed_data"

    def run(self):
        """
        Execute the metadata generation process.

        This method processes the metadata file, generates biosamples (if needed) 
        and metadata, and manages the workflow for generating NOM analysis data.
        """
        file_ext = '.d'

        raw_dir_zip = self.data_dir / Path("raw_zip/")
        raw_dir_zip.mkdir(parents=True, exist_ok=True)

        results_dir = self.data_dir / Path("results/")
        results_dir.mkdir(parents=True, exist_ok=True)

        failed_files = results_dir / "nom_failed_files.json"

        registration_dir = self.data_dir / 'registration'
        registration_dir.mkdir(parents=True, exist_ok=True)
        registration_file = registration_dir / self.database_dump_json_path

        failed_list = []

        nmdc_database = self.start_nmdc_database()

        # Initialize parser
        parser = MetadataParser(metadata_file=self.metadata_file)

        # Load metadata spreadsheet with Biosample metadata into dataframe
        metadata_df = parser.load_metadata_file()

        # Set up logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
        )

        logging.info("Starting metadata processing...")

        # Iterate through each row in df to generate metadata
        for index, row in tqdm(metadata_df.iterrows(), total=metadata_df.shape[0], desc="Processing rows"):
            # Do not generate biosamples if biosample_id exists in spreadsheet
            if parser.check_for_biosamples(row):
                emsl_metadata = parser.parse_no_biosample_metadata(row)
                biosample_id = emsl_metadata.biosample_id
                logging.info(f"Biosample already exists for {emsl_metadata.data_path}, will not generate Biosample...")
            else:
                # Generate biosamples if no biosample_id in spreadsheet
                emsl_metadata = parser.parse_biosample_metadata(row)
                biosample = self.generate_biosample(biosamp_metadata=emsl_metadata)

                # Append biosample instance to nmdc_database
                nmdc_database.biosample_set.append(biosample)
                biosample_id = biosample.id
                logging.info(f"Generating Biosamples for {emsl_metadata.data_path}")

            # Create raw_file_path
            raw_file_path = self.data_dir / emsl_metadata.data_path.with_suffix(file_ext)

            # Run nmdc workflow
            issue, ms = run_nmdc_workflow((raw_file_path, self.ref_calibration_path, self.field_strength))

            try:
                if ms:
                    if raw_file_path.suffix == '.d':
                        raw_file_to_upload_path = Path(raw_dir_zip / raw_file_path.stem)
                        # Create a zip file
                        shutil.make_archive(raw_file_to_upload_path, 'zip', raw_file_path)
                    else:
                        raw_file_to_upload_path = raw_file_path

                    result_file_name = Path(raw_file_path.name)
                    output_file_path = results_dir / result_file_name.with_suffix('.csv')
                    ms.to_csv(output_file_path, write_metadata=False)

                    self.create_nmdc_metadata(raw_data_path=raw_file_to_upload_path.with_suffix('.zip'),
                                               data_product_path=output_file_path,
                                               emsl_metadata=emsl_metadata,
                                               biosample_id=biosample_id,
                                               nom_metadata_db=nmdc_database)
                else:
                    logging.warning(f"Workflow issue for {raw_file_path}: {issue}")
                    failed_list.append(str(raw_file_path))

            except Exception as inst:
                logging.error(f"Error processing {raw_file_path}: {str(inst)}")
                failed_list.append(str(raw_file_path))

        self.dump_nmdc_database(nmdc_database, registration_file)

        with failed_files.open('w+') as json_file:
            json_file.write(json.dumps(failed_list, indent=1))

        logging.info("Metadata processing completed.")
  
    def create_nmdc_metadata(self, raw_data_path: Path, data_product_path: Path,
                              emsl_metadata: object, biosample_id: str,
                              nom_metadata_db: nmdc.Database):
        """
        Create NMDC metadata entries.

        Parameters
        ----------
        raw_data_path : Path
            The path to the raw data file.
        data_product_path : Path
            The path to the processed data product.
        emsl_metadata : object
            The EMSL metadata object containing information about the sample.
        biosample_id : str
            The ID of the biosample.
        nom_metadata_db : nmdc.Database
            The database instance to store the generated metadata.
        """
        # Generate mass spectrometry instance
        mass_spectrometry = self.generate_mass_spectrometry(metadata_obj=emsl_metadata,
                                                            file_path=raw_data_path,
                                                            biosample_id=biosample_id,
                                                            raw_data_id="nmdc:placeholder")

        # Generate raw data object / create a raw data object description.
        eluent_intro_pretty = emsl_metadata.eluent_intro.replace("_", " ")
        raw_data_object_desc = f"Raw {emsl_metadata.instrument_used} {eluent_intro_pretty} data."
        raw_data_object = self.generate_data_object(file_path=raw_data_path,
                                                    data_category=self.raw_data_category,
                                                    data_object_type=self.raw_data_object_type,
                                                    description=raw_data_object_desc,
                                                    was_generated_by=mass_spectrometry.id)

        # Generate nom analysis instance
        nom_analysis = self.generate_nom_analysis(file_path=raw_data_path,
                                                  raw_data_id=raw_data_object.id,
                                                  data_gen_id=mass_spectrometry.id,
                                                  processed_data_id="nmdc:placeholder")

        # Generate processed data object
        processed_data_object_desc = (f"EnviroMS {emsl_metadata.instrument_used} "
                                       "natural organic matter workflow molecular formula assignment output details")
        processed_data_object = self.generate_data_object(file_path=data_product_path,
                                                          data_category=self.processed_data_category,
                                                          data_object_type=self.processed_data_object_type,
                                                          description=processed_data_object_desc,
                                                          was_generated_by=nom_analysis.id)

        # Update the outputs for mass_spectrometry and nom_analysis
        self.update_outputs(mass_spec_obj=mass_spectrometry,
                            analysis_obj=nom_analysis,
                            raw_data_obj=raw_data_object,
                            processed_data_obj=processed_data_object)

        # Add instances to database
        nom_metadata_db.data_object_set.append(raw_data_object)
        nom_metadata_db.workflow_execution_set.append(nom_analysis)
        nom_metadata_db.data_generation_set.append(mass_spectrometry)
        nom_metadata_db.data_object_set.append(processed_data_object)

    def mint_nmdc_id(self, nmdc_type: str) -> list[str]:
        """
        Mint a new NMDC ID.

        Parameters
        ----------
        nmdc_type : str
            The type of NMDC entity to mint an ID for.

        Returns
        -------
        list[str]
            A list containing the minted NMDC ID(s).
        """
        # TODO: Update api_base_url to regular url once Berkeley is integrated

        config = yaml.safe_load(open(self.minting_client_config_path))
        client = oauthlib.oauth2.BackendApplicationClient(
            client_id=config['client_id'])
        oauth = requests_oauthlib.OAuth2Session(client=client)

        api_base_url = 'https://api-berkeley.microbiomedata.org'

        token = oauth.fetch_token(token_url=f'{api_base_url}/token',
                                  client_id=config['client_id'],
                                  client_secret=config['client_secret'])

        nmdc_mint_url = f'{api_base_url}/pids/mint'

        payload = {
            "schema_class": {"id": nmdc_type
                             },
            "how_many": 1
        }

        response = oauth.post(nmdc_mint_url, data=json.dumps(payload))
        list_ids = response.json()

        return list_ids
    
    def generate_biosample(self, biosamp_metadata: BiosampleIncludedMetadata) -> nmdc.Biosample:
        """
        Generate a biosample from the given metadata.

        Parameters
        ----------
        biosamp_metadata : object
            The metadata object containing biosample information.

        Returns
        -------
        object
            The generated biosample instance.
        """

        # Drop non biosample-related keys from EmslMetadata to create Biosample dict 
        biosamp_dict = asdict(biosamp_metadata)
        non_biosamp_keys = {"data_path", "dms_dataset_id", "myemsl_link", "instrument_used", "eluent_intro", "mass_spec_config", "chrom_config_name"}
        for key in non_biosamp_keys:
            biosamp_dict.pop(key, None)

        # Remove keys with NaN values
        biosamp_dict = {k: v for k, v in biosamp_dict.items() if not (isinstance(v, float) and np.isnan(v))}

        # If no biosample id in spreadsheet, mint biosample ids
        if biosamp_dict['biosample_id'] is None:
            biosamp_dict['biosample_id'] = self.mint_nmdc_id(nmdc_type=NmdcTypes.Biosample)[0]

        # Change biosample_id to id and add type slot
        biosamp_dict['id'] = biosamp_dict.pop('biosample_id')
        biosamp_dict['type'] = NmdcTypes.Biosample

        # Filter dictionary to remove any key/value pairs with None as the value
        biosamp_dict = {k: v for k, v in biosamp_dict.items() if v is not None}
    
        biosample_object = nmdc.Biosample(**biosamp_dict)

        return biosample_object
        
    def generate_mass_spectrometry(self, metadata_obj: object, file_path: Path, biosample_id: str, raw_data_id: str) -> nmdc.DataGeneration:
        """
        Generate a mass spectrometry object from the provided metadata.

        Parameters
        ----------
        metadata_obj : object
            The metadata object containing mass spectrometry information.
        file_path : Path
            The file path of the mass spectrometry data file.
        biosample_id : str
            The ID of the associated biosample.
        raw_data_id : str
            The ID of the raw data associated with the mass spectrometry.

        Returns
        -------
        nmdc.DataGeneration
            The generated mass spectrometry object.
        """

        nmdc_id = self.mint_nmdc_id(nmdc_type=NmdcTypes.MassSpectrometry)[0]

        # Look up instrument id by name slot using API
        api_instrument_getter = ApiInfoRetriever(
            collection_name="instrument_set")
        instrument_id = api_instrument_getter.get_id_by_slot_from_collection(slot_name="name", slot_field_value=metadata_obj.instrument_used)

        # Instantiate configuration_set info retriever
        api_config_getter = ApiInfoRetriever(
            collection_name="configuration_set")
    
        # Get the mass spec configuration id based on the mass_spec_config_name
        mass_spec_config_id = api_config_getter.get_id_by_slot_from_collection(slot_name="name", slot_field_value=metadata_obj.mass_spec_config)
        
        data_dict = {
            "id": nmdc_id,
            "name": file_path.stem,
            "instrument_used": instrument_id,
            "description": f"{metadata_obj.eluent_intro} ultra high resolution mass spectrum",
            "add_date": datetime.now().strftime('%Y-%m-%d'),
            "eluent_introduction_category": metadata_obj.eluent_intro,
            "has_mass_spectrometry_configuration": mass_spec_config_id,
            "analyte_category": self.analyte_category,
            "has_input": [biosample_id],
            "has_output": [raw_data_id],
            "associated_studies": metadata_obj.associated_studies,
            "processing_institution": self.processing_institution,
            "type": NmdcTypes.MassSpectrometry
            }

        # Add the chromatography configuration if given.
        if metadata_obj.chrom_config_name:
            config_id = api_config_getter.get_id_by_slot_from_collection(slot_name="name", slot_field_value=metadata_obj.chrom_config_name)
            data_dict["has_chromatography_configuration"] = config_id
        
        massSpectrometry = nmdc.DataGeneration(**data_dict)

        return massSpectrometry

    def generate_data_object(self, file_path: Path, data_category: str, data_object_type: str,
                             description: str, was_generated_by: str = None) -> nmdc.DataObject:
        """
        Generate a data object from the provided file information.

        Parameters
        ----------
        file_path : Path
            The file path of the data object.
        data_category : str
            The category of the data object.
        data_object_type : str
            The type of the data object.
        description : str
            A description of the data object.
        was_generated_by : str, optional
            The ID of the entity that generated the data object.

        Returns
        -------
        nmdc.DataObject
            The generated data object.
        """
        
        nmdc_id = self.mint_nmdc_id(nmdc_type=NmdcTypes.DataObject)[0]

        data_dict = {
            "id": nmdc_id,
            "data_category": data_category,
            "data_object_type": data_object_type,
            "name": str(file_path.name),
            "description": description,
            "file_size_bytes": file_path.stat().st_size,
            "md5_checksum": hashlib.md5(file_path.open('rb').read()).hexdigest(),
            "was_generated_by": was_generated_by,
            "url": self.base_url + str(file_path.name),
            "type": NmdcTypes.DataObject
        }

        dataObject = nmdc.DataObject(**data_dict)

        return dataObject

    def generate_nom_analysis(self, file_path: Path, raw_data_id: str, data_gen_id: str, processed_data_id: str) -> nmdc.MetabolomicsAnalysis:
        """
        Generate a metabolomics analysis object from the provided file information.

        Parameters
        ----------
        file_path : Path
            The file path of the metabolomics analysis data file.
        raw_data_id : str
            The ID of the raw data associated with the analysis.
        data_gen_id : str
            The ID of the data generation process that informed this analysis.
        processed_data_id : str
            The ID of the processed data resulting from this analysis.

        Returns
        -------
        nmdc.MetabolomicsAnalysis
            The generated metabolomics analysis object.
        """

        nmdc_id = self.mint_nmdc_id(
            nmdc_type=NmdcTypes.NomAnalysis)[0]
        
        # Lookup calibration id by md5 checksum of ref_calibration_path file
        calib_md5 = hashlib.md5(self.ref_calibration_path.open('rb').read()).hexdigest()
        api_calib_do_getter = ApiInfoRetriever(
            collection_name="data_object_set")
        
        try:
            calib_do_id = api_calib_do_getter.get_id_by_slot_from_collection(slot_name="md5_checksum", slot_field_value=calib_md5)
            api_calibration_getter = ApiInfoRetriever(collection_name="calibration_set")
            calibration_id = api_calibration_getter.get_id_by_slot_from_collection(slot_name="calibration_object", slot_field_value=calib_do_id)

        except ValueError as e:
            print(f"Calibration object does not exist: {e}")

        except Exception as e:
            print(f"An error occurred: {e}")

        data_dict = {
            'id': nmdc_id,
            'name': f'{self.workflow_analysis_name} for {file_path.name}',
            'description': self.workflow_description,
            'has_calibration': calibration_id,
            'processing_institution': self.processing_institution,
            'execution_resource': self.execution_resource,
            'git_url': self.workflow_git_url,
            'version': self.workflow_version,
            'was_informed_by': data_gen_id,
            'has_input': [raw_data_id],
            'has_output': [processed_data_id],
            'started_at_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'ended_at_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'type': NmdcTypes.NomAnalysis,
        }

        nomAnalysis = nmdc.NomAnalysis(**data_dict)

        return nomAnalysis

    def update_outputs(self, mass_spec_obj: object, analysis_obj: object, raw_data_obj: object,
                       processed_data_obj):
        """
        Update the output references for mass spectrometry and analysis objects.

        Parameters
        ----------
        mass_spec_obj : object
            The mass spectrometry object to update.
        analysis_obj : object
            The analysis object to update.
        raw_data_obj : object
            The raw data object to reference as output.
        processed_data_obj : object
            The processed data object to reference as output.
        """
       

        mass_spec_obj.has_output = [raw_data_obj.id]
        analysis_obj.has_output = [processed_data_obj.id]

    def start_nmdc_database(self) -> nmdc.Database:
       
        return nmdc.Database()

    def dump_nmdc_database(self, nmdc_database: nmdc.Database, output_file_path: str):
        """
        Dump the NMDC database to a JSON file.

        Parameters
        ----------
        nmdc_database : nmdc.Database
            The NMDC database instance to dump.
        output_file_path : str
            The file path where the database will be dumped.
        """

        json_dumper.dump(nmdc_database, output_file_path)
        print(
            f"Database successfully dumped in {output_file_path}")