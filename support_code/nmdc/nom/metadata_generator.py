from api_info_retriever import NmdcApiInfoRetriever
from pathlib import Path
from dataclasses import asdict
from datetime import datetime
from metadata_parser import BiosampleIncludedMetadata, MetadataParser, NmdcTypes
from nom_workflow import run_nmdc_workflow
from linkml_runtime.dumpers import json_dumper
from tqdm import tqdm

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
# TODO: uncomment and comment testing lines (modules), ms, issue, etc. see below in run


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
    config_path : str
        The path to the minting client and BioPortal api_key configuration file.
    """

    def __init__(self, metadata_file: str, data_dir: str, ref_calibration_path: str,
                 raw_data_object_type: str, processed_data_object_type: str,
                 database_dump_json_path: str, execution_resource: str,
                 field_strength: int, workflow_version: str,
                 config_path: str):
        self.metadata_file = metadata_file
        self.data_dir = data_dir
        self.ref_calibration_path = ref_calibration_path
        self.raw_data_object_type = raw_data_object_type
        self.processed_data_object_type = processed_data_object_type
        self.database_dump_json_path = database_dump_json_path
        self.execution_resource = execution_resource
        self.field_strength = field_strength
        self.workflow_version = workflow_version
        self.config_path = config_path
        self.analyte_category = "nom"
        self.processing_institution = "EMSL"
        self.raw_data_category = "instrument_data"
        self.base_url = "https://nmdcdemo.emsl.pnnl.gov/"
        self.workflow_analysis_name = "NOM Analysis"
        self.workflow_description = ("Natural Organic Matter analysis of raw mass "
                                      "spectrometry data.")
        self.workflow_git_url = "https://github.com/microbiomedata/enviroMS"
        self.processed_data_category = "processed_data"
        self.workflow_param_data_category = "workflow_parameter_data"
        self.workflow_param_data_object_type = "Configuration toml"

    def run(self):
        """
        Execute the metadata generation process.

        This method processes the metadata file, generates biosamples (if needed) 
        and metadata, and manages the workflow for generating NOM analysis data.
        """

        file_ext = '.d'
        raw_dir_zip, results_dir, registration_dir = self.setup_directories()
        registration_file = registration_dir / self.database_dump_json_path

        # Dictionary to track failures
        failed_metadata = {
            'validation_errors': [],
            'processing_errors': []
        }

        nmdc_database = self.start_nmdc_database()

        # Initialize parser
        parser = MetadataParser(metadata_file=self.metadata_file, config_path=self.config_path)

        # Load metadata spreadsheet with Biosample metadata into dataframe
        metadata_df = parser.load_metadata_file()

        tqdm.write("\033[92mStarting metadata processing...\033[0m")

        # Iterate through each row in df to generate metadata
        for index, row in tqdm(metadata_df.iterrows(), total=metadata_df.shape[0], desc="\033[95mProcessing rows\033[0m"):
            # Do not generate biosamples if biosample_id exists in spreadsheet
            try:
                
                # Check if biosample_id is in metadata_csv. If no biosample_id, then will generate biosamples,
                # if biosample_id exists, will return None for biosample.
                emsl_metadata, biosample_id, biosample = self.handle_biosample(parser, row)

                # Create raw_file_path
                raw_file_path = self.data_dir / emsl_metadata.data_path.with_suffix(file_ext)

                # Run nmdc workflow
                issue, ms = run_nmdc_workflow((raw_file_path, self.ref_calibration_path, self.field_strength))

                if ms:

                    # Process data files
                    raw_file_to_upload_path, output_file_path, toml_file_path = self.process_data_files(
                        ms, raw_file_path, raw_dir_zip, results_dir
                    )

                    # Generate NMDC metadata
                    self.create_nmdc_metadata(raw_data_path=raw_file_to_upload_path.with_suffix('.zip'),
                                            data_product_path=output_file_path,
                                            emsl_metadata=emsl_metadata,
                                            biosample_id=biosample_id,
                                            toml_workflow_param_path=toml_file_path,
                                            nom_metadata_db=nmdc_database)
                    
                    # Add biosample to database if it was newly generated
                    if biosample:
                        nmdc_database.biosample_set.append(biosample)
                
                else:
                    self.record_processing_error(failed_metadata, index, raw_file_path, f"Workflow issue: {issue}")

            except Exception as e:
                # Record the failed row with its error
                self.record_processing_error(
                    failed_metadata, 
                    index,
                    row.get('LC-MS filename', 'Unknown'),
                    str(e)
                )
                continue

        # At the end of processing, save the failed metadata if there are any errors
        self.save_error_log(failed_metadata, results_dir)

        self.dump_nmdc_database(nmdc_database, registration_file)

        tqdm.write("\033[92mMetadata processing completed.\033[0m")

    def setup_directories(self) -> tuple[Path, Path, Path]:
        """
        Create directory structure for storing raw, processed and registration data.

        Creates three directories:
        - raw_zip: For storing zipped raw data files
        - results: For storing processed output files
        - registration: For storing registration/metadata files

        Returns
        -------
        tuple[Path, Path, Path]
            A tuple containing (raw_dir_zip, results_dir, registration_dir) Paths
        """

        raw_dir_zip = self.data_dir / Path("raw_zip/")
        raw_dir_zip.mkdir(parents=True, exist_ok=True)

        results_dir = self.data_dir / Path("results/")
        results_dir.mkdir(parents=True, exist_ok=True)

        registration_dir = self.data_dir / 'registration'
        registration_dir.mkdir(parents=True, exist_ok=True)

        return raw_dir_zip, results_dir, registration_dir
    
    def handle_biosample(self, parser: MetadataParser, row: pd.Series) -> tuple:
        """
        Process biosample information from metadata row.

        Checks if a biosample ID exists in the row. If it does, returns the existing 
        biosample information. If not, generates a new biosample.

        Parameters
        ----------
        parser : MetadataParser
            Parser instance for processing metadata
        row : pd.Series
            A row from the metadata DataFrame containing biosample information

        Returns
        -------
        tuple
            A tuple containing:
            - emsl_metadata : BiosampleIncludedMetadata or NoBiosampleIncludedMetadata
                The parsed metadata from the row
            - biosample_id : str
                The ID of the biosample (existing or newly generated)
            - biosample : Biosample or None
                The generated biosample object if new, None if existing
        """

        if parser.check_for_biosamples(row):
            emsl_metadata = parser.parse_no_biosample_metadata(row)
            biosample_id = emsl_metadata.biosample_id
            tqdm.write(f"Biosample already exists for {emsl_metadata.data_path}, will not generate Biosample...")
            return emsl_metadata, biosample_id, None
        else:
            # Generate biosamples if no biosample_id in spreadsheet
            emsl_metadata = parser.parse_biosample_metadata(row)
            biosample = self.generate_biosample(biosamp_metadata=emsl_metadata)
            biosample_id = biosample.id
            tqdm.write(f"Generating Biosamples for {emsl_metadata.data_path}")
            return emsl_metadata, biosample_id, biosample

    def process_data_files(self, ms, raw_file_path: Path, raw_dir_zip: Path, results_dir: Path) -> tuple[Path, Path, Path]:
        """
        Process raw data files and generate output files.

        Takes a raw data file, zips it if needed (.d extension), generates CSV output,
        and creates associated TOML metadata file.

        Parameters
        ----------
        ms : object
            Mass spectrometry object with to_csv method
        raw_file_path : Path
            Path to the raw data file to be processed
        raw_dir_zip : Path  
            Directory path for storing zipped raw files
        results_dir : Path
            Directory path for storing output files

        Returns
        -------
        tuple[Path, Path, Path]
            A tuple containing:
            - raw_file_to_upload_path : Path
                Path to the processed raw file (zipped if .d extension)
            - output_file_path : Path
                Path to the generated CSV output file
            - toml_file_path : Path
                Path to the generated TOML metadata file

        Notes
        -----
        The ms.to_csv() call with write_metadata=True generates two files:
        1. A CSV file with the processed data
        2. A TOML file containing metadata configuration
        """

        if raw_file_path.suffix == '.d':
            raw_file_to_upload_path = Path(raw_dir_zip / raw_file_path.stem)
            # Create a zip file
            shutil.make_archive(raw_file_to_upload_path, 'zip', raw_file_path)
        else:
            raw_file_to_upload_path = raw_file_path

        result_file_name = Path(raw_file_path.name)
        output_file_path = results_dir / result_file_name.with_suffix('.csv')
        # to_csv with write_metadata=True will save two files: a csv and a toml file.
        ms.to_csv(output_file_path, write_metadata=True)

        # Get workflow parameter toml path
        toml_file_path = output_file_path.with_suffix('.toml')

        return raw_file_to_upload_path, output_file_path, toml_file_path
    
    def record_processing_error(self, failed_metadata: dict, index: int, filename: str, error: str) -> None:
        """
        Record processing errors in tracking dictionary and display status message.

        Records details about processing errors that occur during metadata generation,
        including row index, filename, and error message. Also displays the error to
        the console using tqdm.write.

        Parameters
        ----------
        failed_metadata : dict
            Dictionary tracking all processing errors that occur during execution
        index : int
            Zero-based row index from the metadata DataFrame
        filename : str
            Name of the file that encountered the error
        error : str
            Description or message explaining the error that occurred

        Returns
        -------
        None
            Updates failed_metadata dictionary in place and prints error message

        Notes
        -----
        Row indices are incremented by 2 in the output to account for:
        1. Zero-based DataFrame indexing
        2. Header row in original spreadsheet
        """
        failed_metadata['processing_errors'].append({
            'row_index': index + 2,
            'filename': str(filename),
            'error': str(error)
        })
        tqdm.write(f"Error processing row {index + 2}: {str(error)}")

    def save_error_log(self, failed_metadata: dict, results_dir: Path) -> None:
        """
        Save processing errors to JSON file and display notification.

        If any errors occurred during metadata processing, saves them to a JSON file
        in the results directory and displays a colored notification message.

        Parameters
        ----------
        failed_metadata : dict
            Dictionary containing lists of validation and processing errors. Expected 
            structure is {'validation_errors': [], 'processing_errors': []}
        results_dir : Path
            Directory path where the error log file should be saved

        Returns
        -------
        None
            Writes error log file if errors exist and displays status message

        Notes
        -----
        - Error log is saved as 'failed_metadata_rows.json' in the results directory
        - Output JSON is formatted with indentation level of 2 for readability
        - Uses ANSI color code 91 (bright red) for the console notification message

        See Also
        --------
        record_processing_error : Method for recording individual processing errors
        """
        if any(failed_metadata.values()):
            error_file = results_dir / 'failed_metadata_rows.json'
            with open(error_file, 'w') as f:
                json.dump(failed_metadata, f, indent=2)
            tqdm.write(f"\n\033[91mSome rows failed processing. See {error_file} for details.\033[0m")

    def create_nmdc_metadata(self, raw_data_path: Path, data_product_path: Path,
                              emsl_metadata: object, biosample_id: str,
                              toml_workflow_param_path: Path,
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
        toml_workflow_param_path: Path
            The path to the workflow parameter metadata toml file.
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
        
        # Generate workflow parameter data object
        workflow_param_data_object_desc = (f"CoreMS processing parameters for natural organic matter analysis "
                                           "used to generate {processed_data_object.id}")
        parameter_data_object = self.generate_data_object(file_path=toml_workflow_param_path,
                                                          data_category=self.workflow_param_data_category,
                                                          data_object_type=self.workflow_param_data_object_type,
                                                          description=workflow_param_data_object_desc)


        # Update the outputs for mass_spectrometry and nom_analysis
        self.update_outputs(mass_spec_obj=mass_spectrometry,
                            analysis_obj=nom_analysis,
                            raw_data_obj=raw_data_object,
                            processed_data_obj=processed_data_object,
                            workflow_param_obj=parameter_data_object)

        # Add instances to database
        nom_metadata_db.data_object_set.append(raw_data_object)
        nom_metadata_db.workflow_execution_set.append(nom_analysis)
        nom_metadata_db.data_generation_set.append(mass_spectrometry)
        nom_metadata_db.data_object_set.append(processed_data_object)
        nom_metadata_db.data_object_set.append(parameter_data_object)

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

        config = yaml.safe_load(open(self.config_path))
        client = oauthlib.oauth2.BackendApplicationClient(
            client_id=config['client_id'])
        oauth = requests_oauthlib.OAuth2Session(client=client)

        api_base_url = 'https://api.microbiomedata.org'

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
        api_instrument_getter = NmdcApiInfoRetriever(
            collection_name="instrument_set")
        instrument_id = api_instrument_getter.get_id_by_slot_from_collection(slot_name="name", slot_field_value=metadata_obj.instrument_used)

        # Instantiate configuration_set info retriever
        api_config_getter = NmdcApiInfoRetriever(
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
        api_calib_do_getter = NmdcApiInfoRetriever(
            collection_name="data_object_set")
        
        try:
            calib_do_id = api_calib_do_getter.get_id_by_slot_from_collection(slot_name="md5_checksum", slot_field_value=calib_md5)
            api_calibration_getter = NmdcApiInfoRetriever(collection_name="calibration_set")
            calibration_id = api_calibration_getter.get_id_by_slot_from_collection(slot_name="calibration_object", slot_field_value=calib_do_id)

        except ValueError as e:
            print(f"Calibration object does not exist: {e}")

        except Exception as e:
            print(f"An error occurred: {e}")

        data_dict = {
            'id': f"{nmdc_id}.1",
            'name': f'{self.workflow_analysis_name} for {file_path.name}',
            'description': self.workflow_description,
            'uses_calibration': calibration_id,
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
                       processed_data_obj: object, workflow_param_obj: object):
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
        analysis_obj.has_input.append(workflow_param_obj.id)

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
        
