from ..nom.metadata_generator import MetadataGenerator
from ..nom.metadata_parser import MetadataParser, NmdcTypes
from ..nom.api_info_retriever import NmdcApiInfoRetriever
from tqdm import tqdm
from datetime import datetime
from pathlib import Path

import nmdc_schema.nmdc as nmdc
import hashlib

#### This child class of MetadataGenerator is for NMDC's metabalomics. I added as many of the TODOs I can see below.
#### The biggest things are adding functionalities for things that are specific to the nmdc metabolomics - like the 
#### calibration files (generates_calibration in the DataGeneration and uses_calibration in the WorkflowExecution), adding
#### the metabolites identified, and adding the call to the metabolomics workflow in the run method. 
#### It assumes that we will use the same shape metadata file so I added the example csv to this folder
#### but perhaps more/different columns would need to be added. It also uses the same config.toml file for the nom metadata generation
#### script which may need to be edited/adjust as well. I left the general classes in the nom folder, but these may want to be 
#### moved since they are shared by both nom and metabolomics (the metadata_generator.py, metadata_parser.py, and 
#### api_info_retriever.py files)
# TODO: Refine hardcoded attributes (analyte_category, workflow_analysis_name, workflow_description, workflow_git_url)
# TODO: fix run method specifically for metabolomics. Try to use same methods/outline as the run method in the 
        # nom_metadata_generator.py. (what is copied below) Try to use setup_directories, handle_biosample, process_data_files, 
        # record_processing_error, save_error_log. If not using these, move these methods under the NOMMetadataGenerator
        # as they won't be shared across child classes. Might be able to just replace the run_nmdc_workflow line with whateveer
        # the metabolomics workflow is
# TODO: Adjust the create_metadata method for metabolomics. It seems like there will be a configuration data object that will
        # be created (generates_calibration slot) for the data generation. May need to add attributes and what not for this.  
# TODO: in create_metadata method: add functionality for has_metabolite_identifications in the metab_analysis bit.
# TODO: in create_metadata_method: Fix mass_spec_description (if different than f"{emsl_metadata.eluent_intro} ultra high resolution mass spectrum")
# TODO: in generate_metab_analysis method: Fix calibration logic if needed. Not sure if it will work the same as NOM.
# TODO: Build out has_metabolite_identifications. Will need to add a step to calculate these and then validate them in the 
        # generate_metab_identifications method.

class MetabolomicsMetadataGenerator(MetadataGenerator):
    analyte_category="metabolome",
    workflow_analysis_name="Metabolomics Analysis",
    workflow_description=("Metabolomics analysis of raw mass "
                                    "spectrometry data."),
    workflow_git_url="https://github.com/microbiomedata/coreMS"

    def __init__(self, metadata_file: str, data_dir: str, ref_calibration_path: str,
                 raw_data_object_type: str, processed_data_object_type: str,
                 database_dump_json_path: str, execution_resource: str,
                 field_strength: int, workflow_version: str,
                 config_path: str):
        super().__init__(metadata_file, data_dir, ref_calibration_path,
                 raw_data_object_type, processed_data_object_type,
                 database_dump_json_path, execution_resource,
                 field_strength, workflow_version,
                 config_path)
        
    def run(self):
        #TODO: fix this method specicially for metablomics. Try to use same methods/outline as the run method in the 
        # nom_metadata_generator.py. (what is copied below) Try to use setup_directories, handle_biosample, process_data_files, 
        # record_processing_error, save_error_log. If not using these, move these methods under the NOMMetadataGenerator
        # as they won't be shared across child classes. Might be able to replace the run_nmdc_workflow_line with whatever the
        # metablomics workflow is.

        """
        Execute the metadata generation process.

        This method processes the metadata file, generates biosamples (if needed) 
        and metadata, and manages the workflow for generating Metabolomics analysis data.
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
                # TODO: replace this line with whatever is called for metabolomics workflow
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
                                            metab_metadata_db=nmdc_database)
                    
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

    def create_nmdc_metadata(self, raw_data_path: Path, data_product_path: Path,
                              emsl_metadata: object, biosample_id: str,
                              toml_workflow_param_path: Path,
                              metab_metadata_db: nmdc.Database):
        
        # TODO: Adjust for metabolomics. It seems like there will be a configuration data object that will be created
        # (generates_calibration) for the data generation. May need to add attributes and what not for this.  
        # TODO: add functionality for has_metabolite_identifications in the metab_analysis bit.
        # TODO: Fix mass_spec_description (if dfferent than f"{emsl_metadata.eluent_intro} ultra high resolution mass spectrum")
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
        metab_metadata_db : nmdc.Database
            The database instance to store the generated metadata.
        """
        # Generate mass spectrometry instance
        mass_spectrometry = self.generate_mass_spectrometry(metadata_obj=emsl_metadata,
                                                            file_path=raw_data_path,
                                                            biosample_id=biosample_id,
                                                            raw_data_id="nmdc:placeholder",
                                                            mass_spec_description=f"{emsl_metadata.eluent_intro} ultra high resolution mass spectrum")

        # Generate raw data object / create a raw data object description.
        eluent_intro_pretty = emsl_metadata.eluent_intro.replace("_", " ")
        raw_data_object_desc = f"Raw {emsl_metadata.instrument_used} {eluent_intro_pretty} data."
        raw_data_object = self.generate_data_object(file_path=raw_data_path,
                                                    data_category=self.raw_data_category,
                                                    data_object_type=self.raw_data_object_type,
                                                    description=raw_data_object_desc,
                                                    was_generated_by=mass_spectrometry.id)

        # TODO: Build out Generate has_metabolite_identifications - I believe it will go here, or add a placeholder
        # to the metabolite_identifications field in the metab_analysis step and add them later. Will need to add a method
        # to get the metabolite_identifcations from the data. This method just creates the metadata for 1 instance
        # of the metabolite identifications, will need to append them to a list and add ot eh metab_analysis step as a list:
        metabolite_identification = self.generate_metab_identifications(highest_similarity_score="",
                                                                        alt_ids="")
        
        # Generate metab analysis instance
        metab_analysis = self.generate_metab_analysis(file_path=raw_data_path,
                                                  raw_data_id=raw_data_object.id,
                                                  data_gen_id=mass_spectrometry.id,
                                                  processed_data_id="nmdc:placeholder",
                                                  metabolite_identifications=metabolite_identifications)

        # Generate processed data object
        processed_data_object_desc = (f"CoreMS {emsl_metadata.instrument_used} "
                                       "metabolomics workflow molecular formula assignment output details")
        processed_data_object = self.generate_data_object(file_path=data_product_path,
                                                          data_category=self.processed_data_category,
                                                          data_object_type=self.processed_data_object_type,
                                                          description=processed_data_object_desc,
                                                          was_generated_by=metab_analysis.id)
        
        # Generate workflow parameter data object
        workflow_param_data_object_desc = (f"CoreMS processing parameters for metabolomics analysis "
                                           "used to generate {processed_data_object.id}")
        parameter_data_object = self.generate_data_object(file_path=toml_workflow_param_path,
                                                          data_category=self.workflow_param_data_category,
                                                          data_object_type=self.workflow_param_data_object_type,
                                                          description=workflow_param_data_object_desc)


        # Update the outputs for mass_spectrometry and metab_analysis
        self.update_outputs(mass_spec_obj=mass_spectrometry,
                            analysis_obj=metab_analysis,
                            raw_data_obj=raw_data_object,
                            processed_data_obj=processed_data_object,
                            workflow_param_obj=parameter_data_object)

        # Add instances to database
        metab_metadata_db.data_object_set.append(raw_data_object)
        metab_metadata_db.workflow_execution_set.append(metab_analysis)
        metab_metadata_db.data_generation_set.append(mass_spectrometry)
        metab_metadata_db.data_object_set.append(processed_data_object)
        metab_metadata_db.data_object_set.append(parameter_data_object)

    def generate_metab_identifications(self, highest_similarity_score: float, alt_ids:list):
        """Generate the metabolite identifications for the workflow"""
        
        nmdc_id = self.mint_nmdc_id(
            nmdc_type=NmdcTypes.MetaboliteIdentification)[0]
        
        data_dict = {
            'id': nmdc_id,
            'highest_similarity_score': highest_similarity_score,
            'alternative_identifiers': alt_ids,
            'type': NmdcTypes.MetaboliteIdentification,
        }

        # Just used as a validation step, will not be returned as this will be an inlined list in the Workflow execution step.
        metabIdentications = nmdc.MetaboliteIdentification(**data_dict)

        return data_dict

    def generate_metab_analysis(self, file_path: Path, raw_data_id: str, data_gen_id: str, 
                                processed_data_id: str, metabolite_identifications: list) -> nmdc.MetabolomicsAnalysis:
        # TODO: Fix calibration logic if needed. Not sure if it will work the same as NOM.
        
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
        metabolite_identifications: list
            A list of the dictionaries for metabolite identifications.

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
            'has_metabolite_identifications': metabolite_identifications,
            'type': NmdcTypes.NomAnalysis,
        }

        metabAnalysis = nmdc.MetabolomicsAnalysis(**data_dict)

        return metabAnalysis