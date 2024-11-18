from metadata_generator import MetadataGenerator
import toml
from pathlib import Path

def main():
    """
    Runs the MetadataGenerator using the configuration provided in a TOML file.

    The TOML configuration file must include the following fields:

    - metadata_file: str
        The path to the .csv, .xlsx, or .xls file where the biosample metadata is stored (include file extension).
        Example: "metadata_file.csv"
    - data_dir: str
        The directory path where the data lives or will be stored.
        Example: "/path/to/data/"
    - ref_calibration_path: str
        The path to the .ref file where the reference calibration data is stored (include .ref).
        Example: "/path/to/calibration_file.ref"
    - raw_data_object_type: str
        The raw data object type. Must match one of the following options:
        - 'Direct Infusion FT ICR-MS Raw Data'
        - 'LC-DDA-MS/MS Raw Data'
    - processed_data_object_type: str
        The processed data object type. Must match one of the following options:
        - 'FT ICR-MS Analysis Results'
        - 'GC-MS Metabolomics Results'
    - registration_file: str
        The desired name of the output file where the data will be dumped (include .json).
        Example: "output_file.json"
    - execution_resource: str
        The execution resource used for the analysis. Must match one of the following options:
        - 'EMSL-RZR'
        - 'EMSL'
    - field_strength: int
        The field strength for the NOM analysis. Must match one of the following values:
        - 7
        - 12
        - 15
        - 21
    - minting_client_config_path: str
        The path to the NMDC minting client configuration file. Defaults to 'enviroMS/nmdc_metadata_generation/config.yaml'.

    Notes
    -----
    This function assumes:
    - If a biosample_id does not exist for the sample in the metadata_file, a biosample will be generated. On the other hand, 
      if a biosample_id does exist in the spreadsheet, no biosample will be generated. And the biosample_id present will be used to 
      generate all other metadata.
    - The metadata_file conforms to a predefined structure. See example spreadsheet: https://docs.google.com/spreadsheets/d/1-xHGkkG5Gpw5Pen1iM_JUP2XmphuF19ps2ZYONvtDUs/edit?gid=1112301083#gid=1112301083
    - Necessary configuration and calibrations are already added to MongoDB. If new configurations are needed, they must be added beforehand.
    """


    # Load arguments from TOML file
    config_data = toml.load('enviroMS/nmdc_metadata_generation/config.toml')

    generator = MetadataGenerator(
        metadata_file=config_data['metadata_file'],
        data_dir=Path(config_data['data_dir']),
        ref_calibration_path=Path(config_data['ref_calibration_path']),
        raw_data_object_type=config_data['raw_data_object_type'],
        processed_data_object_type=config_data['processed_data_object_type'],
        database_dump_json_path=config_data['registration_file'],
        execution_resource=config_data['execution_resource'],
        field_strength=config_data['field_strength'],
        workflow_version=config_data['workflow_version'],
        minting_client_config_path=config_data['minting_client_config_path']
    )
    
    generator.run()

if __name__ == "__main__":
    main()