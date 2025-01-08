import pandas as pd
import numpy as np

from dataclasses import dataclass
from typing import Optional
from pathlib import Path
from api_info_retriever import NmdcApiInfoRetriever, BioOntologyInfoRetriever

@ dataclass
class BiosampleIncludedMetadata:

    data_path: str
    dms_dataset_id: str
    myemsl_link: str
    sample_type: str
    samp_name: str
    env_medium: dict
    name: str
    geo_loc_name: dict
    lat_lon: dict
    env_broad_scale: dict
    env_local_scale: dict
    description: str
    collection_date: dict
    associated_studies: list
    samp_store_temp: dict
    samp_collec_device: str
    elev: float
    instrument_used: str
    eluent_intro: str
    mass_spec_config: str
    chrom_config_name: Optional[str] = None
    size_frac: Optional[dict] = None
    samp_size: Optional[dict] = None
    ecosystem_subtype: Optional[str] = None
    ecosystem: Optional[str] = None
    soil_type: Optional[dict] = None
    light_regm: Optional[dict] = None
    soil_horizon: Optional[str] = None
    depth: Optional[dict] = None
    biosample_categories: Optional[list] = None
    air_temp_regm: Optional[list] = None
    sample_collection_site: Optional[str] = None
    ncbi_taxonomy_name: Optional[str] = None
    ecosystem_type: Optional[str] = None
    location: Optional[str] = None
    img_identifiers: Optional[list] = None
    habitat: Optional[str] = None
    ecosystem_category: Optional[str] = None
    biosample_id: Optional[str] = None

@ dataclass
class NoBiosampleIncludedMetadata:
    biosample_id: str
    data_path: str
    dms_dataset_id: str
    myemsl_link: str
    associated_studies: list
    instrument_used: str
    eluent_intro: str
    mass_spec_config: str
    chrom_config_name: Optional[str] = None

@dataclass
class NmdcTypes:

    Biosample: str = "nmdc:Biosample"
    MassSpectrometry: str = "nmdc:MassSpectrometry"
    NomAnalysis: str = "nmdc:NomAnalysis"
    DataObject: str = "nmdc:DataObject"
    OntologyClass: str = "nmdc:OntologyClass"
    ControlledIdentifiedTermValue: str = "nmdc:ControlledIdentifiedTermValue"
    TextValue: str = "nmdc:TextValue"
    GeolocationValue: str = "nmdc:GeolocationValue"
    TimeStampValue: str = "nmdc:TimestampValue"
    QuantityValue: str = "nmdc:QuantityValue"
    CalibrationInformation: str = "nmdc:CalibrationInformation"
    MetaboliteIdentification: str = "nmdc:MetaboliteIdentification"

class MetadataParser:
    """Parsers metadata from input metadata spreadsheet."""

    def __init__(self, metadata_file, config_path):
        """
        Parameters
        ----------
        metadata_file : str
            Path to the metadata file to be loaded.
        """

        self.metadata_file = metadata_file
        self.config_path = config_path

    def load_metadata_file(self) -> pd.DataFrame:
        """
        Load metadata from a file, either CSV or Excel format.

        Raises
        ------
        ValueError
            If the file extension is not supported.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the loaded metadata.
        """

        metadata_file_path = Path(self.metadata_file)
        if metadata_file_path.suffix == '.csv':
            metadata_df = pd.read_csv(metadata_file_path)
        elif metadata_file_path.suffix in {'.xlsx', '.xls'}:
            metadata_df = pd.read_excel(metadata_file_path)
        else:
            raise ValueError(f'Unsupported file extension: {metadata_file_path.suffix}')
        
        # Check that all configs are valid
        self.check_for_valid_configs(metadata_df)
        # Check that 'chrom__config_name' column value is present if eluent_intro is 'liquid_chromatography)
        self.check_chrom_config(metadata_df)

        return metadata_df
    
    def check_chrom_config(self, df: pd.DataFrame):
        """
        Check the configuration for chromatography.

        This method verifies that the 'chrom_config_name' is provided 
        for liquid chromatography and gas chromatography, and that it 
        is empty for direct infusion methods.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame containing metadata to be checked.

        Raises
        ------
        ValueError
            If 'chrom_config_name' is missing or incorrectly set 
            for the specified 'eluent_intro' types.
        """
       # Ensure chrom_config_name is present for liquid_chromatography and gas_chromatography
        required_mask = ((df['eluent_intro'] == 'liquid_chromatography') | (df['eluent_intro'] == 'gas_chromatography')) \
                        & (df['chrom_config_name'].isna() | (df['chrom_config_name'] == ''))
        
        invalid_required_rows = df[required_mask]

        if not invalid_required_rows.empty:
            raise ValueError(f"Missing 'chrom_config_name' for the following rows where 'eluent_intro' is 'liquid_chromatography' or 'gas_chromatography': {(invalid_required_rows.index+2).tolist()}")

        # Ensure chrom_config_name is empty for direct_infusion_syringe or direct_infusion_autosampler
        no_config_needed_mask = ((df['eluent_intro'] == 'direct_infusion_syringe') | (df['eluent_intro'] == 'direct_infusion_autosampler')) \
                                & ~(df['chrom_config_name'].isna() | (df['chrom_config_name'] == ''))
        
        invalid_no_config_rows = df[no_config_needed_mask]

        if not invalid_no_config_rows.empty:
            raise ValueError(f"'chrom_config_name' should be empty for the following rows where 'eluent_intro' is 'direct_infusion_syringe' or 'direct_infusion_autosampler': {(invalid_no_config_rows.index+2).tolist()}")

    def check_for_valid_configs(sefl, df: pd.DataFrame):
        """
        Get unique values for all columns containing 'config' in their name and verify they exist in API.
        
        Parameters
        ----------
        df : pd.DataFrame
            DataFrame to analyze
            
        Raises
        ------
        ValueError
            If any config value doesn't exist in the API
        """

        # Instantiate configuration_set info retriever
        api_config_getter = NmdcApiInfoRetriever(
            collection_name="configuration_set")
        
        config_columns = ['chrom_config_name', 'mass_spec_config']
        invalid_configs = []

        for col in config_columns:
            unique_vals = [val for val in df[col].unique() if pd.notna(val)]

            for val in unique_vals:
                try:
                    # skip empty values
                    if not val:
                        continue
                    api_config_getter.get_id_by_slot_from_collection(slot_name="name", slot_field_value=val)
                except ValueError:
                    invalid_configs.append((col, val))

        # If any invalid conifgs were found, raise error
        if invalid_configs:
            error_msg = "The following configurations were not found in the API:\n"
            for col, val in invalid_configs:
                error_msg += f"  Column '{col}': '{val}'\n"
            raise ValueError(error_msg)

    def check_for_biosamples(self, row: pd.Series) -> bool:
        """
        Check if the biosample_id is not None, NaN, or empty.

        Parameters
        ----------
        row : pd.Series
            A row from the DataFrame containing metadata.

        Returns
        -------
        bool
            True if biosample_id is valid; False otherwise.
        """

        value = row.get('biosample_id')
        if pd.isna(value) or value == '':
            return False
        return True
    
    # Helper function to handle missing or NaN values
    def get_value(self, row: pd.Series, key: str, default=None):
        """
        Retrieve a value from a row, handling missing or NaN values.

        Parameters
        ----------
        row : pd.Series
            A row from the DataFrame.
        key : str
            The key to retrieve the value for.
        default : optional
            Default value to return if the key does not exist or is NaN.

        Returns
        -------
        The value associated with the key, or default if not found.
        """
        value = row.get(key, default)
        if isinstance(value, float) and np.isnan(value):
            return default
        return value
        
    def parse_no_biosample_metadata(self, row: pd.Series) -> NoBiosampleIncludedMetadata:
        """
        Parse the metadata row if it does not include biosample information.

        Parameters
        ----------
        row : pd.Series
            A row from the DataFrame containing metadata.

        Returns
        -------
        NoBiosampleIncludedMetadata
            An instance of NoBiosampleIncludedMetadata populated with the parsed values.
        """

        #Initialize metadata dictionary
        metadata_dict = {
            'biosample_id': self.get_value(row, 'biosample_id'),
            'data_path': Path(self.get_value(row, 'LC-MS filename')),
            'dms_dataset_id': self.get_value(row, 'DMS Dataset ID'),
            'myemsl_link': self.get_value(row, 'MyEMSL link'),
            'associated_studies': [study.strip() for study in self.get_value(row, 'NMDC Study ID').split(',')] if self.get_value(row, 'NMDC Study ID') else None,
            'instrument_used': self.get_value(row, 'instrument_used'),
            'eluent_intro': self.get_value(row, 'eluent_intro'),
            'mass_spec_config': self.get_value(row, 'mass_spec_config'),
            'chrom_config_name': self.get_value(row, 'chrom_config_name')
        }

        # Create and return the EmslMetadata instance
        metadata = NoBiosampleIncludedMetadata(**metadata_dict)

        return metadata

    def parse_biosample_metadata(self, row: pd.Series) -> BiosampleIncludedMetadata:
        """
        Parse the metadata row if it includes biosample information.

        Parameters
        ----------
        row : pd.Series
            A row from the DataFrame containing metadata.

        Returns
        -------
        BiosampleIncludedMetadata
            An instance of BiosampleIncludedMetadata populated with the parsed values.
        """

        # Initialize BioOntologyInfoRetriever
        envo_retriever = BioOntologyInfoRetriever(config_path=self.config_path)
        
        # Initialize the metadata dictionary
        metadata_dict = {
            'data_path': Path(self.get_value(row, 'LC-MS filename')),
            'dms_dataset_id': self.get_value(row,'DMS Dataset ID'),
            'myemsl_link': self.get_value(row,'MyEMSL link'),
            'sample_type': self.get_value(row,'Sample Type'),
            'samp_name': self.get_value(row,'samp_name'),
            'env_medium': self.create_controlled_identified_term_value(self.get_value(row, 'env_medium'), envo_retriever.get_envo_terms(self.get_value(row, 'env_medium'))) if self.get_value(row, 'env_medium') else None,
            'name': self.get_value(row,'name'),
            'geo_loc_name': self.create_text_value(self.get_value(row, 'geo_loc_name'), is_list=False) if self.get_value(row, 'geo_loc_name') else None,
            'lat_lon': self.create_geo_loc_value(raw_value=self.get_value(row, 'lat_lon: has_raw_value'), lat_value=self.get_value(row,'latitude'), long_value=self.get_value(row, 'longitude'))
                        if self.get_value(row, 'lat_lon: has_raw_value') and self.get_value(row, 'latitude') and self.get_value(row,'longitude') else None,
            'env_broad_scale': self.create_controlled_identified_term_value(self.get_value(row,'env_broad_scale'), envo_retriever.get_envo_terms(self.get_value(row, 'env_broad_scale'))) if self.get_value(row, 'env_broad_scale') else None,
            'env_local_scale': self.create_controlled_identified_term_value(self.get_value(row, 'env_local_scale'), envo_retriever.get_envo_terms(self.get_value(row, 'env_local_scale'))) if self.get_value(row,'env_local_scale') else None,
            'description': self.get_value(row, 'description'),
            'collection_date': self.create_time_stamp_value(row_value=self.get_value(row, 'collection_date:has_raw_value')) if self.get_value(row, 'collection_date:has_raw_value') else None,
            'associated_studies': [study.strip() for study in self.get_value(row, 'NMDC Study ID').split(',')] if self.get_value(row, 'NMDC Study ID') else None,
            'samp_store_temp': {'has_raw_value': self.get_value(row, 'samp_store_temp.has_raw_value'), 'type': NmdcTypes.QuantityValue} if self.get_value(row, 'samp_store_temp.has_raw_value') else None,
            'samp_size': self.create_quantity_value(raw_value=self.get_value(row, 'samp_size.has_raw_value')) if self.get_value(row, 'samp_size.has_raw_value') else None,
            'samp_collec_device': self.get_value(row, 'samp_collec_device'),
            'elev': self.get_value(row, 'elev'),
            'size_frac': self.create_text_value(self.get_value(row, 'size_frac.has_raw_value'), is_list=False) if self.get_value(row, 'size_frac.has_raw_value') else None,
            'air_temp_regm': self.create_text_value(self.get_value(row, 'air_temp_regm.has_raw_value'), is_list=True) if self.get_value(row,'air_temp_regm.has_raw_value') else None,
            'biosample_categories': [category.strip() for category in self.get_value(row,'biosample_categories').split(',')] if self.get_value(row,'biosample_categories') else None,
            'depth': self.create_quantity_value(num_value=self.get_value(row,'depth_has_numeric_value'), min_num_value=self.get_value(row,'depth.has_minimum_numeric_value'), unit_value=self.get_value(row,'depth.has_unit')) if self.get_value(row,'depth_has_numeric_value') or self.get_value(row,'depth.has_minimum_numeric_value') or self.get_value(row, 'depth.has_unit') else None,
            'soil_horizon': self.get_value(row, 'soil_horizon'),
            'light_regm': self.create_text_value(self.get_value(row, 'light_regm.has_raw_value'), is_list=False) if self.get_value(row, 'light_regm.has_raw_value') else None,
            'soil_type': self.create_text_value(self.get_value(row, 'soil_type.has_raw_value'), is_list=False) if self.get_value(row, 'soil_type.has_raw_value') else None,
            'img_identifiers': [identifier.strip() for identifier in self.get_value(row, 'img_identifiers').split(',')] if self.get_value(row, 'img_identifiers') else None,
            'ncbi_taxonomy_name': self.get_value(row, 'ncbi_taxonomy_name'),
            'ecosystem_type': self.get_value(row, 'ecosystem_type'),
            'location': self.get_value(row, 'location'),
            'habitat': self.get_value(row, 'habitat'),
            'ecosystem_category': self.get_value(row, 'ecosystem_category'),
            'ecosystem': self.get_value(row, 'ecosystem'),
            'sample_collection_site': self.get_value(row, 'sample_collection_site'),
            'ecosystem_subtype': self.get_value(row, 'ecosystem_subtype'),
            'biosample_id': self.get_value(row, 'biosample_id'),
            'instrument_used': self.get_value(row, 'instrument_used'),
            'eluent_intro': self.get_value(row, 'eluent_intro'),
            'mass_spec_config': self.get_value(row, 'mass_spec_config'),
            'chrom_config_name': self.get_value(row, 'chrom_config_name')
        }

        # Create and return the EmslMetadata instance
        metadata = BiosampleIncludedMetadata(**metadata_dict)
        
        return metadata
    
    def create_controlled_identified_term_value(self, row_value: str, slot_enum_dict: dict):
        """
        Create a controlled identified term value.

        Parameters
        ----------
        raw_value : str
            The raw value to be converted.
        control_terms : dict
            A mapping of controlled terms.

        Returns
        -------
        dict
            A dictionary representing the controlled identified term.
        """

        nmdc_controlled_term_slot = {
            'has_raw_value': row_value,
            'term': {
                'id': row_value,
                'name': slot_enum_dict.get(row_value),
                'type': NmdcTypes.OntologyClass
            },
            'type': NmdcTypes.ControlledIdentifiedTermValue
        }
        
        return nmdc_controlled_term_slot
    
    def create_text_value(self, row_value: str, is_list: bool):
        """
        Create a text value representation.

        Parameters
        ----------
        row_value : str
            The raw value to convert.
        is_list : bool
            Whether to treat the value as a list.

        Returns
        -------
        dict
            A dictionary representing the text value.
        """
        
        if is_list:
            output_list = []
            values = row_value.split(',')
            for value in values:
                stripped_value = value.strip()
                nmdc_text_value = {
                    'has_raw_value': stripped_value,
                    'type': NmdcTypes.TextValue
                    }
                output_list.append(nmdc_text_value)

            return output_list
            
        else: 
            nmdc_text_value = {
                    'has_raw_value': row_value,
                    'type': NmdcTypes.TextValue
                    }
            
            return nmdc_text_value

    def create_geo_loc_value(self, raw_value: str, lat_value: str, long_value: str):
        """
        Create a geolocation value representation.

        Parameters
        ----------
        raw_value : str
            The raw value associated with geolocation.
        lat_value : str
            The latitude value.
        long_value : str
            The longitude value.

        Returns
        -------
        dict
            A dictionary representing the geolocation value.
        """

        nmdc_geo_loc_value = {
            'has_raw_value': raw_value,
            'latitude': lat_value,
            'longitude': long_value,
            'type': NmdcTypes.GeolocationValue}
        
        return nmdc_geo_loc_value
    
    def create_time_stamp_value(self, row_value:str):
        """
        Create a timestamp value representation.

        Parameters
        ----------
        row_value : str
            The raw value to convert to a timestamp.

        Returns
        -------
        dict
            A dictionary representing the timestamp value.
        """

        nmdc_time_stamp_value = {
            'has_raw_value': row_value,
            'type': NmdcTypes.TimeStampValue}
        
        return nmdc_time_stamp_value
    
    def create_quantity_value(self, num_value:str = None, min_num_value:str = None, max_num_value:str = None, unit_value:str = None, raw_value:str = None):
        """
        Create a quantity value representation.

        Parameters
        ----------
        raw_value : str
            The raw value to convert to a quantity.

        Returns
        -------
        dict
            A dictionary representing the quantity value.
        """
        
        nmdc_quant_value = {
            'type': NmdcTypes.QuantityValue
            }
        
        if num_value:
            nmdc_quant_value['has_numeric_value'] = num_value
        if min_num_value:
            nmdc_quant_value['has_minimum_numeric_value'] = min_num_value
        if max_num_value:
            nmdc_quant_value['has_maximum_numeric_value'] = max_num_value
        if unit_value:
            nmdc_quant_value['has_unit'] = unit_value
        if num_value and unit_value: 
            nmdc_quant_value['has_raw_value'] = str(num_value) + str(unit_value)
        if raw_value:
            nmdc_quant_value['has_raw_value'] = raw_value
        
        return nmdc_quant_value
