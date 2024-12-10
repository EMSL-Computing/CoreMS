import requests
from typing import Optional
import yaml

# TODO: Update og_url to be regular nmdc api url when berkeley is implemented

class NmdcApiInfoRetriever:
    def __init__(self, collection_name: str):
        self.collection_name = collection_name

    def get_id_by_slot_from_collection(self, slot_name: str, slot_field_value: str):
        """
        Retrieve the NMDC identifier from a specified collection based on a slot name and field value.

        Parameters
        ----------
        slot_name : str
            The name of the slot to filter by.
        slot_field_value : str
            The value of the slot field to filter for. Trailing whitespace will be removed.

        Returns
        -------
        str
            The identifier corresponding to the specified slot name and field value.

        Raises
        ------
        ValueError
            If the request to the API fails or if no resources are found for the given slot name and value.
        """
        # trim trailing white spaces
        slot_field_value = slot_field_value.rstrip()

        filter = f'{{"{slot_name}": "{slot_field_value}"}}'
        field = "id"

        og_url = f'https://api.microbiomedata.org/nmdcschema/{self.collection_name}?&filter={filter}&projection={field}'
        resp = requests.get(og_url)

        # Check if the response status is 200
        if resp.status_code != 200:
            raise ValueError(f"Failed to retrieve data from {self.collection_name}, response code: {resp.status_code}")
        
        data = resp.json()

        # Ensure there is at least one resource in the response
        if not data['resources']:
            raise ValueError(f"No resources in Mongo found for '{slot_name}' slot in {self.collection_name} with value {slot_field_value}")
        
        identifier = data['resources'][0]['id']

        return identifier

class BioOntologyInfoRetriever:
    """
    Client for retrieving ENVO term information from BioPortal API.

    A class to handle authentication and retrieval of Environmental Ontology (ENVO)
    terms using the BioPortal REST API service.

    Parameters
    ----------
    config_path : str
        Path to YAML configuration file containing BioPortal API credentials

    Notes
    -----
    The configuration file should contain an 'api_key' field with a valid
    BioPortal API key.

    Examples
    --------
    >>> retriever = BioOntologyInfoRetriever('config.yaml')
    >>> envo_terms = retriever.get_envo_terms('ENVO:00002042')
    >>> print(envo_terms)
    {'ENVO:00002042': 'surface water'}
    """
    def __init__(self, config_path: str):
        self.config = config_path

    def get_envo_terms(self, envo_id: dict):
        """
        Look up an ENVO term label using BioPortal API.

        Parameters
        ----------
        envo_id : str
            The ENVO identifier to look up (e.g., 'ENVO:00002042')

        Returns
        -------
        dict
            Dictionary with envo_id as key and term label as value
            Example: {'ENVO:00002042': 'surface water'}

        Raises
        ------
        requests.exceptions.RequestException
            If the API request fails
        KeyError
            If the response doesn't contain expected data format
        yaml.YAMLError
            If the config file cannot be parsed
        FileNotFoundError
            If the config file is not found

        Notes
        -----
        Makes an authenticated request to BioPortal API to retrieve the
        preferred label (prefLabel) for the given ENVO term.
        """
        
        config = yaml.safe_load(open(self.config))
        api_key = config['api_key']

        url = f"http://data.bioontology.org/ontologies/ENVO/classes/{envo_id}"
        headers = {"Authorization": f"apikey token={api_key}"}

        response = requests.get(url, headers=headers)
        response.raise_for_status()

        data = response.json()
        return {envo_id: data['prefLabel']}