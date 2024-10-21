import requests
from typing import Optional

# TODO: Update og_url to be regular nmdc api url when berkeley is implemented

class ApiInfoRetriever:
    def __init__(self, collection_name: str):
        self.collection_name = collection_name

    def get_id_by_slot_from_collection(self, slot_name: str, slot_field_value: str):
        """
        Retrieve the identifier from a specified collection based on a slot name and field value.

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

        og_url = f'https://api-berkeley.microbiomedata.org/nmdcschema/{self.collection_name}?&filter={filter}&projection={field}'
        resp = requests.get(og_url)

        # Check if the response status is 200
        if resp.status_code != 200:
            raise ValueError(f"Failed to retrieve data from {self.collection_name}, response code: {resp.status_code}")
        
        data = resp.json()

        # Ensure there is at least one resource in the response
        if not data['resources']:
            raise ValueError(f"No resources found for {slot_name} with value {slot_field_value}")
        
        identifier = data['resources'][0]['id']

        return identifier

