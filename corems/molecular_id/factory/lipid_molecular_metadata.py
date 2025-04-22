__author__ = "Katherine R. Heal"
__date__ = "2024-01-24"

from dataclasses import dataclass

from .EI_SQL import MetaboliteMetadata

@dataclass
class LipidMetadata(MetaboliteMetadata):
    def __init__(self, casno: str, structure_level: str, lipid_summed_name: str, lipid_subclass: str, lipid_class: str, lipid_category: str, **kwargs):
        """
        Initialize LipidMetadata with specific attributes and pass additional arguments to the superclass.

        Parameters
        ----------
        casno : str
            The CAS number of the lipid
        structure_level : str
            The structure level of the lipid
        lipid_summed_name : str
            The summed name of the lipid
        lipid_subclass : str
            The subclass of the lipid
        lipid_class : str
            The class of the lipid
        lipid_category : str
            The category of the lipid
        kwargs : dict
            Additional arguments for the superclass
        """
        super().__init__(**kwargs)
        self.casno = casno
        self.structure_level = structure_level
        self.lipid_summed_name = lipid_summed_name
        self.lipid_subclass = lipid_subclass
        self.lipid_class = lipid_class
        self.lipid_category = lipid_category