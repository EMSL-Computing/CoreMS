__author__ = "Katherine R. Heal"
__date__ = "2024-01-24"

from dataclasses import dataclass

from .EI_SQL import MetaboliteMetadata


@dataclass
class LipidMetadata(MetaboliteMetadata):
    """Dataclass for the Lipid Metadata

    Parameters
    ----------
    name : str
        The name of the lipid, using the LIPID MAPS nomenclature
    casno : str
        The CAS number of the lipid
    formula : str
        The molecular formula of the lipid
    pubchem_id : str
        The PubChem ID of the lipid
    structure_level : str
        The structure level of the lipid, following the LIPID MAPS classification
    lipid_summed_name : str
        The summed name of the lipid, aka lipid species,
        following the LIPID MAPS classification
    lipid_subclass : str
        The subclass of the lipid, following the LIPID MAPS classification
    lipid_class : str
        The class of the lipid, following the LIPID MAPS classification
    lipid_category : str
        The category of the lipid, following the LIPID MAPS classification
    """

    name: str
    casno: str
    formula: str
    pubchem_id: str
    structure_level: str

    lipid_summed_name: str
    lipid_subclass: str
    lipid_class: str
    lipid_category: str
