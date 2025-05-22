__author__ = "Katherine R. Heal"
__date__ = "2024-01-24"

from dataclasses import dataclass

from .EI_SQL import MetaboliteMetadata

_no_default = object()
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

    casno: str = _no_default
    pubchem_id: str = _no_default
    structure_level: str = _no_default

    lipid_summed_name: str = _no_default
    lipid_subclass: str = _no_default
    lipid_class: str = _no_default
    lipid_category: str = _no_default

    def __post_init__(self):
        for field in self.__dataclass_fields__:
            if getattr(self, field) is _no_default:
                raise TypeError(
                    f"__init__ missing 1 required argument: '{field}'"
                )