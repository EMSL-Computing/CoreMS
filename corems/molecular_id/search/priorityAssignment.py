"""Legacy oxygen-priority molecular formula assignment (retired).

Canonical formula assignment uses
:class:`~corems.molecular_id.search.molecularFormulaSearch.SearchMolecularFormulas`
(DI / FT-ICR) or
:class:`~corems.molecular_id.search.molecularFormulaSearch.SearchMolecularFormulasLC`
(LC-MS). :class:`FindOxygenPeaks` remains available for oxygen-series utilities
and calibration workflows.
"""


class OxygenPriorityAssignment:
    """Removed oxygen-priority formula assignment strategy.

    This class is no longer supported. Use
    :class:`~corems.molecular_id.search.molecularFormulaSearch.SearchMolecularFormulas`
    for molecular formula assignment. Construction raises
    :class:`NotImplementedError`. The name is retained as a stub until the next
    major release.

    Parameters
    ----------
    *args
        Ignored; retained only for call-site compatibility.
    **kwargs
        Ignored; retained only for call-site compatibility.
    """

    def __init__(self, *args, **kwargs) -> None:
        raise NotImplementedError(
            "OxygenPriorityAssignment is no longer supported. "
            "Use SearchMolecularFormulas for molecular formula assignment. "
            "This stub will be removed in the next major release."
        )
