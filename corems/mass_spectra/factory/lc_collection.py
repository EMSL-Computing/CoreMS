class MassSpectraCollectionBase:
    """Base class for a collection of MassSpectra objects.
    
    Attributes
    ----------
    _mass_spectra : list
        A list of MassSpectraBase objects.
    """
    def __init__(self):
        self._mass_spectra = []