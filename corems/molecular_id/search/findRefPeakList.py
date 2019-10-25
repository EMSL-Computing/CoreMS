__author__ = "Yuri E. Corilo"
__date__ = "Oct 24, 2019"

from threading import Thread
from pathlib import Path
import sys

from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula 
from corems.encapsulation.constant import Labels

class SearchReferencePeaks(Thread):
    
    '''
        Class to walk 14Da units over oxygen space for negative ion mass spectrum of natural organic matter
        Returns a list of MSPeak class cotaining the possibles Molecular Formula class objects.  
        
        Parameters
        ----------
        mass_spectrum_obj : MassSpec class
            This is where we store MassSpec class obj,   
        
        lookupTableSettings:  MoleculaLookupTableSettings class
            This is where we store MoleculaLookupTableSettings class obj
        
        min_O , max_O : int
            minum and maxium of oxigen to allow the software to look for
            it will override the settings at lookupTableSettings.usedAtoms
            default min = 1, max = 30

        Attributes
        ----------
        mass_spectrum_obj : MassSpec class
            This is where we store MassSpec class obj,   
        lookupTableSettings:  MoleculaLookupTableSettings class
            This is where we store MoleculaLookupTableSettings class obj
        
        Methods
        ----------
            run()    
                will be called when the instantiated class method start is called
            get_list_found_peaks()
                returns a list of MSpeaks classes cotaining all the MolecularFormula canditates inside the MSPeak
                for more details of the structure see MSPeak class and MolecularFormula class    
            set_mass_spec_indexes_by_found_peaks()
                set the mass spectrum to interate over only the selected indexes
    '''
    
    