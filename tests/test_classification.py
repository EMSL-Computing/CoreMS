import sys
sys.path.append('.')
from pathlib import Path

import pytest

from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.factory.classification import  HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from test_molecularFormulaSearch import create_mass_spectrum


def test_heteroatoms_classification():

    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error  = -5
    MSParameters.molecular_search.max_ppm_error = 3
    MSParameters.molecular_search.mz_error_range = 1
    MSParameters.molecular_search.isProtonated = True 
    MSParameters.molecular_search.isRadical= False 
    MSParameters.molecular_search.isAdduct= False 

    mass_spec_obj = create_mass_spectrum()
    
    assignOx = OxygenPriorityAssignment(mass_spec_obj) 
    
    assignOx.run() 

    #test classification 
    mass_spec_obj.percentile_assigned()

    mass_spectrum_by_classes = HeteroatomsClassification(mass_spec_obj)

    mass_spectrum_by_classes.plot_ms_assigned_unassigned()
    
    mass_spectrum_by_classes.atoms_ratio_all("H", "C")

    mass_spectrum_by_classes.dbe_all()

    mass_spectrum_by_classes.carbon_number_all()

    mass_spectrum_by_classes.abundance_assigned()

    mass_spectrum_by_classes.mz_exp_assigned()

    mass_spectrum_by_classes.abundance_count_percentile(Labels.unassigned)

    mass_spectrum_by_classes.peaks_count_percentile(Labels.unassigned)

if __name__ == "__main__":
     
     test_heteroatoms_classification()  