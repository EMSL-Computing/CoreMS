import sys
sys.path.append('.')
from pathlib import Path

import pytest

from corems.encapsulation.factory.parameters import MSParameters
from corems.molecular_id.factory.classification import  HeteroatomsClassification, Labels
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from test_molecularFormulaSearch import create_mass_spectrum


def test_heteroatoms_classification():

    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error  = -10
    MSParameters.molecular_search.max_ppm_error = 10
    MSParameters.molecular_search.mz_error_range = 1
    MSParameters.molecular_search.isProtonated = True 
    MSParameters.molecular_search.isRadical= False 
    MSParameters.molecular_search.isAdduct= False 
    
    MSParameters.molecular_search.usedAtoms['C'] = (1, 100)
    MSParameters.molecular_search.usedAtoms['H'] = (4, 200)
    MSParameters.molecular_search.usedAtoms['O'] = (1, 18)
    #MSParameters.molecular_search.usedAtoms = usedatoms
    
    mass_spec_obj = create_mass_spectrum()
    
    assignOx = SearchMolecularFormulas(mass_spec_obj).run_worker_mass_spectrum()
    
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