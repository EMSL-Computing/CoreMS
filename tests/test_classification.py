

import sys
sys.path.append('.')
from pathlib import Path

import pytest

from corems.encapsulation.settings.molecular_id.MolecularIDSettings import  MolecularSearchSettings
from corems.mass_spectrum.factory.classification import  HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from test_MolecularFormulaSearch import create_mass_spectrum


def test_heteroatoms_classification():

    MolecularSearchSettings.error_method = 'None'
    MolecularSearchSettings.min_mz_error = -5
    MolecularSearchSettings.max_mz_error = 3
    MolecularSearchSettings.mz_error_range = 1
    MolecularSearchSettings.isProtonated = True 
    MolecularSearchSettings.isRadical= False 
    MolecularSearchSettings.isAdduct= False 

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