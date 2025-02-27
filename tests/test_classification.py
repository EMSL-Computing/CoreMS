import sys

from corems.molecular_id.factory.classification import  HeteroatomsClassification, Labels
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas


def test_heteroatoms_classification(mass_spectrum_ftms, postgres_database):
    mass_spectrum_ftms.molecular_search_settings.url_database = postgres_database
    mass_spectrum_ftms.molecular_search_settings.error_method = 'None'
    mass_spectrum_ftms.molecular_search_settings.min_ppm_error  = -10
    mass_spectrum_ftms.molecular_search_settings.max_ppm_error = 10
    mass_spectrum_ftms.molecular_search_settings.mz_error_range = 1
    mass_spectrum_ftms.molecular_search_settings.isProtonated = True
    mass_spectrum_ftms.molecular_search_settings.isRadical = False
    mass_spectrum_ftms.molecular_search_settings.isAdduct = False
    usedAtoms = {'C': (1, 100), 'H': (4, 200), 'O': (1, 18)}
    mass_spectrum_ftms.molecular_search_settings.usedAtoms = usedAtoms

    # Check that there are not assigned peaks
    assert mass_spectrum_ftms.percentile_assigned()[2] == 0
    
    SearchMolecularFormulas(mass_spectrum_ftms).run_worker_mass_spectrum()
    
    # Check if search was successful
    assert mass_spectrum_ftms.percentile_assigned()[2] > 0

    mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum_ftms)

    # Check that the plot is created
    mass_spectrum_by_classes.plot_ms_assigned_unassigned()

    # Check that ratios, DBE, carbon number, abundance and mz_exp are calculated
    
    assert len(mass_spectrum_by_classes.atoms_ratio_all("H", "C")) > 0
    assert len(mass_spectrum_by_classes.dbe_all()) > 0
    assert len(mass_spectrum_by_classes.abundance_assigned()) > 0
    assert len(mass_spectrum_by_classes.mz_exp_assigned()) > 0
    assert mass_spectrum_by_classes.abundance_count_percentile(Labels.unassigned) > 0
    assert mass_spectrum_by_classes.peaks_count_percentile(Labels.unassigned) > 0