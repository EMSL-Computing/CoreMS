import sys

sys.path.append(".")

from corems.molecular_id.factory.classification import HeteroatomsClassification
from corems.mass_spectrum.input.numpyArray import ms_from_array_centroid
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment


def test_run_molecular_formula_search(postgres_database):
    """Test for generating accurate molecular formula from mass and isotope using the local sql database"""
    # Generate a mass spectrum object from a list of mz and abundance
    mz = [760.58156938877, 761.58548]
    abundance = [1, 0.4]
    rp, s2n = [[1, 1], [1, 1]]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "single mf search", polarity=1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    # Set the settings for the molecular search on the mass spectrum object
    mass_spectrum_obj.molecular_search_settings.url_database = postgres_database
    mass_spectrum_obj.molecular_search_settings.error_method = "None"
    mass_spectrum_obj.molecular_search_settings.min_ppm_error = -5
    mass_spectrum_obj.molecular_search_settings.max_ppm_error = 5
    mass_spectrum_obj.molecular_search_settings.mz_error_range = 1
    mass_spectrum_obj.molecular_search_settings.isProtonated = True
    mass_spectrum_obj.molecular_search_settings.isRadical = False
    mass_spectrum_obj.molecular_search_settings.isAdduct = False
    mass_spectrum_obj.molecular_search_settings.usedAtoms = {
        "C": (1, 57),
        "H": (4, 200),
        "N": (0, 1),
    }

    # Process the mass spectrum object
    mass_spectrum_obj.process_mass_spec()

    # Run the molecular formula search on the mass spectrum object
    SearchMolecularFormulas(
        mass_spectrum_obj, find_isotopologues=True
    ).run_worker_ms_peaks([mass_spectrum_obj[0]])
    assert mass_spectrum_obj.to_dataframe().shape[0] > 1
    assert mass_spectrum_obj[0][0].string == "C56 H73 N1"
    assert mass_spectrum_obj[1][0].string == "C55 H73 N1 13C1"


def test_run_molecular_formula_search_adduct(postgres_database):
    """Test for generating accurate molecular formula from mass and isotope, for an adduct, using the local sql database"""
    # Generate a mass spectrum object from a list of mz and abundance
    mz = [782.563522, 783.566877]  # Na+ adduct of C56H73N1 and its M+1
    abundance = [1, 0.4]
    rp, s2n = [[1, 1], [1, 1]]
    mass_spectrum_obj = ms_from_array_centroid(
        mz, abundance, rp, s2n, "single mf search", polarity=1, auto_process=False
    )
    mass_spectrum_obj.settings.noise_threshold_method = "absolute_abundance"
    mass_spectrum_obj.settings.noise_threshold_absolute_abundance = 0

    # Set the settings for the molecular search on the mass spectrum object
    mass_spectrum_obj.molecular_search_settings.url_database = postgres_database
    mass_spectrum_obj.molecular_search_settings.error_method = "None"
    mass_spectrum_obj.molecular_search_settings.min_ppm_error = -5
    mass_spectrum_obj.molecular_search_settings.max_ppm_error = 5
    mass_spectrum_obj.molecular_search_settings.mz_error_range = 1
    mass_spectrum_obj.molecular_search_settings.isProtonated = True
    mass_spectrum_obj.molecular_search_settings.isRadical = False
    mass_spectrum_obj.molecular_search_settings.isAdduct = True
    mass_spectrum_obj.molecular_search_settings.usedAtoms = {
        "C": (1, 57),
        "H": (4, 200),
        "N": (0, 1),
    }
    mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False
    mass_spectrum_obj.molecular_search_settings.use_isotopologue_filter = False

    # Process the mass spectrum object
    mass_spectrum_obj.process_mass_spec()

    # Run the molecular formula search on the mass spectrum object
    SearchMolecularFormulas(
        mass_spectrum_obj, find_isotopologues=True
    ).run_worker_ms_peaks([mass_spectrum_obj[0]])
    assert mass_spectrum_obj.to_dataframe().shape[0] > 1
    assert mass_spectrum_obj[0][0].string == "C56 H73 N1"
    assert mass_spectrum_obj[0][0].H_C == 73 / 56
    assert mass_spectrum_obj[1][0].string == "C55 H73 N1 13C1"
    assert mass_spectrum_obj[1][0].H_C == 73 / 56


def test_mspeak_search(mass_spectrum_ftms, postgres_database):
    mass_spectrum_ftms.molecular_search_settings.url_database = postgres_database
    mass_spectrum_ftms.molecular_search_settings.usedAtoms = {
        "C": (1, 100),
        "H": (4, 200),
        "O": (0, 10),
        "N": (0, 1),
        "P": (0, 1),
    }
    mass_spectrum_ftms.molecular_search_settings.isAdduct = False
    mass_spectrum_ftms.molecular_search_settings.isRadical = False
    mspeak_obj = mass_spectrum_ftms.most_abundant_mspeak
    assert round(mass_spectrum_ftms.most_abundant_mspeak.mz_exp, 2) == 421.04
    SearchMolecularFormulas(mass_spectrum_ftms).run_worker_ms_peaks([mspeak_obj])
    assert mspeak_obj.is_assigned
    # Try each of the possible filters
    mspeak_obj.molecular_formula_earth_filter()
    mspeak_obj.molecular_formula_water_filter()
    mspeak_obj.molecular_formula_air_filter()
    mspeak_obj.cia_score_S_P_error()
    mspeak_obj.cia_score_N_S_P_error()
    assert mspeak_obj.best_molecular_formula_candidate.string == "C29 H11 O2 P1"
    mspeak_obj[0].string_formated
    mspeak_obj[0].mz_error


def test_molecular_formula_search_db(mass_spectrum_ftms, postgres_database):
    mass_spectrum_ftms.molecular_search_settings.url_database = postgres_database
    mass_spectrum_ftms.molecular_search_settings.usedAtoms = {
        "C": (1, 100),
        "H": (4, 200),
        "O": (0, 10),
        "N": (0, 1),
        "P": (0, 1),
    }

    SearchMolecularFormulas(mass_spectrum_ftms, first_hit=True).run_worker_mass_spectrum()

    i = 0
    j = 0
    error = list()
    mass = list()
    abundance = list()

    for mspeak in mass_spectrum_ftms.sort_by_abundance():
        if mspeak.is_assigned:
            i += 1
            for mformula in mspeak:
                mass.append(mspeak.mz_exp)
                error.append(mformula.mz_error)
                abundance.append(mspeak.abundance)
        else:
            j += 1
            pass
    fraction_assigned = i / (i + j)
    assert fraction_assigned > 0.7


def test_priorityAssignment(mass_spectrum_ftms, postgres_database):
    mass_spectrum_ftms.molecular_search_settings.url_database = postgres_database
    mass_spectrum_ftms.molecular_search_settings.error_method = "None"
    mass_spectrum_ftms.molecular_search_settings.min_ppm_error = -3
    mass_spectrum_ftms.molecular_search_settings.max_ppm_error = 5
    mass_spectrum_ftms.molecular_search_settings.mz_error_range = 1
    mass_spectrum_ftms.molecular_search_settings.isProtonated = True
    mass_spectrum_ftms.molecular_search_settings.isRadical = True
    mass_spectrum_ftms.molecular_search_settings.isAdduct = False
    usedatoms = {"C": (1, 100), "H": (4, 200), "O": (1, 10)}
    mass_spectrum_ftms.molecular_search_settings.usedAtoms = usedatoms
    mass_spectrum_ftms.process_mass_spec()

    # Run the molecular formula search on the mass spectrum object and check the percentage of assigned peaks
    assignOx = OxygenPriorityAssignment(mass_spectrum_ftms)
    assignOx.run()
    assert mass_spectrum_ftms.percentile_assigned()[0] > 15

    # Test the HeteroatomsClassification class
    mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum_ftms)
    mass_spectrum_by_classes.plot_ms_assigned_unassigned()
    assert mass_spectrum_by_classes.atoms_ratio_all("H", "C")[0] > 0.5
