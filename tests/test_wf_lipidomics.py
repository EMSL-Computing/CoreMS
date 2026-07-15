from pathlib import Path
from unittest.mock import MagicMock

import shutil
import numpy as np
import pandas as pd
import pytest

from corems.encapsulation.constant import Labels
from corems.mass_spectra.input.corems_hdf5 import ReadCoreMSHDFMassSpectra
from corems.mass_spectra.input.mzml import MZMLSpectraParser
from corems.mass_spectra.output.export import LipidomicsExport
from corems.molecular_id.search.database_interfaces import LCLipidLibraryInterface
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulasLC
from corems.encapsulation.factory.parameters import LCMSParameters, reset_lcms_parameters, reset_ms_parameters


def test_import_lcmsobj_mzml():
    # Instantiate parser based on binary file type
    file_mzml = (
        Path.cwd()
        / "tests/tests_data/lcms/"
        / "test_centroid_neg_RP_metab.mzML"
    )

    parser = MZMLSpectraParser(file_mzml)

    # Get the instrument information and creation time
    instrument_info = parser.get_instrument_info()
    assert instrument_info['model'] == "Orbitrap ID-X"
    creation_time = parser.get_creation_time()
    assert creation_time.year == 2022
    
    # Instatiate lc-ms data object using parser and pull in ms1 spectra into dataframe (without storing as MassSpectrum objects to save memory)
    myLCMSobj = parser.get_lcms_obj(spectra="none")
    myLCMSobj = parser.get_lcms_obj(spectra="ms1")
    myLCMSobj.parameters = LCMSParameters(use_defaults=True)

    # Modify parameters to deal with centroid data
    myLCMSobj.parameters.lc_ms.peak_picking_method = "centroided_persistent_homology"
    myLCMSobj.parameters.mass_spectrum[
            "ms1"
        ].mass_spectrum.noise_threshold_method = "relative_abundance"

    myLCMSobj.find_mass_features()
    myLCMSobj.integrate_mass_features()
    myLCMSobj.add_associated_ms1(
        auto_process=True, use_parser=True, spectrum_mode="centroid"
    )
    mass_features_df = myLCMSobj.mass_features_to_df()
    assert mass_features_df.shape[0] == 1183
    assert mass_features_df.shape[1] > 15
    
    # Reset the MSParameters to the original values
    reset_lcms_parameters()
    reset_ms_parameters()


def test_mzml_scan_list_skips_non_ms_controller_duplicates():
    """Multi-controller mzML reuses scan numbers; keep MS spectra, not aux traces."""

    class Spec:
        def __init__(self, scan_id, ms_level, mz, intensity, centroid=True):
            self.ID = scan_id
            self.ms_level = ms_level
            self.mz = np.asarray(mz, dtype=float)
            self.i = np.asarray(intensity, dtype=float)
            self._centroid = centroid

        def get(self, accession):
            if accession == "MS:1000127":
                return True if self._centroid else None
            return None

        def __getitem__(self, key):
            return True if key == "negative scan" else None

    spectra = [
        Spec(2, 2, [100.1, 200.2], [10.0, 20.0]),  # MS
        Spec(2, None, [1.0], [1.0], centroid=False),  # same ID, non-MS controller
        Spec(5, None, [9.0], [9.0], centroid=False),  # non-MS first
        Spec(5, 1, [150.5], [42.0]),  # then MS
    ]
    parser = MZMLSpectraParser.__new__(MZMLSpectraParser)
    parser.file_location = Path("fake_multicontroller.mzML")
    reader = MagicMock()
    reader.__iter__.return_value = iter(spectra)
    parser.load = MagicMock(return_value=reader)

    result = parser.get_mass_spectra_from_scan_list(
        [2, 5], spectrum_mode="centroid", auto_process=False
    )
    assert [list(ms.data_dict[Labels.mz]) for ms in result] == [
        pytest.approx([100.1, 200.2]),
        pytest.approx([150.5]),
    ]


@pytest.mark.lipidomics_db
@pytest.mark.molecular_db
def test_lipidomics_workflow(tmp_path, postgres_database, lcms_obj, lipidomics_sqlite_path):
    # Delete the "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems" directory
    shutil.rmtree(
        "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems",
        ignore_errors=True,
    )

    # Set parmaeters to the defaults for reproducible testing
    lcms_obj.parameters = LCMSParameters(use_defaults=True)

    # Set parameters on the LCMS object that are reasonable for testing
    ## persistent homology parameters
    lcms_obj.parameters.lc_ms.peak_picking_method = "persistent homology"
    lcms_obj.parameters.lc_ms.ph_inten_min_rel = 0.0005
    lcms_obj.parameters.lc_ms.ph_persis_min_rel = 0.05
    lcms_obj.parameters.lc_ms.ph_smooth_it = 0
    lcms_obj.parameters.lc_ms.ms2_min_fe_score = 0.3
    lcms_obj.parameters.lc_ms.ms1_scans_to_average = 5

    ## MSParameters for ms1 mass spectra
    ms1_params = lcms_obj.parameters.mass_spectrum['ms1']
    ms1_params.mass_spectrum.noise_threshold_method = "relative_abundance"
    ms1_params.mass_spectrum.noise_threshold_min_relative_abundance = 0.1
    ms1_params.mass_spectrum.noise_min_mz, ms1_params.mass_spectrum.min_picking_mz = 0, 0
    ms1_params.mass_spectrum.noise_max_mz, ms1_params.mass_spectrum.max_picking_mz = np.inf, np.inf
    ms1_params.ms_peak.legacy_resolving_power = False
    ms1_params.molecular_search.url_database = postgres_database
    ms1_params.molecular_search.usedAtoms = {
        'C': (10, 30),
        'H': (18, 200),
        'O': (1, 23),
        'N': (0, 3),
        'P': (0, 1),
        'S': (0, 1),
    }

    ## settings for ms2 data (HCD scans)
    ms2_params_hcd = ms1_params.copy()
    ms2_params_hcd.molecular_search.ion_types_excluded = ["[M+HCOO]-"]
    lcms_obj.parameters.mass_spectrum['ms2'] = ms2_params_hcd

    ## settings for ms2 data (CID scans)
    ms2_params_cid = ms2_params_hcd.copy()
    ms2_params_cid.molecular_search.max_ppm_error = 200 # wider ppm error for CID scans
    ms2_params_cid.mass_spectrum.noise_threshold_min_relative_abundance = 0.01 # lower noise threshold for CID scans
    lcms_obj.parameters.mass_spectrum['ms2_cid'] = ms2_params_cid

    ## reporting settings
    lcms_obj.parameters.lc_ms.search_as_lipids = True
    lcms_obj.parameters.lc_ms.include_fragment_types = True
    lcms_obj.parameters.lc_ms.export_eics = True
    lcms_obj.parameters.lc_ms.export_profile_spectra = True

    # Use persistent homology to find mass features in the lc-ms data
    # Find mass features, cluster, and integrate them.  Then annotate pairs of mass features that are c13 iso pairs.

    lcms_obj.find_mass_features()
    assert len(lcms_obj.mass_features) == 131
    lcms_obj.integrate_mass_features(drop_if_fail=True)
    lcms_obj.add_peak_metrics()
    lcms_obj.add_associated_ms1(
        auto_process=True, use_parser=False, spectrum_mode="profile"
    )
    lcms_obj.deconvolute_ms1_mass_features()

    mass_spec_decon = lcms_obj.mass_features[1].mass_spectrum_deconvoluted
    assert len(mass_spec_decon.mspeaks) < len(
        lcms_obj.mass_features[1].mass_spectrum.mspeaks
    )
    lcms_obj.find_c13_mass_features()
    assert len(lcms_obj.mass_features) == 128

    # Perform a molecular search on all of the mass features' ms1 peaks
    mol_form_search = SearchMolecularFormulasLC(lcms_obj)
    mol_form_search.run_mass_feature_search()

    # Check results of molecular search
    assert lcms_obj.mass_features[0].ms1_peak[0].string == "C20 H30 O2"
    assert lcms_obj.mass_features_ms1_annot_to_df().shape[0] > 128
    lcms_obj.mass_features[0].mass_spectrum.to_dataframe()

    # Add hcd ms2 data to lcms object, using the ms2 mass spectrum parameters
    og_ms_len = len(lcms_obj._ms)
    lcms_obj.add_associated_ms2_dda(spectrum_mode="centroid", scan_filter="hcd")
    assert len(lcms_obj._ms) > og_ms_len

    # Add cid ms2 data to lcms object, using the ms2_cid mass spectrum parameters
    og_ms_len = len(lcms_obj._ms)
    lcms_obj.add_associated_ms2_dda(spectrum_mode="centroid", ms_params_key="ms2_cid", scan_filter="cid")
    assert len(lcms_obj._ms) > og_ms_len

    lcms_obj.plot_composite_mz_features()

    # Export the mass features to a pandas dataframe
    df = lcms_obj.mass_features_to_df()
    assert df.shape[0] == 128
    assert df.shape[1] > 15

    # Plot a mass feature
    lcms_obj.mass_features[0].plot(return_fig=False)

    mzs = [
        mass_feature.mz
        for mass_feature in lcms_obj.mass_features.values()
        if len(mass_feature.ms2_scan_numbers) > 0 and mass_feature.isotopologue_type is None
    ]
    lipid_library = LCLipidLibraryInterface(db_location=str(lipidomics_sqlite_path))
    spectra_library_fe, lipid_metadata = lipid_library.get_lipid_library(
        mz_list=mzs,
        polarity="negative",
        mz_tol_ppm=5,
        format="flashentropy",
        normalize=True,
        fe_kwargs={
            "normalize_intensity": True,
            "min_ms2_difference_in_da": 0.02,  # for cleaning spectra
            "max_ms2_tolerance_in_da": 0.01,  # for setting search space
            "max_indexed_mz": 3000,
            "precursor_ions_removal_da": None,
            "noise_threshold": 0,
        },
    )

    # Perform a spectral search on the mass features
    hcd_ms2_scan_df = lcms_obj.scan_df[
        lcms_obj.scan_df.scan_text.str.contains("hcd")
        & (lcms_obj.scan_df.ms_level == 2)
    ]
    ms2_scans_oi_hr = [
        x for x in hcd_ms2_scan_df.scan.tolist() if x in lcms_obj._ms.keys()
    ]
    lcms_obj.fe_search(
        scan_list=ms2_scans_oi_hr, fe_lib=spectra_library_fe, peak_sep_da=0.01
    )
    # Export the lcms object to an hdf5 file using the LipidomicsExport class
    export_stem = tmp_path / "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801"
    export_dir = tmp_path / "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.corems"
    exporter = LipidomicsExport(str(export_stem), lcms_obj)
    exporter.to_hdf(overwrite=True)
    exporter.to_parquet(overwrite=True, export_spectra=True)
    exporter.report_to_csv(molecular_metadata=lipid_metadata)
    exporter.report_to_parquet(molecular_metadata=lipid_metadata)
    report = exporter.to_report(molecular_metadata=lipid_metadata)
    assert report['Ion Formula'][1] == 'C24 H47 O2'
    assert report['Lipid Molecular Species'][0] == 'FA 20:5'
    parquet_path = export_dir / "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.parquet"
    assert parquet_path.exists()
    parquet_report = pd.read_parquet(parquet_path)
    assert parquet_report.shape[0] == report.shape[0]
    assert 'Ion Formula' in parquet_report.columns
    assert (export_dir / "scan_info.parquet").exists()
    assert (export_dir / "mass_features.parquet").exists()

    # Import the hdf5 file, assert that its df is same as above and that we can plot a mass feature
    parser = ReadCoreMSHDFMassSpectra(
        export_dir / "Blanch_Nat_Lip_C_12_AB_M_17_NEG_25Jan18_Brandi-WCSH5801.hdf5"
    )
    
    # Check that creation_time was saved and can be retrieved
    creation_time = parser.get_original_creation_time()
    assert creation_time is not None
    assert creation_time.year == 2018  # Based on the filename date
    
    myLCMSobj2 = parser.get_lcms_obj()

    # Check that the parameters match
    assert myLCMSobj2.parameters == lcms_obj.parameters

    # Check that the spectra parser class is the same as the original parser and that we can plot a mass spectrum using the original parser
    assert myLCMSobj2.spectra_parser_class.__name__ == "ImportMassSpectraThermoMSFileReader"
    myLCMSobj2.spectra_parser.get_mass_spectrum_from_scan(1, spectrum_mode="profile").plot_centroid()

    # Check that the mass features dataframe is the same as the original
    df2 = myLCMSobj2.mass_features_to_df()
    assert df2.shape[0] == 128
    assert df2.shape[1] > 15
    myLCMSobj2.mass_features[0].mass_spectrum.to_dataframe()
    assert myLCMSobj2.mass_features[0].ms1_peak[0].string == "C20 H30 O2"
    assert myLCMSobj2.mass_features_ms1_annot_to_df().shape[0] > 130
    myLCMSobj2.mass_features[0].plot(return_fig=False)

    # Reset the MSParameters to the original values
    reset_lcms_parameters()
    reset_ms_parameters()
