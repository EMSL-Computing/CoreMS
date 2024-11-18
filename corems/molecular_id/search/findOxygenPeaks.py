__author__ = "Yuri E. Corilo"
__date__ = "Jul 31, 2019"

from copy import deepcopy
from threading import Thread

from numpy import average, std

from corems.molecular_id.calc.ClusterFilter import ClusteringFilter
from corems.molecular_id.factory.molecularSQL import MolForm_SQL
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas


class FindOxygenPeaks(Thread):
    """Class to find Oxygen peaks in a mass spectrum for formula assignment search

    Class to walk 14Da units over oxygen space for negative ion mass spectrum of natural organic matter
    Returns a list of MSPeak class containing the possible Molecular Formula class objects.

    Parameters
    ----------
    mass_spectrum_obj : MassSpec class
        This is where we store MassSpec class obj,

    lookupTableSettings:  MolecularLookupTableSettings class
        This is where we store MolecularLookupTableSettings class obj

    min_O , max_O : int
        minium and maximum of Oxygen to allow the software to look for
        it will override the settings at lookupTableSettings.usedAtoms
        default min = 1, max = 22

    Attributes
    ----------
    mass_spectrum_obj : MassSpec class
        This is where we store MassSpec class obj,
    lookupTableSettings:  MolecularLookupTableSettings class
        This is where we store MolecularLookupTableSettings class obj

    Methods
    ----------
    * run().
            will be called when the instantiated class method start is called
    * get_list_found_peaks().
            returns a list of MSpeaks classes cotaining all the MolecularFormula candidates inside the MSPeak
            for more details of the structure see MSPeak class and MolecularFormula class
    * set_mass_spec_indexes_by_found_peaks().
            set the mass spectrum to interate over only the selected indexes
    """

    def __init__(
        self, mass_spectrum_obj, sql_db: bool = False, min_O: int = 1, max_O: int = 22
    ):
        Thread.__init__(self)

        self.mass_spectrum_obj = mass_spectrum_obj
        self.min_0 = min_O
        self.max_O = max_O

        if not sql_db:
            self.sql_db = MolForm_SQL(
                mass_spectrum_obj.molecular_search_settings.url_database
            )
        else:
            self.sql_db = sql_db

    def run(self):
        """Run the thread"""
        # save initial settings min peaks per class filter
        initial_min_peak_bool = deepcopy(
            self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter
        )

        # deactivate the usage of min peaks per class filter
        self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = False

        # save initial settings for Ox
        initial_ox = deepcopy(
            self.mass_spectrum_obj.molecular_search_settings.usedAtoms["O"]
        )

        # resets the used atoms to look only for oxygen organic compounds
        self.mass_spectrum_obj.molecular_search_settings.usedAtoms["O"] = (
            self.min_0,
            self.max_O,
        )

        self.list_found_mspeaks = []

        kmd_base = self.mass_spectrum_obj.mspeaks_settings.kendrick_base

        self.mass_spectrum_obj.change_kendrick_base_all_mspeaks(kmd_base)

        # needs to be wrapped inside the mass_spec class
        ClusteringFilter().filter_kendrick(self.mass_spectrum_obj)

        if self.mass_spectrum_obj.parameters.mass_spectrum.verbose_processing:
            print("Start most abundant mass spectral peak search")
        molecular_formula_obj_reference = self.find_most_abundant_formula(
            self.mass_spectrum_obj
        )

        if self.mass_spectrum_obj.parameters.mass_spectrum.verbose_processing:
            print(
                "Select most abundant peak with molecular formula =  %s with a m/z error of %s ppm"
                % (
                    molecular_formula_obj_reference.string,
                    molecular_formula_obj_reference.mz_error,
                )
            )
            print("Started mass spectral peak series search")

        self.list_found_mspeaks = self.find_series_mspeaks(
            self.mass_spectrum_obj, molecular_formula_obj_reference, deltamz=14
        )

        # reset indexes after done with operation that includes a filter (i.e. ClusteringFilter().filter_kendrick())

        self.mass_spectrum_obj.molecular_search_settings.usedAtoms["O"] = initial_ox

        self.mass_spectrum_obj.molecular_search_settings.use_min_peaks_filter = (
            initial_min_peak_bool
        )

        self.mass_spectrum_obj.reset_indexes()

        self.mass_spectrum_obj.filter_by_noise_threshold()
        if self.mass_spectrum_obj.parameters.mass_spectrum.verbose_processing:
            print("Done with mass spectral peak series search")

        self.sql_db.close()

    def find_most_abundant_formula(self, mass_spectrum_obj):
        """Find the most abundant formula in the mass spectrum

        Parameters
        ----------
        mass_spectrum_obj : MassSpec class
            Mass spectrum object

        Returns
        ----------
        MolecularFormula class obj
            most abundant MolecularFormula with the lowest mass error
        """
        # need to find a better way to cut off outliners
        # import matplotlib.pyplot as plt
        # plt.hist(mass_spectrum_obj.abundance, bins=100)
        # plt.show()

        abundances = mass_spectrum_obj.abundance
        abun_mean = average(abundances, axis=0)
        abun_std = std(abundances, axis=0)

        upper_limit = abun_mean + 7 * abun_std
        if mass_spectrum_obj.parameters.mass_spectrum.verbose_processing:
            print(
                "Maximum abundance limit  = %s and max abundance kendrick cluster = %s"
                % (
                    upper_limit,
                    max(mass_spectrum_obj, key=lambda m: m.abundance).abundance,
                )
            )

        mspeak_most_abundant = max(
            mass_spectrum_obj,
            key=lambda m: m.abundance if m.abundance <= upper_limit else 0,
        )

        print("Searching molecular formulas")

        SearchMolecularFormulas(mass_spectrum_obj, self.sql_db).run_worker_ms_peaks(
            [mspeak_most_abundant]
        )

        print("Finished searching molecular formulas")

        if mspeak_most_abundant:
            return mspeak_most_abundant.best_molecular_formula_candidate

        else:
            raise Exception(
                "Could not find a possible molecular formula match for the most abundant peak of m/z %.5f"
                % mspeak_most_abundant.mz_exp
            )

        # return the first option
        # return mspeak_most_abundant[0]

    def find_most_abundant_formula_test(self, mass_spectrum_obj, settings):
        """[Test function] Find the most abundant formula in the mass spectrum

        Parameters
        ----------
        mass_spectrum_obj : MassSpec class
            Mass spectrum object
        settings : MolecularSearchSettings class
            Molecular search settings object

        Returns
        ----------
        MolecularFormula class obj
            most abundant MolecularFormula with the lowest mass error

        """
        # this function is intended for test only.
        # Have to sort by Kendrick to be able to select the most abundant series
        # then select the most abundant peak inside the series
        # or have the user select the reference mspeak on the gui

        mspeak_most_abundant = mass_spectrum_obj.most_abundant_mspeak

        SearchMolecularFormulas(mass_spectrum_obj, self.sql_db).run_worker_ms_peaks(
            [mspeak_most_abundant]
        )

        if mspeak_most_abundant:
            return mspeak_most_abundant.best_molecular_formula_candidate

        else:
            raise Exception(
                "Could not find a possible molecular formula match for the most abundant peak of m/z %.5f"
                % mspeak_most_abundant.mz_exp
            )
        # return the first option
        # return mspeak_most_abundant[0]

    def find_series_mspeaks(
        self, mass_spectrum_obj, molecular_formula_obj_reference, deltamz=14
    ):
        """Find a series of abundant peaks in the mass spectrum for a given molecular formula

        Parameters
        ----------
        mass_spectrum_obj : MassSpec class
            Mass spectrum object
        molecular_formula_obj_reference : MolecularFormula class
            Molecular formula object
        deltamz : float
            delta m/z to look for peaks

        Returns
        ----------
        list
            list of MSpeak class objects
        """
        abundances = mass_spectrum_obj.abundance
        abun_mean = average(abundances, axis=0)
        abun_std = std(abundances, axis=0)
        upper_limit = abun_mean + 7 * abun_std

        list_most_abundant_peaks = list()

        min_mz = mass_spectrum_obj.min_mz_exp

        max_mz = mass_spectrum_obj.max_mz_exp

        initial_nominal_mass = molecular_formula_obj_reference.mz_nominal_calc

        mass = initial_nominal_mass

        nominal_masses = []
        while mass <= max_mz:
            # print "shit 1", mass, min_mz
            mass += deltamz
            nominal_masses.append(mass)

        mass = initial_nominal_mass
        while mass >= min_mz:
            # print "shit 1", mass, min_mz
            mass -= deltamz
            nominal_masses.append(mass)

        nominal_masses = sorted(nominal_masses)

        for nominal_mass in nominal_masses:
            first_index, last_index = (
                mass_spectrum_obj.get_nominal_mz_first_last_indexes(nominal_mass)
            )

            ms_peaks = mass_spectrum_obj[first_index:last_index]

            if ms_peaks:
                #
                # print (nominal_mass, first_index,
                #    last_index,
                #    mass_spectrum_obj[first_index].mz_exp,
                #    mass_spectrum_obj[last_index].mz_exp
                #    )
                #

                mspeak_most_abundant = max(
                    ms_peaks,
                    key=lambda m: m.abundance if m.abundance <= upper_limit else 0,
                )

                # mspeak_most_abundant = max(ms_peaks, key=lambda m: m.abundance)

                list_most_abundant_peaks.append(mspeak_most_abundant)
        if mass_spectrum_obj.parameters.mass_spectrum.verbose_processing:
            print("Start molecular formula search")
        SearchMolecularFormulas(mass_spectrum_obj, self.sql_db).run_worker_ms_peaks(
            list_most_abundant_peaks
        )
        if mass_spectrum_obj.parameters.mass_spectrum.verbose_processing:
            print("Done molecular formula search")
        return [mspeak for mspeak in list_most_abundant_peaks if mspeak]

    def get_list_found_peaks(self):
        """Get the list of found peaks

        Returns
        ----------
        list
            list of MSpeak class objects
        """
        return sorted(self.list_found_mspeaks, key=lambda mp: mp.mz_exp)

    def set_mass_spec_indexes_by_found_peaks(self):
        """Set the mass spectrum to interate over only the selected indexes.

        Notes
        ----------
        Warning!!!!
        set the mass spectrum to interate over only the selected indexes
        don not forget to call mass_spectrum_obj.reset_indexes after the job is done
        """

        indexes = [msp.index for msp in self.list_found_mspeaks]
        self.mass_spectrum_obj.set_indexes(indexes)
