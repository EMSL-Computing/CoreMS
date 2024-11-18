import numpy as np


class FreqDomain_Calibration:
    """Frequency Domain Calibration class for mass spectrum.

    Parameters
    ----------
    mass_spectrum : MassSpectrum
        The mass spectrum object.
    selected_mass_peaks : list
        List of selected mass peaks.
    include_isotopologue : bool, optional
        Flag to include isotopologues, by default False.

    Attributes
    ----------
    mz_exp : ndarray
        Array of experimental m/z values.
    mz_calc : ndarray
        Array of calculated m/z values.
    freq_exp : ndarray
        Array of experimental frequencies.
    mass_spectrum : MassSpectrum
        The mass spectrum object.
    freq_exp_ms : ndarray
        Array of experimental frequencies for mass spectrum.

    Methods
    -------
    * recal_mass_spec(mz_domain, Aterm, Bterm, Cterm).
        Recalibrate the mass spectrum with the given parameters.
    * linear().
        Perform linear calibration.
    * quadratic(iteration=False).
        Perform quadratic calibration.
    * ledford_calibration(iteration=False).
        Perform Ledford calibration.
    * step_fit(steps=4).
        Perform step fit calibration.

    """

    def __init__(self, mass_spectrum, selected_mass_peaks, include_isotopologue=False):
        self.selected_mspeaks = selected_mass_peaks
        error = list()
        freq_exp = list()
        mz_calc = list()
        mz_exp = list()

        for mspeak in selected_mass_peaks:
            if not include_isotopologue:
                molecular_formulas = [
                    formula for formula in mspeak if not formula.is_isotopologue
                ]
            else:
                molecular_formulas = mspeak

            for molecular_formula in molecular_formulas:
                freq_exp.append(mspeak.freq_exp)
                error.append(molecular_formula.mz_error)
                mz_calc.append(molecular_formula.mz_calc)
                mz_exp.append(mspeak.mz_exp)

        self.mz_exp = np.array(mz_exp)
        self.mz_calc = np.array(mz_calc)
        self.freq_exp = np.array(freq_exp)
        self.mass_spectrum = mass_spectrum
        self.freq_exp_ms = np.array([mspeak.freq_exp for mspeak in mass_spectrum])

    def recal_mass_spec(self, mz_domain, Aterm, Bterm, Cterm):
        """Recalibrate the mass spectrum with the given parameters.

        Parameters
        ----------
        mz_domain : ndarray
            Array of m/z values for recalibration.
        Aterm : float
            Aterm parameter for recalibration.
        Bterm : float
            Bterm parameter for recalibration.
        Cterm : float
            Cterm parameter for recalibration.

        """
        self.mass_spectrum._calibration_terms = (Aterm, Bterm, 0)
        self.mass_spectrum.mz_cal = mz_domain

    def linear(self):
        """Perform linear calibration."""
        matrix = np.vstack([1 / self.freq_exp, np.ones(len(self.freq_exp))]).T
        Aterm, Bterm = np.linalg.lstsq(matrix, self.mz_calc, rcond=None)[0]
        if self.mass_spectrum.parameters.mass_spectrum.verbose_processing:
            print("%.2f Aterm,  %.2f Bterm" % (Aterm, Bterm))
            print("Linear Calibration %.2f Aterm,  %.2f Bterm " % (Aterm, Bterm))
        mz_domain = (Aterm / self.freq_exp_ms) + Bterm
        self.recal_mass_spec(mz_domain, Aterm, Bterm, 0)

    def quadratic(self, iteration: bool = False):
        """Perform quadratic calibration.

        Parameters
        ----------
        iteration : bool, optional
            Flag to perform iterative calibration, by default False.

        """
        mz_calc = self.mz_calc
        freq_exp = self.freq_exp
        mz_exp = self.mz_exp

        error = ((mz_exp - mz_calc) / mz_calc) * 1000000
        last_rms = np.sqrt(np.mean(error**2))
        while True:
            matrix = np.vstack(
                [1 / freq_exp, 1 / np.power(freq_exp, 2), np.ones(len(freq_exp))]
            ).T
            Aterm, Bterm, Cterm = np.linalg.lstsq(matrix, self.mz_calc, rcond=None)[0]
            mz_exp = (Aterm / (freq_exp)) + (Bterm / np.power((freq_exp), 2)) + Cterm
            error = ((mz_exp - mz_calc) / mz_calc) * 1000000
            rms = np.sqrt(np.mean(error**2))
            std = np.std(error)
            if self.mass_spectrum.parameters.mass_spectrum.verbose_processing:
                print("%.2f Aterm,  %.2f Bterm" % (Aterm, Bterm))
                print(
                    "Quadratic Calibration %.2f RMS,  %.2f std,  %.2f Aterm,  %.2f Bterm "
                    % (rms, std, Aterm, Bterm)
                )
            if rms < last_rms:
                last_rms = rms
                freq_exp = (
                    Aterm
                    + np.sqrt(np.power(-Aterm, 2) - (4 * Cterm * (mz_exp - Bterm)))
                ) / (2 * mz_exp)

                mz_domain = (
                    (Aterm / (self.freq_exp_ms))
                    + (Bterm / np.power((self.freq_exp_ms), 2))
                    + Cterm
                )
                self.recal_mass_spec(mz_domain, Aterm, Bterm, Cterm)
                if not iteration:
                    break
            else:
                break

    def ledford_calibration(self, iteration: bool = False):
        """Perform Ledford calibration.

        Parameters
        ----------
        iteration : bool, optional
            Flag to perform iterative calibration, by default False.

        """
        mz_calc = self.mz_calc
        freq_exp = self.freq_exp
        mz_exp = self.mz_exp

        error = ((mz_exp - self.mz_calc) / self.mz_calc) * 1000000
        last_rms = np.sqrt(np.mean(error**2))
        while True:
            matrix = np.vstack([1 / freq_exp, 1 / np.power(freq_exp, 2)]).T
            Aterm, Bterm = np.linalg.lstsq(matrix, self.mz_calc, rcond=None)[0]

            mz_exp = (Aterm / (freq_exp)) + (Bterm / np.power((freq_exp), 2))
            error = ((mz_exp - mz_calc) / mz_calc) * 1000000
            rms = np.sqrt(np.mean(error**2))
            std = np.std(error)
            if self.mass_spectrum.parameters.mass_spectrum.verbose_processing:
                print("%.2f Aterm,  %.2f Bterm" % (Aterm, Bterm))
                print(
                    "Ledford Calibration %.2f RMS,  %.2f std,  %.2f Aterm,  %.2f Bterm "
                    % (rms, std, Aterm, Bterm)
                )
            if rms < last_rms:
                last_rms = rms
                freq_exp = (
                    Aterm + np.sqrt(np.power(-Aterm, 2) - (4 * mz_exp - Bterm))
                ) / (2 * mz_exp)
                mz_domain = (Aterm / (self.freq_exp_ms)) + (
                    Bterm / np.power((self.freq_exp_ms), 2)
                )
                self.recal_mass_spec(mz_domain, Aterm, Bterm, 0)
                if not iteration:
                    break
            else:
                break

    def step_fit(self, steps: int = 4):
        """Perform step fit calibration.

        Parameters
        ----------
        steps : int, optional
            Number of steps for step fit calibration, by default 4.

        """

        def f_to_mz(f, A, B, C, a):
            return (A / f) + (B / np.power(f, 2)) + (C * a / np.power(f, 2))

        def mz_to_f(m, A, B, C):
            return -A - m / B

        tuple_indexes = [
            (i, i + steps) for i in range(0, len(self.selected_mspeaks) - steps, steps)
        ]

        for current_index, tuple_index in enumerate(tuple_indexes):
            mspeak_ii, mspeak_fi = tuple_index
            freq_exp = list()
            mz_calc = list()
            mz_exp = list()
            abu = list()

            for i in range(mspeak_ii, mspeak_fi + 1):
                best_formula = self.selected_mspeaks[i].best_molecular_formula_candidate

                freq_exp.append(self.selected_mspeaks[i].freq_exp)
                mz_calc.append(best_formula.mz_calc)
                mz_exp.append(self.selected_mspeaks[i].mz_exp)
                abu.append(self.selected_mspeaks[i].abundance)

            freq_exp = np.array(freq_exp)
            mz_calc = np.array(mz_calc)
            mz_exp = np.array(mz_exp)
            abu = np.array(abu)

            if current_index == len(tuple_indexes) - 1:
                ms_peaks_indexes = (self.selected_mspeaks[mspeak_ii].index, 0)

            elif current_index == 0:
                ms_peaks_indexes = (
                    len(self.mass_spectrum) - 1,
                    self.selected_mspeaks[mspeak_fi].index - 1,
                )
            else:
                ms_peaks_indexes = (
                    self.selected_mspeaks[mspeak_ii].index,
                    self.selected_mspeaks[mspeak_fi].index - 1,
                )

            final_index, start_index = ms_peaks_indexes

            matrix = np.vstack([1 / freq_exp, 1 / np.power(freq_exp, 2)]).T
            A, B = np.linalg.lstsq(matrix, mz_calc, rcond=None)[0]
            C = 0

            for mspeak in self.mass_spectrum[start_index:final_index]:
                mspeak.mz_cal = f_to_mz(mspeak.freq_exp, A, B, C, 0)

        self.mass_spectrum.is_calibrated = True
