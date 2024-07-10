# -*- coding: utf-8 -*-
"""
Created on Wed May 13 02:16:09 2020

@author: Will Kew
"""

# import modules
import pandas as pd
import numpy as np
import os
import csv
from io import BytesIO
from pathlib import Path

from s3path import S3Path
# import corems modules
from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from corems.encapsulation.factory.parameters import MSParameters
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas

# import scipy modules for calibration
from scipy.optimize import minimize

class MzDomainCalibration:
    """ MzDomainCalibration class for recalibrating mass spectra
    
    Parameters
    ----------
    mass_spectrum : CoreMS MassSpectrum Object
        The mass spectrum to be calibrated.
    ref_masslist : str
        The path to a reference mass list.
    mzsegment : tuple of floats, optional
        The mz range to recalibrate, or None. Used for calibration of specific parts of the mz domain at a time.
        Future work - allow multiple mzsegments to be passed.
    
    Attributes
    ----------
    mass_spectrum : CoreMS MassSpectrum Object
        The mass spectrum to be calibrated.
    mzsegment : tuple of floats or None
        The mz range to recalibrate, or None.
    ref_mass_list_path : str or Path
        The path to the reference mass list.
    
    Methods
    -------
    * run(). 
        Main function to run this class.  
    * load_ref_mass_list(). 
        Load reference mass list (Bruker format).  
    * gen_ref_mass_list_from_assigned(min_conf=0.7). 
        Generate reference mass list from assigned masses.  
    * find_calibration_points(df_ref, calib_ppm_error_threshold=(-1, 1), calib_snr_threshold=5). 
        Find calibration points in the mass spectrum based on the reference mass list.  
    * robust_calib(param, cal_peaks_mz, cal_refs_mz, order=1). 
        Recalibration function.  
    * recalibrate_mass_spectrum(cal_peaks_mz, cal_refs_mz, order=1, diagnostic=False). 
        Main recalibration function which uses a robust linear regression.
    

    """

    def __init__(self, mass_spectrum, ref_masslist,mzsegment=None):
        
        self.mass_spectrum = mass_spectrum
        self.mzsegment = mzsegment

        # define reference mass list - bruker .ref format
        self.ref_mass_list_path = ref_masslist
        if self.mass_spectrum.percentile_assigned(verbose=False)[0]!=0:
            print('Warning: calibrating spectra which have already been assigned may yield erroneous results')
        self.mass_spectrum.mz_cal = self.mass_spectrum.mz_exp    
        self.mass_spectrum.mz_cal_profile = self.mass_spectrum._mz_exp  
        
        print("MS Obj loaded - "+str(len(mass_spectrum.mspeaks))+" peaks found.")

        print("MS Obj loaded - " + str(len(mass_spectrum.mspeaks)) + " peaks found.")

    def load_ref_mass_list(self):
        """ Load reference mass list (Bruker format)

        Loads in a reference mass list from a .ref file
        Note that some versions of Bruker's software produce .ref files with a different format. 
        As such, users may need to manually edit the .ref file in a text editor to ensure it is in the correct format. 
        CoreMS includes an example .ref file with the correct format for reference. 

        Returns
        -------
        df_ref : Pandas DataFrame
            reference mass list object.

        """
        refmasslist = Path(self.ref_mass_list_path) if isinstance(self.ref_mass_list_path, str) else self.ref_mass_list_path

        if not refmasslist.exists():
            raise FileExistsError("File does not exist: %s" % refmasslist)

        with refmasslist.open('r') as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read(1024))
            delimiter = dialect.delimiter

        if isinstance(refmasslist, S3Path):
            # data = self.file_location.open('rb').read()
            data = BytesIO(refmasslist.open('rb').read())

        else:
            data = refmasslist

        df_ref = pd.read_csv(data, sep=delimiter, header=None, skiprows=1)

        df_ref = df_ref.rename({0: 'Formula',
                                1: 'm/z',
                                2: 'Charge',
                                3: 'Form2'
                                }, axis=1)

        df_ref.sort_values(by='m/z', ascending=False)
        print("Reference mass list loaded - " + str(len(df_ref)) + " calibration masses loaded.")

        return df_ref

    def gen_ref_mass_list_from_assigned(self, min_conf : float=0.7):
        """ Generate reference mass list from assigned masses

        This function will generate a ref mass dataframe object from an assigned corems mass spec obj
        using assigned masses above a certain minimum confidence threshold.

        This function needs to be retested and check it is covered in the unit tests.

        Parameters
        ----------
        min_conf : float, optional
            minimum confidence score. The default is 0.7.

        Returns
        -------
        df_ref : Pandas DataFrame
            reference mass list - based on calculated masses.

        """
        #TODO this function needs to be retested and check it is covered in the unit tests
        df = self.mass_spectrum.to_dataframe()
        df = df[df['Confidence Score'] > min_conf]
        df_ref = pd.DataFrame(columns=['m/z'])
        df_ref['m/z'] = df['Calculated m/z']
        print("Reference mass list generated - " + str(len(df_ref)) + " calibration masses.")
        return df_ref

    def find_calibration_points(self, df_ref,
                                calib_ppm_error_threshold : tuple[float, float]=(-1, 1),
                                calib_snr_threshold : float=5):
        """Function to find calibration points in the mass spectrum 
        
        Based on the reference mass list.

        Parameters
        ----------
        df_ref : Pandas DataFrame
            reference mass list for recalibration.
        calib_ppm_error_threshold : tuple of floats, optional
            ppm error for finding calibration masses in the spectrum. The default is -1,1.
            Note: This is based on the calculation of ppm = ((mz_measure - mz_theoretical)/mz_theoretical)*1e6. 
                Some software does this the other way around and value signs must be inverted for that to work. 
        calib_snr_threshold : float, optional
            snr threshold for finding calibration masses in the spectrum. The default is 5.

        Returns
        -------
        cal_peaks_mz : list of floats
            masses of measured ions to use in calibration routine
        cal_refs_mz : list of floats
            reference mz values of found calibration points.

        """

        # This approach is much more efficient and expedient than the original implementation.
        peaks_mz = []
        for x in self.mass_spectrum.mspeaks:
            if x.signal_to_noise > calib_snr_threshold:
                if self.mzsegment:
                    if (min(self.mzsegment) <= x.mz_exp <= max(self.mzsegment)):
                        peaks_mz.append(x.mz_exp)
                else:
                    peaks_mz.append(x.mz_exp)
        peaks_mz = np.asarray(peaks_mz)
        
        cal_peaks_mz = []
        cal_refs_mz = []
        for mzref in df_ref['m/z']:
            tmp_peaks_mz = peaks_mz[abs(peaks_mz-mzref)<1]
            for mzmeas in tmp_peaks_mz:
                delta_mass = ((mzmeas-mzref)/mzref)*1e6
                if delta_mass < max(calib_ppm_error_threshold):
                    if delta_mass > min(calib_ppm_error_threshold):
                        cal_peaks_mz.append(mzmeas)
                        cal_refs_mz.append(mzref)

        # To remove entries with duplicated indices (reference masses matching multiple peaks)
        tmpdf = pd.Series(index = cal_refs_mz,data = cal_peaks_mz,dtype=float)
        tmpdf = tmpdf[~tmpdf.index.duplicated(keep=False)]

        cal_peaks_mz = list(tmpdf.values)
        cal_refs_mz = list(tmpdf.index)

        if False:
            min_calib_ppm_error = calib_ppm_error_threshold[0]
            max_calib_ppm_error = calib_ppm_error_threshold[1]
            df_raw = self.mass_spectrum.to_dataframe()

            df_raw = df_raw[df_raw['S/N'] > calib_snr_threshold]
            # optionally further subset that based on minimum S/N, RP, Peak Height
            # to ensure only valid points are utilized
            # in this example, only a S/N threshold is implemented.        
            imzmeas = []
            mzrefs = []

            for mzref in df_ref['m/z']:

                # find all peaks within a defined ppm error threshold
                tmpdf = df_raw[((df_raw['m/z']-mzref)/mzref)*1e6<max_calib_ppm_error]
                # Error is relative to the theoretical, so the divisor should be divisor

                tmpdf = tmpdf[((tmpdf['m/z']-mzref)/mzref)*1e6>min_calib_ppm_error]
                        
                # only use the calibration point if only one peak is within the thresholds
                # This may require some optimization of the threshold tolerances
                if len(tmpdf) == 1:
                    imzmeas.append(int(tmpdf.index.values))
                    mzrefs.append(mzref)

        # it is crucial the mass lists are in same order
        # corems likes to do masses from high to low.
        cal_refs_mz.sort(reverse=False)
        cal_peaks_mz.sort(reverse=False)
        print(str(len(cal_peaks_mz)) + " calibration points matched within thresholds.")
        return cal_peaks_mz, cal_refs_mz

    def robust_calib(self, param : list[float], 
                     cal_peaks_mz : list[float], cal_refs_mz : list[float], 
                     order : int=1):
        """ Recalibration function

        Computes the rms of m/z errors to minimize when calibrating.
        This is adapted from from spike.

        Parameters
        ----------
        param : list of floats
            generated by minimize function from scipy optimize.
        cal_peaks_mz : list of floats
            masses of measured peaks to use in mass calibration.
        cal_peaks_mz : list of floats
            reference mz values of found calibration points.
        order : int, optional
            order of the recalibration function. 1 = linear, 2 = quadratic. The default is 1.

        Returns
        -------
        rmserror : float
            root mean square mass error for calibration points.

        """
        Aterm = param[0]
        Bterm = param[1]
        try:
            Cterm = param[2]
        except IndexError:
            pass

        # get the mspeaks from the mass spectrum object which were calibration points
        #mspeaks = [self.mass_spectrum.mspeaks[x] for x in imzmeas]
        # get their calibrated mass values
        #mspeakmzs = [x.mz_cal for x in mspeaks]
        cal_peaks_mz = np.asarray(cal_peaks_mz)

        # linearz
        if order == 1:
            ref_recal_points = (Aterm * cal_peaks_mz) + Bterm
        # quadratic
        elif order == 2:
            ref_recal_points = (Aterm * (cal_peaks_mz)) + \
                (Bterm * np.power((cal_peaks_mz), 2) + Cterm)

        # sort both the calibration points (measured, recalibrated)
        ref_recal_points.sort()
        # and sort the calibration points (theoretical, predefined)
        cal_refs_mz.sort()

        # calculate the ppm error for each calibration point
        error = ((ref_recal_points - cal_refs_mz) / cal_refs_mz) * 1e6
        # calculate the root mean square error - this is our target to minimize
        rmserror = np.sqrt(np.mean(error**2))
        return rmserror

    def recalibrate_mass_spectrum(self, cal_peaks_mz : list[float], cal_refs_mz : list[float], 
                                        order : int=1, diagnostic : bool=False):

        """ Main recalibration function which uses a robust linear regression

        This function performs the recalibration of the mass spectrum object. 
        It iteratively applies 

        Parameters
        ----------
        cal_peaks_mz : list of float
            masses of measured peaks to use in mass calibration.
        cal_refs_mz : list of float
            reference mz values of found calibration points.
        order : int, optional
            order of the recalibration function. 1 = linear, 2 = quadratic. The default is 1.

        Returns
        -------
        mass_spectrum : CoreMS mass spectrum object
            Calibrated mass spectrum object


        Notes 
        -----
        This function is adapted, in part, from the SPIKE project [1,2] and is based on the robust linear regression method. 

        References
        ----------
        1.  Chiron L., Coutouly M-A., Starck J-P., Rolando C., Delsuc M-A. 
            SPIKE a Processing Software dedicated to Fourier Spectroscopies 
            https://arxiv.org/abs/1608.06777 (2016)
        2.  SPIKE - https://github.com/spike-project/spike 

        """
        # initialise parameters for recalibration
        # these are the 'Aterm, Bterm, Cterm'
        # as spectra are already freq->mz calibrated, these terms are very small
        # may be beneficial to formally separate them from the freq->mz terms
        if order == 1:
            Po = [1, 0]
        elif order == 2:
            Po = [1, 0, 0]

        if len(cal_peaks_mz) >= 2:
            if self.mzsegment: #If only part of the spectrum is to be recalibrated
                mz_exp_peaks = np.array([mspeak.mz_exp for mspeak in self.mass_spectrum])
                # Split the array into two parts - one to recailbrate, one to keep unchanged. 
                mz_exp_peaks_tocal = mz_exp_peaks[(mz_exp_peaks>=min(self.mzsegment)) & (mz_exp_peaks<=max(self.mzsegment))]
                mz_exp_peaks_unchanged = mz_exp_peaks[~(mz_exp_peaks>=min(self.mzsegment)) | ~(mz_exp_peaks<=max(self.mzsegment))]
                # TODO: - segmented calibration needs a way to better track the calibration args/values... 
                if not self.mass_spectrum.is_centroid:
                    mz_exp_profile = np.array(self.mass_spectrum.mz_exp_profile)
                    # Split the array into two parts - one to recailbrate, one to keep unchanged. 
                    mz_exp_profile_tocal = mz_exp_profile[(mz_exp_profile>=min(self.mzsegment)) & (mz_exp_profile<=max(self.mzsegment))]
                    mz_exp_profile_unchanged = mz_exp_profile[~(mz_exp_profile>=min(self.mzsegment)) | ~(mz_exp_profile<=max(self.mzsegment))]
            else: #if just recalibrating the whole spectrum
                mz_exp_peaks_tocal = np.array([mspeak.mz_exp for mspeak in self.mass_spectrum])
                if not self.mass_spectrum.is_centroid:
                    mz_exp_profile_tocal = np.array(self.mass_spectrum.mz_exp_profile)


            minimize_method = self.mass_spectrum.settings.calib_minimize_method
            res = minimize(self.robust_calib, Po, args=(cal_peaks_mz, cal_refs_mz, order), method=minimize_method)
            #TODO add quiet option to these print statements
            print("minimize function completed with RMS error of: {:0.3f} ppm".format(res['fun']))
            print("minimize function performed {:1d} fn evals and {:1d} iterations".format(res['nfev'], res['nit']))
            Pn = res.x

            #mz_exp_ms = np.array([mspeak.mz_exp for mspeak in self.mass_spectrum])

            if order == 1:
                mz_domain = (Pn[0] * mz_exp_peaks_tocal) + Pn[1]
                if not self.mass_spectrum.is_centroid:
                    mz_profile_calc = (Pn[0] * mz_exp_profile_tocal) + Pn[1]

            elif order == 2:
                mz_domain = (Pn[0] * (mz_exp_peaks_tocal)) + \
                    (Pn[1] * np.power((mz_exp_peaks_tocal), 2) + Pn[2])

                if not self.mass_spectrum.is_centroid:
                    mz_profile_calc = (Pn[0] * (mz_exp_profile_tocal)) + \
                        (Pn[1] * np.power((mz_exp_profile_tocal), 2) + Pn[2])

            if self.mzsegment:
                # Recombine the mass domains
                mz_domain = np.concatenate([mz_domain,mz_exp_peaks_unchanged])
                mz_domain.sort()
                if not self.mass_spectrum.is_centroid:
                    mz_profile_calc = np.concatenate([mz_profile_calc,mz_exp_profile_unchanged])
                    mz_profile_calc.sort()
                # Sort them 
                if mz_exp_peaks[0] > mz_exp_peaks[1]: #If originally descending mass order
                    mz_domain = mz_domain[::-1]
                    if not self.mass_spectrum.is_centroid:
                        mz_profile_calc = mz_profile_calc[::-1]

            self.mass_spectrum.mz_cal = mz_domain
            if not self.mass_spectrum.is_centroid:
                self.mass_spectrum.mz_cal_profile = mz_profile_calc

            self.mass_spectrum.calibration_order = order
            self.mass_spectrum.calibration_RMS = float(res['fun'])
            self.mass_spectrum.calibration_points = int(len(cal_refs_mz))
            self.mass_spectrum.calib_ref_mzs = cal_refs_mz
            #self.mass_spectrum.calib_meas_mzs = cal_refs_mz

            self.mass_spectrum.calibration_segment = self.mzsegment

            if diagnostic:
                return self.mass_spectrum,res
            return self.mass_spectrum
        else:
            print("Too few calibration points - aborting.")
            return self.mass_spectrum

    def run(self):
        """ Run the calibration routine
        
        This function runs the calibration routine.
        
        """
        calib_ppm_error_threshold = self.mass_spectrum.settings.calib_sn_threshold
        max_calib_ppm_error = self.mass_spectrum.settings.max_calib_ppm_error
        min_calib_ppm_error = self.mass_spectrum.settings.min_calib_ppm_error
        calib_pol_order = self.mass_spectrum.settings.calib_pol_order

        # load reference mass list
        df_ref = self.load_ref_mass_list()

        # find calibration points
        cal_peaks_mz, cal_refs_mz = self.find_calibration_points(df_ref,
                                                       calib_ppm_error_threshold=(min_calib_ppm_error,
                                                                                  max_calib_ppm_error),
                                                       calib_snr_threshold=calib_ppm_error_threshold)
        if len(cal_peaks_mz)==2:
            self.mass_spectrum.settings.calib_pol_order = 1
            calib_pol_order = 1
            print('Only 2 calibration points found, forcing a linear recalibration')
        elif len(cal_peaks_mz)<2:
            print('Too few calibration points found, function will fail')
        self.recalibrate_mass_spectrum(cal_peaks_mz, cal_refs_mz, order=calib_pol_order)
