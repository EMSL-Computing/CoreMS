# -*- coding: utf-8 -*-
"""
Created on Wed May 13 02:16:09 2020

@author: Will Kew
"""

### import modules
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

    def __init__(self, mass_spectrum, ref_masslist):
        
        self.mass_spectrum = mass_spectrum

        # define reference mass list - bruker .ref format
        self.ref_mass_list_path = ref_masslist

        self.mass_spectrum.mz_cal = self.mass_spectrum.mz_exp    
    
        print("MS Obj loaded - "+str(len(mass_spectrum.mspeaks))+" peaks found.")


    def load_ref_mass_list(self, refmasslist):
        """
        function to load in bruker mass list format into a dataframe


        Parameters
        ----------
        refmasslist : str
            full path to reference mass list (.ref) to import.

        Returns
        -------
        df_ref : pandas dataframe
            reference mass list object.

        """
        refmasslist = Path(refmasslist) if isinstance(refmasslist, str) else refmasslist
        
        if not refmasslist.exists():
            raise FileExistsError("File does not exist: %s" % refmasslist)


        with refmasslist.open('r') as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read(1024))
            delimiter=dialect.delimiter

        if isinstance(refmasslist, S3Path):
            # data = self.file_location.open('rb').read()
            data = BytesIO(refmasslist.open('rb').read())
        
        else:
            data = refmasslist

        df_ref = pd.read_csv(data,sep=delimiter,header=None,skiprows=1)
        
        df_ref = df_ref.rename({0:'Formula',
                            1:'m/z',
                            2:'Charge',
                            3:'Form2'},axis=1)
        
        df_ref.sort_values(by='m/z',ascending=False)
        print("Reference mass list loaded - "+str(len(df_ref))+ " calibration masses loaded.")
        
        return df_ref

    def gen_ref_mass_list_from_assigned(self, mass_spectrum,min_conf=0.7):
        
        """
        This function will generate a ref mass dataframe object 
        from an assigned corems mass spec obj
        using assigned masses above a certain minimum confidence threshold

        Parameters
        ----------
        mass_spectrum : corems mass spec obj
            processed assigned mass spectrum.
        min_conf : float, optional
            minimum confidence score. The default is 0.7.

        Returns
        -------
        df_ref : pandas dataframe
            reference mass list - based on calculated masses.

        """
        df = mass_spectrum.to_dataframe()
        df = df[df['Isotopologue Similarity']>min_conf]
        df_ref = pd.DataFrame(columns=['m/z'])
        df_ref['m/z'] = df['Calculated m/z']
        print("Reference mass list generated - "+str(len(df_ref))+ " calibration masses.")
        return df_ref

    def find_calibration_points(self, mass_spectrum, df_ref,
                                calib_ppm_error_threshold=(-1,1),
                                calib_snr_threshold=5):
        """
        function to find calibration points in the mass spectrum based on the reference
        mass list. 
        
        Parameters
        ----------
        mass_spectrum : corems mass spec obj
            mass spectrum to be calibrated.
        df_ref : pandas dataframe
            reference mass list for recalibration.
        calib_ppm_error_threshold : float, optional
            ppm error for finding calibration masses in the spectrum. The default is 1.
        calib_snr_threshold : float, optional
            snr threshold for finding calibration masses in the spectrum. The default is 5.

        Returns
        -------
        imzmeas : list of ints
            indexes of measured peaks to use in mass calibration.
        mzrefs : list of floats
            reference mz values of found calibration points.

        """

        min_calib_ppm_error = calib_ppm_error_threshold[0]
        max_calib_ppm_error = calib_ppm_error_threshold[1]

        df_raw = mass_spectrum.to_dataframe()
        
        imzmeas = [] 
        mzrefs = [] 
        
        for mzref in df_ref['m/z']:
            
            # find all peaks within a defined ppm error threshold
            tmpdf = df_raw[((df_raw['m/z']-mzref)/df_raw['m/z'])*1e6<max_calib_ppm_error]

            tmpdf = tmpdf[((tmpdf['m/z']-mzref)/tmpdf['m/z'])*1e6>min_calib_ppm_error]
            
            #tmpdf = df_raw[(abs(df_raw['m/z']-mzref)/df_raw['m/z'])*1e6 < calib_ppm_error_threshold]
            
            # optionally further subset that based on minimum S/N, RP, Peak Height
            # to ensure only valid points are utilized
            # in this example, only a S/N threshold is implemented.
            tmpdf = tmpdf[tmpdf['S/N']>calib_snr_threshold]
            
            # only use the calibration point if only one peak is within the thresholds
            # This may require some optimization of the threshold tolerances
            if len(tmpdf)==1: 
                imzmeas.append(int(tmpdf.index.values))
                mzrefs.append(mzref)
                
        # it is crucial the mass lists are in same order
        # corems likes to do masses from high to low.
        mzrefs.sort(reverse=True) 
        print(str(len(imzmeas))+" calibration points matched within thresholds.")
        return imzmeas, mzrefs

    def robust_calib(self, param, imzmeas, mzrefs, mass_spectrum, order=1):
        """
        computes the rms of m/z errors to minimize when calibrating
        This is adapted from from spike

        Parameters
        ----------
        param : list of floats
            generated by minimize function from scipy optimize.
        imzmeas : list of int
            indexes of measured peaks to use in mass calibration.
        mzrefs : list of floats
            reference mz values of found calibration points.
        mass_spectrum : corems mass spec obj
            mass spectrum to be calibrated.
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
        mspeaks = [mass_spectrum.mspeaks[x] for x in imzmeas]
        # get their calibrated mass values
        mspeakmzs = [x.mz_cal for x in mspeaks]
        mspeakmzs = np.array(mspeakmzs)
        
        # linear
        if order == 1:
            ref_recal_points = (Aterm * mspeakmzs) + Bterm
        # quadratic
        elif order == 2:
            ref_recal_points = (Aterm * (mspeakmzs)) + \
            (Bterm * np.power((mspeakmzs), 2) + Cterm)
        
        # sort both the calibration points (measured, recalibrated)
        ref_recal_points.sort()
        # and sort the calibration points (theoretical, predefined)
        mzrefs.sort()
        
        #calculate the ppm error for each calibration point
        error = ((ref_recal_points-mzrefs)/ref_recal_points) * 1000000
        # calculate the root mean square error - this is our target to minimize
        rmserror = np.sqrt(np.mean(error**2))
        return rmserror

    def recalibrate_mass_spectrum(self, mass_spectrum, imzmeas,mzrefs,order=1):
        
        """
        function to recalibrate the mass spectrum object
        
        Parameters
        ----------
        mass_spectrum : corems mass spec obj
            mass spectrum to be calibrated.
        imzmeas : list of int
            indexes of measured peaks to use in mass calibration.
        mzrefs : list of float
            reference mz values of found calibration points.
        order : int, optional
            order of the recalibration function. 1 = linear, 2 = quadratic. The default is 1.

        Returns
        -------
        mass_spectrum : corems mass spec obj
            mass spectrum to be calibrated.

        """
        # initialise parameters for recalibration
        # these are the 'Aterm, Bterm, Cterm'
        # as spectra are already freq->mz calibrated, these terms are very small
        # may be beneficial to formally separate them from the freq->mz terms
        if order == 1:
            Po = [1,0]
        elif order == 2:
            Po = [1,0,0]
        
        if len(imzmeas) > 2:

            minimize_method = mass_spectrum.settings.calib_minimize_method
            res =  minimize(self.robust_calib, Po, args=(imzmeas, mzrefs, mass_spectrum,order), method=minimize_method)
            
            print("minimize function completed with RMS error of: {:0.3f} ppm".format(res['fun']))
            print("minimize function performed {:1d} fn evals and {:1d} iterations".format(res['nfev'],res['nit']))
            Pn = res.x
            
            mz_exp_ms = np.array(
                [mspeak.mz_exp for mspeak in mass_spectrum])  
            
            if order == 1:
                mz_domain = (Pn[0] * mz_exp_ms) + Pn[1]
                if not mass_spectrum.is_centroid:
                    mz_profile_calc = (Pn[0] * mass_spectrum.mz_exp_profile) + Pn[1]
                
            elif order == 2:
                mz_domain = (Pn[0] * (mz_exp_ms)) + \
                    (Pn[1] * np.power((mz_exp_ms), 2) + Pn[2])

                if not mass_spectrum.is_centroid:
                    mz_profile_calc = (Pn[0] * (mass_spectrum.mz_exp_profile)) + \
                            (Pn[1] * np.power((mass_spectrum.mz_exp_profile), 2) + Pn[2])
            
            mass_spectrum.mz_cal = mz_domain
            if not mass_spectrum.is_centroid:
                mass_spectrum.mz_cal_profile = mz_profile_calc

            mass_spectrum.calibration_order = order
            mass_spectrum.calibration_points = len(mzrefs)
            mass_spectrum.calibration_RMS = float(res['fun'])
        
        return mass_spectrum

    def run(self):
        
        calib_ppm_error_threshold = self.mass_spectrum.settings.calib_sn_threshold
        max_calib_ppm_error = self.mass_spectrum.settings.max_calib_ppm_error
        min_calib_ppm_error = self.mass_spectrum.settings.min_calib_ppm_error
        calib_pol_order = self.mass_spectrum.settings.calib_pol_order

        # load reference mass list
        df_ref = self.load_ref_mass_list(self.ref_mass_list_path)
        
        # find calibration points
        imzmeas, mzrefs = self.find_calibration_points(self.mass_spectrum, df_ref,
                                    calib_ppm_error_threshold=(min_calib_ppm_error, max_calib_ppm_error),
                                    calib_snr_threshold=calib_ppm_error_threshold)

        self.recalibrate_mass_spectrum(self.mass_spectrum, imzmeas, mzrefs, order=calib_pol_order)
  
    