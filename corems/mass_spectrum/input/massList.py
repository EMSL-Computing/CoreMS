__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

import numpy as np

from corems.encapsulation.constant import Atoms
from corems.mass_spectrum.input.baseClass import MassListBaseClass
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile, MassSpecCentroid
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula, MolecularFormulaIsotopologue
from corems.encapsulation.constant import Labels, Atoms
from corems.encapsulation.factory.processingSetting  import DataInputSetting

class ReadCoremsMasslist(MassListBaseClass):
    """
    The ReadCoremsMasslist object reads processed mass list data types
    and returns the mass spectrum obj with the molecular formula obj

    **Only available for centroid mass spectrum type:** it will ignore the parameter **isCentroid** 
    Please see MassListBaseClass for more details

    """

    def get_mass_spectrum(self, loadSettings:bool =True) -> MassSpecCentroid:
        """
        Get the mass spectrum object from the processed mass list data.

        Parameters
        ----------
        loadSettings : bool, optional
            Whether to load the settings for the mass spectrum. Default is True.

        Returns
        -------
        MassSpecCentroid
            The mass spectrum object.

        Raises
        ------
        ValueError
            If the input file is not a valid CoreMS file.
        """

        dataframe = self.get_dataframe()

        if not set(['H/C', 'O/C', 'Heteroatom Class', 'Ion Type', 'Is Isotopologue']).issubset(dataframe.columns):
            raise ValueError("%s it is not a valid CoreMS file" % str(self.file_location))

        self.check_columns(dataframe.columns)

        dataframe.rename(columns=self.parameters.header_translate, inplace=True)

        polarity = dataframe['Ion Charge'].values[0]

        output_parameters = self.get_output_parameters(polarity)

        mass_spec_obj = MassSpecCentroid(dataframe.to_dict(orient='list'), output_parameters)

        if loadSettings is True:
            self.load_settings(mass_spec_obj, output_parameters)

        self.add_molecular_formula(mass_spec_obj, dataframe)

        return mass_spec_obj

    def add_molecular_formula(self, mass_spec_obj, dataframe):
        """
        Add molecular formula information to the mass spectrum object.

        Parameters
        ----------
        mass_spec_obj : MassSpecCentroid
            The mass spectrum object to add the molecular formula to.
        dataframe : pandas.DataFrame
            The processed mass list data.
            
        """

        # check if is coreMS file
        if 'Is Isotopologue' in dataframe:

            mz_exp_df = dataframe[Labels.mz].astype(float)
            formula_df = dataframe[dataframe.columns.intersection(Atoms.atoms_order)].copy()
            formula_df.fillna(0, inplace=True)
            formula_df.replace(b'nan', 0, inplace=True)

            ion_type_df = dataframe["Ion Type"]
            ion_charge_df = dataframe["Ion Charge"]
            is_isotopologue_df = dataframe['Is Isotopologue']
            if 'Adduct' in dataframe:
                adduct_df = dataframe['Adduct']
            else:
                adduct_df = None

        mass_spec_mz_exp_list = mass_spec_obj.mz_exp

        for df_index, mz_exp in enumerate(mz_exp_df):

            counts = 0

            ms_peak_index = list(mass_spec_mz_exp_list).index(float(mz_exp))

            if 'Is Isotopologue' in dataframe:

                atoms = list(formula_df.columns.astype(str))
                counts = list(formula_df.iloc[df_index].astype(int))

                formula_dict = dict(zip(atoms, counts))
            if sum(counts) > 0:

                ion_type = str(Labels.ion_type_translate.get(ion_type_df[df_index]))
                if adduct_df is not None:
                    adduct_atom = str(adduct_df[df_index])
                    if adduct_atom == 'None':
                        adduct_atom = None
                else:
                    adduct_atom = None

                # If not isotopologue, cast as MolecularFormula
                if not bool(int(is_isotopologue_df[df_index])):
                    mfobj = MolecularFormula(
                        formula_dict, int(ion_charge_df[df_index]), 
                        mspeak_parent=mass_spec_obj[ms_peak_index] , 
                        ion_type=ion_type, adduct_atom=adduct_atom
                        )
                    
                # if is isotopologue, recast as MolecularFormulaIsotopologue
                if bool(int(is_isotopologue_df[df_index])):

                    # First make a MolecularFormula object for the parent so we can get probabilities etc
                    formula_list_parent = {}
                    for atom in formula_dict:
                        if atom in Atoms.isotopes.keys():
                            formula_list_parent[atom] = formula_dict[atom]
                        else:
                            # remove any numbers from the atom name to cast as a mono-isotopic atom
                            atom_mono = atom.strip('0123456789')
                            if atom_mono in Atoms.isotopes.keys():
                                formula_list_parent[atom_mono] = formula_list_parent[atom_mono]+formula_dict[atom]
                            else:
                                print(f"Atom {atom} not in Atoms.atoms_order")
                    mono_index = int(dataframe.iloc[df_index]['Mono Isotopic Index'])
                    mono_mfobj = MolecularFormula(
                        formula_list_parent, 
                        int(ion_charge_df[df_index]), 
                        mspeak_parent=mass_spec_obj[mono_index], 
                        ion_type=ion_type, 
                        adduct_atom=adduct_atom
                        )
                    
                    # Next, generate isotopologues from the parent
                    isos = list(
                        mono_mfobj.isotopologues(
                        min_abundance = mass_spec_obj[df_index].abundance*0.1, 
                        current_mono_abundance = mass_spec_obj[mono_index].abundance, 
                        dynamic_range = mass_spec_obj.dynamic_range
                         )
                    )

                    # Finally, find the isotopologue that matches the formula_dict
                    matched_isos = isos
                    for iso in isos:
                        if set(iso.atoms) == set(formula_dict.keys()):
                            # Check the values of the atoms match
                            if all([iso[atom] == formula_dict[atom] for atom in formula_dict]):
                                matched_isos = [iso]
                    if len(matched_isos) > 1:
                        raise ValueError("More than one isotopologue matched the formula_dict: {matched_isos}")
                    if len(matched_isos) == 0:
                        raise ValueError("No isotopologue matched the formula_dict")
                    mfobj = matched_isos[0]        

                    # Add the mono isotopic index, confidence score and isotopologue similarity    
                    mfobj.mspeak_index_mono_isotopic = int(dataframe.iloc[df_index]['Mono Isotopic Index'])
                
                # Add the confidence score and isotopologue similarity and average MZ error score
                if 'm/z Error Score' in dataframe:
                    mfobj._mass_error_average_score = float(dataframe.iloc[df_index]['m/z Error Score'])
                if 'Confidence Score' in dataframe:
                    mfobj._confidence_score = float(dataframe.iloc[df_index]['Confidence Score'])
                if 'Isotopologue Similarity' in dataframe:
                    mfobj._isotopologue_similarity = float(dataframe.iloc[df_index]['Isotopologue Similarity'])
                mass_spec_obj[ms_peak_index].add_molecular_formula(mfobj)


class ReadMassList(MassListBaseClass):
    """
    The ReadMassList object reads unprocessed mass list data types
    and returns the mass spectrum object.

    Parameters
    ----------
    MassListBaseClass : class
        The base class for reading mass list data types.

    Methods
    -------
    * get_mass_spectrum(polarity, scan=0, auto_process=True, loadSettings=True). Reads mass list data types and returns the mass spectrum object.

    """

    def get_mass_spectrum(self, polarity:int, scan:int=0, auto_process:bool=True, loadSettings:bool=True):
        """
        Reads mass list data types and returns the mass spectrum object.

        Parameters
        ----------
        polarity : int
            The polarity of the mass spectrum (+1 or -1).
        scan : int, optional
            The scan number of the mass spectrum (default is 0).
        auto_process : bool, optional
            Flag indicating whether to automatically process the mass spectrum (default is True).
        loadSettings : bool, optional
            Flag indicating whether to load settings for the mass spectrum (default is True).

        Returns
        -------
        mass_spec : MassSpecCentroid or MassSpecProfile
            The mass spectrum object.

        """

        # delimiter = "  " or " " or  "," or "\t" etc

        if self.isCentroid:

            dataframe = self.get_dataframe()

            self.check_columns(dataframe.columns)

            self.clean_data_frame(dataframe)

            dataframe.rename(columns=self.parameters.header_translate, inplace=True)

            output_parameters = self.get_output_parameters(polarity)

            mass_spec = MassSpecCentroid(dataframe.to_dict(orient='list'), output_parameters)

            if loadSettings:
                self.load_settings(mass_spec, output_parameters)

            return mass_spec

        else:

            dataframe = self.get_dataframe()

            self.check_columns(dataframe.columns)

            output_parameters = self.get_output_parameters(polarity)

            self.clean_data_frame(dataframe)

            dataframe.rename(columns=self.parameters.header_translate, inplace=True)

            mass_spec = MassSpecProfile(dataframe.to_dict(orient='list'), output_parameters, auto_process=auto_process)

            if loadSettings:
                self.load_settings(mass_spec, output_parameters)

            return mass_spec
    

class ReadBrukerXMLList(MassListBaseClass):
    """
    The ReadBrukerXMLList object reads Bruker XML objects
    and returns the mass spectrum object.
    See MassListBaseClass for details

    Parameters
    ----------
    MassListBaseClass : class
        The base class for reading mass list data types and returning the mass spectrum object.

    Methods
    -------
    * get_mass_spectrum(polarity: bool = None, scan: int = 0, auto_process: bool = True, loadSettings: bool = True). Reads mass list data types and returns the mass spectrum object.

    """

    def get_mass_spectrum(self, polarity: bool = None, scan: int = 0, auto_process: bool = True, loadSettings: bool = True):
        """
        Reads mass list data types and returns the mass spectrum object.

        Parameters
        ----------
        polarity : bool, optional
            The polarity of the mass spectrum. Can be +1 or -1. If not provided, it will be determined from the XML file.
        scan : int, optional
            The scan number of the mass spectrum. Default is 0.
        auto_process : bool, optional
            Whether to automatically process the mass spectrum. Default is True.
        loadSettings : bool, optional
            Whether to load the settings for the mass spectrum. Default is True.

        Returns
        -------
        mass_spec : MassSpecCentroid
            The mass spectrum object representing the centroided mass spectrum.
        """
        # delimiter = "  " or " " or  "," or "\t" etc

        if polarity == None:
            polarity = self.get_xml_polarity()
        dataframe = self.get_dataframe()

        self.check_columns(dataframe.columns)

        self.clean_data_frame(dataframe)

        dataframe.rename(columns=self.parameters.header_translate, inplace=True)

        output_parameters = self.get_output_parameters(polarity)

        mass_spec = MassSpecCentroid(dataframe.to_dict(orient='list'), output_parameters)

        if loadSettings:
            self.load_settings(mass_spec, output_parameters)

        return mass_spec

