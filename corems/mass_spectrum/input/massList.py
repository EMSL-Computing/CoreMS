__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

from corems.mass_spectrum.input.baseClass import MassListBaseClass
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile, MassSpecCentroid
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula
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

                formula_list = [sub[item] for item in range(len(atoms))
                                for sub in [atoms, counts]]
            if sum(counts) > 0:

                ion_type = str(Labels.ion_type_translate.get(ion_type_df[df_index]))
                if adduct_df is not None:
                    adduct_atom = str(adduct_df[df_index])
                else:
                    adduct_atom = None
                mfobj = MolecularFormula(formula_list, int(ion_charge_df[df_index]), mspeak_parent=mass_spec_obj[ms_peak_index] , ion_type=ion_type, adduct_atom=adduct_atom)
                mfobj.is_isotopologue = bool(is_isotopologue_df[df_index])
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

