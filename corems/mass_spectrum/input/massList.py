__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

from corems.mass_spectrum.input.baseClass import MassListBaseClass
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile, MassSpecCentroid
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula
from corems.encapsulation.constant import Labels, Atoms
from corems.encapsulation.factory.processingSetting  import DataInputSetting

class ReadCoremsMasslist(MassListBaseClass):
    '''
    The ReadCoremsMasslist object reads processed mass list data types
    and returns the mass spectrum obj with the molecular formula obj

    **Only available for centroid mass spectrum type: it will ignore the parameter **isCentroid** 
    Please see MassListBaseClass for more details
    
    '''
        
    def get_mass_spectrum(self, scan_number=0, time_index=-1, auto_process=True, loadSettings=True):
        
        dataframe = self.get_dataframe()
       
        if not set(['H/C', 'O/C', 'Heteroatom Class', 'Ion Type', 'Is Isotopologue']).issubset(dataframe.columns):
            raise ValueError("%s it is not a valid CoreMS file" % str(self.file_location))
        
        self.check_columns(dataframe.columns)
        
        dataframe.rename(columns=self.parameters.header_translate, inplace=True)
 
        polarity = dataframe['Ion Charge'].values[0]

        output_parameters = self.get_output_parameters(polarity)

        mass_spec_obj = MassSpecCentroid(dataframe.to_dict(orient='list'), output_parameters)

        if loadSettings: self.load_settings(mass_spec_obj, output_parameters)

        self.add_molecular_formula(mass_spec_obj, dataframe)
        
        return mass_spec_obj

    def add_molecular_formula(self, mass_spec_obj, dataframe):

        # check if is coreMS file
        if 'Is Isotopologue' in dataframe:

            mz_exp_df = dataframe[Labels.mz].astype(float)
            # formula_df = dataframe.loc[:, 'C':].fillna(0)
            # \.replace({b'nan':0})
            formula_df = dataframe[dataframe.columns.intersection(Atoms.atoms_order)]
            formula_df.fillna(0, inplace=True)
            formula_df.replace(b'nan', 0, inplace=True)

            ion_type_df = dataframe["Ion Type"]
            ion_charge_df = dataframe["Ion Charge"]
            is_isotopologue_df = dataframe['Is Isotopologue']

        mass_spec_mz_exp_list = mass_spec_obj.mz_exp

        for df_index, mz_exp in enumerate(mz_exp_df):

            ms_peak_index = list(mass_spec_mz_exp_list).index(float(mz_exp))

            if 'Is Isotopologue' in dataframe:

                atoms = list(formula_df.columns.astype(str))
                counts = list(formula_df.iloc[df_index].astype(int))

                formula_list = [sub[item] for item in range(len(atoms))
                                for sub in [atoms, counts]]
            if sum(counts) > 0:

                ion_type = str(Labels.ion_type_translate.get(ion_type_df[df_index]))
                mfobj = MolecularFormula(formula_list, int(ion_charge_df[df_index]), mspeak_parent=mass_spec_obj[ms_peak_index] , ion_type=ion_type)
                mfobj.is_isotopologue = bool(is_isotopologue_df[df_index])
                mass_spec_obj[ms_peak_index].add_molecular_formula(mfobj)


class ReadMassList(MassListBaseClass):
    
    '''
    The ReadCoremsMasslist object reads unprocessed mass list data types
    and returns the mass spectrum obj 
    See MassListBaseClass for details
    
    '''

    def get_mass_spectrum(self, polarity, scan=0, auto_process=True, loadSettings=True):
        '''
         The MassListBaseClass object reads mass list data types and returns the mass spectrum obj

        Parameters
        ----------
        polarity: int 
            +1 or -1 
        '''
        #delimiter = "  " or " " or  "," or "\t" etc  
        
        if self.isCentroid:

            dataframe = self.get_dataframe()

            self.check_columns(dataframe.columns)
                
            self.clean_data_frame(dataframe)
            
            dataframe.rename(columns=self.parameters.header_translate, inplace=True)
            
            output_parameters = self.get_output_parameters(polarity)

            mass_spec = MassSpecCentroid(dataframe.to_dict(orient='list'), output_parameters)
            
            if loadSettings: self.load_settings(mass_spec, output_parameters)
            
            return mass_spec

        else:

            dataframe = self.get_dataframe()

            self.check_columns(dataframe.columns)
            
            output_parameters = self.get_output_parameters(polarity)

            self.clean_data_frame(dataframe)
            
            dataframe.rename(columns=self.parameters.header_translate, inplace=True)

            mass_spec = MassSpecProfile(dataframe.to_dict(orient='list'), output_parameters, auto_process=auto_process)

            if loadSettings: self.load_settings(mass_spec, output_parameters)

            return mass_spec
    
    