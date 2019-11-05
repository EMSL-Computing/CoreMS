__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

from pandas import read_csv

from corems.mass_spectrum.input.baseClass import MassListBaseClass
from corems.mass_spectrum.factory.MassSpectrumClasses import MassSpecProfile, MassSpecCentroid
from corems.molecular_formula.factory.MolecularFormulaFactory import MolecularFormula
from corems.encapsulation.constant import Labels

class ReadCoremsMasslist(MassListBaseClass):

    def __init__(self, file_location, delimiter=","):
        
        '''
        Constructor
        '''
        
        self.dataframe = read_csv(file_location, delimiter=delimiter, engine='python')

        if not 'Ion Charge' in self.dataframe:
            raise ValueError("%s it is not a valid CoreMS file" % str(file_location))
        
        polarity = self.dataframe['Ion Charge'].values[0]

        super().__init__(file_location, polarity, delimiter=delimiter,isCentroid=True)
        
        
    def get_mass_spectrum(self, auto_process=True):
        
        #delimiter = "  " or " " or  "," or "\t" etc  
        
        self.check_columns(self.dataframe.columns)
        
        self.dataframe.rename(columns=self.name_dict, inplace=True)
 
        output_parameters = self.get_output_parameters(self.polarity)

        mass_spec_obj = MassSpecCentroid(self.dataframe, output_parameters, auto_process=auto_process)

        self.add_molecular_formula(mass_spec_obj)
        
        return mass_spec_obj

    def add_molecular_formula(self, mass_spec_obj):
        
        #check if is coreMS file
        if 'Is Isotopologue' in self.dataframe:
            
            
            mz_exp_df = self.dataframe["m/z"]
            formula_df = self.dataframe.loc[:, 'C':].fillna(0)
            ion_type_df =  self.dataframe["Ion Type"]
            ion_charge_df = self.dataframe["Ion Charge"]
            is_isotopologue_df = self.dataframe['Is Isotopologue']
        
        mass_spec_mz_exp_list = mass_spec_obj.mz_exp
    
        for df_index, mz_exp in enumerate(mz_exp_df):
            
            counts = 0

            ms_peak_index = list(mass_spec_mz_exp_list).index(mz_exp)
            
            if 'Is Isotopologue' in self.dataframe:
                
                atoms = list(formula_df.columns)
                counts = list(formula_df.iloc[df_index])

                formula_list = [sub[item] for item in range(len(atoms)) 
                        for sub in [atoms, counts]] 
            if sum(counts) > 0:
                

                ion_type = Labels.ion_type_translate.get(ion_type_df[df_index])
                mfobj = MolecularFormula(formula_list, ion_charge_df[df_index], mass_spec_obj[ms_peak_index].mz_exp, ion_type=ion_type)
                mfobj.is_isotopologue = is_isotopologue_df[df_index]
                mass_spec_obj[ms_peak_index].add_molecular_formula(mfobj)


class ReadMassList(MassListBaseClass):
    '''
    The ReadMassList object contains lots of MassSpectrum objs

    Parameters
    ----------
    arg : str
        The arg is used for ...
    *args
        The variable arguments are used for ...
    **kwargs
        The keyword arguments are used for ...

    Attributes
    ----------
    arg : str
        This is where we store arg,
    '''

    def get_mass_spectrum(self, auto_process=True):
        
        #delimiter = "  " or " " or  "," or "\t" etc  
        
        dataframe = read_csv(self.file_location, delimiter=self.delimiter, engine='python')
        
        self.check_columns(dataframe.columns)
            
        self.clean_data_frame(dataframe)
        
        dataframe.rename(columns=self.name_dict, inplace=True)
 
        output_parameters = self.get_output_parameters(self.polarity)
            
        if self.isCentroid:
            
            return MassSpecCentroid(dataframe, output_parameters, auto_process=auto_process)

        else:
            
            output_parameters[Labels.label] = Labels.bruker_profile
    
            return MassSpecProfile(dataframe, output_parameters, auto_process=auto_process)
    
    