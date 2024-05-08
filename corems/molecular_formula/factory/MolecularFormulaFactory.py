from corems.molecular_formula.calc.MolecularFormulaCalc import MolecularFormulaCalc
from corems.encapsulation.constant import Atoms, Labels

import re

__author__ = "Yuri E. Corilo"
__date__ = "Jun 24, 2019"

class MolecularFormulaBase(MolecularFormulaCalc):
    """Base class for representing a molecular formula.

    Parameters
    ----------
    molecular_formula : dict, list, str
        The molecular formula.
    ion_charge : int
        The ion charge.
    ion_type : str, optional
        The ion type. Defaults to None.
    adduct_atom : str, optional
        The adduct atom. Defaults to None.
    mspeak_parent : _MSPeak, optional
        The parent mass spectrum peak object instance. Defaults to None.
    external_mz : float, optional
        The external m/z value. Defaults to None.

    Raises
    ------
    TypeError
        If the ion type is not 'DE_OR_PROTONATED', 'RADICAL' or  'ADDUCT'.

    Attributes
    ----------
    isotopologue_count_percentile : float
        The isotopologue count percentile.
    O_C : float
        The O/C ratio.
    H_C : float
        The H/C ratio.
    dbe : float
        The double bond equivalent.
    mz_nominal_calc : int
        The nominal m/z value.
    mz_error : float
        The m/z error.
    mz_calc : float
        The m/z value.
    protonated_mz : float
        The protonated or deprotonated m/z value.
    radical_mz : float
        The radical m/z value.
    neutral_mass : float
        The neutral mass.
    ion_type : str
        The ion type.
    ion_charge : int
        The ion charge.
    atoms : list
        The atoms in the molecular formula.
    confidence_score : float
        The confidence score of the molecular formula identification.
    isotopologue_similarity : float
        The isotopologue similarity score of the molecular formula identification.
    average_mz_error_score : float
        The average m/z error score of the molecular formula identification, including the isotopologues.
    mz_error_score : float
        The m/z error score of the molecular formula identification.
    kmd : float
        The Kendrick mass defect (KMD).
    kendrick_mass : float
        The Kendrick mass.
    knm : float
        The nominal Kendrick mass.
    string : str
        The molecular formula string.
    string_formated : str
        The molecular formula string formated with subscripts and superscripts.
    class_label : str
        The class label.
    class_dict : dict
        The class dictionary.

    Methods
    -------
    * change_kendrick_base(kendrick_dict_base).
        Change the Kendrick base.
    * isotopologues(min_abundance, current_mono_abundance, dynamic_range).
        Calculate the isotopologues.
    * atoms_qnt(atom).
        Get the atom quantity.
    * atoms_symbol(atom).
        Get the atom symbol without the mass number.
    * to_dict().
        Get the molecular formula as a dictionary.
    * to_list().
        Get the molecular formula as a list.
    """    

    def __init__(self, molecular_formula, ion_charge, ion_type=None, 
                adduct_atom=None, mspeak_parent=None, external_mz=None):
        # clear dictionary of atoms with 0 value
        if  type(molecular_formula) is dict:
                self._from_dict(molecular_formula, ion_type, adduct_atom)   
        
        elif type(molecular_formula) is list:
                self._from_list(molecular_formula, ion_type, adduct_atom)   
        
        elif type(molecular_formula) is str:
                self._from_str(molecular_formula, ion_type, adduct_atom)   

        self._ion_charge = ion_charge
        self._external_mz = external_mz
        self._confidence_score = None        
        self._isotopologue_similarity = None
        self._mz_error_score = None
        self._mass_error_average_score = None

        self.is_isotopologue = False
        
        # parent mass spectrum peak obj instance
        self._mspeak_parent = mspeak_parent

        self.expected_isotopologues = []
        self.mspeak_mf_isotopologues_indexes = []
        
        if self._mspeak_parent:
            kendrick_dict_base = self._mspeak_parent._ms_parent.mspeaks_settings.kendrick_base
        else:
            kendrick_dict_base = {'C':1, 'H':2}
        self._kmd, self._kendrick_mass, self._nominal_km = self._calc_kmd(
            kendrick_dict_base)  
        
    def __repr__(self):

        return "MolecularFormula({0},{1},ion type = {2}".format(self._d_molecular_formula, self.ion_charge, self.ion_type)
    
    def __str__(self):

        return "MolecularFormula {0}, ion_charge:{1}, ion type:{2}, m/z:{3} ".format(self.string, self.ion_charge, self.ion_type, self.mz_calc)
    
    def __len__(self):
        
        # crash if keys are not ordered
        return len(self._d_molecular_formula.keys())
        
    def __getitem__(self, atom):
        
            #atom = list(self._d_molecular_formula.keys())[position]
            if atom in self._d_molecular_formula.keys():
                return self._d_molecular_formula[atom]
            else:
                return 0
    def get(self, atom):
        """Get the atom quantity of a specific atom.
        
        Parameters
        ----------
        atom : str
            The atom symbol.
            
        Returns
        -------
        int
            The atom quantity.
        """
        #atom = list(self._d_molecular_formula.keys())[position]
        if atom in self._d_molecular_formula.keys():
            return self._d_molecular_formula[atom]
        else:
            return 0
                
    def _from_dict(self, molecular_formula, ion_type, adduct_atom):
        
        self._d_molecular_formula = {key:val for key, val in molecular_formula.items() if val != 0}
        
        if ion_type is not None:
            self._d_molecular_formula[Labels.ion_type] = ion_type
            
        if adduct_atom:
            if adduct_atom in self._d_molecular_formula:
                self._d_molecular_formula[adduct_atom] += 1 
            else: self._d_molecular_formula[adduct_atom] = 1 
        self.adduct_atom = adduct_atom

    def _from_list(self, molecular_formula_list, ion_type, adduct_atom):
        # list has to be in the format 
        #['C', 10, 'H', 21, '13C', 1, 'Cl', 1, etc]  
        self._d_molecular_formula = {}
        for each in range(0, len(molecular_formula_list),2):
            
            atoms_label =  molecular_formula_list[each]
            atoms_count = int(molecular_formula_list[each+1])
            
            if atoms_count > 0:
                self._d_molecular_formula[atoms_label] = int(atoms_count)
        
        self._d_molecular_formula[Labels.ion_type] = ion_type
        if adduct_atom:
            self.adduct_atom = adduct_atom
            if adduct_atom in self._d_molecular_formula:
                self._d_molecular_formula[adduct_atom] += 1 
            else: self._d_molecular_formula[adduct_atom] = 1 
        else:
            self.adduct_atom = None

    def _from_str(self, molecular_formula_str,  ion_type, adduct_atom):
        # string has to be in the format 
        #'C10 H21 13C1 Cl1 37Cl1 etc'
        molecular_formula_list = molecular_formula_str.split(' ')
        final_formula = []
        for i in molecular_formula_list:
            atoms_count = self.split(Atoms.atoms_order, i)
            final_formula.extend(atoms_count)
        print(final_formula)
        self._from_list(final_formula, ion_type, adduct_atom)

    def split(self, delimiters, string, maxsplit=0): #pragma: no cover
        """Splits the molecular formula string.
        
        Parameters
        ----------
        delimiters : list
            The list of delimiters.
        string : str
            The molecular formula string.
        maxsplit : int, optional
            The maximum number of splits. Defaults to 0.

        Returns
        -------
        list
            The molecular formula list.

        Notes
        -----
        Does not work when formula has atoms with same characters in a row that below to different atoms, i.e. C10H21NNa.
        """
        regexPattern = '|'.join(map(re.escape, delimiters)) #pragma: no cover
        isotopes = re.findall(regexPattern, string) #pragma: no cover
        counts = re.split(regexPattern, string, maxsplit)  #pragma: no cover
       
        return [isotopes[0], int(counts[1])]

    @property
    def isotopologue_count_percentile(self, ):
        if not len(self.expected_isotopologues) == 0:
            return (len(self.mspeak_mf_isotopologues_indexes)/len(self.expected_isotopologues))*100
        else: 
            return 100

    @property
    def O_C(self): 
            
            if 'O' in self._d_molecular_formula.keys():
                return self._d_molecular_formula.get("O")/self._d_molecular_formula.get("C")
            else:
                return 0    
    
    @property
    def H_C(self): return self._d_molecular_formula.get("H")/self._d_molecular_formula.get("C")
    
    @property
    def dbe(self): return self._calc_dbe()
    
    @property
    def mz_nominal_calc(self): return int(self._calc_mz())

    @property    
    def mz_error(self): return self._calc_assignment_mass_error()

    @property
    def mz_calc(self): return self._calc_mz()

    @property
    def protonated_mz(self): return self._protonated_mz(self.ion_charge)
    
    @property
    def radical_mz(self): return self._radical_mz(self.ion_charge)
    
    @property
    def neutral_mass(self): return self._neutral_mass()
    
    def adduct_mz(self, adduct_atom): 
        """Get m/z of an adducted ion version of the molecular formula.
        
        Parameters
        ----------
        adduct_atom : str
            The adduct atom.
            
        Returns
        -------
        float
            The m/z value of the adducted ion version of the molecular formula.
        """
        return self._adduct_mz(adduct_atom, self.ion_charge)

    @property
    def ion_type(self): 
        
        ion_type = self._d_molecular_formula.get(Labels.ion_type)
        if ion_type == Labels.protonated_de_ion:
            if self.ion_charge > 0: 
                return Labels.protonated
            else: 
                return Labels.de_protonated    
        else:
            return ion_type

    @ion_type.setter
    def ion_type(self, ion_type):
        if  ion_type in [Labels.protonated_de_ion, Labels.adduct_ion, Labels.radical_ion]:
            self._d_molecular_formula[Labels.ion_type] = ion_type
        else:
            raise TypeError("Ion type can only be: 'DE_OR_PROTONATED', 'RADICAL' or  'ADDUCT', not %s"%ion_type)   

    @property
    def ion_charge(self): return self._ion_charge
    
    @property
    def atoms(self): 
        """Get the atoms in the molecular formula."""
        # if there is an adduct_atom, them reduce it from the atoms list
        if self.adduct_atom is None:
            return [key for key in self._d_molecular_formula.keys() if key != Labels.ion_type]
        else:
            temp_dict = self._d_molecular_formula.copy()
            temp_dict[self.adduct_atom] -= 1
            return [key for key,val in temp_dict.items() if key != Labels.ion_type and val > 0]

    
    @property
    def confidence_score(self): 
        
        if not self._confidence_score:
            
            self._confidence_score = self._calc_confidence_score()
        
        return self._confidence_score

    @property
    def isotopologue_similarity(self): 
        
        if not self._isotopologue_similarity:
           
           self._isotopologue_similarity = self._calc_isotopologue_confidence()  
       
        return self._isotopologue_similarity  
    
    @property
    def average_mz_error_score(self): 
        
        # includes the isotopologues
        
        if not self._mass_error_average_score:
           
           self._mass_error_average_score = self._calc_average_mz_score()  
        
        return self._mass_error_average_score

    @property
    def mz_error_score(self): 
        if not self._mz_error_score:
           
           self._mz_error_score = self._calc_mz_confidence()  
        
        return self._mz_error_score
    
    @property
    def kmd(self): return self._kmd

    @property
    def kendrick_mass(self): return self._kendrick_mass

    @property
    def knm(self): return self._nominal_km

    def change_kendrick_base(self, kendrick_dict_base):
        """Change the Kendrick base.

        Parameters
        ----------
        kendrick_dict_base : dict
            The Kendrick base dictionary. Ex: {"C": 1, "H": 2}
        """ 
        self._kmd, self._kendrick_mass, self._nominal_km = self._calc_kmd(kendrick_dict_base)
                
    def isotopologues(self, min_abundance, current_mono_abundance, dynamic_range): 
        """Calculate the isotopologues for a given molecular formula.

        Parameters
        ----------
        min_abundance : float
            The minimum abundance.
        current_mono_abundance : float
            The current monoisotopic abundance.
        dynamic_range : float
            The dynamic range.

        Yields
        ------
        MolecularFormulaIsotopologue
            The molecular formula isotopologue.

        Notes
        -----
        This calculation ignores the hydrogen isotopes.
        """
        isotopologues = []
        for mf in self._cal_isotopologues(self._d_molecular_formula, min_abundance, current_mono_abundance, dynamic_range ):
            isotopologues.append(mf)
        
        # To account for differences in how the isotopologue outputs are sorted between IsoSpec versions. 
        sorted_isotopologues = sorted(isotopologues, key=lambda mf: mf[1], reverse=True)

        for mf in sorted_isotopologues:
            yield MolecularFormulaIsotopologue(
                *mf, 
                current_mono_abundance, 
                self.ion_charge, 
                ion_type=self.ion_type, 
                adduct_atom=self.adduct_atom
                )
    
    def atoms_qnt(self,atom): 
        """Get the atom quantity of a specific atom in the molecular formula."""
        if atom in self._d_molecular_formula:
            return self._d_molecular_formula.get(atom)
        else:
            raise Warning('Could not find %s in this Molecular Formula object'%str(atom))
    
    def atoms_symbol(self, atom): 
        """Get the atom symbol without the mass number."""
        return ''.join([i for i in atom if not i.isdigit()])

    @property       
    def string(self):
        """Returns the molecular formula as a string."""
        if self._d_molecular_formula:
            if self.adduct_atom is None:
                mol_form_dict = self._d_molecular_formula
            else:
                mol_form_dict = self._d_molecular_formula.copy()
                if self.adduct_atom not in mol_form_dict.keys():
                    raise Exception("Adduct atom not found in molecular formula dict")
                mol_form_dict[self.adduct_atom] -= 1
                mol_form_dict = {key:val for key, val in mol_form_dict.items() if val != 0}
            formula_srt = ''
            for atom in Atoms.atoms_order:
                if atom in mol_form_dict.keys():
                    formula_srt += atom + str(int(mol_form_dict.get(atom))) + ' '
            return formula_srt.strip()
        
        else:
            raise Exception("Molecular formula identification not performed yet")    
    
    @property
    def string_formated(self):
        
        SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
        
        if self._d_molecular_formula:
            formula_srt = ''
            for atom in Atoms.atoms_order:
                if atom in self.to_dict().keys():
                    formula_srt += atom.translate(SUP) + str(int(self.to_dict().get(atom))).translate(SUB)
            return formula_srt
        
        else:
            raise Exception("Molecular formula identification not performed yet")    

    def to_dict(self):
        """Returns the molecular formula as a dictionary.
        
        Returns
        -------
        dict
            The molecular formula as a dictionary.
        """
        return self._d_molecular_formula

    def to_list(self):
        """Returns the molecular formula as a list.
        
        Returns
        -------
        list
            The molecular formula as a list.
            
        Raises
        ------
        Exception
            If the molecular formula identification was not performed yet.
        """
        #TODO ensure self._d_molecular_formula is a orderedDict
        
        if self._d_molecular_formula:
            formula_list = []    
            
            for atom, atom_number in self._d_molecular_formula.items():
    
                if atom != Labels.ion_type:
                    
                    formula_list.append(atom)
                    formula_list.append(atom_number)
    
            return formula_list
        else:
            raise Exception("Molecular formula identification not performed yet")
    
    @property
    def class_label(self):
        
        if self._d_molecular_formula:
            
            formulalist = self.to_list()
            classstring = '' 
            
            for each in range(0, len(formulalist),2):
                
                if formulalist[each] != 'C' and formulalist[each] != 'H' and formulalist[each] != 'HC':
                     
                    classstring = classstring + str(formulalist[each]) + str(formulalist[each+1]) + ' '    
            
            if classstring == '': classstring = 'HC'
                
            classstring = classstring.strip()
            
            if self._d_molecular_formula.get(Labels.ion_type) == Labels.radical_ion:    
                
                return classstring + ' -R'
            
            #elif self._d_molecular_formula.get(Labels.ion_type) == Labels.adduct_ion:    
                
            #    return classstring + ' -A'

            else: return classstring
            
            #'dict, tuple or string'
        
        else:
            
            raise Exception("Molecular formula identification not performed yet")        
    
    @property
    def class_dict(self):
        
        if self._d_molecular_formula:
            
            class_dict = {}
            
            for atom, qnt in self._d_molecular_formula.items():
    
                if atom != Labels.ion_type and atom !='C' and atom !='H':
                    
                    class_dict[atom] = qnt
                    
            return class_dict
        
        raise Exception("Molecular formula identification not performed yet")           
    

class MolecularFormulaIsotopologue(MolecularFormulaBase):
    """Class for representing a molecular formula isotopologue.
    
    Parameters
    ----------
    _d_molecular_formula : dict
        The molecular formula as a dictionary.
    prob_ratio : float
        The probability ratio.
    mono_abundance : float
        The monoisotopic abundance.
    ion_charge : int
        The ion charge.
    mspeak_parent : object, optional
        The parent mass spectrum peak object instance. Defaults to None.
    ion_type : str, optional
        The ion type. Defaults to None.
    adduct_atom : str, optional
        The adduct atom. Defaults to None.
    
    Attributes
    ----------
    prob_ratio : float
        The probability ratio.
    abundance_calc : float
        The calculated abundance.
    area_error : float
        The area error.
    abundance_error : float
        The abundance error.
    is_isotopologue : bool
        The isotopologue flag. Defaults to True.
    mspeak_index_mono_isotopic : int
        The index of the monoisotopic peak in the mass spectrum peak list. Defaults to None.
    mono_isotopic_formula_index : int
        The index of the monoisotopic formula in the molecular formula list. Defaults to None.
    """
    def __init__(
            self, 
            _d_molecular_formula, 
            prob_ratio, 
            mono_abundance, 
            ion_charge, 
            mspeak_parent=None,
            ion_type = None,
            adduct_atom = None
            ):
        
        if ion_type is None:
            # check if ion type or adduct_atom is in the molecular formula dict
            if Labels.ion_type in _d_molecular_formula:
                ion_type = _d_molecular_formula.get(Labels.ion_type)
            else:
                ion_type = None
        else:
            ion_type = Labels.ion_type_translate.get(ion_type)
        
        if ion_type == Labels.adduct_ion:
            adduct_atom_int = None
            if adduct_atom in _d_molecular_formula.keys():
                adduct_atom_int = adduct_atom
            else:
                # Check to see if adduct_atom should actually be an isotope of the adduct atom
                for adduct_iso in Atoms.isotopes.get(adduct_atom)[1]:
                    if adduct_iso in _d_molecular_formula.keys():
                        adduct_atom_int = adduct_iso
            adduct_atom = adduct_atom_int
            if adduct_atom is None:
                raise Exception("adduct_atom is required for adduct ion")
            _d_molecular_formula[adduct_atom] -= 1
            _d_molecular_formula = {key:val for key, val in _d_molecular_formula.items() if val != 0}

        
        super().__init__(
            molecular_formula =_d_molecular_formula, 
            ion_charge = ion_charge, 
            ion_type=ion_type,
            adduct_atom=adduct_atom
            )
        #prob_ratio is relative to the monoisotopic peak p_isotopologue/p_mono_isotopic
        
        self.prob_ratio = prob_ratio
        
        self.abundance_calc = mono_abundance * prob_ratio

        self.is_isotopologue = True
        
        self.mspeak_index_mono_isotopic = None

        self.mono_isotopic_formula_index = None
        # parent mass spectrum peak obj instance
        self._mspeak_parent = mspeak_parent

    
    @property
    def area_error(self):
        return self._calc_area_error()

    @property
    def abundance_error(self):
        return self._calc_abundance_error()

class LCMSLibRefMolecularFormula(MolecularFormulaBase):
    """Class for representing a molecular formula associated with a molecule in a LCMS library reference.

    Parameters
    ----------
    molecular_formula : dict, list, str
        The molecular formula.
    ion_charge : int
        The ion charge.
    ion_type : str, optional
        The ion type. Defaults to None.
    adduct_atom : str, optional
        The adduct atom. Defaults to None.
    mspeak_parent : object, optional
        The parent mass spectrum peak object instance. Defaults to None.
    name : str, optional
        The name of the reference molecule. Defaults to None.
    kegg_id : str, optional
        The KEGG ID of the reference molecule. Defaults to None.
    cas : str, optional
        The CAS number of the reference molecule. Defaults to None.

    """
    
    def __init__(self, molecular_formula, ion_charge, ion_type=None, 
                    adduct_atom=None, mspeak_parent=None, name=None, kegg_id=None, cas=None) -> None:
        
        super().__init__(molecular_formula, ion_charge, ion_type=ion_type, 
                    adduct_atom=adduct_atom, mspeak_parent=mspeak_parent)

        self._name = name
        self._kegg_id = kegg_id
        self._cas = cas    
    
    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        if isinstance(name, str):
            self._name = name
        else:
            raise TypeError('name: {} should be type string')    

    @property
    def kegg_id(self):
        return self._kegg_id
    
    @kegg_id.setter
    def kegg_id(self, kegg_id):
        self._kegg_id = kegg_id
        #if isinstance(kegg_id, str):
        #    self._kegg_id = kegg_id
        #else:
        #    print(kegg_id)
        #    raise TypeError('name: {} should be type string') 

    @property
    def cas(self):
        return self._cas    
    
    @cas.setter
    def cas(self, cas):
        self._cas = cas
        #if isinstance(cas, str):
        #    self._cas = cas
        #else:
        #    raise TypeError('name: {} should be type string') 
    
class MolecularFormula(MolecularFormulaBase):
    """General class for representing a molecular formula.

    Parameters
    ----------
    molecular_formula : dict, list, str
        The molecular formula.
    ion_charge : int
        The ion charge.
    ion_type : str, optional
        The ion type. Defaults to None.
    adduct_atom : str, optional
        The adduct atom. Defaults to None.
    mspeak_parent : object, optional
        The parent mass spectrum peak object instance. Defaults to None.
    external_mz : float, optional
        The external m/z value. Defaults to False.
    """

    def __init__(self, molecular_formula, ion_charge, ion_type=None, 
                adduct_atom=None, mspeak_parent=None, external_mz=False):
        super().__init__(molecular_formula, ion_charge, ion_type=ion_type, 
                adduct_atom=adduct_atom, mspeak_parent=mspeak_parent, external_mz=external_mz)