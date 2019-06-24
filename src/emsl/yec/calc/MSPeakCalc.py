'''
Created on Jun 14, 2019

@author: eber373
'''
from emsl.yec.structure.static.Constants import Constants

class MassSpecPeakCalculation(object):
    
    '''
    classdocs
    '''
    
    @property.setter
    def kdm(self):
        
        Hmass = Constants.atomic_masses.get('H')
        
        self._kendrick_mass = (14.00000/(12.00000+ (2*Hmass)))*self.exp_mass
        
        self._nominal_km =int(self.kendrick_mass)
        #for ICR
        nominal_mass = int(self.exp_mass) 
        
        kmd = (nominal_mass - self.kendrick_mass) * 1000
        
        #kmd = (nominal_km - km) * 1
        
        self._kdm  = round(kmd,0)
    
    @property.setter
    def __calc_dbe(self, formula_dict):
            
            individual_dbe = 0
            
            for atom in formula_dict.keys():
                if atom != "IonType":
                    n_atom = int(formula_dict.get(atom))
                    valencia = Constants.atoms_valence.get(atom)
                    if valencia > 0 and valencia is not None:
                        #print atom, valencia, n_atom, individual_dbe
                        individual_dbe = individual_dbe + (n_atom * (valencia - 2))
                    else:
                        continue
            
            self._dbe = 1 + (0.5 * individual_dbe)
            
    @property.setter        
    def __calc_theoretical_mz(self, formula_dict):
        
        if self.ion_charge:
            
            mass = 0
            
            for each_atom in formula_dict.keys() :
                
                if each_atom != "IonType" and each_atom != 'HC':
                    
                    mass = mass + Constants.atomic_masses[each_atom]  *  formula_dict.get(each_atom)
                    
                    #print mass
                    
            self._theoretical_mz = mass + ((-self.ion_charge) * Constants.electron_mass)
                    
                    
        else:
            
            raise Exception("Please set ion charge first")
         
        self.__set_mass_error('ppm')
           
    @property.setter
    def __set_mass_error(self, method):
        '''methodo should be ppm or ppb'''
        
        Hum_Milhao = 1000000
        Hum_Bilhao = 1000000000        
        
        if method == 'ppm':
            mult_factor = Hum_Milhao
        
        elif method == 'ppb':
            mult_factor = Hum_Bilhao
        
        else:
            raise Exception("method needs to be ppm or ppb, you have entered %s" % method)
              
        if self.theoretical_mass:
            
            self._assigment_mass_error =  ((self.theoretical_mass - self.exp_mass)/self.theoretical_mass)*mult_factor
        
        else:
            
            raise Exception("Please set theoretical_mass first")
    