'''
Created on Jun 14, 2019

@author: eber373
'''
from emsl.yec.structure.static.Constants import Constants

class MassSpecPeakCalculation(object):
    
    '''
    classdocs
    '''
    def _calc_kdm(self, dict_base):
        '''dict_base = {"C": 1, "H": 2}
        '''
        mass = 0
        for atom in dict_base.keys():
            mass = mass + Constants.atomic_masses.get(atom) * dict_base.get(atom)
        
        kendrick_mass = (int(mass)/mass)*self.exp_mz
        
        nominal_km =int(kendrick_mass)
        #for ICR
        nominal_mass = int(self.exp_mz) 
        
        kmd = (nominal_mass - kendrick_mass) * 1000
        
        #kmd = (nominal_km - km) * 1
        kdm  = round(kmd,0)
        
        return kdm, kendrick_mass, nominal_km
    
        
    def _calc_confidence_score(self): raise Exception("Not implemented yet") 
    