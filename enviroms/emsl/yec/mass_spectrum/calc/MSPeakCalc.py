'''
Created on Jun 14, 2019

@author: eber373
'''
from enviroms.emsl.yec.encapsulation.constant.Constants import Atoms


class MassSpecPeakCalculation(object):
    
    '''
    classdocs
    '''
    def _calc_kdm(self, dict_base):
        '''dict_base = {"C": 1, "H": 2}
        '''
        mass = 0
        for atom in dict_base.keys():
            mass = mass + Atoms.atomic_masses.get(atom) * dict_base.get(atom)
        
        kendrick_mass = (int(mass)/mass)*self.mz_exp
        
        nominal_km =int(kendrick_mass)
        #for ICR
        nominal_mass = int(self.mz_exp) 
        
        kmd = (nominal_mass - kendrick_mass) * 1000
        
        #kmd = (nominal_km - km) * 1
        kdm  = round(kmd,0)
        
        return kdm, kendrick_mass, nominal_km
    
        
    def _calc_confidence_score(self): raise Exception("Not implemented yet") 
    