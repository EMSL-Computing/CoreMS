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
    def _calc_kdm(self):
        
        Hmass = Constants.atomic_masses.get('H')
        
        kendrick_mass = (14.00000/(12.00000+ (2*Hmass)))*self.exp_mass
        
        nominal_km =int(self.kendrick_mass)
        #for ICR
        nominal_mass = int(self.exp_mass) 
        
        kmd = (nominal_mass - self.kendrick_mass) * 1000
        
        #kmd = (nominal_km - km) * 1
        kdm  = round(kmd,0)
        
        return kdm, kendrick_mass, nominal_km
            
    @property.setter
    
        
    @property.setter
    def _calc_confidence_score(self): raise Exception("Not implemented yet") 
    
    @property.setter
    def _calc_resolving_power(self): raise Exception("Not implemented yet") 
        
        