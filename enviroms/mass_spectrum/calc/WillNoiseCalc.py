from enviroms.mass_spectrum.calc.NoiseCalc import NoiseThreshouldCalc

class WillNoiseThreshouldCalc(NoiseThreshouldCalc):
    
    '''Will here is the function that the MSClass call to do the MS calculation
        anyhow don't worry about the intergration part, just do your methods development here
        I will make sure it will get integrated 
        Ideally you should have the ouput to be the baseline noise and the std
        
        'all mass spec objects are available inside this class'
        #print(self.__dict__) #this will give you all the objects availables
    '''
    
    def xrun_noise_threshould_calc(self, auto):
            mz = self.exp_mz
            abun = self.abundance_profile 
            Y_cut = self.cut_mz_domain_noise(auto)
            
            if auto:

                yminima = self.get_abundance_minima_centroide(Y_cut)
                
                #returns baseline_noise and std
                return self.get_noise_average(yminima)

            else:
                
                #returns baseline_noise and std
                return self.get_noise_average(Y_cut)    