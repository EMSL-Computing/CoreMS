__author__ = "Yuri E. Corilo"
__date__ = "Jan 31, 2020"

from collections.abc import Mapping

from corems.encapsulation.constant import Labels

class HeteroatomsClassification(Mapping):
    
    '''Group mass spectrum data by heteroatom classes (Nn, Oo, Ss, NnOo, NnSs, etc..)
        
       class obj behaves as a dictionary of classes and return a list of ms_peak obj

    '''
    
    def __init__(self, mass_spectrum, choose_molecular_formula=True):

        self._ms_grouped_class = dict()
        
        #mapping for ms peaks without any molecular formula associated

        self._ms_grouped_class[Labels.unassigned] = list()

        self.total_peaks = 0

        self.sum_abundance = 0

        check_assign = False

        for ms_peak in mass_spectrum:
            
            self.total_peaks += 1

            self.sum_abundance += ms_peak.abundance

            if not ms_peak.is_assigned:

                self._ms_grouped_class.get(Labels.unassigned).append(ms_peak)
                
            else:    
                
                check_assign = True    

                if choose_molecular_formula:
                    
                    classes =  ms_peak.best_molecular_formula_candidate.class_label
                
                else: 

                    classes = [mf.class_label for mf in ms_peak]

                for classe in classes:
                    
                    if classe in self._ms_grouped_class.keys():

                        self._ms_grouped_class.get(classe).append(ms_peak)
                    
                    else:     

                        self._ms_grouped_class[classe] = [ms_peak]

        if not check_assign:

            raise Exception("No molecular formula associated with any mspeak objects")

    def __len__(self):
        
        return len(self._ms_grouped_class)
        
    def __getitem__(self, classe):
        
        return self._ms_grouped_class.get(classe)

    def __iter__(self):

         return iter(self._ms_grouped_class) 

    def carbon_number(self, classe):

        return [mf.get('C') for mspeak in self[classe] for mf in mspeak ]
    
    def dbe(self, classe):

        return [mf.dbe for mspeak in self[classe] for mf in mspeak]
    
    def atoms_ratio(self, classe, numerator, denominator):

        #return [mspeak.mz_exp for classe in classes for mspeak in self[classe] if classe != Labels.unassigned]

        return [mf.get(numerator)/mf.get(denominator) for mspeak in self[classe] for mf in mspeak ]
    
    def mz_exp(self, classe):

        return [mspeak.mz_exp for mspeak in self[classe]]

    def abundance(self, classe):

        return [mspeak.abundance for mspeak in self[classe]]

    def peaks_count_percentile(self, classe):

        return (len(self[classe])/self.total_peaks)*100

    def abundance_count_percentile(self, classe):

        return (sum([mspeak.abundance for mspeak in self[classe]])/self.sum_abundance)*100

    def mz_exp_assigned(self):

        classes = self.keys()

        return [mspeak.mz_exp for classe in classes for mspeak in self[classe] if classe != Labels.unassigned]
    
    def abundance_assigned(self):

        classes = self.keys()
            
        return [mspeak.abundance for classe in classes for mspeak in self[classe] if classe != Labels.unassigned]

    def carbon_number_all(self):

        classes = self.keys()
            
        return [self.carbon_number(classe) for classe in classes if classe != Labels.unassigned]

    def dbe_all(self):

        classes = self.keys()
            
        return [self.dbe(classe) for classe in classes if classe != Labels.unassigned]

    def atoms_ratio_all(self, numerator, denominator):

        classes = self.keys()
            
        return [self.atoms_ratio(classe, numerator, denominator) for classe in classes if classe != Labels.unassigned]

    def to_dataframe(self,):
        
        from pandas import DataFrame

    def plot_ms_assigned_unassigned(self, assigned_color= 'g', unassigned_color = 'r'):
        
        from matplotlib import pyplot as plt
        
        mz_assigned = self.mz_exp_assigned()
        abundance_assigned = self.abundance_assigned()
    
        mz_not_assigned = self.mz_exp(Labels.unassigned)
        abundance_not_assigned = self.abundance(Labels.unassigned)
        
        ax = plt.gca()

        for plot_obj in ax.stem(mz_assigned,abundance_assigned, linefmt='-',  markerfmt=" ", use_line_collection =True):
        
            plt.setp(plot_obj, 'color', 'g', 'linewidth', 2)
        

        markerline, stemlines, baseline  = ax.stem(mz_not_assigned,abundance_not_assigned, linefmt='-', markerfmt=" ",  use_line_collection =True)
        
        
        plt.setp(markerline, 'color', 'r', 'linewidth', 2)
        plt.setp(stemlines, 'color', 'r', 'linewidth', 2)
        plt.setp(baseline, 'color', 'r', 'linewidth', 2)

        ax.set_xlabel("$\t{m/z}$", fontsize=12)
        ax.set_ylabel('Abundance', fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=12)

        ax.axes.spines['top'].set_visible(False)
        ax.axes.spines['right'].set_visible(False)

        ax.get_yaxis().set_visible(False)
        ax.spines['left'].set_visible(False)
        
        return ax    

