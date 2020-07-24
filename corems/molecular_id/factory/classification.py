__author__ = "Yuri E. Corilo"
__date__ = "Jan 31, 2020"

from collections.abc import Mapping

from matplotlib import pyplot as plt
from numpy import linspace

from corems.encapsulation.constant import Labels
from corems.encapsulation.constant import Atoms

flatten_list = lambda l: [item for sublist in l for item in sublist]

class HeteroatomsClassification(Mapping):
    
    '''Group mass spectrum data by heteroatom classes (Nn, Oo, Ss, NnOo, NnSs, etc..)
        
       class obj behaves as a dictionary of classes and return a list of ms_peak obj

    '''
    
    def __init__(self, mass_spectrum, choose_molecular_formula=True):

        def sort_atoms_method( atom):
            
            return [Atoms.atoms_order.index(atom)]

        self._ms_grouped_class = dict()
        
        self.choose_mf = choose_molecular_formula
        
        #mapping for ms peaks without any molecular formula associated
        self._ms_grouped_class[Labels.unassigned] = list()

        self.total_peaks = 0

        self.sum_abundance = 0

        self.min_max_mz = (mass_spectrum.min_mz_exp, mass_spectrum.max_mz_exp)

        self.min_max_abundance = (mass_spectrum.min_abundance, mass_spectrum.max_abundance)

        self.min_ppm_error = mass_spectrum.molecular_search_settings.min_ppm_error 

        self.max_ppm_error = mass_spectrum.molecular_search_settings.max_ppm_error

        check_assign = False

        all_used_atoms = set()

        for ms_peak in mass_spectrum:
            
            self.total_peaks += 1

            self.sum_abundance += ms_peak.abundance

            if not ms_peak.is_assigned:

                self._ms_grouped_class.get(Labels.unassigned).append(ms_peak)
                
            else:    
                
                check_assign = True    

                if choose_molecular_formula:
                    
                    mf = ms_peak.best_molecular_formula_candidate
                    
                    classes =  [mf.class_label]
                    
                    for atom in mf.atoms:
                        
                        all_used_atoms.add(atom)

                else: 

                    classes = []
                    
                    for mf in ms_peak:
                        
                        classes.append(mf.class_label)
                        
                        for atom in mf.atoms:
                             
                             all_used_atoms.add(atom)

                for classe in classes:
                    
                    if classe in self._ms_grouped_class.keys():

                        self._ms_grouped_class.get(classe).append(ms_peak)
                    
                    else:     

                        self._ms_grouped_class[classe] = [ms_peak]

        self.all_identified_atoms = sorted(all_used_atoms, key=sort_atoms_method)

        if not check_assign:

            raise Exception("No molecular formula associated with any mspeak objects")
    

    def __len__(self):
        
        return len(self._ms_grouped_class)
        
    def __getitem__(self, classe):
        
        return self._ms_grouped_class.get(classe)

    def __iter__(self):

         return iter(self._ms_grouped_class) 

    def get_classes(self, threshold_perc=1, isotopologue=True):
        
        classes = list()
        for classe in self.keys():
            if classe != Labels.unassigned:
                if self.abundance_count_percentile(classe) > threshold_perc:
                    
                    if classe != Labels.unassigned:
                        # access first molecular formula inside the first ms peak and check isotopologue
                        if not isotopologue and self.get(classe)[0][0].is_isotopologue: continue
                    
                    classes.append(classe)
        #TODO sort classes chemically here too
        return classes
        
    def molecular_formula_string(self, classe,):

        if self.choose_mf:
            return [mspeak.best_molecular_formula_candidate for mspeak in self[classe]]
        else:
            return [mf for mspeak in self[classe] for mf in mspeak if mf.class_label == classe]

    def molecular_formula(self, classe,):

        if self.choose_mf:
            return [mspeak.best_molecular_formula_candidate for mspeak in self[classe]]
        else:
            return [mf for mspeak in self[classe] for mf in mspeak if mf.class_label == classe]
            
    def carbon_number(self, classe):

        if self.choose_mf:
            return [mspeak.best_molecular_formula_candidate.get("C") for mspeak in self[classe]]
        else:
            return [mf.get('C') for mspeak in self[classe] for mf in mspeak if mf.class_label == classe]

    def atom_count(self, atom, classe):
        if self.choose_mf:
            return [mspeak.best_molecular_formula_candidate.get(atom) for mspeak in self[classe]]
        else:    
            return [mf.get(atom) for mspeak in self[classe] for mf in mspeak if mf.class_label == classe]

    def dbe(self, classe):
        if self.choose_mf:
            return [mspeak.best_molecular_formula_candidate.dbe for mspeak in self[classe]]
        else:    
            return [mf.dbe for mspeak in self[classe] for mf in mspeak if mf.class_label == classe]
    
    def atoms_ratio(self, classe, numerator, denominator):

        return [mf.get(numerator)/mf.get(denominator) for mf in self.molecular_formula(classe)]
       
    def mz_exp(self, classe):
        
        if self.choose_mf or classe == Labels.unassigned:
            
            return [mspeak.mz_exp for mspeak in self[classe]]
        
        else:
            
            return [mspeak.mz_exp for mspeak in self[classe] for mf in mspeak if mf.class_label == classe]
    
    def abundance(self, classe):

        if self.choose_mf or classe == Labels.unassigned:
            
            return [mspeak.abundance for mspeak in self[classe]]
        
        else:
            
            return [mspeak.abundance for mspeak in self[classe] for mf in mspeak if mf.class_label == classe]

    def mz_error(self, classe):

        if classe != Labels.unassigned:
            
            if self.choose_mf:
                
                return [mspeak.best_molecular_formula_candidate.mz_error for mspeak in self[classe]]
            
            else:
                
                return [mf.mz_error for mspeak in self[classe] for mf in mspeak if mf.class_label == classe]
    
    def mz_calc(self, classe):
        
        if self.choose_mf:
            
            return [mspeak.best_molecular_formula_candidate.mz_calc for mspeak in self[classe]]
        
        else:
            
            return [mf.mz_calc for mspeak in self[classe] for mf in mspeak if mf.class_label == classe] 

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

    def mz_exp_all(self):
        
        classes = self.keys()
        
        return flatten_list([self.mz_exp(classe) for classe in classes if classe != Labels.unassigned])
    
    def mz_error_all(self):
        
        classes = self.keys()
        
        return flatten_list([self.mz_error(classe) for classe in classes if classe != Labels.unassigned])
        
    def carbon_number_all(self):

        classes = self.keys()
            
        return flatten_list([self.carbon_number(classe) for classe in classes if classe != Labels.unassigned])

    def dbe_all(self):

        classes = self.keys()
            
        return flatten_list([self.dbe(classe) for classe in classes if classe != Labels.unassigned])

    def atoms_ratio_all(self, numerator, denominator):

        classes = self.keys()
            
        return flatten_list([self.atoms_ratio(classe, numerator, denominator) for classe in classes if classe != Labels.unassigned])

    def to_dataframe(self, incluse_isotopologue=False, abundance_perc_threshold=5, include_unassigned=False):
        
        from pandas import DataFrame
        
        columns_labels = ['mz', 'calibrated_mz', 'calculated_m_z', 'abundance',
                                'resolving_power', 'sn', 'ion_charge', 'mass_error',
                                'DBE', 'class', 'HC', 'OC', 'ion_type','is_isotopologue',
                                'class_abundance', 'class_count']

        dict_data_list = []

        for classe, list_mspeaks in self.items():

            percent_abundance = self.abundance_count_percentile(classe)
            
            #ignores low abundant classes
            if abundance_perc_threshold < abundance_perc_threshold: continue
                
            peaks_count_percentile = self.peaks_count_percentile(classe)

            for ms_peak in list_mspeaks:
                 
                if ms_peak.is_assigned:
                    
                    for m_formula in ms_peak:
                        
                        #ignores isotopologues
                        if not incluse_isotopologue and m_formula.is_isotopologue: continue
                        
                        formula_dict = m_formula.to_dict()

                        dict_result = {'mz':  ms_peak._mz_exp,
                                'calibrated_mz': ms_peak.mz_exp,
                                'calculated_mz': m_formula.mz_calc,
                                'abundance': ms_peak.abundance,
                                'resolving_power': ms_peak.resolving_power,
                                'sn':  ms_peak.signal_to_noise,
                                'ion_charge': ms_peak.ion_charge,
                                'mass_error': m_formula.mz_error,
                                'DBE':  m_formula.dbe,
                                'class': classe,
                                'HC':  m_formula.H_C,
                                'OC':  m_formula.O_C,
                                'ion_type': str(m_formula.ion_type.lower().encode('utf-8')),
                                'is_isotopologue': int(m_formula.is_isotopologue),
                                'class_abundance': percent_abundance,
                                'class_count': peaks_count_percentile
                                }
                        
                        for atom in formula_dict.keys():
                        
                           dict_result[atom] = formula_dict.get(atom)

                    dict_data_list.append(dict_result)

                else:

                    if not include_unassigned: continue

                    dict_result = {'mz':  ms_peak._mz_exp,
                                'calibrated_mz': ms_peak.mz_exp,
                                'abundance': ms_peak.abundance,
                                'resolving_power': ms_peak.resolving_power,
                                'sn':  ms_peak.signal_to_noise,
                                'ion_charge': ms_peak.ion_charge,
                                'class': classe,
                                'class_abundance': percent_abundance,
                                'class_count': percent_abundance
                                }
                
                    dict_data_list.append(dict_result)                

        columns = columns_labels + self.all_identified_atoms

        return DataFrame(dict_data_list, columns=columns)
     
    def plot_ms_assigned_unassigned(self, assigned_color= 'b', unassigned_color = 'r'):
        
        mz_assigned = self.mz_exp_assigned()
        abundance_assigned = self.abundance_assigned()
    
        mz_not_assigned = self.mz_exp(Labels.unassigned)
        abundance_not_assigned = self.abundance(Labels.unassigned)
        
        ax = plt.gca()

        for plot_obj in ax.stem(mz_assigned,abundance_assigned, linefmt='-',  markerfmt=" ", use_line_collection =True, label="Assigned"):
        
            plt.setp(plot_obj, 'color', assigned_color, 'linewidth', 2)
        
        for plot_obj in ax.stem(mz_not_assigned, abundance_not_assigned, linefmt='-', markerfmt=" ",  use_line_collection =True,label="Unassigned"):
        
            plt.setp(plot_obj, 'color', unassigned_color, 'linewidth', 2)
        
        ax.set_xlabel("$\t{m/z}$", fontsize=12)
        ax.set_ylabel('Abundance', fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=12)

        ax.axes.spines['top'].set_visible(False)
        ax.axes.spines['right'].set_visible(False)

        ax.get_yaxis().set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.legend()
        return ax    

    def plot_mz_error(self, color= 'g'):
        
        ax = plt.gca()

        mz_assigned = self.mz_exp_all()
        mz_error= self.mz_error_all()
        
        ax.scatter( mz_assigned, mz_error, c=color)
        
        ax.set_xlabel("$\t{m/z}$", fontsize=12)
        ax.set_ylabel('Error (ppm)', fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=12)

        ax.axes.spines['top'].set_visible(True)
        ax.axes.spines['right'].set_visible(True)

        ax.get_yaxis().set_visible(True)
        ax.spines['left'].set_visible(True)

        ax.set_xlim(self.min_max_mz)
        ax.set_ylim(self.min_ppm_error , self.max_ppm_error)
    
        return ax

    def plot_mz_error_class(self, classe, color= 'g'):
        
        if classe != Labels.unassigned:
            ax = plt.gca()
            
            abun_perc = self.abundance_count_percentile(classe)
            mz_assigned = self.mz_exp(classe)
            mz_error= self.mz_error(classe)
            
            ax.scatter( mz_assigned, mz_error, c=color)
            
            title = "%s, %.2f %%" % (classe, abun_perc)
            ax.set_title(title)
            ax.set_xlabel("$\t{m/z}$", fontsize=12)
            ax.set_ylabel('Error (ppm)', fontsize=12)
            ax.tick_params(axis='both', which='major', labelsize=12)

            ax.axes.spines['top'].set_visible(True)
            ax.axes.spines['right'].set_visible(True)

            ax.get_yaxis().set_visible(True)
            ax.spines['left'].set_visible(True)

            ax.set_xlim(self.min_max_mz)
            ax.set_ylim(self.min_ppm_error , self.max_ppm_error)
        
            return ax   
            
    def plot_ms_class(self, classe, color= 'g'):
        
        if classe != Labels.unassigned:
            ax = plt.gca()
            
            abun_perc = self.abundance_count_percentile(classe)
            mz_assigned = self.mz_exp(classe)
            abundance_assigned= self.abundance(classe)

            for plot_obj in ax.stem( mz_assigned, abundance_assigned, linefmt='-',  markerfmt=" ", use_line_collection =True):
            
                plt.setp(plot_obj, 'color', color, 'linewidth', 2)
            
            title = "%s, %.2f %%" % (classe, abun_perc)
            ax.set_title(title)
            ax.set_xlabel("$\t{m/z}$", fontsize=12)
            ax.set_ylabel('Abundance', fontsize=12)
            ax.tick_params(axis='both', which='major', labelsize=12)

            ax.axes.spines['top'].set_visible(False)
            ax.axes.spines['right'].set_visible(False)

            ax.get_yaxis().set_visible(False)
            ax.spines['left'].set_visible(False)

            ax.set_xlim(self.min_max_mz)
            ax.set_ylim(self.min_max_abundance)
        
            return ax

    def plot_van_krevelen(self, classe, max_hc=2.5, max_oc=2, ticks_number=5, color="jet"):
        
        if classe != Labels.unassigned:

            # get data 
            abun_perc = self.abundance_count_percentile(classe)
            hc = self.atoms_ratio(classe, "H", "C") 
            oc = self.atoms_ratio(classe, "O", "C") 
            abundance = self.abundance(classe)
            
            #plot data
            ax = plt.gca()

            ax.scatter(oc, hc, c=abundance, alpha=0.5, cmap=color)

            #ax.scatter(carbon_number, dbe, c=color, alpha=0.5)
            
            title = "%s, %.2f %%" % (classe, abun_perc)
            ax.set_title(title)
            ax.set_xlabel("O/C", fontsize=16)
            ax.set_ylabel('H/C', fontsize=16)
            ax.tick_params(axis='both', which='major', labelsize=18)
            ax.set_xticks(linspace(0, max_oc, ticks_number, endpoint=True))
            ax.set_yticks(linspace(0, max_hc, ticks_number, endpoint=True))

            # returns matplot axes obj and the class percentile of the relative abundance 
            
            return ax, abun_perc 
            
    def plot_dbe_vs_carbon_number(self, classe, max_c=50, max_dbe=40, dbe_incr=5, c_incr=10, color="jet"):
        
        if classe != Labels.unassigned:

            # get data 
            abun_perc = self.abundance_count_percentile(classe)
            carbon_number = self.carbon_number(classe)
            dbe = self.dbe(classe)
            abundance = self.abundance(classe)
            
            #plot data
            ax = plt.gca()

            ax.scatter(carbon_number, dbe, c=abundance, alpha=0.5, cmap=color)

            #ax.scatter(carbon_number, dbe, c=color, alpha=0.5)
            
            title = "%s, %.2f %%" % (classe, abun_perc)
            ax.set_title(title)
            ax.set_xlabel("Carbon number", fontsize=16)
            ax.set_ylabel('DBE', fontsize=16)
            ax.tick_params(axis='both', which='major', labelsize=18)
            ax.set_xticks(range(0, max_c, c_incr))
            ax.set_yticks(range(0, max_dbe, dbe_incr))

            # returns matplot axes obj and the class percentile of the relative abundance 
            
            return ax, abun_perc 
