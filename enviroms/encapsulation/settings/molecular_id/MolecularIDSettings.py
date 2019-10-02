class MoleculaLookupDictSettings:
    
    # C, H, N, O, S and P atoms are ALWAYS needed at usedAtoms
    # if you don't want to include one of thoses atoms set the max and min at 0
    # you can include any atom listed at Atoms class inside constants module
    # make sure to include the selected covalence at the used_atoms_valences when adding atoms 
    # to the usedAtoms dicts 
    # NOTE : Adducts atoms have zero covalence
    used_atom_valences = {'C': 4,
                            '13C': 4,
                            'H': 1,
                            'O': 2,
                            '18O': 2,
                            'N': 3,
                            'S': 2,
                            '34S': 2,
                            'P': 3,
                            'Cl': 0,
                            '37Cl': 0,
                            'Br': 0,
                            'Na': 1,
                            'F': 0,
                            }

    def __init__(self):
        
        self.usedAtoms = {'C': (1, 100),
                    'H': (4, 200),
                    'O': (0, 4),
                    'N': (0, 0),
                    'S': (0, 0),
                    'P': (0, 0),
                    'Cl': (0, 0),
                    }
        
        #min_mz changes automatically with mass spectrum
        self.min_mz = 100

        #max_mz changes automatically with mass spectrum
        self.max_mz = 1200

        self.min_dbe = 0

        self.max_dbe = 100

        #overwrites the dbe limits above to DBE = (C + heteroatoms) * 0.9
        self.use_pah_line_rule = False
        
        self.isRadical = True

        self.isProtonated = True

        self.ionization_type = "ESI"

        #ionCharge changes automatically with mass spectrum
        self.ionCharge = -1

        self.hc_filter = 0.3

        self.oc_filter = 1.2

class MoleculaSearchSettings:
    
    # look for close shell ions [M + Adduct]+ only considers metal set in the list adduct_atoms  
    isAdduct = False
    
    adduct_atoms_neg = ['Cl', 'Br', 'F']
    
    adduct_atoms_pos = ['Na', 'K']

    # depending on the polarity mode it looks for [M].+ , [M].-
    # needs to be enabled at the class MoleculaLookupTableSettings
    isRadical = False
    
    # depending on the polarity mode it looks for [M + H]+ , [M - H]+
    # needs to be enabled at the class MoleculaLookupTableSettings
    isProtonated = True

    # empirically set / needs optimization
    min_mz_error = -5 #ppm

    # empirically set / needs optimization    
    max_mz_error = 5 #ppm

    # empirically set / needs optimization
    min_abun_error = -30 # percentage 
    
    # empirically set / needs optimization
    max_abun_error = 70 # percentage 

    # empirically set / needs optimization
    mz_error_range = 1.5

    # 'distance', 'lowest', 'symmetrical','average' 'None'
    error_method = 'average'

    mz_error_average = 0

    min_dbe = 0 

    max_dbe = 50 

