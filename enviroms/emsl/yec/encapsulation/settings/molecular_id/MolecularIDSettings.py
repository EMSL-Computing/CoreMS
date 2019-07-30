

class MoleculaLookupTableSettings:
    
    # C, H, N, O, S and P atoms are ALWAYS needed at usedAtoms
    # if you don't want to include one of thoses atoms set the max and min at 0
    # you can include any atom listed at Atoms class inside Constants module
    # make sure to include the selected valence at the used_atoms_valences when adding atoms 
    # to the usedAtoms dicts 
    # NOTE : Adducts atoms have zero valence

    usedAtoms = {'C': (1, 100),
                 'H': (4, 200),
                 'O': (0, 30),
                 'N': (0, 1),
                 'S': (0, 1),
                 'P': (0, 0),
                 }
    
    used_atom_valences = {'C': 4,
                        'H': 1,
                        'O': 2,
                        'N': 3,
                        'S': 2,
                        'P': 3,
                        }

    #min_mz changes automatically with mass spectrum
    min_mz = 200

    #max_mz changes automatically with mass spectrum
    max_mz = 900

    min_dbe = 0

    max_dbe = 50

    #overwrites the dbe limits above to DBE = (C + heteroatoms) * 0.9
    use_pah_line_rule = True
    
    isRadical = True

    isProtonated = True

    ionization_type = "ESI"

    #ionCharge changes automatically with mass spectrum
    ionCharge = -1

    hc_filter = 0.3

    oc_filter = 1.2

class MoleculaSearchSettings:
    
    #needs to be enabled at the class MoleculaLookupTableSettings
    isRadical = False
    
    #needs to be enabled at the class MoleculaLookupTableSettings
    isProtonated = True

    #empirically set / needs optimization
    min_mz_error = -5 #ppm

    #empirically set / needs optimization    
    max_mz_error = 5 #ppm

    #empirically set / needs optimization
    min_abun_error = -30 # percentage 
    
    #empirically set / needs optimization
    max_abun_error = 70 # percentage 
