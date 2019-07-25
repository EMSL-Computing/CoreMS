

class MolecularSpaceTableSetting:
    used_atom_valences = {'C': 4,
                        'H': 1,
                        'O': 2,
                        'N': 3,
                        'S': 2,
                        'P': 3,
                        }

    # C, H, N, O, S and P atoms are ALWAYS needed here
    # if you don't want to include one of thoses atoms set the max and min at 0
    usedAtoms = {'C': (1, 100),
                 'H': (4, 200),
                 'O': (0, 10),
                 'N': (0, 2),
                 'S': (0, 2),
                 'P': (0, 0),
                 }
    min_dbe = 0

    max_dbe = 50

    use_pah_line_rule = True

    isRadical = True

    isProtonated = True

    ionization_type = "ESI"

    ionCharge = -1
