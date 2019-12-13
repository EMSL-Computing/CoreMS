import pprint, json
from collections import namedtuple, OrderedDict
from pathlib import Path
from pprint import pprint

path = str(Path(__file__).parent.absolute())

filepath = str(path)+'/AtomicWeightsAndIsotopicCompNIST2019.txt'

Atom = namedtuple('Atom', 'atomic_number symbol mass_number atomic_mass natural_abundance')
atoms_dict = dict()

with open(filepath, 'r') as f:
    
    read_lines = f.readlines()
    
    
    replace_dict = {"12C": "C", "2D": "D", "3T": "T", "4He": "He", "14N": "N",
                    "16O": "O", "19F": "F", "31P": "P", "32S": "S", "39K": "K",
                     "40Ca": "Ca", "51V": "V", "56Fe": "Fe", "1H": "H",
                    }
    
    for index in range(0,len(read_lines),8):
        
        atomic_number =  int(read_lines[index].split("=")[1].rstrip())
        
        symbol = str(read_lines[index+1].split("=")[1].rstrip()).strip()
        
        mass_number = str(read_lines[index+2].split("=")[1].rstrip()).strip()
        
        atomic_mass = float("".join(s for s in read_lines[index+3].split("=")[1].rstrip() if '(' not in s and ')' not in s and '#' not in s))
        
        if read_lines[index+4].split("=")[1].rstrip():
        
            natural_abundance = float("".join(s for s in read_lines[index+4].split("=")[1].rstrip() if '(' not in s and ')' not in s and '#' not in s))
        else: 
            natural_abundance = None
        
        chem_symbol = mass_number+symbol
        if chem_symbol in replace_dict.keys():
            chem_symbol = replace_dict.get(chem_symbol)

        #print(chem_symbol, atomic_number,symbol,mass_number,atomic_mass,natural_abundance)        

        atoms_dict[chem_symbol] = Atom(atomic_number= atomic_number,
                                        symbol = symbol,
                                        mass_number=mass_number, 
                                        atomic_mass=atomic_mass, 
                                        natural_abundance=natural_abundance)
    #print(atoms_dict)
    #pprint(atoms_dict, width=30)  
    #print(json.dumps(atoms_dict, indent=4))                                      
    #print(atoms_dict)
    for key, item in atoms_dict.items():
        print ('"'+key+'"', ":", item.natural_abundance, ",")
    #print(10)