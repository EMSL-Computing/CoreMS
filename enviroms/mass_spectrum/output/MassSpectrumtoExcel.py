__author__ = "Yuri E. Corilo"
__date__ = "Jun 12, 2019"

        
class MassSpectoExcel(object):
    '''
    classdocs
    '''
    def __init__(self, out_file_path, mass_spectrum):
        
        '''
        Constructor
        '''
        #change this dict VALUES to match your labels, THE ORDER WON'T MATTER
        #self.header =  

        self.file_location = out_file_path
        
        self.mass_spectrum  = mass_spectrum

    def run(self):
        
        print(self.get_hetero_classes_by_chem_group())

    
    def get_hetero_classes_by_chem_group(self):
        
        import re
        
        classes = set()
        for mspeak in self.mass_spectrum:
            if mspeak:
                for molecular_formula in mspeak:
                    if not molecular_formula.is_isotopologue:
                       classes.add(molecular_formula.class_label)

        classes = {molecular_formula.class_label for molecular_formula in mspeak if not molecular_formula.is_isotopologue for mspeak in self.mass_spectrum if mspeak }
        
        classes_tuple_lista = list()
        
        for classe_str in classes:
            
            classe_list = classe_str.split()
            new_classe_list = list()
            new_classe_dict = {}
            
            for atoms in classe_list:
                
                lista =  re.findall(r'\d+|\D+', atoms)
                print(lista)
                if len(lista1) == 3:
                    lista1 = [lista1[0]+ lista1[1], lista1[2]]
                #if lista1[0] != "13C":
                new_classe_list = new_classe_list + lista1
            
            
            if new_classe_list[-1] == "-R" or new_classe_list[-1] == "-H": 
                
                new_classe_list = new_classe_list[0:-1]
                
            if len(new_classe_list) == 1:
                
                new_classe_list = new_classe_list + [1]
                
            new_classe_list_int = []
            
            for cada in range(0, len(new_classe_list), 2):
                    
                    new_classe_dict[new_classe_list[cada]] = int(new_classe_list[cada+1])
                    new_classe_list_int.append(new_classe_list[cada])
                    new_classe_list_int.append(int(new_classe_list[cada+1]))
                    
            classes_tuple_lista.append((new_classe_list_int, new_classe_dict, classe_str))       