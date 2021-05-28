import json

from pymongo import MongoClient
from corems.encapsulation.factory.parameters import MSParameters

class MolForm_Mongo:
    
    '''this class is outdated
    'exists only for reference
    '''
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        # make sure the dbconnection gets closed
        self.client.close()

    def __enter__(self):
        
        self.client = MongoClient("mongodb://corems-client:esmlpnnl2019@localhost:27017/corems")
        
        db = self.client.corems
        self.molform_collection = db.molform
        return self
       
    def add_all(self, dic_molform_list):
        
        try:
            self.molform_collection.insert_many(dic_molform_list, ordered=False)
        except Exception as e:
            #duplicate_error
            try:
                if e.details["writeErrors"][0]["code"] == 11000:
                    print('ohh nooo')
            except:
                raise  Exception(e)        
    
    def check_entry(self,classe, ion_type):
       
        formulas = self.molform_collection.find({'classe': classe, 'ion_type': ion_type}).limit(1)
        
        if list(formulas):
            return True
        else:
            return False    
        
    def add_entry(self, dic_molform): 

        try:
            self.molform_collection.insert(dic_molform)
        except Exception as e:
            #duplicate_error
            try:
                if e.details["writeErrors"][0]["code"] == 11000:
                    pass
            except:
                raise  Exception(e)   
  
    def get_all(self):

        return self.molform_collection.find({})

    def get_dict_entries(self, classe, ion_type, nominal_mzs):

        dict_res = {}

        formulas = self.molform_collection.find( {'classe': {"$in": classe}, 
                                                'ion_type': ion_type,
                                                'nominal_mz':{"$in": nominal_mzs},  
                                                'O_C' : { '$lt': MSParameters.molecular_search.oc_filter }, 
                                                'H_C' : { '$gt': MSParameters.molecular_search.min_hc_filter},
                                                'DBE' : { '$gt': MSParameters.molecular_search.min_dbe},
                                                'DBE' : { '$lt': MSParameters.molecular_search.max_dbe},
                                                })
        for formula in formulas:
            
            if formula["classe"] in dict_res.keys():
                
                if formula['nominal_mz'] in dict_res[formula['classe']].keys():
                    
                    dict_res.get(formula['classe']).get(formula['nominal_mz']).append(json.loads(formula['mol_formula']) )
                
                else:

                    dict_res.get(formula['classe'])[formula['nominal_mz']] = [json.loads(formula['mol_formula']) ]  
        
            else:
                
                dict_res[formula['classe']] = {formula['nominal_mz']: [json.loads(formula['mol_formula'])] }     
        
        return dict_res
    
    def get_entries(self,classe, ion_type, nominal_mz):
        
        #print (classe, ion_type, nominal_mz)
        
        formulas = self.molform_collection.find( {'classe': classe, 
                                                'ion_type': ion_type,
                                                'nominal_mz': nominal_mz,  
                                                'O_C' : { '$lt': MSParameters.molecular_search.oc_filter }, 
                                                'H_C' : { '$gt': MSParameters.molecular_search.min_hc_filter},
                                                'DBE' : { '$gt': MSParameters.molecular_search.min_dbe},
                                                'DBE' : { '$lt': MSParameters.molecular_search.max_dbe},
                                                })
        

        return [json.loads(formula['mol_formula']) for formula in formulas]
        
   
    def update_entry(self, entry):
        
        raise NotImplementedError 

    def delete_entry(self, entry):
        
        raise NotImplementedError 
