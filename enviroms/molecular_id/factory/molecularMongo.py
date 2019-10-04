import pickle

from pymongo import MongoClient

class MolForm_Mongo:

    def __exit__(self, exc_type, exc_val, exc_tb):
        # make sure the dbconnection gets closed
        self.client.close()

    def __enter__(self):
        
        self.client = MongoClient("mongodb://enviroms-client:esmlpnnl2019@localhost:27017/enviroms")
        
        db = self.client.enviroms
        
        self.molform_collection = db.molform
    
    def add_all(self, dic_molform_list):
        
        try:
            self.molform_collection.insert_many(dic_molform_list, ordered=False)
        except Exception as e:
            #duplicate_error
            try:
                if e.details["writeErrors"][0]["code"] == 11000:
                    pass
            except:
                raise  Exception(e)        
    
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
  
    def read_entry(self,):
        
        formulas = self.molform_collection.find({'classe': "O2"})

        print(len(list(formulas)))
        
        #print(formulas)
        for formula in formulas:
        #    print()
            print(pickle.loads(formula['mol_formula']).to_dict)    

    def update_entry(self, entry):
        
        raise NotImplementedError 

    def delete_entry(self, entry):
        
        raise NotImplementedError 
