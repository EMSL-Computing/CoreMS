import pickle

from pymongo import MongoClient
from enviroms.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaSearchSettings

class MolForm_Mongo:

    def __exit__(self, exc_type, exc_val, exc_tb):
        # make sure the dbconnection gets closed
        self.client.close()

    def __enter__(self):
        
        self.client = MongoClient("mongodb://enviroms-client:esmlpnnl2019@localhost:27017/enviroms")
        
        db = self.client.enviroms
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
        # this is way too slow, create a pos and neg table
        #try:
        #yes = self.session.query(MolecularFormulaTable.id).filter(MolecularFormulaTable.classe==classe).filter(MolecularFormulaTable.ion_charge == MoleculaSearchSettings.ion_charge).scalar() is not None
        
        #except MultipleResultsFound as e:
        #    yes = True
        #except MultipleResultsFound as e:
        #    yes = True
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

    def get_dict_entries(self, classes, ion_type, nominal_mzs):

        dict_res = {}

        formulas = self.molform_collection.find( {'classe': classes, 
                                                'ion_type': ion_type,
                                                'nominal_mz': nominal_mzs,  
                                                'O_C' : { '$gt': MoleculaSearchSettings.oc_filter }, 
                                                'H_C' : { '$gt': MoleculaSearchSettings.hc_filter},
                                                'DBE' : { '$gt': MoleculaSearchSettings.min_dbe},
                                                'DBE' : { '$lt': MoleculaSearchSettings.max_dbe},
                                                })

        for formula in formulas:
            
            if formula['classes'] in dict_res.keys():
                
                if formula['nominal_mz'] in dict_res[formula['classes']].keys():
                    
                    dict_res.get(formula['classes']).get(formula['nominal_mz']).append(formula['mol_formula'])
                
                else:

                    dict_res.get(formula['classes'])[formula['nominal_mz']] = [formula['mol_formula']]  
        
            else:
                
                dict_res[formula['classes']] = {formula['nominal_mz']: [formula['mol_formula']]}     
        
        return dict_res
    
    def get_entries(self,classe, ion_type, nominal_mz):
        
        #print (classe, ion_type, nominal_mz)
        
        formulas = self.molform_collection.find( {'classe': classe, 
                                                'ion_type': ion_type,
                                                'nominal_mz': nominal_mz,  
                                                'O_C' : { '$lt': MoleculaSearchSettings.oc_filter }, 
                                                'H_C' : { '$gt': MoleculaSearchSettings.hc_filter},
                                                'DBE' : { '$gt': MoleculaSearchSettings.min_dbe},
                                                'DBE' : { '$lt': MoleculaSearchSettings.max_dbe},
                                                })
        

        return [pickle.loads(formula['mol_formula']) for formula in formulas]
        
   
    def update_entry(self, entry):
        
        raise NotImplementedError 

    def delete_entry(self, entry):
        
        raise NotImplementedError 
