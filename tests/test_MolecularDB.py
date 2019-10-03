__author__ = "Yuri E. Corilo"
__date__ = "Jul 22, 2019"


from bson.binary import Binary
import pickle

from pymongo import MongoClient
import pymongo
import time, sys, os, pytest
sys.path.append(".")

from enviroms.encapsulation.constant import Atoms, Labels
from enviroms.molecular_id.factory.MolecularLookupTableDB import  MolecularCombinations
from enviroms.molecular_id.factory.MolecularSQLBaseClass import MolForm_SQL
from enviroms.molecular_id.output.export import  MolecularLookUpDictExport
from enviroms.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupDictSettings

def create_lookup_dict(LookupTableSettings):
    
    MolecularCombinations().runworker(LookupTableSettings)

def xtest_query_mongo():

    client = MongoClient("mongodb://enviroms-client:esmlpnnl2019@localhost:27017/enviroms")
    db = client.enviroms
    molform_collection = db.molform
    
    #molform_collection.drop() 
    
    #molform_collection.drop_indexes()
    
    #molform_collection.create_index("mol_formula", unique= True )
    
    formulas = molform_collection.find({'classe': "O2"})

    print(len(list(formulas)))
    
    #print(formulas)
    for formula in formulas:
    #    print()
        print(pickle.loads(formula['mol_formula']).to_dict)    
    client.close()

def xtest_query_sql():

    with MolForm_SQL() as sqldb:

        sqldb.read_entry()
    

def xtest_molecular_lookup_db():    
    
    LookupDictSettings = MoleculaLookupDictSettings()
    #margin_error needs to be optimized by the data rp and sn
    #min_mz,max_mz  needs to be optimized by the data
    LookupDictSettings.min_mz = 100
    LookupDictSettings.max_mz = 1200
    # C, H, N, O, S and P atoms are ALWAYS needed in the dictionary
    #the defaults values are defined at the encapsulation MolecularSpaceTableSetting    
    LookupDictSettings.usedAtoms['C'] = (1,90)
    LookupDictSettings.usedAtoms['H'] = (4,200)
    LookupDictSettings.usedAtoms['O'] = (0,3)
    LookupDictSettings.usedAtoms['N'] = (0,3)
    LookupDictSettings.usedAtoms['S'] = (0,3)

    LookupDictSettings.isRadical = True
    #some atoms has more than one covalence state and the most commun will be used
    # adduct atoms needs covalence 0
    LookupDictSettings.usedAtoms['Cl'] = (0,0)
    possible_valences = Atoms.atoms_covalence.get('Cl')
    valence_one = possible_valences[0]
    # if you want to specify it in needs to be changed here
    # otherwise it will use the lowest covalence, PS needs insure propagation to isotopologues
    MoleculaLookupDictSettings.used_atom_valences['Cl'] =  valence_one

    time0 = time.time()
    create_lookup_dict(LookupDictSettings)
    time1 = time.time()
    print("create the molecular lookup table took %.2f seconds", time1-time0)
    
if __name__ == '__main__':
    
    xtest_molecular_lookup_db()
    #xtest_query_sql()
    #xtest_query()

