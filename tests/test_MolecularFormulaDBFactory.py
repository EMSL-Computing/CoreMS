__author__ = "Yuri E. Corilo"
__date__ = "Jul 22, 2019"


from bson.binary import Binary
import pickle

from pymongo import MongoClient
import pymongo
import time, sys, os, pytest
sys.path.append(".")

from corems.encapsulation.constant import Atoms, Labels
from corems.molecular_id.factory.MolecularLookupTable import  MolecularCombinations
from corems.molecular_id.factory.molecularSQL import MolForm_SQL
from corems.molecular_id.factory.molecularMongo import MolForm_Mongo
from corems.molecular_id.output.export import  MolecularLookUpDictExport
from corems.encapsulation.settings.molecular_id.MolecularIDSettings import MoleculaLookupDictSettings, MoleculaSearchSettings

def create_lookup_dict():
    
    MolecularCombinations().runworker()

def xtest_query_mongo():
    #client = MongoClient("mongodb://corems-client:esmlpnnl2019@localhost:27017/corems")
        
    #db = client.corems.drop_collection('molform')
    with MolForm_Mongo() as mongo_db:
       formula = mongo_db.get_all()
       print(formula[0])
       print(pickle.loads(formula[0]['mol_formula']).mz_theor)

def test_query_sql():

    with MolForm_SQL() as sqldb:
        #sqldb.get_all()

        ion_type = Labels.protonated_de_ion
        print(ion_type)
        classe = 'O8'
        nominal_mz = 501
        print(len(sqldb.get_entries(classe, ion_type, nominal_mz)))

def test_molecular_lookup_db():    
    
    #margin_error needs to be optimized by the data rp and sn
    #min_mz,max_mz  needs to be optimized by the data
    MoleculaSearchSettings.min_mz = 100
    MoleculaSearchSettings.max_mz = 800
    # C, H, N, O, S and P atoms are ALWAYS needed in the dictionary
    #the defaults values are defined at the encapsulation MolecularSpaceTableSetting    
    MoleculaSearchSettings.usedAtoms['C'] = (1,90)
    MoleculaSearchSettings.usedAtoms['H'] = (4,200)
    MoleculaSearchSettings.usedAtoms['O'] = (1,8)
    MoleculaSearchSettings.usedAtoms['N'] = (0,0)
    MoleculaSearchSettings.usedAtoms['S'] = (0,0)

    MoleculaSearchSettings.isRadical = True
    MoleculaSearchSettings.isProtonated = False
    #some atoms has more than one covalence state and the most commun will be used
    # adduct atoms needs covalence 0
    MoleculaSearchSettings.usedAtoms['Cl'] = (0,0)
    possible_valences = Atoms.atoms_covalence.get('Cl')
    valence_one = possible_valences[0]
    
    # if you want to specify it in needs to be changed here
    # otherwise it will use the lowest covalence, PS needs insure propagation to isotopologues
    MoleculaSearchSettings.used_atom_valences['Cl'] =  valence_one
    
    time0 = time.time()
    create_lookup_dict()
    time1 = time.time()
    print("create the molecular lookup table took %.2f seconds", time1-time0)
    
if __name__ == '__main__':
    
    test_molecular_lookup_db()
    test_query_sql()
    #xtest_query_mongo()
    