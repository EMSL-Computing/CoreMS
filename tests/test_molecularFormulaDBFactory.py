__author__ = "Yuri E. Corilo"
__date__ = "Jul 22, 2019"


import pickle

from pathlib import Path
import time, sys, os, pytest
sys.path.append(".")

from corems.encapsulation.constant import Labels
from corems.molecular_id.factory.MolecularLookupTable import  MolecularCombinations
from corems.molecular_id.factory.molecularSQL import MolForm_SQL
from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.encapsulation.factory.processingSetting  import MolecularFormulaSearchSettings

def test_nist_to_sql():

    file_location = Path.cwd() / "tests/tests_data/gcms/" / "PNNLMetV20191015.MSL"

    sqlLite_obj = ReadNistMSI(file_location).get_sqlLite_obj()

    sqlLite_obj.query_min_max_ri((1637.30, 1638.30)) 
    sqlLite_obj.query_min_max_rt((17.111, 18.111))            
    sqlLite_obj.query_min_max_ri_and_rt((1637.30, 1638.30),(17.111, 18.111)) 

def test_query_sql():

    with MolForm_SQL() as sqldb:
        #sqldb.clear_data()

        ion_type = Labels.protonated_de_ion
        print('ion_type', ion_type)
        classe = ['{"O": 2}']
        nominal_mz = [301]
        results = sqldb.get_dict_by_classes(classe, ion_type, nominal_mz, +1, MolecularFormulaSearchSettings())
        
        #print('total mol formulas found: ', len(list( results.get(classe[0]).get(301))))

def generate_database():
    
    '''corems_parameters_file: Path for CoreMS JSON Parameters file
       --jobs: Number of processes to run   
    '''
    
    #url = "postgresql://postgres:labthomson0102@172.22.113.27:5432/"
    #url = "postgresql://doadmin:rn9fenbsdbwqis9v@db-postgresql-corems-do-user-7454084-0.a.db.ondigitalocean.com:25060/defaultdb?sslmode=require"
    #url = "postgresql://postgres:qqmica@34.71.74.212/defaultdb?sslmode=require"
    #molecular_search_settings.url_database = url
    #molecular_search_settings.db_jobs = jobs
    
    molecular_search_settings = MolecularFormulaSearchSettings()
    #molecular_search_settings.usedAtoms['C'] = (1,5)
    #molecular_search_settings.usedAtoms['O'] = (1,20)
    MolecularCombinations().runworker(molecular_search_settings)

if __name__ == '__main__':
    

    test_query_sql()
    #settings_parsers.load_search_setting_yaml()
    #settings_parsers.load_search_setting_json()
    #test_nist_to_sql()    
    #generate_database()
    
   