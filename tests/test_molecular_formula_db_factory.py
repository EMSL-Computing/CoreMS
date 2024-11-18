from pathlib import Path

from corems.encapsulation.constant import Labels
from corems.molecular_id.factory.molecularSQL import MolForm_SQL
from corems.molecular_id.input.nistMSI import ReadNistMSI
from corems.encapsulation.factory.processingSetting  import MolecularFormulaSearchSettings

def test_nist_to_sql():

    file_location = Path.cwd() / "tests/tests_data/gcms/" / "PNNLMetV20191015.MSL"

    sqlLite_obj = ReadNistMSI(file_location).get_sqlLite_obj()

    response = sqlLite_obj.query_min_max_ri((1637.30, 1638.30)) 
    assert len(response) == 6

    response = sqlLite_obj.query_min_max_rt((17.111, 18.111))          
    assert len(response) == 137

    response = sqlLite_obj.query_min_max_ri_and_rt((1637.30, 1638.30),(17.111, 18.111)) 
    assert len(response) == 6

def test_query_sql():

    sqldb = MolForm_SQL()

    ion_type = Labels.protonated_de_ion
    classe = ['{"O": 2}']
    nominal_mz = [301]
    results = sqldb.get_dict_by_classes(classe, ion_type, nominal_mz, +1, MolecularFormulaSearchSettings())
    assert len(results.get(classe[0]).get(301)) == 3
   