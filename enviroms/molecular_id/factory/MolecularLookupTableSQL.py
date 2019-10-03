from sqlalchemy import create_engine, Column, Integer, Binary, String
from sqlalchemy.ext.declarative import declarative_base  
from sqlalchemy.orm import sessionmaker



base = declarative_base()
class MolecularFormulaTable(base):  
    __tablename__ = 'films'

    Column('mol_formula', Binary, primary_key=True),
    Column('nominal_mass', Integer),
    Column('ion_type', String),
    Column('ion_charge', String),
    Column('classe', String),

class MolForm_SQL:

    def __init__(self,):
        
        url = 'sqlite:///{DB}'.format(DB='mydb.sqlite')
        db = create_engine(url)  
        Session = sessionmaker(db)
        self.session = Session()
        base.metadata.create_all(db)
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        # make sure the dbconnection gets closed
        self.session.close()

    def add_all(self, sql_molform_list):

        self.session.add_all( [MolecularFormulaTable(**sql_molform_dict)  for sql_molform_dict in sql_molform_list] )
    
    def add_entry(self,sql_molform): 

        one_formula = MolecularFormulaTable(**sql_molform)  
        self.session.add(one_formula)  

    def commit(self):
        try:
            self.session.commit()  
        except:
            self.session.rollback()
            raise
            
    def read_entry(self,):
        mol_formulas = self.session.query(MolecularFormulaTable)  
        for mol_formula in mol_formulas:  
            print(mol_formula.nominal_mass)

    def update_entry(self, entry):
        entry.title = "Some2016Film"  
        self.session.commit()

    def delete_entry(self, entry):
        try:
            self.session.delete(entry)  
            self.session.commit()  
        except:
            self.session.rollback()
            raise

        