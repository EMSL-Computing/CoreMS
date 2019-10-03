from sqlalchemy import create_engine, Column, Integer, Binary, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker

Base = declarative_base()

class MolecularFormulaTable(Base):  
    __tablename__ = 'molform'

    id = Column( Binary, primary_key=True)
    nominal_mass= Column(Integer, nullable=False)
    ion_type = Column(String, nullable=False)
    ion_charge = Column(Integer, nullable=False)
    classe = Column(String, nullable=False)

    def __init__(self, kargs): 
        
        self.id = kargs['mol_formula']
        self.nominal_mass =kargs['nominal_mass']
        self.ion_type = kargs['ion_type']
        self.ion_charge = kargs['ion_charge']
        self.classe = kargs['classe']
       
    def __repr__(self):
        return "<MolecularFormulaTable(classe='%s', nominal_mass='%i', ion_type='%s', ion_charge='%i')>" % (
                                    self.classe, self.nominal_mass, self.ion_type, self.ion_charge)

class MolForm_SQL:
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        # make sure the dbconnection gets closed
        self.commit()
        self.session.close()
        self.engine.dispose()

    def __enter__(self):
        
        self.engine = create_engine('sqlite:///{DB}'.format(DB='molformulas.sqlite'), connect_args={'timeout': 15})
        
        Base.metadata.create_all(self.engine)

        Session = sessionmaker(bind=self.engine)
        
        self.session = Session()

        return self
    
    def add_all(self, sql_molform_list):
        
        self.session.add_all( [MolecularFormulaTable(sql_molform_dict)  for sql_molform_dict in sql_molform_list] )
    
    def add_entry(self,sql_molform): 

        one_formula = MolecularFormulaTable(sql_molform)  
        self.session.add(one_formula)  

    def commit(self):
        
        try:
            self.session.commit()  
        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))
            
    def read_entry(self,):
        
        mol_formulas = self.session.query(MolecularFormulaTable).filter(MolecularFormulaTable.classe == 'O2')
        
        mol_formulas = mol_formulas.filter(MolecularFormulaTable.ion_charge == -1)
        
        for mol_formula in mol_formulas:  
            print(mol_formula.nominal_mass)

    def update_entry(self, entry):
        
        entry.title = "Some2016Film"  
        self.session.commit()

    def delete_entry(self, entry):
        
        try:
            self.session.delete(entry)  
            self.session.commit()  
        
        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))

        