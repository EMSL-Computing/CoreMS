import pickle, os

from sqlalchemy import create_engine, Column, Integer, LargeBinary, String, Float, exists
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy.pool import QueuePool

Base = declarative_base()

class MolecularFormulaTable(Base):  
    
    __tablename__ = 'molform'

    id = Column( LargeBinary, primary_key=True)
    nominal_mz= Column(Integer, nullable=False)
    ion_type = Column(String, nullable=False)
    ion_charge = Column(Integer, nullable=False)
    classe = Column(String, nullable=False)
    
    C = Column(Integer, nullable=False)
    H = Column(Integer, nullable=True)
    N = Column(Integer, nullable=True)
    O = Column(Integer, nullable=True)
    S = Column(Integer, nullable=True)
    P = Column(Integer, nullable=True)
    DBE = Column(Float, nullable=False)
    O_C = Column(Float, nullable=True)
    H_C = Column(Float, nullable=True)
    
    def __init__(self, kargs): 
        
        self.id = kargs['mol_formula']
        self.nominal_mz =kargs['nominal_mz']
        self.ion_type = kargs['ion_type']
        self.ion_charge = kargs['ion_charge']
        self.classe = kargs['classe']
        self.C = kargs['C']
        self.H = kargs['H']
        self.N = kargs['N']
        self.O = kargs['O']
        self.S = kargs['S']
        self.P = kargs['P']
        self.H_C = kargs['H_C']
        self.O_C = kargs['O_C']
        self.DBE = kargs['DBE']
       
    def __repr__(self):
        return "<MolecularFormulaTable(classe='%s', nominal_mass='%i', ion_type='%s', ion_charge='%i')>" % (
                                    self.classe, self.nominal_mz, self.ion_type, self.ion_charge)

class MolForm_SQL:
    
    def __init__(self, url='sqlite://', echo=False):

        self.engine = create_engine(url)       
        #self.engine = create_engine('postgresql://postgres:docker@localhost:5432/')
        
        Base.metadata.create_all(self.engine)

        Session = sessionmaker(bind=self.engine)
        
        self.session = Session()

    def __exit__(self, exc_type, exc_val, exc_tb):
        # make sure the dbconnection gets closed
        self.commit()
        self.session.close()
        self.engine.dispose()

    def __enter__(self):
        
        directory = os.getcwd()
       
        if not os.path.isdir(directory+'/db'):
                
            os.mkdir(directory+'/db')    
            
        self.engine = create_engine('sqlite:///{DB}'.format(DB=directory+'/db'+'/molformulas.sqlite'), 
                    poolclass=QueuePool)
        
        #self.engine = create_engine('postgresql://postgres:docker@localhost:5432/')
        
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
            
    def get_dict_entries(self, classes, ion_type, nominal_mzs, molecular_search_settings):
        
        '''Known issue, when using SQLite:
         if the number of classes and nominal_m/zs are higher than 1,998 the query will fail
         Solution: use postgres or split query''' 
        
        print("Started molecular formula table generation")
        
        def query(class_list):
            
            return self.session.query(MolecularFormulaTable).filter(
                MolecularFormulaTable.nominal_mz.in_(nominal_mzs),
                MolecularFormulaTable.classe.in_(class_list), 
                MolecularFormulaTable.ion_type == ion_type,
                MolecularFormulaTable.DBE >= molecular_search_settings.min_dbe, 
                MolecularFormulaTable.DBE <= molecular_search_settings.max_dbe, 
                MolecularFormulaTable.ion_charge == molecular_search_settings.ion_charge,
                MolecularFormulaTable.O_C <= molecular_search_settings.oc_filter,
                MolecularFormulaTable.H_C >= molecular_search_settings.hc_filter)
                
            
        def add_dict_formula(formulas):
            
            for formula in formulas:
                
                if formula.O and formula.P:

                    if  not (formula.O -2)/ formula.P >= molecular_search_settings.op_filter:
                        continue

                if formula.classe in dict_res.keys():
                    
                    if formula.nominal_mz in dict_res[formula.classe].keys():
                        
                        dict_res.get(formula.classe).get(formula.nominal_mz).append(pickle.loads(formula.id) )
                    
                    else:

                        dict_res.get(formula.classe)[formula.nominal_mz] = [pickle.loads(formula.id) ]  
            
                else:
                    
                    dict_res[formula.classe] = {formula.nominal_mz: [pickle.loads(formula.id)] }     

        
        dict_res = {}

        if (len(classes) + len(nominal_mzs)) >= 993:
            
            classes_chunk = int(len(classes)/2)
            
            formulas = query(classes[:classes_chunk])
            add_dict_formula(formulas)
            
            formulas = query(classes[classes_chunk:])
            add_dict_formula(formulas)
            
        else:

            formulas = query(classes)
            add_dict_formula(formulas)
            
            print("Finished molecular formula table generation")
            print()

        return dict_res

    def check_entry(self,classe, ion_type, molecular_search_settings):
        #TODO
        #  get all classes, ion_type, ion charge as str add to a dict or list
        #  then check if class in database
        # 
        #  
        # this is way too slow, create a pos and neg table
        #try:
        #yes = self.session.query(MolecularFormulaTable.id).filter(MolecularFormulaTable.classe==classe).filter(MolecularFormulaTable.ion_charge == molecular_search_settings.ion_charge).scalar() is not None
        
        #except MultipleResultsFound as e:
        #    yes = True
        #except MultipleResultsFound as e:
        #    yes = True
        yes = self.session.query(exists().where(
            (MolecularFormulaTable.classe == classe) &
            (MolecularFormulaTable.ion_type == ion_type) &
            (MolecularFormulaTable.ion_charge == molecular_search_settings.ion_charge))).scalar()
        return yes
    
    
    def get_all_classes(self, ion_type, molecular_search_settings):
        
        query = self.session.query(MolecularFormulaTable.classe.distinct().label("classe"))
        
        classes = [row.classe for row in query.filter(MolecularFormulaTable.ion_type == ion_type,
             MolecularFormulaTable.ion_charge == molecular_search_settings.ion_charge)]
        
        return classes  
    
    def get_all(self,):
        
        mol_formulas = self.session.query(MolecularFormulaTable).all()
        
        #mol_formulas = mol_formulas.filter(ion_type = ion_type)

        #mol_formulas = mol_formulas.filter(ion_charge = molecular_search_settings.ion_charge)
        
        return [pickle.loads(formula.id) for formula in mol_formulas]

    def get_entries(self,classe, ion_type, nominal_mz, molecular_search_settings):
        
        mol_formulas = self.session.query(MolecularFormulaTable).filter(
            MolecularFormulaTable.nominal_mz == nominal_mz,
            MolecularFormulaTable.classe == classe, 
            MolecularFormulaTable.ion_type == ion_type,
            MolecularFormulaTable.DBE >= molecular_search_settings.min_dbe, 
            MolecularFormulaTable.DBE <= molecular_search_settings.max_dbe, 
            MolecularFormulaTable.ion_charge == molecular_search_settings.ion_charge,
            MolecularFormulaTable.O_C <= molecular_search_settings.oc_filter,
            MolecularFormulaTable.H_C >= molecular_search_settings.hc_filter,
            )
        
        #mol_formulas = mol_formulas.filter(ion_type = ion_type)

        #mol_formulas = mol_formulas.filter(ion_charge = molecular_search_settings.ion_charge)
        return [pickle.loads(formula.id) for formula in mol_formulas]
       

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

    def purge(self):
        '''Carefull, this will delete the entire database table'''
        self.session.query(MolecularFormulaTable).delete()
        self.session.commit()  

    def clear_data(self):
        '''clear tables'''
        meta = Base.metadata
        for table in reversed(meta.sorted_tables):
            print ('Clear table %s' % table)
            self.session.execute(table.delete())
        self.session.commit()
   