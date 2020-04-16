import pickle, os

from sqlalchemy import create_engine, Column, Integer, Binary, String, Float, exists
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy.pool import QueuePool


Base = declarative_base()

class MolecularFormulaTable(Base):  
    
    id = Column( String, primary_key=True)
    mz = Column(Float, nullable=False)
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
    
    __tablename__ = 'molform'

    def __init__(self, kargs): 
        
        self.id = kargs['mol_formula']
        self.mz = kargs['mz']
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
        return "<MolecularFormulaTable(classe='%s', mz='%i', ion_type='%s', ion_charge='%i')>" % (
                                    self.classe, self.mz, self.ion_type, self.ion_charge)

class MolForm_SQL:
    
    def __init__(self, polarity, url=None, echo=False):

        self.engine = self.init_engine(polarity, url)

        #self.engine = create_engine('postgresql://postgres:docker@localhost:5432/')
        
        Base.metadata.create_all(self.engine)

        Session = sessionmaker(bind=self.engine)
        
        self.session = Session()

    def __exit__(self, exc_type, exc_val, exc_tb):
        # make sure the dbconnection gets closed
        
        self.commit()
        self.session.close()
        self.engine.dispose()

    def init_engine(self, polarity, url):
        directory = os.getcwd()
        
        if not url:
            
            if not os.path.isdir(directory+'/db'):
                    
                os.mkdir(directory+'/db')    

            polarity_label = 'pos' if polarity > 0 else 'neg'

            url = 'sqlite:///{DB}/db/molformulas_{charge}.sqlite'.format(DB=directory, charge=polarity_label)

        return create_engine(url, poolclass=QueuePool)

    def __enter__(self):
        
        return self
    
    def add_all(self, sql_molform_list):
        
        self.session.add_all(sql_molform_list)
        self.commit()

    def add_entry(self,sql_molform): 
        
        self.session.add(sql_molform)  
        self.session.commit()

    def commit(self):
        
        try:
            self.session.commit()  
        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))
            
    def get_dict_entries(self, classes, ion_type, nominal_mzs, molecular_search_settings):
        
        '''Known issue, when using SQLite:
         if the number of classes and nominal_m/zs are higher than 999 the query will fail
         Solution: use postgres or split query''' 
        
        def query_no_nominal(class_list):
            
            return self.session.query(MolecularFormulaTable).filter(
                MolecularFormulaTable.classe.in_(class_list), 
                MolecularFormulaTable.ion_type == ion_type,
                MolecularFormulaTable.DBE >= molecular_search_settings.min_dbe, 
                MolecularFormulaTable.DBE <= molecular_search_settings.max_dbe, 
                MolecularFormulaTable.ion_charge == molecular_search_settings.ion_charge,
                MolecularFormulaTable.O_C <= molecular_search_settings.oc_filter,
                MolecularFormulaTable.H_C >= molecular_search_settings.hc_filter)

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
                
            
        def add_dict_formula(formulas, check_nominal=False):
            "organize data by heteroatom classes"
            for formula in formulas:
                
                if formula.O and formula.P:

                    if  not (formula.O -2)/ formula.P >= molecular_search_settings.op_filter:
                        continue

                if check_nominal:
                    if formula.nominal_mz not in nominal_mzs:
                        continue

                if formula.classe in dict_res.keys():
                    
                    if formula.nominal_mz in dict_res[formula.classe].keys():
                        
                        dict_res.get(formula.classe).get(formula.nominal_mz).append(formula )
                    
                    else:

                        dict_res.get(formula.classe)[formula.nominal_mz] = [formula ]  
            
                else:
                    
                    dict_res[formula.classe] = {formula.nominal_mz: [formula] }     

        
        dict_res = {}

        if (len(classes) + len(nominal_mzs)) >= 998:
            
            formulas = query_no_nominal(classes)
            add_dict_formula(formulas, check_nominal=True)
            
                        
        else:

            formulas = query(classes)
            add_dict_formula(formulas)
            
        return dict_res

    def check_entry(self,classe, ion_type, molecular_search_settings):
        #TODO
        #  get all classes, ion_type, ion charge as str add to a dict or list
        #  then check if class in database
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
        
        return mol_formulas

    def get_entries(self,classe, ion_type, nominal_mz, molecular_search_settings):
        
        mol_formulas = self.session.query(MolecularFormulaTable).filter(
            #MolecularFormulaTable.nominal_mz == nominal_mz,
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
        return mol_formulas
       

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
   