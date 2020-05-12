import pickle, os

from sqlalchemy import text
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy import create_engine, Column, Integer, Binary, String, Float, exists
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy.pool import QueuePool, NullPool
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import event
from sqlalchemy import exc


from corems import chunks
import warnings
Base = declarative_base()

class MolFormMixin(object):

    id =  Column(Integer, primary_key=True)
    
    mol_formula = Column( String, unique=True)
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
    
class MolecularFormulaTablePos(MolFormMixin, Base):  

    __tablename__ = 'molformulas_pos'
    
class MolecularFormulaTableNeg(MolFormMixin, Base):  
    
    __tablename__ = 'molformulas_neg'

class MolForm_SQL:
    
    def __init__(self, polarity, url=None, echo=False):
        self.molform_model = MolecularFormulaTablePos if polarity > 0 else MolecularFormulaTableNeg
        
        self.engine = self.init_engine(polarity, url)

        self.add_engine_pidguard(self.engine)
        #self.engine = create_engine('postgresql://postgres:docker@localhost:5432/')
        
        Session = sessionmaker(bind=self.engine)

        self.session = Session()

        Base.metadata.create_all(self.engine)

    def __exit__(self, exc_type, exc_val, exc_tb):
        # make sure the dbconnection gets closed
        
        self.commit()
        self.session.close()
        self.engine.dispose()

    def init_engine(self, polarity, url):
        directory = os.getcwd()
        
        if not url:
            
            if not os.path.isdir(directory+'/db'):  os.mkdir(directory+'/db')
            
            url = 'sqlite:///{DB}/db/molformulas.sqlite'.format(DB=directory)

        if url[0:6] == 'sqlite':
            
            engine = create_engine(url, echo = False)
            self.chunks_count = 50
        
        else:
            
            engine = create_engine(url, echo = False, isolation_level="AUTOCOMMIT")
            self.chunks_count = 30000
        
        return engine# poolclass=NullPool

    def __enter__(self):
        
        return self
    
    def add_all_core(self, molform_list):

        #objs = [self.molform_model(**molform) for molform in molform_list]
        #print("ok")
        #self.session.bulk_save_objects(objs)
        #print("ok2")
        self.engine.execute(self.molform_model.__table__.insert(),
               [data for data in molform_list],autocommit=False)
 
    def add_entry(self, molform): 
        
        self.session.add(self.molform_model(**molform))  
        
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
            
            return self.session.query(self.molform_model).filter(
                self.molform_model.classe.in_(class_list), 
                self.molform_model.ion_type == ion_type,
                self.molform_model.DBE >= molecular_search_settings.min_dbe, 
                self.molform_model.DBE <= molecular_search_settings.max_dbe, 
                self.molform_model.ion_charge == molecular_search_settings.ion_charge,
                self.molform_model.O_C <= molecular_search_settings.oc_filter,
                self.molform_model.H_C >= molecular_search_settings.hc_filter)

        def query(class_list):
            
            return self.session.query(self.molform_model).filter(
                self.molform_model.nominal_mz.in_(nominal_mzs),
                self.molform_model.classe.in_(class_list), 
                self.molform_model.ion_type == ion_type,
                self.molform_model.DBE >= molecular_search_settings.min_dbe, 
                self.molform_model.DBE <= molecular_search_settings.max_dbe, 
                self.molform_model.ion_charge == molecular_search_settings.ion_charge,
                self.molform_model.O_C <= molecular_search_settings.oc_filter,
                self.molform_model.H_C >= molecular_search_settings.hc_filter)
                
        def add_dict_formula(formulas, check_nominal=False):
            "organize data by heteroatom classes"
            
            for formula in formulas:
                
                if formula.O and formula.P:

                    if  not (formula.O -2)/ formula.P >= molecular_search_settings.op_filter:
                        continue

                if check_nominal:
                    if len(nominal_mzs) > 999:            
                        #filter manually
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

        if len(classes) > 900 or (len(classes) + len(nominal_mzs)) >= 900:
            
            for class_chunk in chunks(classes, 500):
                
                formulas = query_no_nominal(class_chunk)
                
                if len(nominal_mzs) < 100:

                    formulas.filter(self.molform_model.nominal_mz.in_(nominal_mzs))
                
                add_dict_formula(formulas, check_nominal=True)
        
        else:

            formulas = query(classes)
            add_dict_formula(formulas)

        # dump all objs to memory
        self.session.expunge_all()
        
        return dict_res

    def check_entry(self,classe, ion_type, molecular_search_settings):
        #TODO
        #  get all classes, ion_type, ion charge as str add to a dict or list
        #  then check if class in database
        yes = self.session.query(exists().where(
            (self.molform_model.classe == classe) &
            (self.molform_model.ion_type == ion_type) &
            (self.molform_model.ion_charge == molecular_search_settings.ion_charge))).scalar()
        
        return yes
    
    
    def get_all_classes(self, ion_type, molecular_search_settings):
        
        query = self.session.query(self.molform_model.classe.distinct().label("classe"))
        
        classes = [row.classe for row in query.filter(self.molform_model.ion_type == ion_type,
             self.molform_model.ion_charge == molecular_search_settings.ion_charge)]
        
        return classes  
    
    def get_all(self,):
        
        mol_formulas = self.session.query(self.molform_model).all()
        
        #mol_formulas = mol_formulas.filter(ion_type = ion_type)

        #mol_formulas = mol_formulas.filter(ion_charge = molecular_search_settings.ion_charge)
        
        return mol_formulas

    def get_entries(self,classe, ion_type, nominal_mz, molecular_search_settings):
        
        mol_formulas = self.session.query(self.molform_model).filter(
            #self.molform_model.nominal_mz == nominal_mz,
            self.molform_model.classe == classe, 
            self.molform_model.ion_type == ion_type,
            self.molform_model.DBE >= molecular_search_settings.min_dbe, 
            self.molform_model.DBE <= molecular_search_settings.max_dbe, 
            self.molform_model.ion_charge == molecular_search_settings.ion_charge,
            self.molform_model.O_C <= molecular_search_settings.oc_filter,
            self.molform_model.H_C >= molecular_search_settings.hc_filter,
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
        self.session.query(self.molform_model).delete()
        self.session.commit()  

    def clear_data(self):
        '''clear tables'''
        meta = Base.metadata
        for table in reversed(meta.sorted_tables):
            print ('Clear table %s' % table)
            self.session.execute(table.delete())
        self.session.commit()

    def close(self):
        # make sure the dbconnection gets closed
        
        self.commit()
        self.session.close()
        self.engine.dispose()    
   
    def add_engine_pidguard(self, engine):
        """Add multiprocessing guards.

        Forces a connection to be reconnected if it is detected
        as having been shared to a sub-process.

        """

        @event.listens_for(engine, "connect")
        def connect(dbapi_connection, connection_record):
            connection_record.info['pid'] = os.getpid()

        @event.listens_for(engine, "checkout")
        def checkout(dbapi_connection, connection_record, connection_proxy):
            pid = os.getpid()
            if connection_record.info['pid'] != pid:
                # substitute log.debug() or similar here as desired
                warnings.warn(
                    "Parent process %(orig)s forked (%(newproc)s) with an open "
                    "database connection, "
                    "which is being discarded and recreated." %
                    {"newproc": pid, "orig": connection_record.info['pid']})
                connection_record.connection = connection_proxy.connection = None
                raise exc.DisconnectionError(
                    "Connection record belongs to pid %s, "
                    "attempting to check out in pid %s" %
                    (connection_record.info['pid'], pid)
                )    
   