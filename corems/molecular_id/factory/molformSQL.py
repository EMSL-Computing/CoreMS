import sys
sys.path.append(".")
import os

from sqlalchemy import create_engine, ForeignKey, Column, Integer, String, Float, SMALLINT
from sqlalchemy.orm import backref, column_property, relationship
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.sql.schema import UniqueConstraint
from sqlalchemy import exc

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.sql.operators import exists

from corems.encapsulation.constant import Atoms
Base = declarative_base()

class HeteroAtoms(Base):
   
    __tablename__ = 'heteroAtoms'
    
    id = Column(Integer, primary_key=True,
                        unique=True,
                        nullable=False)

    name = Column(String, unique=True,  nullable=False)
    
    carbonHydrogen = relationship('CarbonHydrogen', secondary = 'molecularformula')

    def __repr__(self):
        return '<HeteroAtoms Model {} class {}>'.format(self.id, self.name)      
    
class CarbonHydrogen(Base):
   
    __tablename__ = 'carbonHydrogen'
    __table_args__ = ( UniqueConstraint('C', 'H', name='unique_c_h'), )

    id = Column(Integer, primary_key=True,
                        unique=True,
                        nullable=False)

    C = Column(Integer, nullable=False)
    
    H = Column(Integer,  nullable=False)
    
    heteroAtoms = relationship("HeteroAtoms",
                        secondary="molecularformula",
                        backref=backref('carbon_hydrogen'))

    def __repr__(self):
        return '<CarbonHydrogen Model {} C{} H{}>'.format(self.id, self.C, self.H)                     

    @property
    def mass(self):
        return (self.C * Atoms.atomic_masses.get('C')) + (self.H * Atoms.atomic_masses.get('H'))

    @property
    def dbe(self):
        return self.C - (self.H/2) + 1

#264888.88 ms
class MolecularFormulaLink(Base):
    
    __tablename__ = 'molecularformula'
    __table_args__ = ( UniqueConstraint('heteroAtoms_id', 'carbonHydrogen_id', name='unique_molform'), )
    
    #id = Column(Integer, primary_key=True,
    #                    unique=True,
    #                    nullable=False)

    heteroAtoms_id = Column(
        Integer, 
        ForeignKey('heteroAtoms.id'), 
        primary_key=True)

    carbonHydrogen_id = Column(
        Integer, 
        ForeignKey('carbonHydrogen.id'), 
        primary_key=True)
    
    mass = Column(Float)
    
    DBE = Column(Float)
    
    C = association_proxy('carbonHydrogen', 'C')

    H = association_proxy('carbonHydrogen', 'H')

    classe = association_proxy('heteroAtoms', 'name')

    carbonHydrogen = relationship(CarbonHydrogen, backref=backref("heteroAtoms_assoc"))
    heteroAtoms = relationship(HeteroAtoms, backref=backref("carbonHydrogen_assoc"))

    def __repr__(self):
        
        return '<MolecularFormulaLink Model C{} H{} {}>'.format(self.C, self.H, self.classe)       

class MolForm_SQL:
    
    def __init__(self, polarity, url=None, echo=False):
        
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

    def commit(self):
        
        try:
            self.session.commit()  
        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))

    def init_engine(self, polarity, url):
        
        directory = os.getcwd()
        
        if not url:
            
            if not os.path.isdir(directory+'/db'):  os.mkdir(directory+'/db')
            
            url = 'sqlite:///{DB}/db/molformulas.sqlite'.format(DB=directory)

        return create_engine(url)# poolclass=NullPool

    def __enter__(self):
        
        return self
    
    def add_all_core(self, molform_list):

        all_data = [data for data in molform_list]
        
        self.engine.execute(MolecularFormulaLink.__table__.insert(),
               [data for data in molform_list],autocommit=False)

    def add_all(self, molform_list): 
        
        self.session.bulk_save_objects(molform_list)  

    def add_entry(self, molform): 
        
        self.session.add(molform)  
        
    def get_dict_entries(self, classes, ion_type, nominal_mzs, molecular_search_settings):
        
        '''Known issue, when using SQLite:
         if the number of classes and nominal_m/zs are higher than 999 the query will fail
         Solution: use postgres or split query''' 
        
        def query_no_nominal(class_list):
            
            return self.session.query(MolecularFormulaLink).filter(
                MolecularFormulaLink.classe.in_(class_list), 
                MolecularFormulaLink.ion_type == ion_type,
                MolecularFormulaLink.DBE >= molecular_search_settings.min_dbe, 
                MolecularFormulaLink.DBE <= molecular_search_settings.max_dbe, 
                MolecularFormulaLink.O_C <= molecular_search_settings.oc_filter,
                MolecularFormulaLink.H_C >= molecular_search_settings.hc_filter)

        def query(class_list):
            
            return self.session.query(MolecularFormulaLink).filter(
                MolecularFormulaLink.nominal_mz.in_(nominal_mzs),
                MolecularFormulaLink.classe.in_(class_list), 
                MolecularFormulaLink.DBE >= molecular_search_settings.min_dbe, 
                MolecularFormulaLink.DBE <= molecular_search_settings.max_dbe, 
                MolecularFormulaLink.O_C <= molecular_search_settings.oc_filter,
                MolecularFormulaLink.H_C >= molecular_search_settings.hc_filter)
        
    def check_entry(self,classe, ion_type, molecular_search_settings):
        #  get all classes, ion_type, ion charge as str add to a dict or list
        #  then check if class in database
        has_class = self.session.query(exists().where(
            (MolecularFormulaLink.classe == classe)))
        
        return has_class
    
    
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

    def get_entries(self,classe, nominal_mz, molecular_search_settings):
        
        mol_formulas = self.session.query(self.molform_model).filter(
            #self.molform_model.nominal_mz == nominal_mz,
            self.molform_model.classe == classe, 
            self.molform_model.DBE >= molecular_search_settings.min_dbe, 
            self.molform_model.DBE <= molecular_search_settings.max_dbe, 
            self.molform_model.O_C <= molecular_search_settings.oc_filter,
            self.molform_model.H_C >= molecular_search_settings.hc_filter,
            )
        
        return mol_formulas
    
    def delete_entry(self, row):
        
        try:
            self.session.delete(row)  
            self.session.commit()  
        
        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))

    def purge(self, cls):
        '''Carefull, this will delete the entire database table'''
        self.session.query(cls).delete()
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
        import os, warnings
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

if __name__ == "__main__":
    
    #from sqlalchemy.orm import sessionmaker
    
    #session = corems_session()
    
    #CH = CarbonHydrogen(C=10, H=21)
    
    
    #classe = HeteroAtoms(name = "O2")
    #classe.carbonHydrogen.append(ch)
    #ch = CarbonHydrogen.query().filter(CarbonHydrogen.C == 10, CarbonHydrogen.H == 22).first()
    
    #classe.mf_property.append(mf)
    #classe.carbonHydrogen.append(ch)
    #session.add(classe)
    
    #for molecularFormulas in HeteroAtoms.query().filter(HeteroAtoms.name == 'O2',  HeteroAtoms.C == 10):
    #    print(molecularFormulas.name, molecularFormulas.carbon_hydrogen.H, molecularFormulas.carbon_hydrogen.C)
    
    #for each in classe.query().all():
    #    print(each.name)
    
    session.commit()

    CH = CarbonHydrogen(C=10, H=23)
    HA = HeteroAtoms(name = "O4")
    molecularFormula = MolecularFormulaLink(heteroAtoms=HA, carbonHydrogen=CH, mass=200, dbe=1)
    
    session.add(molecularFormula)
    session.commit()
    #for x in session.query( CarbonHydrogen, HeteroAtoms).filter(MolecularFormulaLink.heteroAtoms_id == HeteroAtoms.id, 
    #    MolecularFormulaLink.carbonHydrogen_id == CarbonHydrogen.id).order_by(MolecularFormulaLink.carbonHydrogen_id).all():
    #    print ("MolecularFormulaLink: {} {}".format(x.CarbonHydrogen.name, x.HeteroAtoms.name))
    #    print ("MolecularFormulaLink: {}".format(x.molecularFormula()))
    
    #session.close()

def join_example():
    records = session.query(CarbonHydrogen).\
    	join(HeteroAtoms, HeteroAtoms.id == CarbonHydrogen.heteroatom_id).all()
    for record in records:
        recordObject = {'name': record.name,
                        'position': record.position,
                        'team_name': record.team.name,
                        'team_city': record.team.city}
        print(recordObject)   