import sys
sys.path.append(".")
import os

from sqlalchemy import create_engine, ForeignKey, Column, Integer, String, Float, SMALLINT
from sqlalchemy.orm import backref, column_property, relationship
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.sql.schema import UniqueConstraint
from sqlalchemy import exc

from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.ext.hybrid import hybrid_property, hybrid_method
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.sql.operators import exists
from sqlalchemy import event, and_
from sqlalchemy import func

from corems.encapsulation.constant import Atoms, Labels
import json
from corems.encapsulation.factory.processingSetting import MolecularFormulaSearchSettings
from sqlalchemy.orm.scoping import scoped_session
from corems import chunks
import tqdm

Base = declarative_base()

class HeteroAtoms(Base):

    __tablename__ = 'heteroAtoms'

    id = Column(Integer, primary_key=True,
                unique=True,
                nullable=False)

    name = Column(String, unique=True, nullable=False)

    halogensCount = Column(Integer, unique=False, nullable=False)

    carbonHydrogen = relationship('CarbonHydrogen', secondary='molecularformula')

    def __repr__(self):
        return '<HeteroAtoms Model {} class {}>'.format(self.id, self.name)      

    @hybrid_property
    def halogens_count(cls):
        return cls.halogensCount.cast(Float)

    def to_dict(self):

        return json.loads(self.name)


class CarbonHydrogen(Base):

    __tablename__ = 'carbonHydrogen'
    __table_args__ = (UniqueConstraint('C', 'H', name='unique_c_h'), )

    id = Column(Integer, primary_key=True,
                unique=True,
                nullable=False)

    C = Column(Integer, nullable=False)

    H = Column(Integer, nullable=False)

    heteroAtoms = relationship("HeteroAtoms",
                               secondary="molecularformula",
                               )

    def __repr__(self):
        return '<CarbonHydrogen Model {} C{} H{}>'.format(self.id, self.C, self.H)                     

    @property
    def mass(self):
        return (self.C * Atoms.atomic_masses.get('C')) + (self.H * Atoms.atomic_masses.get('H'))

    @hybrid_property
    def c(cls):
        return cls.C.cast(Float)

    @hybrid_property
    def h(cls):
        return cls.H.cast(Float)

    @hybrid_property
    def dbe(cls):
        # return cls.C.cast(Float) - (cls.H.cast(Float) / 2) + 1
        return float(cls.C) - float(cls.H / 2) + 1

# 264888.88 ms
class MolecularFormulaLink(Base):

    __tablename__ = 'molecularformula'
    __table_args__ = (UniqueConstraint('heteroAtoms_id', 'carbonHydrogen_id', name='unique_molform'), )

    # id = Column(Integer, primary_key=True,
    #                    unique=True,
    #                    nullable=False)

    heteroAtoms_id = Column(Integer,
                            ForeignKey('heteroAtoms.id'),
                            primary_key=True)

    carbonHydrogen_id = Column(Integer,
                               ForeignKey('carbonHydrogen.id'), 
                               primary_key=True)

    mass = Column(Float)

    DBE = Column(Float)

    carbonHydrogen = relationship(CarbonHydrogen, backref=backref("heteroAtoms_assoc"))

    heteroAtoms = relationship(HeteroAtoms, backref=backref("carbonHydrogen_assoc"))

    C = association_proxy('carbonHydrogen', 'C')

    H = association_proxy('carbonHydrogen', 'H')

    classe = association_proxy('heteroAtoms', 'name')

    @property
    def formula_dict(self):

        carbon = {'C': self.C, 'H': self.H}
        classe = json.loads(self.classe)
        if self.classe == '{"HC": ""}':
            return {**carbon}
        else:
            return {**carbon, **classe}

    @property
    def formula_string(self):
        class_dict = self.formula_dict
        class_str = ' '.join([atom + str(class_dict[atom]) for atom in class_dict.keys()])
        return class_str.strip()

    @property
    def classe_string(self):
        class_dict = json.loads(self.classe)
        class_str = ' '.join([atom + str(class_dict[atom]) for atom in class_dict.keys()])
        return class_str.strip()

    @hybrid_method
    def adduct_mass(self, ion_charge, adduct_atom):
        return (self.mass + (Atoms.atomic_masses.get(adduct_atom)) + (ion_charge * -1 * Atoms.electron_mass))/ abs(ion_charge)

    @hybrid_method
    def protonated_mass(self, ion_charge):
        return (self.mass + (ion_charge * Atoms.atomic_masses.get("H")) + (ion_charge * -1 * Atoms.electron_mass))/abs(ion_charge)

    @hybrid_method
    def radical_mass(self, ion_charge):
        return (self.mass + (ion_charge * -1 * Atoms.electron_mass))/ abs(ion_charge)

    def __repr__(self):

        return '<MolecularFormulaLink Model {}>'.format(self.formula_string)       


class MolForm_SQL:

    def __init__(self, url=None, echo=False):

        self.engine = self.init_engine(url)

        self.add_engine_pidguard(self.engine)

        session_factory = sessionmaker(bind=self.engine)

        Session = scoped_session(session_factory)

        self.session = session_factory()

        Base.metadata.create_all(self.engine)

        self.session.commit()

    def __exit__(self, exc_type, exc_val, exc_tb):
        # make sure the dbconnection gets closed

        self.commit()
        self.session.close()
        self.engine.dispose()

    def initiate_database(self, url, database_name):  #CREATION

        engine = create_engine(url)
        conn = engine.connect()
        conn.execute("commit")
        conn.execute("create database " + database_name)
        conn.close()

    def commit(self):

        try:
            self.session.commit()  
        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))

    def init_engine(self, url):

        if not url or url == 'None' or url == 'False':
            directory = os.getcwd()

            if not os.path.isdir(directory + '/db'):
                os.mkdir(directory + '/db')

            url = 'sqlite:///{DB}/db/molformulas.sqlite'.format(DB=directory)

        if url[0:6] == 'sqlite':
            self.type = 'sqlite'
        else:
            self.type = 'normal'
            
        if url[0:6] == 'sqlite':
            engine = create_engine(url, echo = False)
            self.chunks_count = 50
        
        elif url[0:10] == 'postgresql' or url[0:8] == 'postgres':
            #postgresql
            self.chunks_count = 50000
            engine = create_engine(url, echo = False, isolation_level="AUTOCOMMIT")
        
        return engine# poolclass=NullPool

    def __enter__(self):
        
        return self
    
    def get_dict_by_classes(self, classes, ion_type, nominal_mzs, ion_charge, molecular_search_settings, adducts=None):
                                    
        '''Known issue, when using SQLite:
         if the number of classes and nominal_m/zs are higher than 999 the query will fail
         Solution: use postgres or split query''' 
        def query_normal(class_list, len_adduct):
            
            base_query = self.session.query(MolecularFormulaLink, CarbonHydrogen, HeteroAtoms)\
                                .filter(MolecularFormulaLink.carbonHydrogen_id == CarbonHydrogen.id)\
                                .filter(MolecularFormulaLink.heteroAtoms_id == HeteroAtoms.id)
            
            return base_query.filter(
                and_(
                    HeteroAtoms.name.in_(class_list), 
                    and_(
                        MolecularFormulaLink.DBE >= molecular_search_settings.min_dbe, 
                        MolecularFormulaLink.DBE <= molecular_search_settings.max_dbe, 
                        and_(
                            ((CarbonHydrogen.h + HeteroAtoms.halogens_count - len_adduct) / CarbonHydrogen.c) >= molecular_search_settings.min_hc_filter,
                            ((CarbonHydrogen.h + HeteroAtoms.halogens_count - len_adduct) / CarbonHydrogen.c) <= molecular_search_settings.max_hc_filter,
                            CarbonHydrogen.C >= molecular_search_settings.usedAtoms.get("C")[0],
                            CarbonHydrogen.c <= molecular_search_settings.usedAtoms.get("C")[1],
                            CarbonHydrogen.h >= molecular_search_settings.usedAtoms.get("H")[0],
                            CarbonHydrogen.h <= molecular_search_settings.usedAtoms.get("H")[1],
                        )
                    )
                )
            )

        def add_dict_formula(formulas, ion_type, ion_charge, adduct_atom=None):
            "organize data by heteroatom classes"
            dict_res = {}

            def nominal_mass_by_ion_type(formula_obj):
                
                if ion_type == Labels.protonated_de_ion:
                
                    return int(formula_obj.protonated_mass(ion_charge))
                
                elif ion_type == Labels.radical_ion:
                    
                    return int(formula_obj.radical_mass(ion_charge))

                elif ion_type == Labels.adduct_ion and adduct_atom:
                    
                    return int(formula_obj.adduct_mass(ion_charge, adduct_atom))
            
            for formula_obj, ch_obj, classe_obj in tqdm.tqdm(formulas, desc="Loading molecular formula database"):
                
                nominal_mz = nominal_mass_by_ion_type(formula_obj)
                
                if self.type != 'normal':
                    if not nominal_mz in nominal_mzs:
                        continue
                classe = classe_obj.name

                # classe_str = formula.classe_string
                
                # pbar.set_description_str(desc="Loading molecular formula database for class %s " % classe_str)
                
                formula_dict = formula_obj.formula_dict

                if formula_dict.get("O"):
                    
                    if formula_dict.get("O") / formula_dict.get("C") >= molecular_search_settings.max_oc_filter:
                        # print(formula_dict.get("O") / formula_dict.get("C"), molecular_search_settings.max_oc_filter)
                        continue
                    elif formula_dict.get("O") / formula_dict.get("C") <= molecular_search_settings.min_oc_filter:
                        # print(formula_dict.get("O") / formula_dict.get("C"), molecular_search_settings.min_oc_filter)
                        continue
                    #if formula_dict.get("P"):

                    #    if  not (formula_dict.get("O") -2)/ formula_dict.get("P") >= molecular_search_settings.min_op_filter:
                            
                    #        continue
        
                if classe in dict_res.keys():
                    
                    if nominal_mz in dict_res[classe].keys():
                        
                        dict_res.get(classe).get(nominal_mz).append(formula_obj)
                    
                    else:

                        dict_res.get(classe)[nominal_mz] = [formula_obj ]  
            
                else:
                    
                    dict_res[classe] = {nominal_mz: [formula_obj] }     
            
            return dict_res
        
        
        len_adducts = 0
        if ion_type == Labels.adduct_ion:
            len_adducts = 1
        
        query = query_normal(classes, len_adducts)

        if ion_type == Labels.protonated_de_ion:
            if self.type == 'normal':
                query = query.filter(and_(
                                MolecularFormulaLink.protonated_mass(ion_charge).cast(Integer).in_(nominal_mzs)
                                ))
            return add_dict_formula(query, ion_type, ion_charge)
        
        if ion_type == Labels.radical_ion:
            if self.type == 'normal':
                query = query.filter(MolecularFormulaLink.radical_mass(ion_charge).cast(Integer).in_(nominal_mzs))    
            return add_dict_formula(query, ion_type, ion_charge)
        
        if ion_type == Labels.adduct_ion:
            dict_res = {}
            if adducts: 
                for atom in adducts:
                    if self.type == 'normal':
                        query = query.filter(MolecularFormulaLink.adduct_mass(ion_charge, atom).cast(Integer).in_(nominal_mzs))    
                    dict_res[atom] = add_dict_formula(query, ion_type, ion_charge, adduct_atom=atom)
                return dict_res
        # dump all objs to memory
        self.session.expunge_all()
        
    def check_entry(self,classe, ion_type, molecular_search_settings):
        #  get all classes, ion_type, ion charge as str add to a dict or list
        #  then check if class in database
        has_class = self.session.query(exists().where(
            (MolecularFormulaLink.classe == classe)))
        
        return has_class
    
    def get_all_classes(self):
        
        query = self.session.query(MolecularFormulaLink.classe.distinct().label("classe"))
        
        return query.all()  
    
    def get_all(self,):
        
        mol_formulas = self.session.query(MolecularFormulaLink).all()
        
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

    def close(self, commit=True):
        # make sure the dbconnection gets closed
        
        if commit: self.commit()
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
    
    sql = MolForm_SQL(url='sqlite:///')
    
    dict_data = {"name": '{"O": 12}'}
    dict_data2 = {"name": '{"O": 13}'}
    hetero_obj = HeteroAtoms(**dict_data)
    hetero_obj2 = HeteroAtoms(**dict_data2)
    sql.session.add(hetero_obj)
    sql.session.add(hetero_obj2)

    print(sql.session.query(HeteroAtoms).all())
    #molecular_search_settings = MolecularFormulaSearchSettings()
    #sql = MolForm_SQL()
    #query = sql.session.query(MolecularFormulaLink).filter_by(classe = '{"O": 12}').filter(MolecularFormulaLink.adduct_mass(+2, "Na") < 250)
    #query = sql.get_by_classe('{"O": 12}', molecular_search_settings).filter(MolecularFormulaLink.adduct_mass(+2, "Na") < 250)
    #classes = ['{"O": 12}']*1
    #for i, classe in enumerate(classes):
        #query = sql.get_by_classe(classe, molecular_search_settings)
        #query = sql.session.query(MolecularFormulaLink).filter_by(classe = '{"O": 12}').filter(MolecularFormulaLink.adduct_mass(+2, "Na") < 250)
        #for i in query.filter(MolecularFormulaLink.mass < 250):
            
        #    print(i.radical_mass(-1), i.protonated_mass(-1), i.adduct_mass(+2, "Na"), i.mass, i.formula_dict, i.formula_string)
    #
 