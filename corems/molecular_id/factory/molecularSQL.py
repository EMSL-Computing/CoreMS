import sys
sys.path.append(".")
import os

from sqlalchemy import Numeric, create_engine, ForeignKey, Column, Integer, String, Float, func
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
    """ HeteroAtoms class for the heteroAtoms table in the SQLite database.
    
    Attributes
    ----------
    id : int
        The primary key for the table.
    name : str
        The name of the heteroAtoms class.
    halogensCount : int
        The number of halogens in the heteroAtoms class.
    carbonHydrogen : relationship
        The relationship to the carbonHydrogen table.
    
    Methods
    -------
    * __repr__()
        Returns the string representation of the object.
    * to_dict()
        Returns the heteroAtoms class as a dictionary.
    * halogens_count()
        Returns the number of halogens as a float.
    

    """
    __tablename__ = 'heteroAtoms'

    id = Column(Integer, primary_key=True,
                         unique = True,
                         nullable = False)

    name = Column(String, unique=True, nullable=False)

    halogensCount = Column(Integer, unique=False, nullable=False)

    carbonHydrogen = relationship('CarbonHydrogen', secondary='molecularformula',  viewonly=True)

    def __repr__(self):
        return '<HeteroAtoms Model {} class {}>'.format(self.id, self.name)      

    @hybrid_property
    def halogens_count(cls):
        """ Returns the number of halogens as a float."""
        return cls.halogensCount.cast(Float)

    def to_dict(self):
        """ Returns the heteroAtoms class as a dictionary."""
        return json.loads(self.name)


class CarbonHydrogen(Base):
    """ CarbonHydrogen class for the carbonHydrogen table in the SQLite database.
    
    Attributes
    ----------
    id : int
        The primary key for the table.
    C : int
        The number of carbon atoms.
    H : int
        The number of hydrogen atoms.
    heteroAtoms : relationship
        The relationship to the heteroAtoms table.
    
    Methods
    -------
    * __repr__()
        Returns the string representation of the object.
    * mass()
        Returns the mass of the carbonHydrogen class as a float.
    * c()
        Returns the number of carbon atoms as a float.
    * h()
        Returns the number of hydrogen atoms as a float.
    * dbe()
        Returns the double bond equivalent as a float.
    
    """

    __tablename__ = 'carbonHydrogen'
    __table_args__ = (UniqueConstraint('C', 'H', name='unique_c_h'), )

    id = Column(Integer, primary_key=True,
                unique=True,
                nullable=False)

    C = Column(Integer, nullable=False)

    H = Column(Integer, nullable=False)

    heteroAtoms = relationship("HeteroAtoms",
                               secondary="molecularformula",
                               viewonly=True
                               )

    def __repr__(self):
        """ Returns the string representation of the object."""
        return '<CarbonHydrogen Model {} C{} H{}>'.format(self.id, self.C, self.H)                     

    @property
    def mass(self):
        """ Returns the mass of the carbonHydrogen class as a float."""
        return (self.C * Atoms.atomic_masses.get('C')) + (self.H * Atoms.atomic_masses.get('H'))

    @hybrid_property
    def c(cls):
        """ Returns the number of carbon atoms as a float."""
        return cls.C.cast(Float)

    @hybrid_property
    def h(cls):
        """ Returns the number of hydrogen atoms as a float."""
        return cls.H.cast(Float)

    @hybrid_property
    def dbe(cls):
        """ Returns the double bond equivalent as a float."""
        # return cls.C.cast(Float) - (cls.H.cast(Float) / 2) + 1
        return float(cls.C) - float(cls.H / 2) + 1

# 264888.88 ms
class MolecularFormulaLink(Base):
    """ MolecularFormulaLink class for the molecularformula table in the SQLite database.

    Attributes
    ----------
    heteroAtoms_id : int
        The foreign key for the heteroAtoms table.
    carbonHydrogen_id : int
        The foreign key for the carbonHydrogen table.
    mass : float
        The mass of the molecular formula.
    DBE : float
        The double bond equivalent of the molecular formula.
    carbonHydrogen : relationship
        The relationship to the carbonHydrogen table.
    heteroAtoms : relationship
        The relationship to the heteroAtoms table.
    C : association_proxy
        The association proxy for the carbonHydrogen table.
    H : association_proxy
        The association proxy for the carbonHydrogen table.
    classe : association_proxy
        The association proxy for the heteroAtoms table.
    
    Methods
    -------
    * __repr__()
        Returns the string representation of the object.
    * to_dict()
        Returns the molecular formula as a dictionary.
    * formula_string()
        Returns the molecular formula as a string.
    * classe_string()
        Returns the heteroAtoms class as a string.
    * _adduct_mz(ion_charge, adduct_atom)
        Returns the m/z of the adduct ion as a float.
    * _protonated_mz(ion_charge)
        Returns the m/z of the protonated ion as a float.
    * _radical_mz(ion_charge)
        Returns the m/z of the radical ion as a float.
    


    """
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

    def to_dict(self):
        """ Returns the molecular formula as a dictionary.
        
        Returns
        -------
        dict
            The molecular formula as a dictionary.  
        """
        carbon = {'C': self.C, 'H': self.H}
        classe = json.loads(self.classe)
        if self.classe == '{"HC": ""}':
            return {**carbon}
        else:
            return {**carbon, **classe}

    @property
    def formula_string(self):
        """ Returns the molecular formula as a string."""
        class_dict = self.to_dict()
        class_str = ' '.join([atom + str(class_dict[atom]) for atom in class_dict.keys()])
        return class_str.strip()

    @property
    def classe_string(self):
        """ Returns the heteroAtoms class as a string."""
        class_dict = json.loads(self.classe)
        class_str = ' '.join([atom + str(class_dict[atom]) for atom in class_dict.keys()])
        return class_str.strip()

    @hybrid_method
    def _adduct_mz(self, ion_charge, adduct_atom):
        """ Returns the m/z of the adduct ion as a float."""
        return (self.mass + (Atoms.atomic_masses.get(adduct_atom)) + (ion_charge * -1 * Atoms.electron_mass))/ abs(ion_charge)

    @hybrid_method
    def _protonated_mz(self, ion_charge):
        """ Returns the m/z of the protonated ion as a float."""
        return (self.mass + (ion_charge * Atoms.atomic_masses.get("H")) + (ion_charge * -1 * Atoms.electron_mass))/abs(ion_charge)

    @hybrid_method
    def _radical_mz(self, ion_charge):
        """ Returns the m/z of the radical ion as a float."""
        return (self.mass + (ion_charge * -1 * Atoms.electron_mass))/ abs(ion_charge)

    def __repr__(self):
        """ Returns the string representation of the object."""
        return '<MolecularFormulaLink Model {}>'.format(self.formula_string)       


class MolForm_SQL:
    """ MolForm_SQL class for the SQLite database.
    
    Attributes
    ----------
    engine : sqlalchemy.engine.base.Engine
        The SQLAlchemy engine.
    session : sqlalchemy.orm.session.Session
        The SQLAlchemy session.
    type : str
        The type of database.
    chunks_count : int
        The number of chunks to use when querying the database.
    
    Methods
    -------
    * __init__(url=None, echo=False)
        Initializes the database.
    * __exit__(exc_type, exc_val, exc_tb)
        Closes the database.
    * initiate_database(url, database_name)
        Creates the database.
    * commit()
        Commits the session.
    * init_engine(url)
        Initializes the SQLAlchemy engine.
    * __enter__()

    * get_dict_by_classes(classes, ion_type, nominal_mzs, ion_charge, molecular_search_settings, adducts=None)
        Returns a dictionary of molecular formulas.
    * check_entry(classe, ion_type, molecular_search_settings)
        Checks if a molecular formula is in the database.
    * get_all_classes()
        Returns a list of all classes in the database.
    * get_all()
        Returns a list of all molecular formulas in the database.
    * delete_entry(row)
        Deletes a molecular formula from the database.
    * purge(cls)
        Deletes all molecular formulas from the database.
    * clear_data()
        Clears the database.
    * close(commit=True)
        Closes the database.
    * add_engine_pidguard(engine)
        Adds multiprocessing guards.
    
    """
    def __init__(self, url=None, echo=False):

        self.engine = self.init_engine(url)

        self.add_engine_pidguard(self.engine)

        session_factory = sessionmaker(bind=self.engine)

        Session = scoped_session(session_factory)

        self.session = session_factory()

        Base.metadata.create_all(self.engine)

        self.session.commit()

    def __exit__(self, exc_type, exc_val, exc_tb):
        """ Closes the database.
        
        Parameters
        ----------
        exc_type : str
            The exception type.
        exc_val : str
            The exception value.
        exc_tb : str
            The exception traceback.
        """
        # make sure the dbconnection gets closed

        self.commit()
        self.session.close()
        self.engine.dispose()

    def initiate_database(self, url, database_name):  #CREATION
        """ Creates the database.
        
        Parameters
        ----------
        url : str
            The URL for the database.
        database_name : str
            The name of the database.
        """
        engine = create_engine(url)
        conn = engine.connect()
        conn.execute("commit")
        conn.execute("create database " + database_name)
        conn.close()

    def commit(self):
        """ Commits the session.
        """
        try:
            self.session.commit()  
        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))

    def init_engine(self, url):
        """ Initializes the SQLAlchemy engine.
        
        Parameters
        ----------
        url : str
            The URL for the database.
        
        Returns
        -------
        sqlalchemy.engine.base.Engine
            The SQLAlchemy engine.
        
        """
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
        """ Returns the object.
        """
        return self
    
    def get_dict_by_classes(self, classes, ion_type, nominal_mzs, ion_charge, molecular_search_settings, adducts=None):
        """ Returns a dictionary of molecular formulas.
        
        Parameters
        ----------
        classes : list
            The list of classes.
        ion_type : str
            The ion type.
        nominal_mzs : list
            The list of nominal m/z values.
        ion_charge : int
            The ion charge.
        molecular_search_settings : MolecularFormulaSearchSettings
            The molecular formula search settings.
        adducts : list, optional
            The list of adducts. Default is None.
        
        Returns
        -------
        dict
            The dictionary of molecular formulas.
        
        Notes
        -----
        Known issue, when using SQLite:
        if the number of classes and nominal_m/zs are higher than 999 the query will fail
        Solution: use postgres or split query
        """                     
         
        def query_normal(class_list, len_adduct):
            """ query for normal database
            
            Parameters
            ----------
            class_list : list
                The list of classes.
            len_adduct : int
                The length of the adduct.
            
            Returns
            -------
            sqlalchemy.orm.query.Query
                The query.
            """
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
            """ add molecular formula to dict
            
            Parameters
            ----------
            formulas : sqlalchemy.orm.query.Query
                The query.
            ion_type : str
                The ion type.
            ion_charge : int
                The ion charge.
            adduct_atom : str, optional
                The adduct atom. Default is None.
            
            Returns
            -------
            dict
                The dictionary of molecular formulas.
            
            """
            "organize data by heteroatom classes"
            dict_res = {}

            def nominal_mass_by_ion_type(formula_obj):
                
                if ion_type == Labels.protonated_de_ion:
                
                    return int(formula_obj._protonated_mz(ion_charge))
                
                elif ion_type == Labels.radical_ion:
                    
                    return int(formula_obj._radical_mz(ion_charge))

                elif ion_type == Labels.adduct_ion and adduct_atom:
                    
                    return int(formula_obj._adduct_mz(ion_charge, adduct_atom))
            
            for formula_obj, ch_obj, classe_obj in tqdm.tqdm(formulas, desc="Loading molecular formula database"):
                
                nominal_mz = nominal_mass_by_ion_type(formula_obj)
                
                if self.type != 'normal':
                    if not nominal_mz in nominal_mzs:
                        continue
                classe = classe_obj.name

                # classe_str = formula.classe_string
                
                # pbar.set_description_str(desc="Loading molecular formula database for class %s " % classe_str)
                
                formula_dict = formula_obj.to_dict()

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
                
                query = query.filter(
                                func.floor(MolecularFormulaLink._protonated_mz(ion_charge)).in_(nominal_mzs)
                                )
                                
                                
            return add_dict_formula(query, ion_type, ion_charge)
        
        if ion_type == Labels.radical_ion:
            if self.type == 'normal':
                query = query.filter(func.floor(MolecularFormulaLink._radical_mz(ion_charge)).in_(nominal_mzs))
            return add_dict_formula(query, ion_type, ion_charge)
        
        if ion_type == Labels.adduct_ion:
            dict_res = {}
            if adducts: 
                for atom in adducts:
                    if self.type == 'normal':
                        query = query.filter(func.floor(MolecularFormulaLink._adduct_mz(ion_charge, atom)).in_(nominal_mzs))    
                    dict_res[atom] = add_dict_formula(query, ion_type, ion_charge, adduct_atom=atom)
                return dict_res
        # dump all objs to memory
        self.session.expunge_all()
        
    def check_entry(self,classe, ion_type, molecular_search_settings):
        """ Checks if a molecular formula is in the database.

        Parameters
        ----------
        classe : str
            The class of the molecular formula.
        ion_type : str
            The ion type.
        molecular_search_settings : MolecularFormulaSearchSettings
            The molecular formula search settings.
        
        Returns
        -------
        sqlalchemy.orm.query.Query
            The query.
        """
        #  get all classes, ion_type, ion charge as str add to a dict or list
        #  then check if class in database
        has_class = self.session.query(exists().where(
            (MolecularFormulaLink.classe == classe)))
        
        return has_class
    
    def get_all_classes(self):
        """ Returns a list of all classes in the database."""
        query = self.session.query(MolecularFormulaLink.classe.distinct().label("classe"))
        
        return query.all()  
    
    def get_all(self,):
        """ Returns a list of all molecular formulas in the database."""
        mol_formulas = self.session.query(MolecularFormulaLink).all()
        
        return mol_formulas

    def delete_entry(self, row):
        """ Deletes a molecular formula from the database."""
        try:
            self.session.delete(row)  
            self.session.commit()  
        
        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))

    def purge(self, cls):
        """ Deletes all molecular formulas from the database.
        
        Notes 
        -------
        Careful, this will delete the entire database table

        """
        self.session.query(cls).delete()
        self.session.commit()  

    def clear_data(self):
        """ Clears the database.
        """
        meta = Base.metadata
        for table in reversed(meta.sorted_tables):
            print ('Clear table %s' % table)
            self.session.execute(table.delete())
        self.session.commit()

    def close(self, commit=True):
        """ Closes the database.
        
        Parameters
        ----------
        commit : bool, optional
            Whether to commit the session. Default is True.
        """
        # make sure the dbconnection gets closed
        
        if commit: self.commit()
        self.session.close()
        self.engine.dispose()    
   
    def add_engine_pidguard(self, engine):
        """ Adds multiprocessing guards.
        
        Forces a connection to be reconnected if it is detected
        as having been shared to a sub-process.

        Parameters
        ----------
        engine : sqlalchemy.engine.base.Engine
            The SQLAlchemy engine.
        
        """
        import os, warnings
     

        @event.listens_for(engine, "connect")
        def connect(dbapi_connection, connection_record):
            """ Forces a connection to be reconnected if it is detected
            
            Parameters
            ----------
            dbapi_connection : sqlalchemy.engine.base.Engine
                The SQLAlchemy engine.
            connection_record : sqlalchemy.engine.base.Engine
                The SQLAlchemy engine.
            """
            connection_record.info['pid'] = os.getpid()

        @event.listens_for(engine, "checkout")
        def checkout(dbapi_connection, connection_record, connection_proxy):
            """ Forces a connection to be reconnected if it is detected
            
            Parameters
            ----------
            dbapi_connection : sqlalchemy.engine.base.Engine
                The SQLAlchemy engine.
            connection_record : sqlalchemy.engine.base.Engine
                The SQLAlchemy engine.
            connection_proxy : sqlalchemy.engine.base.Engine
                The SQLAlchemy engine.
            
            Raises
            ------
            exc.DisconnectionError
                If the connection record belongs to a different process.
            
            """
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
    #query = sql.session.query(MolecularFormulaLink).filter_by(classe = '{"O": 12}').filter(MolecularFormulaLink._adduct_mz(+2, "Na") < 250)
    #query = sql.get_by_classe('{"O": 12}', molecular_search_settings).filter(MolecularFormulaLink._adduct_mz(+2, "Na") < 250)
    #classes = ['{"O": 12}']*1
    #for i, classe in enumerate(classes):
        #query = sql.get_by_classe(classe, molecular_search_settings)
        #query = sql.session.query(MolecularFormulaLink).filter_by(classe = '{"O": 12}').filter(MolecularFormulaLink._adduct_mz(+2, "Na") < 250)
        #for i in query.filter(MolecularFormulaLink.mass < 250):
            
        #    print(i._radical_mz(-1), i._protonated_mz(-1), i._adduct_mz(+2, "Na"), i.mass, i.to_dict(), i.formula_string)
    #
 