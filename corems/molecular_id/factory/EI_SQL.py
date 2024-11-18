__author__ = "Yuri E. Corilo"
__date__ = "Feb 12, 2020"

import os
from dataclasses import dataclass

from numpy import array, frombuffer
from sqlalchemy import (
    Column,
    Float,
    ForeignKey,
    Integer,
    LargeBinary,
    String,
    create_engine,
)
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.pool import QueuePool

Base = declarative_base()


class Metadatar(Base):
    """This class is used to store the metadata of the compounds in the database

    Attributes
    -----------
    id : int
        The id of the compound.
    cas : str
        The CAS number of the compound.
    inchikey : str
        The InChiKey of the compound.
    inchi : str
        The InChi of the compound.
    chebi : str
        The ChEBI ID of the compound.
    smiles : str
        The SMILES of the compound.
    kegg : str
        The KEGG ID of the compound.
    iupac_name : str
        The IUPAC name of the compound.
    traditional_name : str
        The traditional name of the compound.
    common_name : str
        The common name of the compound.
    data_id : int
        The id of the compound in the molecularData table.
    data : LowResolutionEICompound
        The compound object.
    """

    __tablename__ = "metaDataR"

    id = Column(Integer, primary_key=True)
    cas = Column(String, nullable=True)
    inchikey = Column(String, nullable=False)
    inchi = Column(String, nullable=False)
    chebi = Column(String, nullable=True)
    smiles = Column(String, nullable=True)
    kegg = Column(String, nullable=True)
    iupac_name = Column(String, nullable=True)
    traditional_name = Column(String, nullable=True)
    common_name = Column(String, nullable=True)

    data_id = Column(Integer, ForeignKey("molecularData.id"))
    data = relationship("LowResolutionEICompound", back_populates="metadatar")


class LowResolutionEICompound(Base):
    """This class is used to store the molecular and spectral data of the compounds in the low res EI database

    Attributes
    -----------
    id : int
        The id of the compound.
    name : str
        The name of the compound.
    classify : str
        The classification of the compound.
    formula : str
        The formula of the compound.
    ri : float
        The retention index of the compound.
    retention_time : float
        The retention time of the compound.
    source : str
        The source of the compound.
    casno : str
        The CAS number of the compound.
    comment : str
        The comment of the compound.
    source_temp_c : float
        The source temperature of the spectra.
    ev : float
        The electron volts of the spectra.
    peaks_count : int
        The number of peaks in the spectra.
    mz : numpy.ndarray
        The m/z values of the spectra.
    abundance : numpy.ndarray
        The abundance values of the spectra.
    metadatar : Metadatar
        The metadata object.
    """

    __tablename__ = "molecularData"

    id = Column(Integer, primary_key=True)

    name = Column(String, nullable=False)
    classify = Column(String, nullable=True)
    formula = Column(String, nullable=True)
    ri = Column(Float, nullable=False)
    retention_time = Column(Float, nullable=False)

    source = Column(String, nullable=True)
    casno = Column(String, nullable=False)
    comment = Column(String, nullable=True)

    derivativenum = Column(String, nullable=True)
    derivatization = Column(String, nullable=True)

    source_temp_c = Column(Float, nullable=True)
    ev = Column(Float, nullable=True)

    peaks_count = Column(Integer, nullable=False)

    mz = Column(LargeBinary, nullable=False)
    abundance = Column(LargeBinary, nullable=False)

    metadatar = relationship("Metadatar", uselist=False, back_populates="data")

    # metadatar = relationship('Metadatar', backref='smile', lazy='dynamic')

    def __init__(self, **dict_data):
        self.id = dict_data.get("id")

        self.name = dict_data.get("NAME")
        self.classify = dict_data.get("classify")
        self.formula = dict_data.get("FORM")
        self.ri = dict_data.get("RI")
        self.retention_time = dict_data.get("RT")

        self.source = dict_data.get("SOURCE")
        self.casno = dict_data.get("CASNO")
        self.comment = dict_data.get("COMMENT")

        self.derivativenum = dict_data.get("derivativenum")
        self.derivatization = dict_data.get("derivatization")

        self.peaks_count = dict_data.get("NUM PEAKS")

        # mz and abun are numpy arrays of 64 bits integer
        # when using postgres array might be a better option

        self.mz = array(dict_data.get("mz"), dtype="int32").tobytes()
        self.abundance = array(dict_data.get("abundance"), dtype="int32").tobytes()

        self.metadatar = dict_data.get("metadatar", None)

    def __repr__(self):
        return (
            "<LowResolutionEICompound(name= %s , cas number = %s, formula = %s, Retention index= %.1f, Retention time= %.1f comment='%s')>"
            % (
                self.name,
                self.casno,
                self.formula,
                self.ri,
                self.retention_time,
                self.comment,
            )
        )


@dataclass
class MetaboliteMetadata:
    """Dataclass for the Metabolite Metadata

    Attributes
    -----------
    id : int
        The id of the compound.
    cas : str
        The CAS number of the compound.
    inchikey : str
        The InChiKey of the compound.
    inchi : str
        The InChi of the compound.
    chebi : str
        The ChEBI ID of the compound.
    smiles : str
        The SMILES of the compound.
    kegg : str
        The KEGG ID of the compound.
    iupac_name : str
        The IUPAC name of the compound.
    traditional_name : str
        The traditional name of the compound.
    common_name : str
        The common name of the compound.
    data_id : int
        The id of the compound in the molecularData table.

    """

    id: int
    cas: str
    inchikey: str
    inchi: str
    chebi: str
    smiles: str
    kegg: str
    data_id: int
    iupac_name: str
    traditional_name: str
    common_name: str


@dataclass
class LowResCompoundRef:
    """Dataclass for the Low Resolution Compound Reference

    This class is used to store the molecular and spectral data of the compounds in the low res EI database

    Parameters
    -----------
    compounds_dict : dict
        A dictionary representing the compound.

    Attributes
    -----------
    id : int
        The id of the compound.
    name : str
        The name of the compound.
    ri : str
        The retention index of the compound.
    retention_time : str
        The retention time of the compound.
    casno : str
        The CAS number of the compound.
    comment : str
        The comment of the compound.
    peaks_count : int
        The number of peaks in the spectra.
    classify : str
        The classification of the compound.
    derivativenum : str
        The derivative number of the compound.
    derivatization : str
        The derivatization applied to the compound.
    mz : numpy.ndarray
        The m/z values of the spectra.
    abundance : numpy.ndarray
        The abundance values of the spectra.
    source_temp_c : float
        The source temperature of the spectra.
    ev : float
        The electron volts of the spectra.
    formula : str
        The formula of the compound.
    source : str
        The source of the spectra data.
    classify : str
        The classification of the compound.
    metadata : MetaboliteMetadata
        The metadata object.
    similarity_score : float
        The similarity score of the compound.
    ri_score : float
        The RI score of the compound.
    spectral_similarity_score : float
        The spectral similarity score of the compound.
    spectral_similarity_scores : dict
        The spectral similarity scores of the compound.

    """

    # this class is use to store the results inside the GCPeak class
    def __init__(self, compounds_dict):
        self.id = compounds_dict.get("id")
        self.name = compounds_dict.get("name")
        self.ri = compounds_dict.get("ri")
        self.retention_time = compounds_dict.get("rt")
        self.casno = compounds_dict.get("casno")
        self.comment = compounds_dict.get("comment")
        self.peaks_count = compounds_dict.get("peaks_count")

        self.classify = compounds_dict.get("classify")
        self.derivativenum = compounds_dict.get("derivativenum")
        self.derivatization = compounds_dict.get("derivatization")

        self.mz = compounds_dict.get("mz")
        self.abundance = compounds_dict.get("abundance")

        self.source_temp_c = compounds_dict.get("source_temp_c")
        self.ev = compounds_dict.get("ev")
        self.formula = compounds_dict.get("formula")
        self.source = compounds_dict.get("source")

        self.classify = compounds_dict.get("classify")

        if compounds_dict.get("metadata"):
            self.metadata = MetaboliteMetadata(**compounds_dict.get("metadata"))

        else:
            self.metadata = None

        self.similarity_score = None
        self.ri_score = None
        self.spectral_similarity_score = None
        self.spectral_similarity_scores = {}


class EI_LowRes_SQLite:
    """
    A class for interacting with a SQLite database for low-resolution EI compounds.

    Parameters
    -----------
    url : str, optional
        The URL of the SQLite database. Default is 'sqlite://'.

    Attributes
    -----------
    engine : sqlalchemy.engine.Engine
        The SQLAlchemy engine for connecting to the database.
    session : sqlalchemy.orm.Session
        The SQLAlchemy session for executing database operations.

    Methods
    --------
    * __init__(self, url='sqlite://').
        Initializes the EI_LowRes_SQLite object.
    * __exit__(self, exc_type, exc_val, exc_tb).
        Closes the database connection.
    * init_engine(self, url).
        Initializes the SQLAlchemy engine.
    * __enter__(self).
        Returns the EI_LowRes_SQLite object.
    * add_compound_list(self, data_dict_list).
        Adds a list of compounds to the database.
    * add_compound(self, data_dict).
        Adds a single compound to the database.
    * commit(self).
        Commits the changes to the database.
    * row_to_dict(self, row).
        Converts a database row to a dictionary.
    * get_all(self).
        Retrieves all compounds from the database.
    * query_min_max_rt(self, min_max_rt).
        Queries compounds based on retention time range.
    * query_min_max_ri(self, min_max_ri).
        Queries compounds based on RI range.
    * query_names_and_rt(self, min_max_rt, compound_names).
        Queries compounds based on compound names and retention time range.
    * query_min_max_ri_and_rt(self, min_max_ri, min_max_rt).
        Queries compounds based on RI range and retention time range.
    * delete_compound(self, compound).
        Deletes a compound from the database.
    * purge(self).
        Deletes all compounds from the database table.
    * clear_data(self).
        Clears all tables in the database.
    """

    def __init__(self, url="sqlite://"):
        self.engine = self.init_engine(url)

        Base.metadata.create_all(self.engine)

        Session = sessionmaker(bind=self.engine)

        self.session = Session()

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Closes the database connection."""
        self.commit()
        self.session.close()
        self.engine.dispose()

    def init_engine(self, url):
        """Initializes the SQLAlchemy engine.

        Parameters
        -----------
        url : str
            The URL of the SQLite database.

        Returns
        --------
        sqlalchemy.engine.Engine
            The SQLAlchemy engine for connecting to the database.
        """
        directory = os.getcwd()
        if not url:
            if not os.path.isdir(directory + "/db"):
                os.mkdir(directory + "/db")
            url = "sqlite:///{DB}/db/pnnl_lowres_gcms_compounds.sqlite".format(
                DB=directory
            )
        return create_engine(url, poolclass=QueuePool)

    def __enter__(self):
        """Returns the EI_LowRes_SQLite object."""
        return self

    def add_compound_list(self, data_dict_list):
        """Adds a list of compounds to the database.

        Parameters
        -----------
        data_dict_list : list of dict
            A list of dictionaries representing the compounds.
        """
        for data_dict in data_dict_list:
            # print(data_dict.get('NUM PEAKS'))
            if not data_dict.get("NUM PEAKS"):
                data_dict["NUM PEAKS"] = len(data_dict.get("mz"))
            if not data_dict.get("CASNO"):
                data_dict["CASNO"] = data_dict.get("CAS")

        self.session.add_all(
            [LowResolutionEICompound(**data_dict) for data_dict in data_dict_list]
        )

    def add_compound(self, data_dict):
        """Adds a single compound to the database.

        Parameters
        -----------
        data_dict : dict
            A dictionary representing the compound.

        """
        one_compound = LowResolutionEICompound(**data_dict)
        self.session.add(one_compound)
        self.commit()

    def commit(self):
        """Commits the changes to the database."""
        try:
            self.session.commit()
        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))

    def row_to_dict(self, row):
        """Converts a database row to a dictionary.

        Parameters
        -----------
        row : sqlalchemy.engine.row.Row
            A row from the database.

        Returns
        --------
        dict
            A dictionary representing the compound.
        """
        data_dict = {c.name: getattr(row, c.name) for c in row.__table__.columns}

        data_dict["mz"] = frombuffer(data_dict.get("mz"), dtype="int32")
        data_dict["abundance"] = frombuffer(data_dict.get("abundance"), dtype="int32")

        if row.metadatar:
            data_dict["metadata"] = {
                c.name: getattr(row.metadatar, c.name)
                for c in row.metadatar.__table__.columns
            }

        else:
            data_dict["metadata"] = None

        return data_dict

    def get_all(
        self,
    ):
        """Retrieves all compounds from the database.

        Returns
        --------
        list
            A list of dictionaries representing the compounds.
        """
        compounds = self.session.query(LowResolutionEICompound).all()

        return [self.row_to_dict(compound) for compound in compounds]

    def query_min_max_rt(
        self,
        min_max_rt,
    ):
        """Queries compounds based on retention time range.

        Parameters
        -----------
        min_max_rt : tuple
            A tuple containing the minimum and maximum retention time values.

        Returns
        --------
        list
            A list of dictionaries representing the compounds.
        """
        min_rt, max_rt = min_max_rt

        compounds = self.session.query(LowResolutionEICompound).filter(
            LowResolutionEICompound.retention_time.between(min_rt, max_rt)
        )

        return [self.row_to_dict(compound) for compound in compounds]

    def query_min_max_ri(self, min_max_ri):
        """Queries compounds based on RI range.

        Parameters
        -----------
        min_max_ri : tuple
            A tuple containing the minimum and maximum RI values.
        """
        min_ri, max_ri = min_max_ri

        compounds = (
            self.session.query(LowResolutionEICompound)
            .filter(LowResolutionEICompound.ri.between(min_ri, max_ri))
            .all()
        )

        return [self.row_to_dict(compound) for compound in compounds]

    def query_names_and_rt(self, min_max_rt, compound_names):
        """Queries compounds based on compound names and retention time range.

        Parameters
        -----------
        min_max_rt : tuple
            A tuple containing the minimum and maximum retention time values.
        compound_names : list
            A list of compound names.

        Returns
        --------
        list
            A list of dictionaries representing the compounds.

        """
        min_rt, max_rt = min_max_rt

        compounds = (
            self.session.query(LowResolutionEICompound)
            .filter(LowResolutionEICompound.name.in_(compound_names))
            .filter(
                LowResolutionEICompound.retention_time >= min_rt,
                LowResolutionEICompound.retention_time <= max_rt,
            )
        )

        # self.session.query.select(LowResolutionEICompound).where(between(LowResolutionEICompound.ri, min_ri, max_ri))
        # x = [self.row_to_dict(compound) for compound in compounds]

        return [self.row_to_dict(compound) for compound in compounds]

    def query_min_max_ri_and_rt(
        self,
        min_max_ri,
        min_max_rt,
    ):
        """Queries compounds based on RI range and retention time range.

        Parameters
        -----------
        min_max_ri : tuple
            A tuple containing the minimum and maximum RI values.
        min_max_rt : tuple
            A tuple containing the minimum and maximum retention time values.

        Returns
        --------
        list
            A list of dictionaries representing the compounds.

        """
        min_ri, max_ri = min_max_ri

        min_rt, max_rt = min_max_rt

        compounds = self.session.query(LowResolutionEICompound).filter(
            LowResolutionEICompound.ri <= max_ri,
            LowResolutionEICompound.ri >= min_ri,
            LowResolutionEICompound.ri >= min_rt,
            LowResolutionEICompound.ri >= max_rt,
        )

        # self.session.query.select(LowResolutionEICompound).where(between(LowResolutionEICompound.ri, min_ri, max_ri))

        return [self.row_to_dict(compound) for compound in compounds]

    def delete_compound(self, compound):
        """Deletes a compound from the database.

        Parameters
        -----------
        compound : LowResolutionEICompound
            A compound object.

        """
        try:
            self.session.delete(compound)
            self.session.commit()

        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))

    def purge(self):
        """Deletes all compounds from the database table.

        Notes
        ------
        Careful, this will delete the entire database table.
        """
        self.session.query(LowResolutionEICompound).delete()
        self.session.commit()

    def clear_data(self):
        """Clears all tables in the database."""
        meta = Base.metadata
        for table in reversed(meta.sorted_tables):
            print("Clear table %s" % table)
            self.session.execute(table.delete())
        self.session.commit()
