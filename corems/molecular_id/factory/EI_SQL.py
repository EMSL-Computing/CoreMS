

__author__ = "Yuri E. Corilo"
__date__ = "Feb 12, 2020"

import os 
from dataclasses import dataclass

from sqlalchemy import create_engine, Column, Integer, String, Float, Binary, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound
from sqlalchemy.pool import QueuePool
from sqlalchemy import between

from numpy import array, frombuffer

Base = declarative_base()

class Metadatar(Base):
    __tablename__ = 'metaDataR'

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

    data_id = Column(Integer, ForeignKey('molecularData.id'))
    data = relationship("LowResolutionEICompound", back_populates="metadatar")

class LowResolutionEICompound(Base):

    __tablename__ = 'molecularData'

    id = Column(Integer, primary_key=True)

    name = Column(String, nullable=False)
    classify = Column(String, nullable=True)

    formula = Column(String, nullable=True)
    ri = Column(Float, nullable=False)
    rt = Column(Float, nullable=False)

    classify = Column(String, nullable=True)
    derivativenum = Column(String, nullable=True)
    derivatization = Column(String, nullable=True)

    source = Column(String, nullable=True)
    casno = Column(String, nullable=False)
    comment = Column(String, nullable=True)

    source_temp_c = Column(Float, nullable=True)
    ev = Column(Float, nullable=True)

    peaks_count = Column(Integer, nullable=False)
    mz = Column(Binary, nullable=False)
    abundance = Column(Binary, nullable=False)

    metadatar = relationship("Metadatar", uselist=False, back_populates="data")

    # metadatar = relationship('Metadatar', backref='smile', lazy='dynamic')

    def __init__(self, dict_data):

        self.id = dict_data.get('id')

        self.name = dict_data.get('NAME')
        self.formula = dict_data.get('FORM')
        self.ri = dict_data.get('RI')
        self.rt = dict_data.get('RT')

        self.source = dict_data.get('SOURCE')
        self.casno = dict_data.get('CASNO')
        self.comment = dict_data.get('COMMENT')

        self.classify = dict_data.get('classify')
        self.derivativenum = dict_data.get('derivativenum')
        self.derivatization = dict_data.get('derivatization')

        self.peaks_count = dict_data.get('NUM PEAKS')

        # mz and abun are numpy arrays of 64 bits integer
        # when using postgres array might be a better option

        self.mz = array(dict_data.get('mz'), dtype='int32').tobytes()
        self.abundance = array(dict_data.get('abundance'), dtype="int32").tobytes()

    def __repr__(self):
        return "<LowResolutionEICompound(name= %s , cas number = %s, formula = %s, Retention index= %.1f, Retention time= %.1f comment='%i')>" % (
                                    self.name, self.casno, self.formula, self.ri, self.rt, self.comment)
@dataclass
class MetaboliteMetadata:

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
    "this class is use to store the results inside the GCPeak class"
    def __init__(self, compounds_dict):

        self.id = compounds_dict.get("id")
        self.name = compounds_dict.get("name")
        self.ri = compounds_dict.get("ri")
        self.rt = compounds_dict.get("rt")
        self.casno = compounds_dict.get("casno")
        self.comment = compounds_dict.get("comment")
        self.peaks_count = compounds_dict.get("peaks_count")

        self.classify = compounds_dict.get('classify')
        self.derivativenum = compounds_dict.get('derivativenum')
        self.derivatization = compounds_dict.get('derivatization')

        self.mz = compounds_dict.get('mz')
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

    def __init__(self, url='sqlite://'):

        self.engine = self.init_engine(url)

        Base.metadata.create_all(self.engine)

        Session = sessionmaker(bind=self.engine)

        self.session = Session()

    def __exit__(self, exc_type, exc_val, exc_tb):
        # make sure the db connection gets closed
        #
        self.commit()
        self.session.close()
        self.engine.dispose()

    def init_engine(self, url):

        directory = os.getcwd()

        if not url:
            
            if not os.path.isdir(directory + '/db'):
                    
                os.mkdir(directory + '/db')    
            
            url = 'sqlite:///{DB}/db/pnnl_lowres_gcms_compounds.sqlite'.format(DB=directory)

        return create_engine(url, poolclass=QueuePool)

    def __enter__(self):

        return self

    def add_compound_list(self, data_dict_list):

        # print(data_dict_list[0])
        for data_dict in data_dict_list:
            # print(data_dict.get('NUM PEAKS'))
            if not data_dict.get('NUM PEAKS'):
                data_dict['NUM PEAKS'] = len(data_dict.get('mz'))
            if not data_dict.get('CASNO'):
                data_dict['CASNO'] = data_dict.get('CAS')

        self.session.add_all([LowResolutionEICompound(data_dict) for data_dict in data_dict_list] )

    def add_compound(self, data_dict):

        one_compound = LowResolutionEICompound(data_dict)
        self.session.add(one_compound)
        self.commit()

    def commit(self):

        try:
            self.session.commit()
        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))

    def row_to_dict(self, row):

        data_dict = {c.name: getattr(row, c.name) for c in row.__table__.columns}        

        data_dict['mz'] = frombuffer(data_dict.get('mz'), dtype="int32")
        data_dict['abundance'] = frombuffer(data_dict.get('abundance'), dtype="int32")

        if row.metadatar:
            data_dict['metadata'] = {c.name: getattr(row.metadatar, c.name) for c in row.metadatar.__table__.columns}

        else:
            data_dict['metadata'] = None

        return data_dict

    def get_all(self,):

        compounds = self.session.query(LowResolutionEICompound).all()

        return [self.row_to_dict(compound) for compound in compounds]

    def query_min_max_rt(self, min_max_rt, ):

        min_rt, max_rt = min_max_rt

        compounds = self.session.query(LowResolutionEICompound).filter(LowResolutionEICompound.rt.between(min_rt, max_rt))    

        return [self.row_to_dict(compound) for compound in compounds]

    def query_min_max_ri(self, min_max_ri):

        min_ri, max_ri = min_max_ri

        compounds = self.session.query(LowResolutionEICompound).filter(LowResolutionEICompound.ri.between(min_ri, max_ri)).all()

        return [self.row_to_dict(compound) for compound in compounds]

    def query_names_and_rt(self, min_max_rt, compound_names):

        min_rt, max_rt = min_max_rt

        compounds = self.session.query(LowResolutionEICompound).filter(LowResolutionEICompound.name.in_(compound_names)).filter(
                                        LowResolutionEICompound.rt >= min_rt,
                                        LowResolutionEICompound.rt <= max_rt,
                                        )
        
        #self.session.query.select(LowResolutionEICompound).where(between(LowResolutionEICompound.ri, min_ri, max_ri))    
        # x = [self.row_to_dict(compound) for compound in compounds]
        
        return [self.row_to_dict(compound) for compound in compounds]

    def query_min_max_ri_and_rt(self, min_max_ri, min_max_rt, ):
        
        min_ri, max_ri = min_max_ri
        
        min_rt, max_rt = min_max_rt

        compounds = self.session.query(LowResolutionEICompound).filter(
            LowResolutionEICompound.ri <= max_ri,
            LowResolutionEICompound.ri >= min_ri,
            LowResolutionEICompound.ri >= min_rt,
            LowResolutionEICompound.ri >= max_rt,
            )
        
        
        #self.session.query.select(LowResolutionEICompound).where(between(LowResolutionEICompound.ri, min_ri, max_ri))    
        
        return [self.row_to_dict(compound) for compound in compounds]

    def delete_compound(self, compound):
        
        try:
            self.session.delete(compound)  
            self.session.commit()  
        
        except SQLAlchemyError as e:
            self.session.rollback()
            print(str(e))

    def purge(self):
        '''Carefull, this will delete the entire database table'''
        self.session.query(LowResolutionEICompound).delete()
        self.session.commit()  

    def clear_data(self):
        '''clear tables'''
        meta = Base.metadata
        for table in reversed(meta.sorted_tables):
            print ('Clear table %s' % table)
            self.session.execute(table.delete())
        self.session.commit()
   