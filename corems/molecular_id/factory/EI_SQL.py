

__author__ = "Yuri E. Corilo"
__date__ = "Feb 12, 2020"

import os 

from sqlalchemy import create_engine, Column, Integer, Binary, LargeBinary, String, Float,  exists
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.exc import MultipleResultsFound
from sqlalchemy.pool import QueuePool
from sqlalchemy import between

from numpy import array

Base = declarative_base()

class LowResolutionEICompound(Base):  
    
    __tablename__ = 'low_res_ei_lib'

    id = Column( String, primary_key=True)
    
    name = Column(String, nullable=False)
    formula = Column(String, nullable=True)
    ri = Column(Float, nullable=False)
    rt = Column(Float, nullable=False)
    
    source = Column(String, nullable=True)
    casno = Column(String, nullable=False)
    comment = Column(String, nullable=True)
    
    source_temp_c =  Column(Float, nullable=True)
    ev = Column(Float, nullable=True)

    peaks_count = Column(Integer, nullable=False)
    mz = Column(Binary,  nullable=False)
    abundance = Column(Binary,  nullable=False)

    def __init__(self, dict_data): 
        
        self.id = dict_data.get('NAME')
        
        self.name =dict_data.get('NAME')
        self.formula = dict_data.get('FORM')
        self.ri =dict_data.get('RI')
        self.rt = dict_data.get('RT')

        self.source = dict_data.get('SOURCE')
        self.casno =dict_data.get('CASNO')
        self.comment = dict_data.get('COMMENT')

        self.peaks_count = dict_data.get('NUM PEAKS')
        
        # mz and abun are numpy arrays of 64 bits integer
        # when using postgres array might be a better option
        
        self.mz = array(dict_data.get('mz')).tobytes()
        self.abundance = array(dict_data.get('abundance')).tobytes()
        
    def __repr__(self):
        return "<LowResolutionEICompound(name= %s , cas number = %s, formula = %s, Retention index= %.1f, Retention time= %.1f comment='%i')>" % (
                                    self.name, self.casno, self.formula, self.ri, self.rt, self.comment)

class EI_LowRes_SQLite:
    
    def __init__(self, url='sqlite://'):

        self.engine = create_engine(url,echo=False)

        # to use a online database: 
        # comment previous statement and uncomment next 
        # add database address    
        #self.engine = create_engine('postgresql://postgres:docker@localhost:5432/')
        
        Base.metadata.create_all(self.engine)

        Session = sessionmaker(bind=self.engine)
        
        self.session = Session()
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        # make sure the db connection gets closed
        # 
        self.commit()
        self.session.close()
        self.engine.dispose()

    def __enter__(self):
        
        
        #local engine     
        #directory = os.getcwd()
        #if not os.path.isdir(directory+'/db'): os.mkdir(directory+'/db')    
        #self.engine = create_engine('sqlite:///{DB}'.format(DB=directory+'/db'+'/ei_low_res.sqlite'), poolclass=QueuePool) 
        
        # in-memory
        self.engine = create_engine('sqlite://',echo=False)

        # to use a online database: 
        # comment previous statement and uncomment next 
        # add database address    
        #self.engine = create_engine('postgresql://postgres:docker@localhost:5432/')
        
        Base.metadata.create_all(self.engine)

        Session = sessionmaker(bind=self.engine)
        
        self.session = Session()

        return self
    
    def add_compound_list(self, data_dict_list):
        
        #print(data_dict_list[0])
        for data_dict in data_dict_list:
            #print(data_dict.get('NUM PEAKS'))
            if not data_dict.get('NUM PEAKS'):
                data_dict['NUM PEAKS'] = len(data_dict.get('mz'))
            if not data_dict.get('CASNO'):
                data_dict['CASNO'] = data_dict.get('CAS')
                
        self.session.add_all( [LowResolutionEICompound(data_dict) for data_dict in data_dict_list] )
    
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

          return {c.name: str(getattr(row, c.name)) for c in row.__table__.columns}

    def get_all(self,):
        
        compounds = self.session.query(LowResolutionEICompound).all()
       
        return [self.row_to_dict(compound) for compound in compounds]

    def query_min_max_rt(self, min_max_rt, ):
        
        min_rt, max_rt = min_max_rt

        compounds = self.session.query(LowResolutionEICompound).filter(LowResolutionEICompound.rt.between(min_rt, max_rt))    
        
        return [self.row_to_dict(compound) for compound in compounds]

    def query_min_max_ri(self, min_max_ri):
        
        min_ri, max_ri = min_max_ri

        compounds = self.session.query(LowResolutionEICompound).filter(LowResolutionEICompound.ri.between(min_ri, max_ri))    
        
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
   