import sys
sys.path.append(".")

from sqlalchemy import create_engine, ForeignKey, Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

engine = create_engine('sqlite:///db/molecularformulas.db', echo = True)
Base = declarative_base()

class HeteroAtoms(Base):
   
   __tablename__ = 'heteroAtoms'
   id = Column(Integer, primary_key = True)
   name = Column(String, unique=True)
   carbonHydrogen = relationship('CarbonHydrogen', secondary = 'molecularformula')
   
class CarbonHydrogen(Base):
   
   __tablename__ = 'carbonHydrogen'
   id = Column(Integer, primary_key = True)
   name = Column(String, unique=True)
   heteroAtoms = relationship(HeteroAtoms,secondary='molecularformula')

class MolecularFormula(Base):
    
    __tablename__ = 'molecularformula'
    
    heteroAtoms_id = Column(
        Integer, 
        ForeignKey('heteroAtoms.id'), 
        primary_key = True)

    carbonHydrogen_id = Column(
        Integer, 
        ForeignKey('carbonHydrogen.id'), 
        primary_key = True)

if __name__ == "__main__":
    
    from sqlalchemy.orm import sessionmaker
    Session = sessionmaker(bind = engine)
    session = Session()
    Base.metadata.create_all(engine)

    for x in session.query( CarbonHydrogen, HeteroAtoms).filter(MolecularFormula.heteroAtoms_id == HeteroAtoms.id, 
        MolecularFormula.carbonHydrogen_id == CarbonHydrogen.id).order_by(MolecularFormula.carbonHydrogen_id).all():
        print ("MolecularFormula: {} {}".format(x.CarbonHydrogen.name, x.HeteroAtoms.name))
    
    session.close()