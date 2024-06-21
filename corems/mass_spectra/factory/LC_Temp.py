from dataclasses import dataclass, field
from typing import List


@dataclass
class TIC_Data:
     '''
    Scans: [int]
        original thermo scan numbers
    Time: [floats]
        list of retention times
    TIC: [floats]
        total ion current [chromatogram]
    BPC: [floats]
        base peak [chromatogram]
    Apexes: [int]    
        original thermo apex scan number after peak picking 
     '''
     
     scans : List[int] = field(default_factory=list)
     time : List[float] = field(default_factory=list)
     tic : List[float] = field(default_factory=list)
     bpc : List[float] = field(default_factory=list)
     apexes : List[int] = field(default_factory=list)

@dataclass
class EIC_Data:
     '''
    Scans: [int]
        original thermo scan numbers
    Time: [floats]
        list of retention times
    EIC: [floats]
        extracted ion chromatogram
    Apexes: [int]    
        original thermo apex scan number after peak picking 
    
     '''
     
     scans : List[int] = field(default_factory=list)
     time : List[float] = field(default_factory=list)
     eic : List[float] = field(default_factory=list)
     apexes : List[int] = field(default_factory=list)