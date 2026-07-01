from dataclasses import dataclass, field
from typing import List


@dataclass
class TIC_Data:
    """A class to represent total ion chromatogram data.

    scans: [int]
        original scan numbers
    time: [floats]
        list of retention times
    tic: [floats]
        total ion current [chromatogram]
    bpc: [floats]
        base peak [chromatogram]
    Apexes: [int]
        original thermo apex scan number after peak picking
    """

    scans: List[int] = field(default_factory=list)
    time: List[float] = field(default_factory=list)
    tic: List[float] = field(default_factory=list)
    bpc: List[float] = field(default_factory=list)
    apexes: List[int] = field(default_factory=list)


@dataclass
class EIC_Data:
    """A class to represent extracted ion chromatogram data.

    scans: [int]
        original scan numbers
    time: [floats]
        list of retention times
    eic: [floats]
        extracted ion chromatogram
    eic_smoothed: [floats]
        extracted ion chromatogram smoothed
    apexes: [int]
        original apex scan number after peak picking
    areas:  [floats]
        area under the curve for each apex
    """

    scans: List[int] = field(default_factory=list)
    time: List[float] = field(default_factory=list)
    eic: List[float] = field(default_factory=list)
    eic_smoothed: List[float] = field(default_factory=list)
    apexes: List[int] = field(default_factory=list)
    areas: List[float] = field(default_factory=list)
