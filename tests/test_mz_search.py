__author__ = "Yuri E. Corilo"
__date__ = "Jun 11, 2021"

import pytest
import sys
sys.path.append(".")

from corems.mass_spectra.calc.MZSearch import MZSearch
from corems.transient.input.brukerSolarix import ReadBrukerSolarix

def test_mzsearch():

    exp = [212.121, 232.121, 232.122, 234.123, 123.34]
    calculated = [212.12101, 232.1213]

    searchmz = MZSearch(exp, calculated, 1000)
    searchmz.start()
    searchmz.join()

    print(searchmz.results)

if __name__ == '__main__':

    test_mzsearch()
