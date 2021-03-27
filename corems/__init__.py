__author__ = 'Yuri E. Corilo'
__version__ = '1.0.0'
__doc__ = '''
CoreMS - a powerful framework for mass spectrometry data processing and analysis of small molecules
=====================================================================

    Data handling and software development for modern mass spectrometry (MS) is an interdisciplinary endeavor 
    requiring skills in computational science and a deep understanding of MS. To enable scientific software development to 
    keep pace with fast improvements in MS technology, we have developed a Python software framework named CoreMS. 
    The goal of the framework is to provide a fundamental, high-level basis for working with all mass spectrometry 
    data types, allowing custom workflows for data signal processing, annotation, and curation. The data structures were 
    designed with an intuitive, mass spectrometric hierarchical structure, thus allowing organized and easy access to the
    data and calculations. Moreover, CoreMS supports direct access for almost all vendors data formats, allowing for 
    the centralization and automation of all data processing workflows from the raw signal to data annotation and 
    curation. 

Main Features
-------------
Here are just a few of the things that CoreMS does well:

    Data input formats

    - Bruker Solarix CompassXtract
    - Bruker Solarix transients, ser and fid (FT and magnitude mode)
    - ThermoFisher raw file
    - Spectroswiss Signal booster data-acquisition station HFD 5 
    - Midas (.dat) from MagLab ICR data-acquisition station (FT-magnitude mode)
    - Mass list in Profile and Centroid Mode (include all delimiters types and Excel)
    - CoreMS exported processed mass list files (Excel, csv, txt, etc)
    - Panda dataframe, processed or unprocessed

    Data output formats

    - Text Files (csv, tab separated txt, etc)
    - Microsoft Excel (.xlsx)
    - Hierarchical Data Format (.h5) (`TODO`)
    - Pandas data frame (can be saved using pickle, h5, etc)

    Data structure type

    - GC-MS     (`TODO`)
    - LC-MS
    - IMS-MS    (`TODO`)
    - LC-IMS-MS (`TODO`)
    - Transient
    - Mass Spectra
    - Mass Spectrum

    Available features

        Signal Processing

        - Apodization, Zerofilling, and Magnitude mode FT
        - Manual and automatic noise threshold calculation
        - Peak picking apex quadratic fitting
        - Resolving Power calculation

        Calibration

        - Frequency and m/z domain calibration functions:
        - ledford equation [ref]
        - linear equation
        - quadratic equation
        - Automatic search most abundant **Ox** homologue series
        - step fit ('walking calibration") based on the ledford equation [ref]

        Molecular formulae search and assignment

        - Automatic local or external database generation
        - Automatic molecular formulae assignments for ESI(-) MS for natural organic matter analysis
        - Automatic fine isotopic structure calculation and search
        - Flexible Kendrick normalization base
        - Kendrick filter using density-based clustering

        Mass spectrum simulations

        - Peak shape (Lorentz and Gaussian)
        - m/z error distribution
        - ICR Resolving Power based on magnetic field (B), and transient time(T)
'''
import time
import os
import sys
import hashlib

		
def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print( "%r  %2.2f ms" % (method.__name__, (te - ts) * 1000))
        return result
    return timed


class SuppressPrints:
    
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout  

def get_filename(app=None):
    
    from PySide2.QtCore import Qt, QCoreApplication
    from PySide2.QtWidgets import QApplication, QFileDialog 
    from pathlib import Path

    app = QApplication(sys.argv)
    file_dialog = QFileDialog()
    file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    file_location, _ = file_dialog.getOpenFileName()
    
    if file_location:
        QCoreApplication.processEvents()
        return Path(file_location)
    
    else:
        
        return None

def get_dirname(app=None):
    
    from PySide2.QtCore import Qt, QCoreApplication
    from PySide2.QtWidgets import QApplication, QFileDialog 
    from pathlib import Path

    app = QApplication(sys.argv)
    file_dialog = QFileDialog()
    file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    file_location = file_dialog.getExistingDirectory()
    
    if file_location:
        QCoreApplication.processEvents()
        return Path(file_location)
    
    else:
        
        return None

def get_dirnames(app=None):   
    
    from PySide2.QtCore import Qt, QCoreApplication
    from PySide2.QtWidgets import QApplication, QFileDialog, QTreeView, QListView, QAbstractItemView
    from pathlib import Path
    
    if not app:
        app = QApplication(sys.argv)
    #file_dialog = QFileDialog()
    #file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    #file_location = file_dialog.getOpenFileNames()
    
    file_dialog = QFileDialog()
    file_dialog.setFileMode(QFileDialog.DirectoryOnly)
    file_dialog.setOption(QFileDialog.DontUseNativeDialog, True)
    file_view = file_dialog.findChild(QListView, 'listView')

    # to make it possible to select multiple directories:
    if file_view:
        file_view.setSelectionMode(QAbstractItemView.MultiSelection)
    f_tree_view = file_dialog.findChild(QTreeView)
    if f_tree_view:
        f_tree_view.setSelectionMode(QAbstractItemView.MultiSelection)

    if file_dialog.exec():
        paths = file_dialog.selectedFiles()
        QCoreApplication.processEvents()
        for path in paths:
            yield Path(path)

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def corems_md5(fname):

    bytes_io = fname.open('rb').read()

    md5_returned = hashlib.sha256(bytes_io).hexdigest()
    
    return "{}:{}".format("sha256", md5_returned)
        