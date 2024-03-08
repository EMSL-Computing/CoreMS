"""
.. include:: ../README.md
"""

__author__ = 'Yuri E. Corilo'
__version__ = '2.0.1'

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
            print("%r  %2.2f ms" % (method.__name__, (te - ts) * 1000))
        return result
    return timed


class SuppressPrints:

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def get_filenames(app=None):

    from PySide2.QtCore import Qt, QCoreApplication
    from PySide2.QtWidgets import QApplication, QFileDialog
    from pathlib import Path

    app = QApplication(sys.argv)
    file_dialog = QFileDialog()
    file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    file_location, _ = file_dialog.getOpenFileNames()

    if file_location:
        QCoreApplication.processEvents()
        return file_location

    else:

        return None

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
    # file_dialog = QFileDialog()
    # file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
    # file_location = file_dialog.getOpenFileNames()

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
