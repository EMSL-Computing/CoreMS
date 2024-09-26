import sys
import hashlib

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
