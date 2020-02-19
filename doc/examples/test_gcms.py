import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("C:\\Users\\eber373\\Desenvolvimento\\Projects-Python\\CoreMS")

from pathlib import Path

from PySide2.QtWidgets import QFileDialog, QApplication
from PySide2.QtCore import Qt

from corems.molecular_id.input.nistMSI import ReadNistMSI


app = QApplication(sys.argv)
file_dialog = QFileDialog()
file_dialog.setWindowFlags(Qt.WindowStaysOnTopHint)
file_location = file_dialog.getOpenFileName()[0]
app.quit()

print(file_location)
sqlLite_obj = ReadNistMSI(file_location).get_sqlLite_obj()

min_max_rt = (18.037, 19.037)

min_max_ri = (1637.30, 1737.30)

sqlLite_obj.query_min_max_ri((1637.30, 1638.30)) 
sqlLite_obj.query_min_max_rt((17.111, 18.111))            
sqlLite_obj.query_min_max_ri_and_rt((1637.30, 1638.30),(17.111, 18.111)) 