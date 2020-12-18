#Qt stuff:
from gisansexplorer.utils import Frozen
from gisansexplorer.GUI import MyTabs
import PyQt5.QtWidgets as qtw
from PyQt5.QtGui import QIcon
from gisansexplorer.utils import __DEPLOYED__

from pathlib import Path
import sys

def icon_file():
    if not __DEPLOYED__:
            filepath = Path('gisansexplorer/resources/Icon.png')
    else:
            filepath = Path(sys.prefix+'/resources/Icon.png')
    return str(filepath)

class App(qtw.QMainWindow,Frozen):
    def __init__(self):
        super().__init__()
        self.title = 'gisansexplorer'
        self.myTabs = MyTabs()
        self.setCentralWidget(self.myTabs)
        self._freeze()
        self.setWindowTitle(self.title)
        icon = QIcon(icon_file())
        self.setWindowIcon(icon)
        self.show()

    def closeEvent(self, event):

        quit_msg = "Are you sure you want to exit?"
        reply = qtw.QMessageBox.question(self, 'Message',
                     quit_msg, qtw.QMessageBox.Yes, qtw.QMessageBox.No)

        if reply == qtw.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


    def addTab(self):
        self.myTabs.addTab()
