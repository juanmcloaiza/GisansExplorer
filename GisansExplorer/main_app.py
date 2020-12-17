#Qt stuff:
from GisansExplorer.utils import Frozen
from GisansExplorer.GUI import MyTabs
import PyQt5.QtWidgets as qtw
from PyQt5.QtGui import QIcon

import traceback

class App(qtw.QMainWindow,Frozen):
    def __init__(self):
        super().__init__()
        self.title = 'GisansExplorer'
        self.myTabs = MyTabs()
        self.setCentralWidget(self.myTabs)
        self._freeze()
        self.setWindowTitle(self.title)
        icon = QIcon('./Icon.png')
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
