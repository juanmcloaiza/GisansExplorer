#!/bin/env python

#Qt stuff:
from PyQt5.QtWidgets import QMainWindow, QFrame
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QDoubleSpinBox, QPushButton, QFormLayout, QMessageBox, QListWidget
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog, QLabel, QTableWidget, QTabWidget
from PyQt5.QtCore import Qt, pyqtSlot, QTimer

#plot stuff:
import pyqtgraph as pg

import numpy as np
import sys
import os
import traceback
import gzip

# Modules for profiling:
import cProfile, pstats, io

def profile_dec(fnc):
    """
    A decorator that uses cProfile to profile a function
    """

    def inner(*args, **kwargs):
        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval

    return inner

@profile_dec
def profile_function_with_arguments(function, *args, **kwargs):
    return function(*args, **kwargs)



class FrozenClass(object):
    __isfrozen = False

    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("%r is a frozen class" % self)
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True


class Canvas(pg.GraphicsLayoutWidget):

    def __init__(self):
        pg.GraphicsLayoutWidget.__init__(self)
        self.w2 = self.addViewBox()
        self.w2.setAspectLocked(True)

        #self.test()

        return


    def test(self):
        n = 1000
        x = np.random.normal(size=n)
        y = np.random.normal(size=n)
        X,Y = np.meshgrid(x,y)
        c = x**2 + y**2
        self.scatter(x,y,c)

        return True


    def scatter(self, xvals, yvals, cvals):
        print("Scatter plot...")
        n = len(xvals)
        curve = pg.ScatterPlotItem(x=xvals, y=yvals, size=5, pen=pg.mkPen(None), pxMode=True)
        spots = [{'pos': (xvals[i],yvals[i]), 'brush':pg.mkColor(cvals[i])} for i in range(n)]
        #curve.addPoints(spots)
        self.w2.addItem(curve)
        return True


class Experiment(FrozenClass):
    def __init__(self):
        self.selector_lambda = None
        self.angle_of_incidence = 0

        # Define the position of the direct beam on the detector
        # and the sensitivity map file
        self.qyc = 528 #530
        self.qzc = 211 #220 #zc=512
        self.sens = None
        self.meansens = None
        self.monitor_counts = None
        self.qy = []
        self.qz = []
        self.I = []
        self.inputd = []

        self.I2d = []
        self.qyrange = []
        self.qzrange = []

        self.cut_qy_xaxis = []
        self.cut_qy = []
        self.cut_qz_xaxis = []
        self.cut_qz = []

        self._freeze()

        return


    @property
    def two_pi_over_lambda(self):
        return 2*np.pi/float(self.selector_lambda)


    @property
    def sin_alpha_i(self):
        return np.sin(np.pi*float(self.angle_of_incidence)/180.0)


    @property
    def cos_alpha_i(self):
        return np.cos(np.pi*float(self.angle_of_incidence)/180.0)


    def sin_2theta_f(self, pixel_i):
        return 0.5755*(self.qyc-pixel_i)/np.sqrt(1990*1990+(0.5755*(self.qyc-pixel_i))*(0.5755*(self.qyc-pixel_i)))


    def sin_alpha_f(self, pixel_j):
        return 0.5755*(pixel_j-self.qzc)/np.sqrt(1990*1990+(0.5755*(self.qzc-pixel_j))*(0.5755*(self.qzc-pixel_j)))


    def cos_alpha_f(self, pixel_j):
        return np.sqrt(1 - self.sin_alpha_f(pixel_j)**2)


class Settings(FrozenClass):
    def __init__(self):
        self.dataPath = None
        self.datFileName = None
        self.yamlFileName = None
        self.gzFileName = None
        self.sensFileName = "sensitivity_map"
        self.cbar_min = None
        self.cbar_max = None

        self._freeze()
        return


    def datFilePath(self):
        return os.path.join(self.dataPath,self.datFileName)


    def yamlFilePath(self):
        return os.path.join(self.dataPath,self.yamlFileName)


    def gzFilePath(self):
        return os.path.join(self.dataPath,self.gzFileName)


    def sensFilePath(self):
        return os.path.join(self.dataPath, self.sensFileName)


    def basename(self):
        return os.path.splitext(self.datFilePath())[0]


    def gisans_map_filepath(self):
        return self.basename()+"_GISANS.map"


    def gisans_cut_filepath(self, y_or_z = "z"):
        return os.path.join(self.basename(),f"_line_cut_q{y_or_z}.out")


class App(QMainWindow,FrozenClass):
    def __init__(self):
        super().__init__()
        self.title = 'Alexandros GISANS Viewer'
        self.myTabs = MyTabs()
        self.setCentralWidget(self.myTabs)
        self._freeze()
        self.setWindowTitle(self.title)
        self.show()


    def addTab(self):
        self.myTabs.addTab()


class MyTabs(QTabWidget,FrozenClass):
    def __init__(self):
        super().__init__()
        self.addTab()
        q = QTabWidget()


    def addTab(self):
        frame = MyFrame()
        super().addTab(frame, "tab" + str(self.count()))
        return



class MyFrame(QFrame,FrozenClass):

    def __init__(self):
        super().__init__()
        self.layout = QHBoxLayout()
        self.centralpanel = QVBoxLayout()
        self.leftpanel = QVBoxLayout()
        self.rightpanel = QVBoxLayout()
        self.minSpinBox = QDoubleSpinBox()
        self.maxSpinBox = QDoubleSpinBox()
        self.settings = Settings()
        self.experiment = Experiment()
        self.canvas = Canvas()
        self.infoTable = QTableWidget()
        self.fileList = QListWidget()
        self.tabs = QTabWidget()
        self.initFrame()
        self._freeze()


    def initFrame(self):
        self.layout.setAlignment(Qt.AlignCenter)
        self.addExperimentInfo()
        self.addFileList()
        self.addWelcomeMessage()
        self.addCanvas()
        self.addMinMaxSpinBoxes()
        self.addFunctionalityButtons()
        self.addPanels()
        self.setLayout(self.layout)


    def addFileList(self):
        self.fileList.setMaximumWidth(self.fileList.width()/2.)
        self.leftpanel.addWidget(QLabel("File:"))
        self.leftpanel.addWidget(self.fileList)
        return


    def addPanels(self):
        self.layout.addLayout(self.leftpanel)
        self.layout.addLayout(self.centralpanel)
        self.layout.addLayout(self.rightpanel)
        return


    def addExperimentInfo(self):
        self.infoTable.setMaximumWidth(self.infoTable.width()/2.)
        self.rightpanel.addWidget(QLabel("Info:"))
        self.rightpanel.addWidget(self.infoTable)
        return


    def addFunctionalityButtons(self):
        buttonOpenDialog = QPushButton("Press here")
        buttonOpenDialog.clicked.connect(self.on_click_open_file)
        self.leftpanel.addWidget(buttonOpenDialog)

        buttonLogLinear = QPushButton("Log / Linear")
        buttonLogLinear .clicked.connect(self.on_click)
        self.rightpanel.addWidget(buttonLogLinear)

        buttonSavePng = QPushButton("Save png")
        buttonSavePng.clicked.connect(self.on_click)
        self.rightpanel.addWidget(buttonSavePng)

        buttonSaveAscii = QPushButton("Save ascii")
        buttonSaveAscii.clicked.connect(self.on_click)
        self.rightpanel.addWidget(buttonSaveAscii)

        return



    @staticmethod
    def init_spinbox(spinbox, slot):
        spinbox.setMinimum(-1.0)
        spinbox.setMaximum(1000000)
        spinbox.setValue(-1.0)
        spinbox.valueChanged.connect(slot)
        return


    def addMinMaxSpinBoxes(self):
        self.init_spinbox(self.minSpinBox, self.on_value_change)
        self.init_spinbox(self.maxSpinBox, self.on_value_change)
        formLayout = QFormLayout()
        formLayout.addRow(self.tr("&Min Intensity"), self.minSpinBox)
        formLayout.addRow(self.tr("&Max Intensity"), self.maxSpinBox)
        formLayout.setFormAlignment(Qt.AlignBottom)
        self.rightpanel.addLayout(formLayout)
        return


    def addWelcomeMessage(self):
        message = QLabel(
            "" + 
            "GISANS viewer for MARIA data - Jan 2019\n"+
            "for questions contact a.koutsioumpas@fz_juelich.de\n"+
            "\n"
        )
        message.setTextInteractionFlags(Qt.TextSelectableByMouse)
        message.setAlignment(Qt.AlignCenter)
        self.centralpanel.addWidget(message)
        return


    def addCanvas(self):
        self.centralpanel.addWidget(self.canvas)
        return


    @pyqtSlot()
    def on_click(self):
        print('PyQt5 button click')


    @pyqtSlot()
    def on_click_open_file(self):
        self.settings = Settings()
        self.doStuff()


    @pyqtSlot()
    def on_value_change(self):
        pass
        return


    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Measurement dat file (*.dat)", options=options)
        # self.openFileNamesDialog()
        # self.saveFileDialog()
        if fileName:
            return fileName 
        return None


    def doStuff(self):
        #self.canvas.test()
        #return
        if not self.read_cbar_min_max():
            return
        if not self.read_dat_file():
            return
        if not self.read_yaml_file():
            return
        if not self.read_sensitivity_file():
            return
        if not self.read_intensity_file():
            return
        if not self.line_cut_y():
            return
        if not self.line_cut_z():
            return
        if not self.show_gisans_map():
            return
        return


    def read_cbar_min_max(self):
        self.settings.cbar_min = self.minSpinBox.value()
        self.settings.cbar_max = self.maxSpinBox.value()
        return True


    def safe_parse(self, parse_func, file_path):
        try:
            with open(file_path, 'r') as fp:
                tf = parse_func(fp)
            print(f" Parsed {file_path}")
            return tf
        except Exception as e:
            self.handle_exception(e)
            return False


    def safe_parse_numpy(self, parse_func, file_path, dtype='i', delimiter=' '):
        try:
            nparray = np.loadtxt(file_path, dtype=dtype, delimiter=delimiter)
            print(f"{file_path}:")
            print(f"Loaded array with shape: {nparray.shape}")
            tf = parse_func(nparray)
            print(f" Parsed {file_path}")
            return tf
        except Exception as e:
            self.handle_exception(e)
            return False


    @staticmethod
    def handle_exception(e):
        msg = (f"Exception: {e}\n")
        msg += ("-"*60+"\n")
        msg += traceback.format_exc()
        msg += ("-"*60+"\n")
        print(msg)
        # traceback.#.print_exc(file=sys.stdout)
        pop_up = QMessageBox()
        pop_up.setWindowTitle(f"Exception: {e}\n")   
        pop_up.setText(msg)
        pop_up.setIcon(QMessageBox.Critical)
        x = pop_up.exec_()
        return


    def read_dat_file(self):
        # Open and read the dat file
        #datFilePath = self.openFileNameDialog()
        datFilePath = "/home/juan/Development/git/imagesNICOS/datafiles/p15347_00001341.dat"
        if datFilePath:
            path, filename = os.path.split(datFilePath)
            self.settings.datFileName = filename
            self.settings.dataPath = path
            return self.safe_parse(self.parse_dat, self.settings.datFilePath())
        return False


    def read_yaml_file(self):
        # Open and read the yaml file
        yamlFileName = self.settings.yamlFileName
        if yamlFileName:
            return self.safe_parse(self.parse_yaml, self.settings.yamlFilePath())
        return False


    def read_sensitivity_file(self):
        if self.settings.sensFileName:
            fpath = self.settings.sensFilePath()
            func = self.parse_sensitivity_map
            return self.safe_parse_numpy(func, fpath, dtype=float, delimiter=' ')
        return False


    def read_intensity_file(self):
        if self.settings.gzFileName:
            fpath = self.settings.gzFilePath()
            func = self.parse_intensity_map
            return self.safe_parse_numpy(func, fpath, dtype=float, delimiter=' ')
        return False


    def parse_dat(self, file):
        for line in file:
            if line.find('omega_value')>0:
                omega_line_list = line.split()
                omega=omega_line_list[3]
                self.experiment.angle_of_incidence = omega
                print('Angle of incidence (degrees): '+str(omega))
                print(f"Original qzc = {self.experiment.qzc}")
                self.experiment.qzc += int( ( 1990.0 * np.tan( np.pi * float(omega) / 180.0 ) ) / 0.5755 )
                print(f"Corrected qzc = {self.experiment.qzc}")
            if line.find('selector_lambda_value')>0:
                lambda_line_list = line.split()
                selector_lambda=lambda_line_list[3]
                self.experiment.selector_lambda = selector_lambda
                print('Neutron wavelength (Angtrom): '+str(selector_lambda))
            if line.find('.gz')>0:
                line_list = line.split()
                for word in line_list:
                    if word.find('.gz')>0:
                        gzFileName = word
                        self.settings.gzFileName = gzFileName
                        print(f"Found gz path: {self.settings.gzFilePath()}")
            if line.find('.yaml')>0:
                line_list = line.split()
                for word in line_list:
                    if word.find('.yaml')>0:
                        yamlFileName = word
                        self.settings.yamlFileName = yamlFileName
                        print(f"Found yaml path: {self.settings.yamlFilePath()}")
        return True


    def parse_yaml(self, fp):
        # Open and read the yaml file
        line1=''
        line2=''
        line3=''
        line4=''
        for line in fp:
            if line4.find('name: mon1')>0:
                line_list=line1.split()
                monitor=line_list[1]
                print('Monitor conuts: '+str(monitor)+'\n')
            line4=line3
            line3=line2
            line2=line1
            line1=line
        self.experiment.monitor_counts = monitor
        return True


    def parse_sensitivity_map(self, sens):
        num=0.0
        meansens=0.0
        pix0 = 162
        pixf = 862
        ipix_range = np.asarray(range(pix0, pixf+1))
        jpix_range = np.asarray(range(pix0, pixf+1))

        for i in ipix_range:
            for j in jpix_range:
                meansens += float(sens[i,j])
                num += 1.0
        meansens=meansens/num
        self.experiment.sens = sens 
        self.experiment.meansens = meansens
        return True


    def parse_intensity_map(self, inputd):
        pix0 = 162
        pixf = 862
        ipix_range = np.asarray(range(pix0, pixf+1))
        jpix_range = np.asarray(range(pix0, pixf+1))

        sens = self.experiment.sens
        meansens = self.experiment.meansens
        monitor = self.experiment.monitor_counts
        two_pi_over_lambda = self.experiment.two_pi_over_lambda
        sin_alpha_i = self.experiment.sin_alpha_i
        sin_2theta_f = self.experiment.sin_2theta_f(ipix_range)
        cos_alpha_f = self.experiment.cos_alpha_f(jpix_range)
        sin_alpha_f = self.experiment.sin_alpha_f(jpix_range)

        qy = np.zeros(len(ipix_range)*len(jpix_range))
        qz = np.zeros(len(ipix_range)*len(jpix_range))
        I = np.zeros(len(ipix_range)*len(jpix_range))
        idx = 0

        with open(self.settings.gisans_map_filepath(), "w") as fp:
            for i in ipix_range:
                for j in jpix_range: 
                    if float(sens[i,j]) > 0.0:
                        I[idx] = ((meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor))
                        qy[idx] = (two_pi_over_lambda *  cos_alpha_f[j-pix0] * sin_2theta_f[i-pix0])
                        qz[idx] = (two_pi_over_lambda * (sin_alpha_f[j-pix0] + sin_alpha_i))
                        fp.write(str(qy[idx])+' '+str(qz[idx])+' '+str(I[idx])+'\n')
                        idx += 1

        self.experiment.inputd = inputd
        self.experiment.qy = qy 
        self.experiment.qz = qz
        self.experiment.I = I 
        self.experiment.I2d = np.nan_to_num(float(meansens) * inputd[pix0:pixf+1,pix0:pixf+1] / sens[pix0:pixf+1,pix0:pixf+1] / float(monitor))
        self.experiment.qyrange = np.flip(two_pi_over_lambda *  cos_alpha_f[jpix_range-pix0] * sin_2theta_f[ipix_range-pix0])
        self.experiment.qzrange = (two_pi_over_lambda * (sin_alpha_f[jpix_range-pix0] + sin_alpha_i))

        return True


    def line_cut_y(self, dqzvalue=None, qzvalue=None):
        print("Performing line cut y")

        self.experiment.cut_qy_xaxis = np.linspace(-1,1)
        self.experiment.cut_qy = 10 + np.sin(self.experiment.cut_qy_xaxis)
        return True
#####TODO
        sens = self.experiment.sens
        meansens = self.experiment.meansens
        monitor = self.experiment.monitor_counts
        inputd = self.experiment.inputd
        two_pi_over_lambda = self.experiment.two_pi_over_lambda
        sin_alpha_i = self.experiment.sin_alpha_i
        cos_alpha_f = self.experiment.cos_alpha_f
        sin_alpha_f = self.experiment.sin_alpha_f
        sin_2theta_f = self.experiment.sin_2theta_f

        cut_qy       = np.zeros(len(self.experiment.qz))
        cut_qy_xaxis = np.zeros(len(self.experiment.qz))
        
        dqzvalue =       0.05         if dqzvalue is None else dqzvalue
        qzvalue  = self.experiment.qzc if  qzvalue is None else  qzvalue

        for i in range(162, 863):
            for j in range(162,863):
                cut_qy_xaxis[i] = two_pi_over_lambda * cos_alpha_f(j) * sin_2theta_f(i)
                cqz = two_pi_over_lambda * (sin_alpha_f(j) + sin_alpha_i)
#                if cqz >= (qzvalue - dqzvalue) and cqz <= (qzvalue + dqzvalue):
                cut_qy[i] += (meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor)

        self.experiment.cut_qy_xaxis = cut_qy_xaxis
        self.experiment.cut_qy = cut_qy
        
        return True
######


    def line_cut_z(self, dqyvalue=None, qyvalue=None):
        print("Performing line cut z")

        self.experiment.cut_qz_xaxis = np.linspace(-1,1)
        self.experiment.cut_qz = 10 + np.cos(self.experiment.cut_qz_xaxis)
        return True
###TODO
        sens = self.experiment.sens
        meansens = self.experiment.meansens
        monitor = self.experiment.monitor_counts
        inputd = self.experiment.inputd
        two_pi_over_lambda = self.experiment.two_pi_over_lambda
        sin_alpha_i = self.experiment.sin_alpha_i
        cos_alpha_f = self.experiment.cos_alpha_f
        sin_alpha_f = self.experiment.sin_alpha_f
        sin_2theta_f = self.experiment.sin_2theta_f

        cut_qz       = np.zeros(len(self.experiment.qy))
        cut_qz_xaxis = np.zeros(len(self.experiment.qy))
        
        dqyvalue =       0.05         if dqyvalue is None else dqyvalue
        qyvalue  = self.experiment.qyc if  qyvalue is None else  qyvalue

        for i in range(162, 863):
            for j in range(162,863):
                cut_qz_xaxis[j] = two_pi_over_lambda * (sin_alpha_f(j) + sin_alpha_i)
                cqy = two_pi_over_lambda * cos_alpha_f(j) * sin_2theta_f(i)
                #if cqy >= (qyvalue - dqyvalue) and cqy <= (qyvalue + dqyvalue):
                cut_qz[j] += (meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor)
        
        self.experiment.cut_qz_xaxis = cut_qz_xaxis
        self.experiment.cut_qz = cut_qz

        return True
###

    def show_gisans_map(self):
        try:
            if self.settings.cbar_min < 0:
                min_value = min(self.experiment.I)+float(1)/float(self.experiment.monitor_counts)
            else:
                min_value = self.settings.cbar_min

            if self.settings.cbar_max < 0:
                max_value=max(self.experiment.I)
            else:
                max_value = self.settings.cbar_max

            self.canvas.scatter(self.experiment.qy,
                            self.experiment.qz,
                            self.experiment.I,
                            #self.experiment.cut_qy_xaxis,
                            #self.experiment.cut_qz_xaxis,
                            #self.experiment.cut_qy,
                            #self.experiment.cut_qz,
                            #vmin=min_value,
                            #vmax=max_value
                            )
        except Exception as e:
            App.handle_exception(e)
            return False

        return True

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    #sys.exit(app.exec_())
    x = profile_function_with_arguments(app.exec_)
    sys.exit(x)