#!/bin/env python

#Qt stuff:
from PyQt5.QtWidgets import QMainWindow, QFrame, QToolButton, QTableWidgetItem
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QDoubleSpinBox, QPushButton, QFormLayout, QMessageBox, QListWidget
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog, QLabel, QTableWidget, QTabWidget
from PyQt5.QtCore import Qt, pyqtSlot, pyqtSignal

#plot stuff:
from matplotlib.ticker import NullFormatter
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

import numpy as np
import sys
import os
import traceback
import gzip

# Modules for profiling:
import cProfile, pstats, io

class Canvas(FigureCanvas):
    def __init__(self, parent=None, width=1000, height=28):
        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        spacing = 0.005
        height_hist = 0.2
        width_cbar = 0.05

        rect_histx = [left + height_hist + spacing, bottom, width, height_hist]
        rect_histy = [left, bottom + height_hist + spacing, height_hist, width]
        rect_scatter = [left + height_hist + spacing, bottom + height_hist + spacing, width, height]
        rect_cbar = [left + height_hist + spacing + width + spacing, bottom + height_hist + spacing, width_cbar, height]

        # start with a rectangular Figure
        fig = plt.figure(figsize=(10, 10))
        
        ax_center = plt.axes(rect_scatter)
#       ax_center.tick_params(direction='in', top=True, right=True)

        ax_histx = plt.axes(rect_histx)
#       ax_histx.tick_params(direction='in', labelbottom=False)

        ax_histy = plt.axes(rect_histy)
#       ax_histy.tick_params(direction='in', labelleft=False)

        ax_cbar = plt.axes(rect_cbar)
        
       
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        self.ax_center = ax_center
        self.ax_histx = ax_histx
        self.ax_histy = ax_histy
        self.ax_cbar = ax_cbar
        self.colorbar = None
        #self.test()
        return


    def test(self):
        x = np.linspace(1,10,1025)
        y = np.linspace(3,7,1025)
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)
        T = np.arctan(Y/X)
        Z = (10 + 1. / (1 + (X - x[512])**2 + (Y - y[512])**2))
        xcut = 10**np.cos(x)
        ycut = 10**np.sin(y)
        self.my_plot(R, T, Z, x, xcut, y, ycut)
        return True


    def my_plot(self, x_grid, y_grid, z_vals, x_range, zx_vals, y_range, zy_vals, use_logscale="true"):
        self.ax_center.cla()
        self.ax_histx.cla()
        self.ax_histy.cla()
        self.ax_cbar.cla()
        print("generating figure...")

        fig = self.figure
        ax_center = self.ax_center
        ax_histx = self.ax_histx
        ax_histy = self.ax_histy

        ax_center.set_xticks([])
        ax_center.set_yticks([])
        z_vals += 1e-10

        vmin = np.min(z_vals)
        vmax = np.max(z_vals)


        if use_logscale:
            zcontours = ax_center.imshow(z_vals, norm = LogNorm())
            ax_histx.set_yscale('log')
            ax_histy.set_xscale('log')
        else:
            zcontours = ax_center.imshow(z_vals)

        xcontours = ax_center.contour(x_grid.T, 10, colors='black', linewidth=.5, )
        ycontours = ax_center.contour(y_grid.T, 10, colors='black', linewidth=.5)
        ax_center.clabel(xcontours, fontsize=20, inline=1)
        ax_center.clabel(ycontours, fontsize=20, inline=1)

        ax_histx.plot(x_range, zx_vals)
        ax_histx.set_xlim(x_range.min(), x_range.max())
        
        ax_histy.plot(zy_vals, y_range)
        ax_histy.set_ylim(y_range.min(), y_range.max())
 
        ax_center.set_title("GISANS map")

        # Now adding the colorbar
        self.colorbar = self.figure.colorbar(zcontours, cax = self.ax_cbar)
        self.colorbar.set_label("Intensity")
        ax_histx.set_xlabel(r'$Q_{y}(\AA^{-1})$', fontsize = 20)
        ax_histy.set_ylabel(r'$Q_{z}(\AA^{-1})$', fontsize = 20) 
        self.draw()
        print("Figure Generated")

        return True


    def scatter(self, xarr, yarr, Iarr, xcutaxis, ycutaxis, xcut, ycut, vmin = None, vmax = None):
        print("generating figure...")
        if vmin is None:
            vmin = min(Iarr)

        if vmax is None:
            vmax = min(Iarr)

        fig = self.figure
        ax_center = self.ax_center
        ax_histx = self.ax_histx
        ax_histy = self.ax_histy
        #choice_size = 1024
        #idx = np.random.choice(range(len(xarr)), size=choice_size, replace=False)
        x = xarr
        y = yarr
        I = Iarr

        # the scatter plot:
        im = ax_center.scatter(x, y, edgecolors='none', c=I, 
                            norm = LogNorm(vmin=vmin, vmax=vmax),
                            marker = '.')

        # now determine nice limits by hand:
        ax_center.set_xlim((x.min(), x.max()))
        ax_center.set_ylim((y.min(), y.max()))

        ax_histx.scatter(xcutaxis, xcut)
        ax_histx.set_yscale('log')
        ax_histx.set_xlim(ax_center.get_xlim())
        
        ax_histy.scatter(ycut, ycutaxis)
        ax_histy.set_xscale('log')
        ax_histy.set_ylim(ax_center.get_ylim())
 
        ax_center.set_title("GISANS map")
        self.colorbar = self.figure.colorbar(im)
        self.colorbar.set_label("Intensity")
        ax_center.set_xlabel(r'$Q_{y}(\AA^{-1})$')
        ax_center.set_ylabel(r'$Q_{z}(\AA^{-1})$') 
        print("Figure Generated")

        return True

    def save_png(self,filepath):
        self.figure.savefig(filepath)



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



class Experiment(FrozenClass):
    def __init__(self):
        self.selector_lambda = None
        self.angle_of_incidence = 0

        # Define the position of the direct beam on the detector
        # and the sensitivity map file
        self.qyc = 528 
        self.qzc = 211 
        self.sens = None
        self.meansens = None
        self.monitor_counts = None
        self.qy = []
        self.qz = []
        self.I = []
        self.inputd = []

        self.Imatrix = []
        self.qymatrix = []
        self.qzmatrix = []

        self.cut_Iz = []
        self.cut_Iy = []

        self._freeze()

        return


        #sin_alpha_f(j)        = px_size*(j-zc)/np.sqrt(sdd*sdd+(px_size*(zc-j))*(px_size*(zc-j)))
        #sin_2theta_f(i)       = px_size*(i-yc)/np.sqrt(sdd*sdd+(px_size*(yc-i))*(px_size*(yc-i)))
        #cos_alpha_f(j)        = sdd/np.sqrt(sdd*sdd+(px_size*(zc-j))*(px_size*(zc-j)))


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
        return 0.5755*(pixel_i-self.qyc)/np.sqrt(1990*1990+(0.5755*(self.qyc-pixel_i))*(0.5755*(self.qyc-pixel_i)))


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
        self.use_logscale = True

        self._freeze()
        return


    def datFilePath(self):
        if self.dataPath is None:
            return None
        return os.path.join(self.dataPath,self.datFileName)


    def yamlFilePath(self):
        if self.dataPath is None:
            return None
        return os.path.join(self.dataPath,self.yamlFileName)


    def gzFilePath(self):
        if self.dataPath is None:
            return None
        return os.path.join(self.dataPath,self.gzFileName)


    def sensFilePath(self):
        if self.dataPath is None:
            return None
        return os.path.join(self.dataPath, self.sensFileName)


    def basename(self):
        if self.dataPath is None:
            return None
        return os.path.splitext(self.datFilePath())[0]


    def gisans_map_filepath(self):
        if self.dataPath is None:
            return None
        return self.basename()+"_GISANS.map"


    def gisans_cut_filepath(self, y_or_z = "z"):
        if self.dataPath is None:
            return None
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

    @staticmethod
    def handle_exception(e):
        msg = (f"Exception: {e}\n")
        msg += ("-"*60+"\n")
        msg += traceback.format_exc()
        msg += ("-"*60+"\n")
        print(msg)
        pop_up = QMessageBox()
        pop_up.setWindowTitle(f"Exception: {e}\n")
        pop_up.setText(msg)
        pop_up.setIcon(QMessageBox.Critical)
        x = pop_up.exec_()
        return



class MyTabs(QTabWidget,FrozenClass):
    def __init__(self):
        super().__init__()
        self.tabButton_add = QToolButton()
        self.tabButton_rmv = QToolButton()
        self.frameList =[]
        self.last_num = 0
        self.addTab()
        self.initCornerButton()
        self._freeze()

    def initCornerButton(self):
        self.setCornerWidget(self.tabButton_add,corner=Qt.TopLeftCorner)
        self.tabButton_add.setText('+')
        font = self.tabButton_add.font()
        font.setBold(True)
        self.tabButton_add.setFont(font)
        self.tabButton_add.clicked.connect(self.addTab)

        self.setCornerWidget(self.tabButton_rmv,corner=Qt.TopRightCorner)
        self.tabButton_rmv.setText('-')
        font = self.tabButton_rmv.font()
        font.setBold(True)
        self.tabButton_rmv.setFont(font)
        self.tabButton_rmv.clicked.connect(self.removeTab)
        return



    def addTab(self):
        frame = MyFrame()
        super().addTab(frame, "New Experiment " + str(1 + self.last_num))
        self.setCurrentIndex(self.last_num)
        self.last_num += 1
        self.frameList.append(frame)

        for i, f in enumerate(self.frameList):
            name = f.settings.datFileName
            if name is not None:
                self.setTabText(i,name)


        return


    def removeTab(self):
        idx = self.currentIndex()
        del self.frameList[idx]
        super().removeTab(idx)
        if len(self.frameList) == 0:
            self.last_num = 0
            self.addTab()

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
        self.finishedDoingStuff = pyqtSignal()
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
        self.infoTable.setColumnCount(1)
        self.infoTable.horizontalHeader().hide()
        self.rightpanel.addWidget(QLabel("Info:"))
        self.rightpanel.addWidget(self.infoTable)
        return


    def addFunctionalityButtons(self):
        buttonOpenDialog = QPushButton("Press here")
        buttonOpenDialog.clicked.connect(self.on_click_open_file)
        self.leftpanel.addWidget(buttonOpenDialog)

        buttonLogLinear = QPushButton("Log / Linear")
        buttonLogLinear .clicked.connect(self.on_click_loglinear)
        self.rightpanel.addWidget(buttonLogLinear)

        buttonSavePng = QPushButton("Save png or pdf")
        buttonSavePng.clicked.connect(self.on_click_save_png)
        self.rightpanel.addWidget(buttonSavePng)

        buttonSaveAscii = QPushButton("Save ascii")
        buttonSaveAscii.clicked.connect(self.on_click_save_ascii)
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
    def on_click_loglinear(self):
        self.settings.use_logscale = not self.settings.use_logscale
        self.show_gisans_map()


    @pyqtSlot()
    def on_click_save_png(self):
        try:
            filepath = self.saveFileNameDialog()
            if filepath is None:
                return

            extension = os.path.splitext(filepath)[-1]
            if extension != ".pdf" and extension != ".png":
                filepath += ".png"

            self.canvas.save_png(filepath)
            print(f"Figure saved: {filepath}")
        except Exception as e:
            App.handle_exception(e)


    @pyqtSlot()
    def on_click_save_ascii(self):
        try:
            filepath = self.saveFileNameDialog()
            if filepath is None:
                return

            noextpath, ext = os.path.splitext(filepath)
            if ext == '': ext = ".txt"
            file_qz = noextpath+"_qz_"+ext
            file_qy = noextpath+"_qy_"+ext
            file_Iyz = noextpath+"_Iyz_"+ext
            np.savetxt(file_qz, self.experiment.qzmatrix)
            np.savetxt(file_qy, self.experiment.qymatrix)
            np.savetxt(file_Iyz,  self.experiment.Imatrix)

            print(f"Arrays saved:\n {file_qz}\n {file_qy}\n {file_Iyz}\n")
        except Exception as e:
            App.handle_exception(e)


    @pyqtSlot()
    def on_click_open_file(self):
        self.settings = Settings()
        self.doStuff()
        self.populateWidgets()
        sum1 = self.experiment.Imatrix.sum()
        sum2 = self.experiment.cut_Iy.sum()
        sum3 = self.experiment.cut_Iz.sum()
        print("If the following three sums are not equal, there's an error:")
        print(sum1, sum2, sum3)



    @pyqtSlot()
    def on_value_change(self):
        pass
        return


    def saveFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Save File", "","All Files (*);;png file (*.png);;pdf file (*.pdf)", options=options)
        if fileName:
            return fileName
        return None

    def openFileNameDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Measurement dat file (*.dat)", options=options)
        # self.openFileNamesDialog()
        if fileName:
            return fileName
        return None




    def populateWidgets(self):
        x = self.settings.datFileName
        y = self.experiment.__dict__

        for i, k in enumerate(y.keys()):
            if k[0] == "_":
                continue
            item_k = QTableWidgetItem(str(k))
            item_v = QTableWidgetItem(str(y[k]))
            self.infoTable.insertRow(i)
            self.infoTable.setVerticalHeaderItem(i,item_k) 
            self.infoTable.setItem(i,0,item_v)

        

        return True


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
        #if not self.line_cut_at_constant_y():
        #    return
        #if not self.line_cut_at_constant_z():
        #    return
        if not self.show_gisans_map():
            return

        self.finishedDoingStuff
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
            App.handle_exception(e)
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
            App.handle_exception(e)
            return False



    def read_dat_file(self):
        # Open and read the dat file
        datFilePath = self.openFileNameDialog()
        #datFilePath = "/home/juan/Development/git/imagesNICOS/datafiles/p15347_00001341.dat"
        if datFilePath:
            path, filename = os.path.split(datFilePath)
            self.settings.datFileName = filename
            self.settings.dataPath = path
            self.fileList.addItem(path)
            self.fileList.addItem(filename)
            return self.safe_parse(self.parse_dat, self.settings.datFilePath())
        return False


    def read_yaml_file(self):
        # Open and read the yaml file
        yamlFileName = self.settings.yamlFileName
        if yamlFileName:
            self.fileList.addItem(yamlFileName)
            return self.safe_parse(self.parse_yaml, self.settings.yamlFilePath())
        return False


    def read_sensitivity_file(self):
        if self.settings.sensFileName:
            self.fileList.addItem(self.settings.sensFileName)
            fpath = self.settings.sensFilePath()
            func = self.parse_sensitivity_map
            return self.safe_parse_numpy(func, fpath, dtype=float, delimiter=' ')
        return False


    def read_intensity_file(self):
        if self.settings.gzFileName:
            self.fileList.addItem(self.settings.gzFileName)
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

        self.experiment.Imatrix = np.nan_to_num(float(meansens) * inputd[pix0:pixf+1,pix0:pixf+1] / sens[pix0:pixf+1,pix0:pixf+1] / float(monitor))
        self.experiment.cut_Iz = self.experiment.Imatrix.sum(axis=1)
        self.experiment.cut_Iy = self.experiment.Imatrix.sum(axis=0)

        QY_i, QY_j = sin_2theta_f[ipix_range-pix0], cos_alpha_f[jpix_range-pix0]
        QYmatrix = np.einsum('i,j->ij', QY_i, QY_j)
        self.experiment.qymatrix = two_pi_over_lambda * QYmatrix
        
        QZ_j = sin_alpha_f[jpix_range-pix0]
        QZ_i = np.ones(ipix_range.shape)
        QZmatrix = np.einsum('i,j->ij', QZ_i, QZ_j)
        self.experiment.qzmatrix = two_pi_over_lambda * (QZmatrix + sin_alpha_i)

        return True


    def line_cut_at_constant_y(self, dqyvalue=None, qyvalue=None):
        print("Performing line cut at constant qy")
        qzmatrix = self.experiment.qzmatrix
        Imatrix = self.experiment.Imatrix
        qz = self.experiment.qz
        qzunique = np.unique(qzmatrix)
        Iz = np.zeros(qzunique.shape)


        for idx, qzbin in enumerate(qzunique):
            indices_for_which_qz_is_equal_to_something = np.argwhere(qzmatrix == qzbin)
            Iz[idx] = np.take(Imatrix, indices_for_which_qz_is_equal_to_something).sum()

        self.experiment.cut_Iz = Iz
        return True


    def line_cut_at_constant_z(self, dqzvalue=None, qzvalue=None):
        print("Performing line cut at constant qz")
        qymatrix = self.experiment.qymatrix
        Imatrix = self.experiment.Imatrix
        qy = self.experiment.qy
        qyunique = np.unique(qymatrix)
        Iy = np.zeros(qyunique.shape)

#        for idx, qybin in enumerate(np.unique(qymatrix)):
#            indices_for_which_qy_is_equal_to_something = np.argwhere(qy == qybin)
#            Iy[idx] = np.take(Imatrix, indices_for_which_qy_is_equal_to_something).sum()

        self.experiment.cut_Iy = Iy
        return True


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

            self.canvas.my_plot(
                            self.experiment.qymatrix,
                            self.experiment.qzmatrix,
                            self.experiment.Imatrix,
                            self.experiment.qymatrix.T[0],
                            self.experiment.cut_Iy,
                            self.experiment.qzmatrix[0],
                            self.experiment.cut_Iz,
                            use_logscale=self.settings.use_logscale 
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