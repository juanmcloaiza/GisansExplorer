#!/bin/env python

#Qt stuff:
import PyQt5.QtWidgets as qtw
from PyQt5.QtCore import Qt, pyqtSlot, pyqtSignal
from PyQt5.QtGui import QValidator

#plot stuff:
from matplotlib.ticker import NullFormatter
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.widgets import RectangleSelector
from matplotlib.figure import Figure

import numpy as np
import sys
import os
import re
import traceback
import gzip

# Modules for profiling:
import cProfile, pstats, io

_DEBUG_ = False


def enable_high_dpi_scaling():
    if hasattr(Qt, 'AA_EnableHighDpiScaling'):
        qtw.QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)

    if hasattr(Qt, 'AA_UseHighDpiPixmaps'):
        qtw.QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

    return


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


#### sciSpinBox

# Regular expression to find floats. Match groups are the whole string, the
# whole coefficient, the decimal part of the coefficient, and the exponent
# part.
_float_re = re.compile(r'(([+-]?\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)')

def valid_float_string(string):
    match = _float_re.search(string)
    return match.groups()[0] == string if match else False


class FloatValidator(QValidator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def validate(self, string, position):
        if valid_float_string(string):
            return QValidator.Acceptable, string, position
        if string == "" or string[position-1] in 'e.-+':
            return QValidator.Intermediate, string, position
        return QValidator.Invalid, string, position

    def fixup(self, text):
        match = _float_re.search(text)
        return match.groups()[0] if match else ""


class mySciSpinBox(qtw.QDoubleSpinBox):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setMinimum(-np.inf)
        self.setMaximum(np.inf)
        self.validator = FloatValidator()
        self.setDecimals(10)

    def validate(self, text, position):
        try:
            asdf = self.validator.validate(text, position)
        except Exception as e:
            App.handle_exception(e)
        return asdf

    def fixup(self, text):
        return self.validator.fixup(text)

    def valueFromText(self, text):
        return float(text)

    def textFromValue(self, value):
        return format_float(value)

    def stepBy(self, steps):
        text = self.cleanText()
        groups = _float_re.search(text).groups()
        decimal = float(groups[1])
        decimal += steps
        new_string = "{:g}".format(decimal) + (groups[3] if groups[3] else "")
        self.lineEdit().setText(new_string)


def format_float(value):
    """Modified form of the 'g' format specifier."""
    string = "{:g}".format(value).replace("e+", "e")
    string = re.sub("e(-?)0*(\d+)", r"e\1\2", string)
    return string

#### sciSpinBox end

class FrozenClass(object):
    __isfrozen = False

    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("%r is a frozen class" % self)
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True




class AreaSelector(FrozenClass):
    def __init__(self, ax, line_select_callback):
        self.ax = ax
        self.rs = RectangleSelector(ax, line_select_callback,
                                                drawtype='box' , useblit=False,
                                                button=[1, 3],  # don't use middle button
                                                minspanx=0, minspany=0,
                                                spancoords='pixels',
                                                interactive=True)


    def __call__(self, event):
        self.rs.update()
        if self.ax == event.inaxes:
            if event.key in ['Q', 'q']:
                self.rs.to_draw.set_visible(False)
                self.rs.set_active(False)
            if event.key in ['A', 'a']:
                self.rs.to_draw.set_visible(True)
                self.rs.set_active(True)

        return #__call__


class PlotData(FrozenClass):
    def __init__(self):
        self.X = np.zeros((10,10))
        self.Y = np.zeros((10,10))
        self.Z = np.zeros((10,10))

        self.Xzoom = np.zeros((10,10))
        self.Yzoom = np.zeros((10,10))
        self.Zzoom = np.zeros((10,10))
        return


class PlotParams(FrozenClass):
    def __init__(self):
        self.x1 = 0
        self.x2 = 10
        self.y1 = 0
        self.y2 = 10

        self.xmin = -1
        self.xmax = 1
        self.ymin = -1
        self.ymax = 1

        self.zmin = -1
        self.zmax = 1

        self.log_scale = False
        self.reset_limits_required = True

        self.title = ""
        return


class MyGraphView(qtw.QWidget):
    finishedUpdating = pyqtSignal()
    def __init__(self, graph_title, parent = None):
        super(MyGraphView, self).__init__(parent)

        self.graph_title = graph_title

        self.dpi = 100
        self.fig = Figure((10.0, 5.0), dpi = self.dpi, facecolor = (1,1,1), edgecolor = (0,0,0))
        self.canvas = FigureCanvas(self.fig)
        self.define_axes()

        self.init_data_and_parameters()
        self.init_xyzLabel()
        #self.commands = MyGraphCommands(self.update_graph)
        self.define_layout()
        self.init_canvas_connections()
        self.canvas.draw()
        self.test_show()
        return #__init__()


    def init_xyzLabel(self):
        self.xyzLabel = qtw.QLabel(self)
        self.xyzLabel.setText("")
        return #init_xyzLabel


    def define_layout(self):
        self.layout = qtw.QVBoxLayout()
        self.layout.addWidget(self.xyzLabel)
        self.layout.addWidget(self.canvas)
        self.layout.setStretchFactor(self.canvas, 1)
        self.setLayout(self.layout)
        return #define_layout


    def define_axes(self):
        self.ax =  self.canvas.figure.add_axes([0.1,0.1,0.30,0.7])
        self.cax = self.canvas.figure.add_axes([0.505,0.1,0.025,0.7])
        self.zoom_ax = self.canvas.figure.add_axes([0.55,0.25,0.25,0.5])
        self.xax = self.canvas.figure.add_axes([0.55,0.1,0.25,0.1])
        self.yax = self.canvas.figure.add_axes([0.85,0.25,0.05,0.5])
        self.area_selector = AreaSelector(self.ax, self.line_select_callback)
        return #define_axes


    def init_canvas_connections(self):
        # connect mouse events to canvas
        self.canvas.mpl_connect('scroll_event', self.on_mouse_wheel)
        self.canvas.mpl_connect('key_press_event', self.area_selector)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        #self.canvas.mpl_connect('draw_event', self.area_selector.mycallback)
        self.canvas.setFocusPolicy( Qt.ClickFocus)
        self.canvas.setFocus()

        return #init_canvas_connections


    def init_data_and_parameters(self):
        self.data = PlotData()
        self.params = PlotParams()
        return #init_data_and_parameters


    def on_mouse_wheel(self, event):
        if self.cax == event.inaxes:
            if event.button == 'up':
                func = lambda c: (np.sign(c) * abs(c)**0.9) if self.params.log_scale else 0.9*c
            elif event.button == 'down':
                func = lambda c: (np.sign(c) * abs(c)**1.1) if self.params.log_scale else 1.1*c
            else:
                return

            vmin, vmax = (func(c) for c in self.cbar.get_clim())
            print("Rescaling colorbar:", vmin,vmax)
            self.update_graph(zmin=vmin, zmax=vmax)
        return #on_mouse_wheel


    def on_mouse_move(self, event):
        x, y = event.xdata,event.ydata
        if x is None or y is None:
            return
        z = self.data.Z[int(y), int(x)]
        self.xyzLabel.setText(f"(x, y; z) = ({x:3.2g}, {y:3.2g}; {z:3.2g})")
        return #on_mouse_move


    def line_select_callback(self, eclick, erelease):
        x1, y1 = int(eclick.xdata), int(eclick.ydata)
        x2, y2 = int(erelease.xdata), int(erelease.ydata)
        print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
        print(" The button you used were: %s %s" % (eclick.button, erelease.button))
        self.update_graph(x1=x1, x2=x2, y1=y1, y2=y2)
        return #line_select_callback


    def update_data(self, **kwargs):
        for k in self.data.__dict__.keys():
            if k in kwargs.keys():
                self.data.__setattr__(k,kwargs[k])
        return #update_data


    def update_params(self, **kwargs):
        for k in self.params.__dict__.keys():
            if k in kwargs.keys():
                self.params.__setattr__(k,kwargs[k])

        if self.params.reset_limits_required:
            self.params.zmin = self.data.Z.min()
            self.params.zmax = self.data.Z.max()
            self.params.reset_limits_required = False
        return #update_params


    def update_axes(self, **kwargs):
        self.update_ax(**kwargs)
        self.update_cax()
        self.update_zoom_ax()
        self.update_xax()
        self.update_yax()
        return #update_axes


    def update_graph(self, **kwargs):
        self.canvas.figure.clear()
        self.define_axes()
        self.update_data(**kwargs)
        self.update_params(**kwargs)
        self.build_norm(**kwargs)
        self.update_axes(**kwargs)
        self.update_area_selector(**kwargs)
        self.canvas.figure.suptitle(self.params.title)
        self.canvas.draw()
        self.save(**kwargs)
        print(self.params.__dict__)
        self.finishedUpdating.emit()
        return


    def save(self, **kwargs):
        savekey = "save_to_file"
        if savekey in kwargs.keys():
            filePath = kwargs[savekey]
            extension = os.path.splitext(filePath)[-1]
            if extension in [".png", ".pdf"]:
                self.canvas.figure.savefig(filePath)
            elif extension in [".txt"]:
                np.savetxt(filePath,self.data.Z)
            else:
                raise NotImplementedError
        return


    def build_norm(self, **kwargs):
        if self.params.log_scale:
            #self.take_care_of_negative_values()
            #self.norm = mpl.colors.LogNorm(vmin=self.params.zmin, vmax=self.params.zmax)
            thres = np.abs(self.data.Z.std()/1e8)
            self.norm = mpl.colors.SymLogNorm(vmin=self.params.zmin, vmax=self.params.zmax, linthresh=thres)
        else:
            self.norm = mpl.colors.Normalize(vmin=self.params.zmin, vmax=self.params.zmax)
        return #build_norm


    def take_care_of_negative_values(self):
        if self.params.zmin <= 0:
            self.params.zmin = np.abs(self.data.Z.mean()/1e3)
        if self.params.zmax <= 0:
            self.params.zmax = 1+np.abs(self.data.Z.mean()/1e3)
        return #take_care_of_negative_values


    def update_ax(self, **kwargs):
        self.ax_imshow = self.ax.imshow(self.data.Z, norm=self.norm, vmin=self.norm.vmin, vmax=self.norm.vmax)
        self.cont_x = self.ax.contour(self.data.X,colors='k', linestyles='solid')
        self.cont_y = self.ax.contour(self.data.Y,colors='k', linestyles='solid')
        return #update_ax


    def update_area_selector(self, **kwargs):
        self.area_selector.rs.to_draw.set_visible(True)
        self.area_selector.rs.extents = (self.params.x1, self.params.x2, self.params.y1, self.params.y2)
        #self.area_selector.rs.update()
        return #update_area_selector


    def update_cax(self):
        self.build_cbar()
        return #update_cax


    def update_zoom_ax(self):
        x1, x2 = self.params.x1, self.params.x2
        y1, y2 = self.params.y1, self.params.y2
        self.data.Zzoom = self.data.Z[y1:y2,x1:x2]

        self.params.zoom_extent = (self.data.X[y1,x1], self.data.X[y2,x2], self.data.Y[y1,x1], self.data.Y[y2,x2])
        self.zoom_ax_imshow = self.zoom_ax.imshow(self.data.Zzoom, norm=self.norm, vmin=self.norm.vmin, vmax=self.norm.vmax,
                                            extent=self.params.zoom_extent)
        self.zoom_ax.set_aspect("auto")
        self.zoom_ax.set_xticks([])
        self.zoom_ax.set_yticks([])
        return #update_zoom_ax


    def update_xax(self):
        if self.params.log_scale:
            self.xax.set_yscale('log')

        integration_x = self.data.Zzoom.sum(axis=0)
        x0, xf = self.params.zoom_extent[0:2]
        rangex = np.linspace(x0, xf, len(integration_x))
        self.xax_line = self.xax.plot(rangex, integration_x)
        self.xax.set_xlim((x0, xf))

        self.xax.xaxis.set_ticks(np.linspace(x0, xf, 5))

        zero = integration_x.min()
        mu =  integration_x.mean()
        sig = integration_x.std()
        self.xax.set_yticks([zero, mu, mu+2*sig])
        self.xax.yaxis.tick_right()
        self.xax.grid(which='both', axis='both')
        return #update_xax


    def update_yax(self):
        if self.params.log_scale:
            self.yax.set_xscale('log')

        integration_y =  self.data.Zzoom.sum(axis=1)
        y0, yf = self.params.zoom_extent[2:4]
        rangey = np.linspace(y0, yf, len(integration_y))
        self.yax_line = self.yax.plot(integration_y, rangey)

        self.yax.set_ylim((yf, y0))
        self.yax.set_yticks(np.linspace(yf,y0,5))
        zero = integration_y.min()
        mu =  integration_y.mean()
        sig = integration_y.std()
        self.yax.set_xticks([zero, mu, mu+2*sig])
        self.yax.yaxis.tick_right()
        self.yax.xaxis.tick_top()
        self.yax.tick_params(axis='x', labelrotation=270)
        self.yax.grid(which='both', axis='both')
        return #update_yax


    def build_cbar(self):
        self.cbar = self.canvas.figure.colorbar(self.ax_imshow, cax=self.cax, orientation='vertical', norm = self.norm)
        self.cax.tick_params(axis='y', direction='in')
        self.cax.yaxis.tick_left()
        self.cbar.mappable.set_clim(self.norm.vmin, self.norm.vmax)
        return #build_cbar


    def test_show(self):
        t = np.linspace(-np.pi,np.pi, 1025)
        y = t#np.sin(t)
        x = t#np.cos(t)
        X, Y = np.meshgrid(x,y)
        Z = np.sin(Y) * np.cos(X)
        self.update_graph(X = X, Y = Y, Z = Z)
        if _DEBUG_:
            np.save("./myNumpyArray.npy", 3 + 10*np.sin(np.sqrt(X**2 + Y**2)))
            np.savetxt("./myNumpyArray.txt", 3 + 10*np.sin(np.sqrt(X**2 + Y**2)))
        return

class Experiment(FrozenClass):
    def __init__(self):
        self.selector_lambda = None
        self.angle_of_incidence = 0

        # Define the position of the direct beam on the detector
        # and the sensitivity map file
        self.qyc = 528
        self.qzc = 211
        self.pix0 = 162
        self.pixf = 862
        self.sens = None
        self.meansens = None
        self.monitor_counts = None
        self.qy = np.asarray([])
        self.qz = np.asarray([])
        self.I = np.asarray([])
        self.inputd = np.asarray([])

        self.Imatrix = np.asarray([])
        self.qymatrix = np.asarray([])
        self.qzmatrix = np.asarray([])

        self.cut_Iz = np.asarray([])
        self.cut_Iy = np.asarray([])

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
        self.dataDirPath = None
        self.datFileName = None
        self.yamlFileName = None
        self.gzFileName = None
        self.sensFileName = "sensitivity_map"
        self._freeze()
        return


    def datFilePath(self):
        if self.dataDirPath is None:
            return None
        return os.path.join(self.dataDirPath,self.datFileName)


    def yamlFilePath(self):
        if self.dataDirPath is None:
            return None
        return os.path.join(self.dataDirPath,self.yamlFileName)


    def gzFilePath(self):
        if self.dataDirPath is None:
            return None
        return os.path.join(self.dataDirPath,self.gzFileName)


    def sensFilePath(self):
        if self.dataDirPath is None:
            return None
        return os.path.join(self.dataDirPath, self.sensFileName)


    def basename(self):
        if self.dataDirPath is None:
            return None
        return os.path.splitext(self.datFilePath())[0]


    def gisans_map_filepath(self):
        if self.dataDirPath is None:
            return None
        return self.basename()+"_GISANS.map"


    def gisans_cut_filepath(self, y_or_z = "z"):
        if self.dataDirPath is None:
            return None
        return os.path.join(self.basename(),f"_line_cut_q{y_or_z}.out")


class App(qtw.QMainWindow,FrozenClass):
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
        pop_up = qtw.QMessageBox()
        pop_up.setWindowTitle(f"Exception: {e}\n")
        pop_up.setText(msg)
        pop_up.setIcon(qtw.QMessageBox.Critical)
        x = pop_up.exec_()
        return



class MyTabs(qtw.QTabWidget,FrozenClass):
    def __init__(self):
        super().__init__()
        self.tabButton_add = qtw.QToolButton()
        self.tabButton_rmv = qtw.QToolButton()
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



    @pyqtSlot()
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


    @pyqtSlot()
    def removeTab(self):
        idx = self.currentIndex()
        del self.frameList[idx]
        super().removeTab(idx)
        if len(self.frameList) == 0:
            self.last_num = 0
            self.addTab()

        return






class MyFrame(qtw.QFrame,FrozenClass):

    def __init__(self):
        super().__init__()
        self.layout = qtw.QHBoxLayout()
        self.centralpanel = qtw.QVBoxLayout()
        self.leftpanel = qtw.QVBoxLayout()
        self.rightpanel = qtw.QVBoxLayout()
        self.minSpinBox = mySciSpinBox()
        self.maxSpinBox = mySciSpinBox()
        self.settings = Settings()
        self.experiment = Experiment()
        self.graphView = MyGraphView("Title")
        self.infoTable = qtw.QTableWidget()
        self.fileList = qtw.QListWidget()
        self.tabs = qtw.QTabWidget()
        self.initFrame()
        self._freeze()


    def initFrame(self):
        self.layout.setAlignment(Qt.AlignCenter)
        self.addExperimentInfo()
        self.addFileList()
        self.addWelcomeMessage()
        self.addMinMaxSpinBoxes()
        self.addFunctionalityButtons()
        self.addCanvas()
        self.addPanels()
        self.setLayout(self.layout)


    def addFileList(self):
        self.fileList.setMaximumWidth(self.fileList.width()/2.)
        self.leftpanel.addWidget(qtw.QLabel("File:"))
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
        self.rightpanel.addWidget(qtw.QLabel("Info:"))
        self.rightpanel.addWidget(self.infoTable)
        return


    def addFunctionalityButtons(self):
        buttonOpenDialog = qtw.QPushButton("Load data")
        buttonOpenDialog.clicked.connect(self.on_click_open_file)
        self.leftpanel.addWidget(buttonOpenDialog)

        buttonLogLinear = qtw.QPushButton("Log / Linear")
        buttonLogLinear .clicked.connect(self.on_click_loglinear)
        self.rightpanel.addWidget(buttonLogLinear)

        buttonSavePng = qtw.QPushButton("Save png or pdf")
        buttonSavePng.clicked.connect(self.on_click_save_png)
        self.rightpanel.addWidget(buttonSavePng)

        buttonSaveAscii = qtw.QPushButton("Save ascii")
        buttonSaveAscii.clicked.connect(self.on_click_save_ascii)
        self.rightpanel.addWidget(buttonSaveAscii)

        return


    @pyqtSlot()
    def on_click_open_file(self):
        try:
            self.settings = Settings()
            self.doStuff()
            sum1 = self.experiment.Imatrix.sum()
            sum2 = self.experiment.cut_Iy.sum()
            sum3 = self.experiment.cut_Iz.sum()
            print("If the following three sums are not equal, there's an error:")
            print(sum1, sum2, sum3)
        except Exception as e:
            App.handle_exception(e)



    @pyqtSlot()
    def on_click_loglinear(self):
        try:
            self.graphView.update_graph(log_scale = not self.graphView.params.log_scale, reset_limits_required=True)
        except Exception as e:
            App.handle_exception(e)
        return


    @pyqtSlot()
    def on_click_save_png(self):
        try:

            fmt_choices = {"All Files(*)":".png", #default
                            "png (*.png)":".png",
                            "pdf (*.pdf)": ".pdf",
                            "ascii (*.txt)": ".txt"}
            choices_str = ";;".join([]+[k for k in fmt_choices.keys()])
            options = qtw.QFileDialog.Options()
            options |= qtw.QFileDialog.DontUseNativeDialog
            filePath, fmtChoice = qtw.QFileDialog.getSaveFileName(self,"Save File", "",
                                                          choices_str, options=options)
            if not filePath:
                return None

            extension = os.path.splitext(filePath)[-1]
            if extension not in fmt_choices.values():
                extension = fmt_choices[fmtChoice]
                filePath+=extension

            self.graphView.update_graph(save_to_file=filePath)
            print(f"Figure saved: {filePath}")
        except Exception as e:
            App.handle_exception(e)
        return #on_click_save_png


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
    def on_spinbox_edit(self):
        try:
            self.graphView.update_graph(zmax=self.maxSpinBox.value(), zmin=self.minSpinBox.value())
        except Exception as e:
            App.handle_exception(e)
        return


    @pyqtSlot()
    def on_graph_updated(self):
        try:
            self.minSpinBox.setValue(self.graphView.params.zmin)
            self.maxSpinBox.setValue(self.graphView.params.zmax)
        except Exception as e:
            App.handle_exception(e)
        return


    @staticmethod
    def init_spinbox(spinbox, slot):
        spinbox.editingFinished.connect(slot)
        return


    def addMinMaxSpinBoxes(self):
        self.init_spinbox(self.minSpinBox, self.on_spinbox_edit)
        self.init_spinbox(self.maxSpinBox, self.on_spinbox_edit)
        formLayout = qtw.QFormLayout()
        formLayout.addRow(self.tr("&Min Intensity"), self.minSpinBox)
        formLayout.addRow(self.tr("&Max Intensity"), self.maxSpinBox)
        formLayout.setFormAlignment(Qt.AlignBottom)
        self.rightpanel.addLayout(formLayout)
        return


    def addWelcomeMessage(self):
        message = qtw.QLabel(
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
        self.centralpanel.addWidget(self.graphView)
        self.graphView.finishedUpdating.connect(self.on_graph_updated)
        self.graphView.update_graph()
        return

    def saveFileNameDialog(self):
        try:
            options = qtw.QFileDialog.Options()
            options |= qtw.QFileDialog.DontUseNativeDialog
            fileName, _ = qtw.QFileDialog.getSaveFileName(self,"Save File", "","All Files (*);;png file (*.png);;pdf file (*.pdf)", options=options)
            if fileName:
                return fileName
        except Exception as e:
            App.handle_exception(e)
        return None

    def openFileNameDialog(self):
        try:
            options = qtw.QFileDialog.Options()
            options |= qtw.QFileDialog.DontUseNativeDialog
            fileName, _ = qtw.QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "",
            "Measurement dat file (*.dat);;All Files (*)",
            options=options)
            # self.openFileNamesDialog()
            if fileName:
                return fileName
        except Exception as e:
            App.handle_exception(e)
        return None





    def update_widgets(self):
        x = self.settings.datFileName
        y = self.experiment.__dict__

        for i, k in enumerate(y.keys()):
            if k[0] == "_":
                continue
            item_k = qtw.QTableWidgetItem(str(k))
            item_v = qtw.QTableWidgetItem(str(y[k]))
            self.infoTable.insertRow(i)
            self.infoTable.setVerticalHeaderItem(i,item_k)
            self.infoTable.setItem(i,0,item_v)

        self.minSpinBox.setValue(self.graphView.params.zmin)
        self.maxSpinBox.setValue(self.graphView.params.zmax)

        return True


    def doStuff(self):
        #self.graphView.test()
        #return
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
        if not self.update_gui():
            return
        return


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
        if _DEBUG_: 
            datFilePath = os.path.join(".","notToVersion","Archive","p15347_00001341.dat")
        else:
            datFilePath = self.openFileNameDialog()

        if datFilePath:
            path, filename = os.path.split(datFilePath)
            self.settings.datFileName = filename
            self.settings.dataDirPath = path
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
        pix0 = self.experiment.pix0
        pixf = self.experiment.pixf
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
        pix0 = self.experiment.pix0
        pixf = self.experiment.pixf
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


    def update_gui(self):
        try:
            self.graphView.update_graph(
                            Y = self.experiment.qymatrix,
                            X = self.experiment.qzmatrix,
                            Z = self.experiment.Imatrix,
                            reset_limits_required=True
                            )
            self.update_widgets()
        except Exception as e:
            App.handle_exception(e)
            return False

        return True


if __name__ == '__main__':
    enable_high_dpi_scaling()
    app = qtw.QApplication(sys.argv)
    ex = App()
    #sys.exit(app.exec_())
    x = profile_function_with_arguments(app.exec_)
    sys.exit(x)