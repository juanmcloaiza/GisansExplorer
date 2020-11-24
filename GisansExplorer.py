#!/bin/env python

#Qt stuff:
import PyQt5.QtWidgets as qtw
from PyQt5.QtCore import Qt, pyqtSlot, pyqtSignal, QThread
from PyQt5.QtGui import QValidator, QColor, QIcon

#plot stuff:
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.widgets import RectangleSelector
from matplotlib.figure import Figure

import numpy as np
import copy
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

class MyThread(QThread):
    progress_signal = pyqtSignal(int)
    def __init__(self, myframe):
        super().__init__()
        self.frame = myframe
        self.retval = None

    def run(self):
        myframe = self.frame
        myframe.settings = Settings()
        myframe.experiment = Experiment()

        self.progress_signal.emit(0)
        if not myframe.read_dat_file():
            self.retval = False, "dat file not read"
            return

        self.progress_signal.emit(10)
        if not myframe.read_yaml_file():
            self.retval = False, "yaml file not read"
            return

        self.progress_signal.emit(25)
        if not myframe.read_sensitivity_file():
            self.retval = False, "Sensitivity file not read"
            return

        self.progress_signal.emit(75)
        if not myframe.read_intensity_file():
            self.retval = False, "Intensity file not read"
            return

        self.progress_signal.emit(90)
        self.retval = True, ""
        return


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
        self.zoom_extent = (0., 1., 0., 1.)

        self.Xc = 0
        self.Xs = 0
        self.Yc = 0

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

        self.log_scale = True
        self.reset_limits_required = True

        self.title = ""
        self._freeze()
        return


class MyGraphView(qtw.QWidget):
    finishedUpdating = pyqtSignal()
    def __init__(self, graph_title, parent = None):
        super(MyGraphView, self).__init__(parent)

        self.graph_title = graph_title

        self.dpi = 100
        self.fig = Figure((10.0, 5.0), dpi = self.dpi, facecolor = (1,1,1), edgecolor = (0,0,0), linewidth=1)
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
        self.layout.addWidget(self.canvas)
        self.layout.addWidget(self.xyzLabel)
        self.layout.setStretchFactor(self.canvas, 1)
        self.setLayout(self.layout)
        return #define_layout


    def define_axes(self):
        #[left, bottom, width, height]
        self.ax =  self.canvas.figure.add_axes([0.1,0.2,0.3,0.5])
        self.cax = self.canvas.figure.add_axes([0.1,0.71,0.3,0.03])
        #self.cax = self.canvas.figure.add_axes([0.55,0.76,0.25,0.025])
        self.zoom_ax = self.canvas.figure.add_axes([0.5,0.2,0.3,0.5])
        self.xax = self.canvas.figure.add_axes([0.5,0.71,0.3,0.1])
        self.yax = self.canvas.figure.add_axes([0.81,0.2,0.05,0.5])
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
        x = np.random.normal(1024,1024)
        self.data.X = x
        self.data.Y = x
        self.data.x1 = 256
        self.data.x2 = 1024-256
        self.data.y1 = 256
        self.data.y2 = 1024-256
        return #init_data_and_parameters


    def on_mouse_wheel(self, event):
        if self.cax == event.inaxes:
            if event.button == 'up':
                func = lambda c: (np.sign(c) * abs(c)**0.9) if self.data.log_scale else 0.9*c
            elif event.button == 'down':
                func = lambda c: (np.sign(c) * abs(c)**1.1) if self.data.log_scale else 1.1*c
            else:
                return

            vmin, vmax = (func(c) for c in self.cbar.get_clim())
            #print("Rescaling colorbar:", vmin,vmax)
            self.update_graph(zmin=vmin, zmax=vmax)
        return #on_mouse_wheel


    def on_mouse_move(self, event):
        if not event.inaxes:
            return
        xd, yd = event.xdata, event.ydata
        xarr = self.data.X[0,:]
        yarr = self.data.Y[:,0]
        if event.inaxes == self.zoom_ax:
            col = np.searchsorted(xarr, xd)-1
            row = np.searchsorted(yarr, yd)-1
        else:
            row, col = int(yd + 0.5), int(xd + 0.5)

        zd = self.data.Z[row, col]
        coord_text = f'x={xd:1.4g}, y={yd:1.4g}, I={zd:1.4e}   [{row},{col}]'
        self.xyzLabel.setText(coord_text)
        #self.xyzLabel.setText(f"(x, y; z) = ({xd:3.2g}, {yd:3.2g}; {z:3.2g})")

        return #on_mouse_move


    def line_select_callback(self, eclick, erelease):
        x1, y1 = int(eclick.xdata), int(eclick.ydata)
        x2, y2 = int(erelease.xdata), int(erelease.ydata)
        #print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
        #print(" The button you used were: %s %s" % (eclick.button, erelease.button))
        self.update_graph(x1=x1, x2=x2, y1=y1, y2=y2)
        return #line_select_callback


    def update_data(self, **kwargs):
        for k in self.data.__dict__.keys():
            if k in kwargs.keys():
                self.data.__setattr__(k,kwargs[k])

        if self.data.reset_limits_required:
            self.data.zmin = 0
            self.data.zmax = 1e-3
            if self.data.log_scale:
                self.data.zmin = 1e-6
            self.data.reset_limits_required = False
        return #update_data


    def update_axes(self, **kwargs):
        self.update_ax(**kwargs)
        self.update_cax()
        self.update_zoom_ax()
        self.update_xax()
        self.update_yax()
        return #update_axes


    def update_graph(self, **kwargs):
        try:
            self.canvas.figure.clear()
            self.define_axes()
            self.update_data(**kwargs)
            self.build_norm(**kwargs)
            self.update_axes(**kwargs)
            self.update_area_selector(**kwargs)
            self.canvas.figure.suptitle(self.data.title, fontsize=5)
            self.canvas.draw()
            self.save(**kwargs)
            self.finishedUpdating.emit()
        except Exception as e:
            App.handle_exception(e)
        return


    def save(self, **kwargs):
        savekey = "save_to_file"
        if savekey in kwargs.keys():
            filePath = kwargs[savekey]
            extension = os.path.splitext(filePath)[-1]
            if extension in [".png", ".pdf"]:
                self.canvas.figure.savefig(filePath)
                self.save_figures(filePath)
            elif extension in [".txt"]:
                header = kwargs["header"]
                np.savetxt(filePath,self.data.Z, header=header,fmt='%.3e')
            else:
                raise NotImplementedError
            print(f"Figure saved: {filePath}")
        return


    def build_norm(self, **kwargs):
        if self.data.log_scale:
            thres = np.abs(self.data.Z.std()/1e8)
            self.norm = mpl.colors.SymLogNorm(vmin=self.data.zmin, vmax=self.data.zmax, linthresh=thres, base=10)
            #self.norm = mpl.colors.LogNorm(vmin=self.data.zmin, vmax=self.data.zmax)
        else:
            self.norm = mpl.colors.Normalize(vmin=self.data.zmin, vmax=self.data.zmax)
        return #build_norm


    def take_care_of_negative_values(self):
        if self.data.zmin <= 0:
            self.data.zmin = np.abs(self.data.Z.mean()/1e3)
        if self.data.zmax <= 0:
            self.data.zmax = 1+np.abs(self.data.Z.mean()/1e3)
        return #take_care_of_negative_values


    def update_ax(self, **kwargs):
        #self.cont_x = self.ax.contour(self.data.X, [0.], colors='k', linestyles='solid', label = "Detector center")#, linewidths=0.5)
        #self.cont_y = self.ax.contour(self.data.Y, [0.], colors='k', linestyles='solid')#, linewidths=0.5)
        self.ax_imshow = self.ax.pcolorfast(self.data.Z, norm=self.norm, vmin=self.norm.vmin, vmax=self.norm.vmax, cmap='jet')
        X_is_zero_at_idx = np.abs(self.data.X[0]).argmin()
        Y_is_zero_at_idx = np.abs(self.data.Y).T[0].argmin()
        self.ax.axvline(x=self.data.Xc, c='r', label = "Detector center (corrected)")
        self.ax.axvline(x=self.data.Xs, c='g', label = "Specular beam")
        self.ax.axvline(x=X_is_zero_at_idx, c='k', label = "Detector center")
        self.ax.axhline(y=self.data.Yc, c='r')
        self.ax.axhline(y=Y_is_zero_at_idx, c='k')
        self.ax.set_aspect("auto")
        self.ax.set_title("Detector View", pad=50)
        self.ax.set_xlim(150,870)
        self.ax.set_ylim(150,870)
        self.ax.legend(loc='upper left', bbox_to_anchor= (-0.3, -0.3), ncol=3,
            borderaxespad=0, frameon=False, fontsize = 6)
        return #update_ax


    def update_area_selector(self, **kwargs):
        self.area_selector.rs.to_draw.set_visible(True)
        self.area_selector.rs.extents = (self.data.x1, self.data.x2, self.data.y1, self.data.y2)
        #self.area_selector.rs.update()
        return #update_area_selector


    def update_cax(self):
        self.build_cbar()
        return #update_cax


    def update_zoom_ax(self):
        x1, x2 = self.data.x1, self.data.x2
        y1, y2 = self.data.y1, self.data.y2
        xc, yc, xs = self.data.Xc, self.data.Yc, self.data.Xs
        x0 = np.abs(self.data.X[0]).argmin()
        y0 = np.abs(self.data.Y).T[0].argmin()
        self.data.Xzoom = self.data.X[y1:y2+1,x1:x2+1]
        self.data.Yzoom = self.data.Y[y1:y2+1,x1:x2+1]
        self.data.Zzoom = self.data.Z[y1:y2,x1:x2]

        self.data.zoom_extent = (self.data.X[y1,x1], self.data.X[y2,x2], self.data.Y[y1,x1], self.data.Y[y2,x2])
        self.zoom_ax_imshow = self.zoom_ax.pcolorfast(self.data.Xzoom, self.data.Yzoom, self.data.Zzoom, norm=self.norm, vmin=self.norm.vmin, vmax=self.norm.vmax, cmap='jet')

        if x0 > x1 and x0 < x2:
            self.zoom_ax.axvline(x=0, c='k', ls='solid')#, lw=0.5)
        if xc > x1 and xc < x2:
            self.zoom_ax.axvline(x=self.data.X[yc,xc], c='r')
        if xs > x1 and xs < x2:
            self.zoom_ax.axvline(x=self.data.X[yc,xs], c='g')

        if yc > y1 and yc < y2:
            self.zoom_ax.axhline(y=self.data.Y[yc,xc], c='r')
        if y0 > y1 and y0 < y2:
            self.zoom_ax.axhline(y=0, c='k', ls='solid')#, lw=0.5)


        self.zoom_ax.set_aspect("auto")
        self.zoom_ax.set_xticks([])
        self.zoom_ax.set_yticks([])
        self.zoom_ax.set_xlabel("$Q_{z}$")
        self.zoom_ax.set_ylabel("$Q_{y}$")

        return #update_zoom_ax


    def update_xax(self):
        if self.data.log_scale:
            try:
                self.xax.set_yscale('log')
            except Exception:
                pass

        integration_x = self.data.Zzoom.sum(axis=0)
        x0, xf = self.data.zoom_extent[0:2]
        rangex = np.linspace(x0, xf, len(integration_x))
        self.xax_line = self.xax.plot(rangex, integration_x)
        self.xax.set_xlim((x0, xf))


        zero = integration_x.min()
        mu =  integration_x.mean()
        sig = integration_x.std()
        self.xax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        self.xax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        self.xax.set_xticks(np.linspace(x0, xf, 5))
        self.xax.set_yticks([zero, mu, mu+2*sig])
        self.xax.yaxis.tick_right()
        self.xax.grid(which='both', axis='both')
        self.xax.tick_params(axis='x', labelrotation=270)
        self.xax.xaxis.tick_top()
        return #update_xax


    def update_yax(self):
        if self.data.log_scale:
            try:
                self.yax.set_xscale('log')
            except Exception:
                pass

        integration_y =  self.data.Zzoom.sum(axis=1)
        y0, yf = self.data.zoom_extent[2:4]
        rangey = np.linspace(y0, yf, len(integration_y))
        self.yax_line = self.yax.plot(integration_y, rangey)

        self.yax.set_ylim((y0, yf))
        zero = integration_y.min()
        mu =  integration_y.mean()
        sig = integration_y.std()
        self.yax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        self.yax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        self.yax.set_xticks([zero, mu, mu+2*sig])
        self.yax.set_yticks(np.linspace(y0,yf,5))
        self.yax.yaxis.tick_right()
        self.yax.tick_params(axis='x', labelrotation=270)
        self.yax.grid(which='both', axis='both')
        return #update_yax


    def build_cbar(self):
        self.cax.tick_params(axis='x', direction='in', labeltop=True, top=True)
        self.cbar = self.canvas.figure.colorbar(self.ax_imshow, cax=self.cax, orientation='horizontal', norm = self.norm)
        self.cbar.mappable.set_clim(self.norm.vmin, self.norm.vmax)
        self.cax.xaxis.tick_top()
        return #build_cbar

    def save_figures(self,filePath):
        self.save_gisans_map(filePath)
        self.save_qz_integration(filePath)
        self.save_qy_integration(filePath)
        return

    def save_gisans_map(self,filePath):
        print("Saving gisans map....")
        x1, x2 = self.data.x1, self.data.x2
        y1, y2 = self.data.y1, self.data.y2
        x0 = np.abs(self.data.X[0]).argmin()
        y0 = np.abs(self.data.Y).T[0].argmin()
        xc, yc, xs = self.data.Xc, self.data.Yc, self.data.Xs
        self.data.Xzoom = self.data.X[y1:y2+1,x1:x2+1]
        self.data.Yzoom = self.data.Y[y1:y2+1,x1:x2+1]
        self.data.Zzoom = self.data.Z[y1:y2,x1:x2]

        self.data.zoom_extent = (self.data.X[y1,x1], self.data.X[y2,x2], self.data.Y[y1,x1], self.data.Y[y2,x2])

        new_fig = plt.figure(figsize=(10,10))
        new_ax = new_fig.add_subplot(111)
        cs = new_ax.pcolorfast(self.data.Xzoom, self.data.Yzoom, self.data.Zzoom, norm=self.norm, vmin=self.norm.vmin, vmax=self.norm.vmax, cmap='jet')

        has_legend = False
        if x0 > x1 and x0 < x2:
            new_ax.axvline(x=0, c='k', ls='solid', label = "Detector center")
        if xc > x1 and xc < x2:
            new_ax.axvline(x=self.data.X[yc,xc], c='r', label = "Detector center (corrected)")
        if xs > x1 and xs < x2:
            new_ax.axvline(x=self.data.X[yc,xs], c='g', label = "Specular beam")

        if yc > y1 and yc < y2:
            new_ax.axhline(y=self.data.Y[yc,xc], c='r')
        if y0 > y1 and y0 < y2:
            new_ax.axhline(y=0, c='k', ls='solid')#, lw=0.5)

        new_ax.set_aspect("auto")
        new_ax.set_xlabel("$Q_{z}$")
        new_ax.set_ylabel("$Q_{y}$")
        cbar = new_fig.colorbar(cs)
        new_ax.set_title(self.data.title, fontsize=10)

        new_ax.legend()#loc='upper left', bbox_to_anchor= (-0.3, -0.3), ncol=3,
            #borderaxespad=0, frameon=False, fontsize = 6)

        no_ext, ext = os.path.splitext(filePath)
        new_fig.savefig(f"{no_ext}-gisans_map{ext}")
        print(f"gisans map saved.")
        return #save_gisans_map#


    def save_qz_integration(self, filePath):
        new_fig = plt.figure(figsize=(10,10))
        new_ax = new_fig.add_subplot(111)
        print("Saving qz integration....")
        if self.data.log_scale:
            try:
                new_ax.set_yscale('log')
            except Exception:
                pass

        integration_x = self.data.Zzoom.sum(axis=0)
        x0, xf = self.data.zoom_extent[0:2]
        rangex = np.linspace(x0, xf, len(integration_x))
        xax_line = new_ax.plot(rangex, integration_x)
        new_ax.set_xlim((x0, xf))

        new_ax.xaxis.set_ticks(np.linspace(x0, xf, 5))

        zero = integration_x.min()
        mu =  integration_x.mean()
        sig = integration_x.std()
        #new_ax.set_yticks([zero, mu, mu+2*sig])
        new_ax.grid(which='both', axis='both')
        new_ax.set_xlabel("$Q_{z}$")
        new_ax.set_ylabel("I($Q_{z})$")
        new_ax.set_title(self.data.title, fontsize=10)

        no_ext, ext = os.path.splitext(filePath)
        new_fig.savefig(f"{no_ext}-integration_qz{ext}")
        print("qz integration saved.")
        return #save_xax


    def save_qy_integration(self, filePath):
        new_fig = plt.figure(figsize=(10,10))
        new_ax = new_fig.add_subplot(111)
        if self.data.log_scale:
            try:
                new_ax.set_yscale('log')
            except Exception:
                pass

        integration_y =  self.data.Zzoom.sum(axis=1)
        y0, yf = self.data.zoom_extent[2:4]
        rangey = np.linspace(y0, yf, len(integration_y))
        new_ax_line = new_ax.plot(rangey, integration_y)

        new_ax.set_xlim((y0, yf))
        new_ax.set_xticks(np.linspace(y0,yf,5))
        zero = integration_y.min()
        mu =  integration_y.mean()
        sig = integration_y.std()
        #new_ax.set_yticks([zero, mu, mu+2*sig])
        new_ax.grid(which='both', axis='both')
        new_ax.set_xlabel("$Q_{y}$")
        new_ax.set_ylabel("I($Q_{y})$")
        new_ax.set_title(self.data.title, fontsize=10)

        no_ext, ext = os.path.splitext(filePath)
        new_fig.savefig(f"{no_ext}-integration_qy{ext}")
        print("qy integration saved.")
        return #update_yax


    def test_show(self):
        #Z = np.loadtxt('./notToVersion/WelcomePlot.txt')
        #Z = (Z - Z.min()) / (Z.max() - Z.min()) * (1e-3 - 1e-6) + 1e-6
        #np.save("./show_test_map.npy", Z)
        Z = np.load('./show_test_map.npy')

        #This creates a cross to test the correct alignment between the
        #different plot axes:
        Z[:,:] = 1e-5
        Z[500:511,:] = 1e-4
        Z[:,500:511] = 1e-4

        print(Z.min())
        print(Z.max())
        print(Z.shape)

        x = np.linspace(-1,1, Z.shape[0])
        y = np.linspace(-1,1, Z.shape[1])
        X, Y = np.meshgrid(x,y)
        Xc = len(x) // 2
        Yc = len(y) //2
        self.update_graph(X = X, Y = Y, Z = Z, Xc = 0, Yc = 0)
        #if _DEBUG_:
        #    np.savetxt("./myNumpyArray.txt", 3 + 10*np.sin(np.sqrt(X**2 + Y**2)))
        return

class Experiment(FrozenClass):
    def __init__(self):
        self.selector_lambda = None
        self.angle_of_incidence = 0
        self.sens = None
        self.meansens = None
        self.monitor_counts = None

        # Define the position of the direct beam on the detector
        # and the sensitivity map file
        self.qyc = 0
        self.qzc = 0
        self.qzc_corr = 0
        self.qzc_spec = 0
        self.x0 = 128
        self.y0 = 128
        self.xf = 256
        self.yf = 256
        self.min_intensity = 1e-6
        self.max_intensity = 1e-3
        #self.qy = np.asarray([])
        #self.qz = np.asarray([])
        #self.I = np.asarray([])
        #self.inputd = np.asarray([])

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
        return 0.5755*(pixel_j-self.qzc_corr)/np.sqrt(1990*1990+(0.5755*(self.qzc_corr-pixel_j))*(0.5755*(self.qzc_corr-pixel_j)))


    def cos_alpha_f(self, pixel_j):
        return np.sqrt(1 - self.sin_alpha_f(pixel_j)**2)



class Settings(FrozenClass):
    def __init__(self):
        self.dataDirPath = None
        self.datFileName = None
        self.yamlFileNames = []
        self.gzFileNames = []
        self.sensFileName = "sensitivity_map"
        self._freeze()
        return


    def datFilePath(self):
        if self.dataDirPath is None:
            return None
        return os.path.join(self.dataDirPath,self.datFileName)


    def yamlFilePaths(self):
        if self.dataDirPath is None:
            return None
        return [os.path.join(self.dataDirPath,f) for f in self.yamlFileNames]


    def gzFilePaths(self):
        if self.dataDirPath is None:
            return None
        return [os.path.join(self.dataDirPath,f) for f in self.gzFileNames]


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

    @staticmethod
    def handle_exception(e):
        msg0= (f"Exception: {e}\n")
        msg = ("-"*60+"\n")
        msg += traceback.format_exc()
        msg += ("-"*60+"\n")
        print(msg0, msg)
        pop_up = qtw.QMessageBox()
        pop_up.setWindowTitle(f"Exception: {e}\n")
        pop_up.setText(msg0)
        pop_up.setInformativeText(msg)
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
        self.splitter = qtw.QSplitter()
        self.centralpanel = qtw.QVBoxLayout()
        self.botpanel = qtw.QHBoxLayout()
        self.rightpanel = qtw.QVBoxLayout()
        #self.minSpinBox = mySciSpinBox()
        #self.maxSpinBox = mySciSpinBox()
        self.settings = Settings()
        self.settings_dict = {}
        self.experiment = Experiment()
        self.experiment_dict = {}
        self.graphView = MyGraphView("Title")
        self.infoTable = qtw.QTableWidget()
        self.fileList = qtw.QListWidget()
        self.subtractCheckBox = qtw.QCheckBox("Subtract Intensities")
        self.dirtree = qtw.QTreeView()
        self.tabs = qtw.QTabWidget()
        self.progress_bar =qtw.QProgressBar()
        self.thread = MyThread(self)

        self.initFrame()
        self._freeze()


    def initFrame(self):
        self.layout.setAlignment(Qt.AlignCenter)
        self.addExperimentInfo()
        #self.addMinMaxSpinBoxes()
        self.addFunctionalityButtons()
        self.addPanels()
        self.setLayout(self.layout)
        self.progress_bar.setMaximum(100)
        self.thread.progress_signal.connect(self.on_progress_emited)
        self.thread.finished.connect(self.update_gui)


    def addFileTreeAndList(self, botSplitter):
        model = qtw.QFileSystemModel()
        model.setRootPath('')
        filters = ["*.dat"]
        model.setNameFilters(filters)
        model.setNameFilterDisables(False)

        self.dirtree.setModel(model)
        self.dirtree.setRootIndex(model.index('./'))

        self.dirtree.setAnimated(True)
        self.dirtree.setIndentation(20)
        self.dirtree.setSortingEnabled(True)
        self.dirtree.doubleClicked.connect(self.on_click_open_file)
        self.fileList.itemSelectionChanged.connect(self.on_file_selection_changed)
        self.fileList.setSelectionMode(3) #https://doc.qt.io/archives/qt-5.11/qabstractitemview.html#SelectionMode-enum

        botSplitter.setOrientation(Qt.Horizontal)
        botSplitter.addWidget(self.dirtree)

        fileListSplitter=qtw.QSplitter(Qt.Vertical)
        fileListSplitter.addWidget(self.fileList)
        fileListSplitter.addWidget(self.subtractCheckBox)
        self.subtractCheckBox.setChecked(False)
        self.subtractCheckBox.setEnabled(False)
        self.subtractCheckBox.stateChanged.connect(self.on_subtract_checkbox_changed)
        botSplitter.addWidget(fileListSplitter)
        return


    def addPanels(self):

        bottom_splitter = qtw.QSplitter()
        left_splitter = qtw.QSplitter(orientation=Qt.Vertical)
        self.addFileTreeAndList(bottom_splitter)
        left_splitter.addWidget(self.graphView)
        left_splitter.addWidget(bottom_splitter)
        left_splitter.addWidget(self.progress_bar)
        self.splitter.addWidget(left_splitter)
        rightlayoutwidget = qtw.QWidget()
        rightlayoutwidget.setLayout(self.rightpanel)
        self.splitter.addWidget(rightlayoutwidget)
        self.layout.addWidget(self.splitter)
        return


    def addExperimentInfo(self):
        #self.infoTable.setMaximumWidth(self.infoTable.width()/2.)
        self.infoTable.setColumnCount(1)
        self.infoTable.setRowCount(0)
        self.infoTable.horizontalHeader().hide()
        self.infoTable.horizontalHeader().setStretchLastSection(True)
        self.infoTable.cellChanged.connect(self.on_cell_changed)
        self.rightpanel.addWidget(qtw.QLabel("Info:"))
        self.rightpanel.addWidget(self.infoTable)
        return


    def addFunctionalityButtons(self):
        buttonUpdateFromTable = qtw.QPushButton("Update")
        buttonUpdateFromTable.clicked.connect(self.on_click_update)
        self.rightpanel.addWidget(buttonUpdateFromTable)

        buttonLogLinear = qtw.QPushButton("Log / Linear")
        buttonLogLinear .clicked.connect(self.on_click_loglinear)
        self.rightpanel.addWidget(buttonLogLinear)


        buttonSavePng = qtw.QPushButton("Save png or pdf")
        buttonSavePng.clicked.connect(self.on_click_save_png)
        self.rightpanel.addWidget(buttonSavePng)

        buttonSaveAscii = qtw.QPushButton("Save ROI as ascii columns")
        buttonSaveAscii.clicked.connect(self.on_click_save_ascii)
        self.rightpanel.addWidget(buttonSaveAscii)

        buttonOpenDialog = qtw.QPushButton("Open...")
        buttonOpenDialog.clicked.connect(self.on_click_open_file)
        self.rightpanel.addWidget(buttonOpenDialog)

        self.graphView.finishedUpdating.connect(self.on_graph_updated)

        return

    @pyqtSlot(int)
    def on_progress_emited(self, value):
        self.progress_bar.setValue(value)

    @pyqtSlot()
    def on_click_open_file(self):
        try:
            #self.settings = Settings()
            #self.experiment = Experiment()
            self.doStuff()
        except Exception as e:
            App.handle_exception(e)

    @pyqtSlot()
    def on_subtract_checkbox_changed(self):
        self.update_from_selection_list()

    @pyqtSlot()
    def on_file_selection_changed(self):
        selectedListEntries = self.fileList.selectedItems()

        if len(selectedListEntries) < 1:
            self.fileList.setCurrentRow( self.fileList.count() - 1 )

        self.subtractCheckBox.setChecked(False)
        if len(selectedListEntries) == 2:
            self.subtractCheckBox.setEnabled(True)
        else:
            self.subtractCheckBox.setEnabled(False)
        self.update_from_selection_list()

    @pyqtSlot()
    def on_click_update(self):
        self.update_from_info_table()

    @pyqtSlot()
    def on_cell_changed(self):
        self.color_outdated()


    @pyqtSlot()
    def on_click_loglinear(self):
        try:
            self.graphView.update_graph(log_scale = not self.graphView.data.log_scale, reset_limits_required=False)
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

            self.graphView.update_graph(save_to_file=filePath, header=self.build_ascii_header(), title=self.build_ascii_header())

        except Exception as e:
            App.handle_exception(e)
        return #on_click_save_png


    def build_ascii_header(self):
        gz_filenames = [s.text() for s in self.fileList.selectedItems()]
        dat_filenames= [self.settings_dict[s].datFileName for s in gz_filenames]
        op_str = ""
        if self.subtractCheckBox.isChecked():
            op_str = "Subtraction of:\n"
        else:
            op_str = "Addition of:\n"
        header = op_str + "\n".join([d + ": " + gz for (d,gz) in zip(dat_filenames, gz_filenames)])
        return header


    @pyqtSlot()
    def on_click_save_ascii(self):
        try:
            filepath = self.saveFileNameDialog()
            if filepath is None:
                return

            noextpath, ext = os.path.splitext(filepath)
            if ext == '': ext = ".txt"

            x1, x2 = self.graphView.data.x1, self.graphView.data.x2
            y1, y2 = self.graphView.data.y1, self.graphView.data.y2

            header = self.build_ascii_header()

            filename = noextpath+"_xyI_"+ext
            original_shape = self.graphView.data.X[y1:y2,x1:x2].shape
            x_to_save = self.graphView.data.X[y1:y2,x1:x2].flatten()
            y_to_save = self.graphView.data.Y[y1:y2,x1:x2].flatten()
            z_to_save = self.graphView.data.Z[y1:y2,x1:x2].flatten()
            columns = np.vstack((x_to_save, y_to_save, z_to_save)).T
            np.savetxt(filename, columns, header=header+f"\n#Original shape: {original_shape}",fmt='%.3e')

            filename = noextpath+"_xI_"+ext
            original_shape = self.graphView.data.X[0,x1:x2].shape
            x_to_save = self.graphView.data.X[0,x1:x2].flatten()
            z_to_save = self.graphView.data.Z[y1:y2,x1:x2].sum(axis=0).flatten()
            columns = np.vstack((x_to_save, z_to_save)).T
            np.savetxt(filename, columns, header=header+f"\n#Original shape: {original_shape}",fmt='%.3e')

            filename = noextpath+"_yI_"+ext
            original_shape = self.graphView.data.X[y1:y2,0].shape
            y_to_save = self.graphView.data.Y[y1:y2,0].flatten()
            z_to_save = self.graphView.data.Z[y1:y2,x1:x2].sum(axis=1).flatten()
            columns = np.vstack((y_to_save, z_to_save)).T
            np.savetxt(filename, columns, header=header+f"\n#Original shape: {original_shape}",fmt='%.3e')

            print(f"Arrays saved:\n {filename}\n")
        except Exception as e:
            App.handle_exception(e)


#     @pyqtSlot()
#     def on_spinbox_edit(self):
#         try:
#             self.graphView.update_graph(zmax=self.maxSpinBox.value(), zmin=self.minSpinBox.value())
#         except Exception as e:
#             App.handle_exception(e)
#         return


    @pyqtSlot()
    def on_graph_updated(self):
        try:
            selectedListEntries = self.fileList.selectedItems()

            for currentListItem in selectedListEntries:
                currentListEntry = currentListItem.text()
                self.experiment = self.experiment_dict[currentListEntry]
                self.experiment.min_intensity = self.graphView.data.zmin
                self.experiment.max_intensity = self.graphView.data.zmax
                self.experiment.x0 = self.graphView.data.x1
                self.experiment.y0 = self.graphView.data.y1
                self.experiment.xf = self.graphView.data.x2
                self.experiment.yf = self.graphView.data.y2

            self.update_table()
        except Exception as e:
            App.handle_exception(e)
        return


#    @staticmethod
#    def init_spinbox(spinbox, slot):
#        spinbox.editingFinished.connect(slot)
#        return


#     def addMinMaxSpinBoxes(self):
#         self.init_spinbox(self.minSpinBox, self.on_spinbox_edit)
#         self.init_spinbox(self.maxSpinBox, self.on_spinbox_edit)
#         formLayout = qtw.QFormLayout()
#         formLayout.addRow(self.tr("&Min Intensity"), self.minSpinBox)
#         formLayout.addRow(self.tr("&Max Intensity"), self.maxSpinBox)
#         formLayout.setFormAlignment(Qt.AlignBottom)
#         self.rightpanel.addLayout(formLayout)
#         return

    def saveFileNameDialog(self):
        try:
            fmt_choices = {"All Files(*)":".png", #default
                            "png (*.png)":".png",
                            "pdf (*.pdf)": ".pdf",
                            "ascii (*.txt)": ".txt"}
            choices_str = ";;".join([]+[k for k in fmt_choices.keys()])
            options = qtw.QFileDialog.Options()
            options |= qtw.QFileDialog.DontUseNativeDialog
            fileName, fmtChoice = qtw.QFileDialog.getSaveFileName(self,"Save File", "",
                                                          choices_str, options=options)
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


    @staticmethod
    def color_validate(table_item, value_a, value_b):
        green = QColor(000, 255, 0, 127)
        if value_a != value_b:
            table_item.setBackground(green)
            table_item.setSelected(False)
        else:
            table_item.setBackground(Qt.white)
        return True

    def update_single_experiment_values(self, experiment):
        expdict = experiment.__dict__
        for i in range(self.infoTable.rowCount()):
            key = self.infoTable.verticalHeaderItem(i).text()
            current_item = self.infoTable.item(i,0)
            if key in ['qyc', 'qzc', 'x0', 'y0', 'xf', 'yf']:
                value = int(current_item.text())
            elif key in ["min_intensity", "max_intensity"]:
                value = float(current_item.text())
            else:
                continue
            expdict[key] = value
        return True

    def update_multi_experiment_values(self):
        selectedListEntries = self.fileList.selectedItems()
        for currentListItem in selectedListEntries:
            currentListEntry = currentListItem.text()
            self.settings = self.settings_dict[currentListEntry]
            self.experiment = self.experiment_dict[currentListEntry]
            self.update_single_experiment_values(self.experiment)
        return True


    def color_outdated(self):
        expdict = self.experiment.__dict__

        try:
            for i in range(self.infoTable.rowCount()):
                key = self.infoTable.verticalHeaderItem(i).text()
                current_item = self.infoTable.item(i,0)
                try:
                    if key in ['qyc', 'qzc', 'x0', 'y0', 'xf', 'yf']:
                        value = int(current_item.text())
                    elif key in ["min_intensity", "max_intensity"]:
                        value = float(current_item.text())
                    else:
                        continue

                    self.color_validate(current_item, value, expdict[key])
                except Exception as e:
                    red = QColor(255, 000, 0, 127)
                    current_item.setBackground(red)
                    current_item.setSelected(False)
                    return False

        except Exception as e:
            App.handle_exception(e)
            return False
        return True


    def update_table(self):

        expdict = self.experiment.__dict__

        self.infoTable.setRowCount(0)
        for k in expdict.keys():
            if k[0] == "_":
                continue

            if k in ["selector_lambda", "qyc", "qzc", "x0", "y0", "xf", "yf",
                     "min_intensity", "max_intensity",
                     "meansens", "monitor_counts",
                     "angle_of_incidence", "qzc_corr", "qzc_spec"]:
                i = self.infoTable.rowCount()
                item_k = qtw.QTableWidgetItem(str(k))
                item_v = qtw.QTableWidgetItem(str(expdict[k]))
                self.infoTable.insertRow(i)
                self.infoTable.setVerticalHeaderItem(i,item_k)
                self.infoTable.setItem(i,0,item_v)

        return True


    def doStuff(self):
        self.thread.start()
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
            datFilePath = os.path.join(".","notToVersion","Archive","p15496_00000992.dat")
        else:
            if len(self.dirtree.selectedIndexes()) < 1:
                datFilePath = self.openFileNameDialog()
            else:
                datFilePath = self.dirtree.model().filePath(self.dirtree.currentIndex())

        if datFilePath:
            path, filename = os.path.split(datFilePath)
            self.settings.datFileName = filename
            self.settings.dataDirPath = path
            return self.safe_parse(self.parse_dat, self.settings.datFilePath())
        return False


    def read_sensitivity_file(self):
        func = self.parse_sensitivity_map
        fpath = self.settings.sensFilePath()
        for key in self.settings.gzFileNames:
            self.experiment = self.experiment_dict[key]
            tf = self.safe_parse_numpy(func, fpath, dtype=float, delimiter=' ')
            if not tf:
                return False
        return True


    def read_yaml_file(self):
        # Open and read the yaml file
        for yamlFileName, key in zip(self.settings.yamlFilePaths(), self.settings.gzFileNames):
            self.experiment = self.experiment_dict[key]
            tf = self.safe_parse(self.parse_yaml, yamlFileName)
            if not tf:
                return False
        return True


    def read_intensity_file(self):
        func = self.parse_intensity_map
        for fpath, key in zip(self.settings.gzFilePaths(),self.settings.gzFileNames):
            self.experiment = self.experiment_dict[key]
            tf = self.safe_parse_numpy(func, fpath, dtype=float, delimiter=' ')
            if not tf:
                return False
        return True


    def parse_dat(self, file):
        for line in file:
            if line.find('omega_value')>0:
                omega_line_list = line.split()
                omega = omega_line_list[3]
                self.experiment.angle_of_incidence = omega
                print('Angle of incidence (degrees): '+str(omega))
                print(f"Original qzc = {self.experiment.qzc}")
            if line.find('selector_lambda_value')>0:
                lambda_line_list = line.split()
                selector_lambda=lambda_line_list[3]
                self.experiment.selector_lambda = selector_lambda
                print('Neutron wavelength (Angtrom): '+str(selector_lambda))
            if line.find('.yaml')>0:
                line_list = line.split()
                for word in line_list:
                    if word.find('.yaml')>0:
                        yamlFileName = word
                        self.settings.yamlFileNames.append(yamlFileName)
                        print(f"Found yaml: {yamlFileName}")
            if line.find('.gz')>0:
                line_list = line.split()
                for word in line_list:
                    if word.find('.gz')>0:
                        gzFileName = word
                        self.settings.gzFileNames.append(gzFileName)
                        self.experiment_dict[gzFileName] = copy.deepcopy(self.experiment)
                        self.settings_dict[gzFileName] = copy.deepcopy(self.settings)
                        print(f"Found gz: {gzFileName}")
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
        meansens = sens.astype(float).mean()
        self.experiment.sens = sens
        self.experiment.meansens = meansens
        return True


    def parse_intensity_map(self, inputd):
        sens = self.experiment.sens
        meansens = self.experiment.meansens
        monitor = self.experiment.monitor_counts
        quotient = np.divide(inputd, sens, out=np.zeros_like(inputd), where=sens!=0)
        self.experiment.Imatrix = float(meansens) * quotient / float(monitor)
        self.experiment.cut_Iz = self.experiment.Imatrix.sum(axis=1)
        self.experiment.cut_Iy = self.experiment.Imatrix.sum(axis=0)
        return True

    def compute_Q(self):
        experiment = self.experiment
        experiment.qzc_corr = experiment.qzc + int( ( 1990.0 * np.tan( np.pi * float(experiment.angle_of_incidence) / 180.0 ) ) / 0.5755 )
        experiment.qzc_spec = experiment.qzc - int( ( 1990.0 * np.tan( np.pi * float(experiment.angle_of_incidence) / 180.0 ) ) / 0.5755 )
        #print(f"Corrected qzc = {experiment.qzc_corr}")
        #print(f"Specular  qzc = {experiment.qzc_spec}")

        Imatrix = experiment.Imatrix
        ipix_range = np.asarray(range(Imatrix.shape[0]))
        jpix_range = np.asarray(range(Imatrix.shape[1]))

        two_pi_over_lambda = experiment.two_pi_over_lambda
        sin_alpha_i = experiment.sin_alpha_i
        sin_2theta_f = experiment.sin_2theta_f(ipix_range)
        cos_alpha_f = experiment.cos_alpha_f(jpix_range)
        sin_alpha_f = experiment.sin_alpha_f(jpix_range)


        QY_i, QY_j = sin_2theta_f[ipix_range], cos_alpha_f[jpix_range]
        QYmatrix = np.einsum('i,j->ij', QY_i, QY_j)
        experiment.qymatrix = two_pi_over_lambda * QYmatrix

        QZ_j = sin_alpha_f[jpix_range]
        QZ_i = np.ones(ipix_range.shape)
        QZmatrix = np.einsum('i,j->ij', QZ_i, QZ_j)
        experiment.qzmatrix = two_pi_over_lambda * (QZmatrix + sin_alpha_i)

        return True


    def save_gisans_map_filepath(self, inputd):
        raise NotImplementedError
#        qy = np.zeros(len(ipix_range)*len(jpix_range))
#        qz = np.zeros(len(ipix_range)*len(jpix_range))
#        I = np.zeros(len(ipix_range)*len(jpix_range))
#        idx = 0
#        with open(self.settings.gisans_map_filepath(), "w") as fp:
#            for i in ipix_range:
#                for j in jpix_range:
#                    if float(sens[i,j]) > 0.0:
#                        I[idx] = ((meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor))
#                        qy[idx] = (two_pi_over_lambda *  cos_alpha_f[j-pix0] * sin_2theta_f[i-pix0])
#                        qz[idx] = (two_pi_over_lambda * (sin_alpha_f[j-pix0] + sin_alpha_i))
#                        fp.write(str(qy[idx])+' '+str(qz[idx])+' '+str(I[idx])+'\n')
#                        idx += 1
#
#        self.experiment.inputd = inputd
#        self.experiment.qy = qy
#        self.experiment.qz = qz
#        self.experiment.I = I
#
##
#
#        I  = ((meansens*float(inputd[pix0:pixf+1,pix0:pixf+1])/float(sens[pix0:pixf+1,pix0:pixf+1]))/float(monitor))
#        qy = (two_pi_over_lambda *  cos_alpha_f[0:pixf+1-pix0] * sin_2theta_f[0:pixf+1-pix0])
#        qz = (two_pi_over_lambda * (sin_alpha_f[0:pixf+1-pix0] + sin_alpha_i))
        return


    def update_from_info_table(self):
        try:
            self.update_multi_experiment_values()
            self.color_outdated()
            self.update_from_selection_list()
        except Exception as e:
            App.handle_exception(e)
            return False
        return True


    def update_from_selection_list(self):
        selectedListEntries = self.fileList.selectedItems()
        if len(selectedListEntries) < 1:
            self.fileList.setCurrentRow( self.fileList.count() - 1 )
        try:
            if self.subtractCheckBox.isChecked():
                Imap = self.subtract_intensities_from_selected_files()
            else:
                Imap = self.sum_intensities_from_selected_files()

            self.compute_Q()
            self.graphView.update_graph(
                            Y = self.experiment.qymatrix,
                            X = self.experiment.qzmatrix,
                            Z = Imap,
                            Xc = self.experiment.qzc_corr,
                            Xs = self.experiment.qzc_spec,
                            Yc = self.experiment.qyc,
                            zmin = self.experiment.min_intensity,
                            zmax = self.experiment.max_intensity,
                            reset_limits_required=False,
                            x1=self.experiment.x0,
                            y1=self.experiment.y0,
                            x2=self.experiment.xf,
                            y2=self.experiment.yf,
                            title=self.build_ascii_header()
                        )
        except Exception as e:
            App.handle_exception(e)


    def sum_intensities_from_selected_files(self):
        selectedListEntries = self.fileList.selectedItems()
        Isum = []
        for currentListItem in selectedListEntries:
            currentListEntry = currentListItem.text()
            self.settings = self.settings_dict[currentListEntry]
            self.experiment = self.experiment_dict[currentListEntry]
            print(f"\n-- Calculating map for {currentListEntry} --")
            Isum += [self.experiment.Imatrix]

        if len(Isum) > 1:
            Isum = np.asarray(Isum).sum(axis=0)
        else:
            Isum = Isum[0]

        return Isum


    def subtract_intensities_from_selected_files(self):
        selectedListEntries = self.fileList.selectedItems()
        if len(selectedListEntries) != 2:
            raise ValueError("Substraction is only possible between two intensity maps")

        Isum = []
        for currentListItem in selectedListEntries:
            currentListEntry = currentListItem.text()
            self.settings = self.settings_dict[currentListEntry]
            self.experiment = self.experiment_dict[currentListEntry]
            print(f"\n-- Calculating map for {currentListEntry} --")
            Isum += [self.experiment.Imatrix]

        return np.abs(Isum[1] - Isum[0])


    @pyqtSlot()
    def update_gui(self):
        try:
            did_stuff, why_not = self.thread.retval
            if not did_stuff:
                raise Exception(f"Did not complete data processing: {why_not}")
            for gzfn in self.settings.gzFileNames:
                self.fileList.addItem(gzfn)
                sum1 = self.experiment_dict[gzfn].Imatrix.sum()
                sum2 = self.experiment_dict[gzfn].cut_Iy.sum()
                sum3 = self.experiment_dict[gzfn].cut_Iz.sum()
                print("If the following three sums are not equal, there's an error:")
                print(sum1, sum2, sum3)
            self.update_table()
            self.update_from_info_table()
            self.progress_bar.setValue(100)
        except Exception as e:
            App.handle_exception(e)
            return False
        finally:
            self.progress_bar.setValue(0)

        return True

if __name__ == '__main__':
    enable_high_dpi_scaling()
    app = qtw.QApplication(sys.argv)
    app.setStyleSheet("QMessageBox { messagebox-text-interaction-flags: 5;}")

    ex = App()
    #sys.exit(app.exec_())
    #x = profile_function_with_arguments(app.exec_)
    x = app.exec()
    sys.exit(x)