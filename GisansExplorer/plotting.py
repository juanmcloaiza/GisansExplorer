#!/bin/env python

#Qt stuff:
import PyQt5.QtWidgets as qtw
from PyQt5.QtCore import Qt, pyqtSignal


import numpy as np
import os
import sys
from sys import platform
from pathlib import Path

#plot stuff:
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter

if platform == "linux" or platform == "linux2":
    from matplotlib import pyplot as plt
elif platform == "darwin":
    gui_env = ['Qt5Agg','TKAgg','GTKAgg','WXAgg']
    for gui in gui_env:
        try:
            print("testing", gui)
            mpl.use(gui,warn=False, force=True)
            from matplotlib import pyplot as plt
            break
        except:
            continue
    print("Using:", mpl.get_backend())
elif platform == "win32":
    from matplotlib import pyplot as plt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.widgets import RectangleSelector
from matplotlib.figure import Figure



from gisansexplorer.utils import Frozen, handle_exception, __DEPLOYED__


class AreaSelector(Frozen):
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


class PlotData(Frozen):
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
###
class PlotStyle(object):
    fontSize = 20
    borderWidth = 3
    figSize = (12,9)
    lineWidth = 0.2 * fontSize
    majorTickLength = fontSize
    minorTickLength = 0.3*fontSize
    axisLabelFontSize = 1.5*fontSize
    titleFontSize = 0.5*fontSize

    @classmethod
    def apply_style(cls,ax):
        ax.tick_params(labelsize=0.7*cls.fontSize)
        ax.tick_params(width=cls.borderWidth, length=0.5*cls.fontSize, which='major')
        ax.tick_params(width=cls.borderWidth, length=0.2*cls.fontSize, which='minor')
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(cls.lineWidth)

        return


def create_gisans_figure(data, cnorm):
        ps = PlotStyle()
        x1, x2 = data.x1, data.x2
        y1, y2 = data.y1, data.y2
        x0 = np.abs(data.X[0]).argmin()
        y0 = np.abs(data.Y).T[0].argmin()
        xc, yc, xs = data.Xc, data.Yc, data.Xs
        data.Xzoom = data.X[y1:y2+1,x1:x2+1]
        data.Yzoom = data.Y[y1:y2+1,x1:x2+1]
        data.Zzoom = data.Z[y1:y2,x1:x2]

        data.zoom_extent = (data.X[y1,x1], data.X[y2,x2], data.Y[y1,x1], data.Y[y2,x2])

        new_fig = plt.figure(figsize=ps.figSize)
        new_ax = new_fig.add_subplot(111)
        cs = new_ax.pcolorfast(data.Xzoom, data.Yzoom, data.Zzoom, norm=cnorm, vmin=cnorm.vmin, vmax=cnorm.vmax, cmap='jet')

        has_legend = False
        if x0 > x1 and x0 < x2:
            new_ax.axvline(x=0, c='k', ls='solid', lw=ps.lineWidth, label = "Detector center")
        if xc > x1 and xc < x2:
            new_ax.axvline(x=data.X[yc,xc], c='r', lw=ps.lineWidth , label = "Detector center (corrected)")
        if xs > x1 and xs < x2:
            new_ax.axvline(x=data.X[yc,xs], c='g', lw=ps.lineWidth , label = "Specular beam")

        if yc > y1 and yc < y2:
            new_ax.axhline(y=data.Y[yc,xc], c='r', lw=ps.lineWidth)
        if y0 > y1 and y0 < y2:
            new_ax.axhline(y=0, c='k', ls='solid', lw=ps.lineWidth)

        new_ax.set_aspect("auto")
        new_ax.set_xlabel("$Q_{z}$", fontsize=ps.axisLabelFontSize)
        new_ax.set_ylabel("$Q_{y}$", fontsize=ps.axisLabelFontSize)
        cbar = new_fig.colorbar(cs)
        cbar.outline.set_linewidth(ps.lineWidth)
        new_ax.set_title(data.title, fontsize=ps.titleFontSize)

        new_ax.legend(fontsize=ps.fontSize)
            #loc='upper left', bbox_to_anchor= (0,1.06),
        #ncol=3, borderaxespad=0, frameon=False, fontsize = ps.fontSize)

        for ax in new_fig.axes:
            ps.apply_style(ax)

        return new_fig

def create_qz_integration_figure(data):
        ps = PlotStyle()
        new_fig = plt.figure(figsize=ps.figSize)
        new_ax = new_fig.add_subplot(111)
        if data.log_scale:
            try:
                new_ax.set_yscale('log')
            except Exception:
                pass

        integration_x = data.Zzoom.sum(axis=0)
        x0, xf = data.zoom_extent[0:2]
        rangex = np.linspace(x0, xf, len(integration_x))
        #xax_line = new_ax.plot(rangex, integration_x, lw=ps.lineWidth)
        new_ax.set_xlim((x0, xf))

        new_ax.xaxis.set_ticks(np.linspace(x0, xf, 5))

        zero = integration_x.min()
        mu =  integration_x.mean()
        sig = integration_x.std()
        #new_ax.set_yticks([zero, mu, mu+2*sig])
        new_ax.grid(which='both', axis='both', lw=ps.lineWidth)
        new_ax.fill_between(rangex, 0, integration_x,
                            facecolor='lightsalmon', edgecolor='orangered', alpha = 0.5,
                            linewidth=ps.lineWidth, zorder=10)
        new_ax.set_xlabel("$Q_{z}$", fontsize=ps.axisLabelFontSize)
        new_ax.set_ylabel("I($Q_{z})$", fontsize=ps.axisLabelFontSize)
        new_ax.set_title(data.title, fontsize=ps.titleFontSize)
        ps.apply_style(new_ax)

        return new_fig

def create_qy_integration_figure(data):
        ps = PlotStyle()
        new_fig = plt.figure(figsize=ps.figSize)
        new_ax = new_fig.add_subplot(111)
        if data.log_scale:
            try:
                new_ax.set_yscale('log')
            except Exception:
                pass

        integration_y =  data.Zzoom.sum(axis=1)
        y0, yf = data.zoom_extent[2:4]
        rangey = np.linspace(y0, yf, len(integration_y))
        #new_ax_line = new_ax.plot(rangey, integration_y, lw=ps.lineWidth)

        new_ax.set_xlim((y0, yf))
        new_ax.set_xticks(np.linspace(y0,yf,5))
        zero = integration_y.min()
        mu =  integration_y.mean()
        sig = integration_y.std()
        #new_ax.set_yticks([zero, mu, mu+2*sig])
        new_ax.grid(which='both', axis='both', lw=ps.lineWidth)
        new_ax.fill_between(rangey, 0, integration_y,
                            facecolor='cornflowerblue', edgecolor='royalblue', alpha=0.5,
                            linewidth=ps.lineWidth, zorder=10)
        new_ax.set_xlabel("$Q_{y}$", fontsize=ps.axisLabelFontSize)
        new_ax.set_ylabel("I($Q_{y})$", fontsize=ps.axisLabelFontSize)
        new_ax.set_title(data.title, fontsize=ps.titleFontSize)
        ps.apply_style(new_ax)

        return new_fig
###

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
        self.canvas.mpl_connect('button_press_event', self.on_mouse_click)
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

    def on_mouse_click(self, event):
        if event.dblclick:
            try:
                self.show_figures()
            except Exception as e:
                handle_exception(e)

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
            self.canvas.figure.suptitle("", fontsize=5)
            self.canvas.draw()
            self.save(**kwargs)
            self.finishedUpdating.emit()
        except Exception as e:
            handle_exception(e)
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
            self.norm = mpl.colors.SymLogNorm(vmin=self.data.zmin, vmax=self.data.zmax, linthresh=thres)#, base=10)
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
        self.ax.set_title(self.data.title + "\n - Detector View - ", fontsize = 6, pad=35)
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
        #self.xax_line = self.xax.plot(rangex, integration_x)

        self.xax_lin = self.xax.fill_between(rangex, 0, integration_x,
                            facecolor='lightsalmon', edgecolor='orangered', alpha = 0.5,
                            zorder=10)

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
        #self.yax_line = self.yax.plot(integration_y, rangey)

        self.yax_line = self.yax.fill_betweenx(rangey, integration_y,
                            facecolor='cornflowerblue', edgecolor='royalblue', alpha=0.5,
                            zorder=10)

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
        #self.cbar.mappable.set_clim(self.norm.vmin, self.norm.vmax)
        self.cbar.set_clim(self.norm.vmin, self.norm.vmax)
        self.cax.xaxis.tick_top()
        return #build_cbar

    def save_figures(self,filePath):
        self.save_gisans_map(filePath)
        self.save_qz_integration(filePath)
        self.save_qy_integration(filePath)
        return

    def show_figures(self):
        plt.close('all')
        create_gisans_figure(self.data, self.norm)
        create_qz_integration_figure(self.data)
        create_qy_integration_figure(self.data)
        plt.show()
        return

    def save_gisans_map(self,filePath=None):
        print("Saving gisans map....")
        new_fig = create_gisans_figure(self.data, self.norm)

        no_ext, ext = os.path.splitext(filePath)
        new_fig.savefig(f"{no_ext}-gisans_map{ext}")
        print(f"gisans map saved.")
        return #save_gisans_map#


    def save_qz_integration(self, filePath=None):
        print("Saving qz integration...")
        new_fig = create_qz_integration_figure(self.data)

        no_ext, ext = os.path.splitext(filePath)
        new_fig.savefig(f"{no_ext}-integration_qz{ext}")
        print("qz integration saved.")
        return #save_xax


    def save_qy_integration(self, filePath=None):
        print("Saving qy integration...")
        new_fig = create_qy_integration_figure(self.data)
        no_ext, ext = os.path.splitext(filePath)
        new_fig.savefig(f"{no_ext}-integration_qy{ext}")
        print("qy integration saved.")
        return #update_yax


    def test_show(self):
        #Z = np.loadtxt('./notToVersion/WelcomePlot.txt')
        #Z = (Z - Z.min()) / (Z.max() - Z.min()) * (1e-3 - 1e-6) + 1e-6
        #np.save("./show_test_map.npy", Z)
        if not __DEPLOYED__:
            filepath = Path('gisansexplorer/resources/show_test_map.npy')
        else:
            filepath = Path(sys.prefix+'/resources/show_test_map.npy')
        try:
            Z = np.load(filepath)
        except Exception as e:
            print(f"Failed to load: {filepath}")
            Z = np.abs(np.random.normal(0,1e-5,(1025,1025)))

        #This creates a cross to test the correct alignment between the
        #different plot axes:
        #Z[:,:] = 1e-5
        #Z[500:511,:] = 1e-4
        #Z[:,500:511] = 1e-4

        #print(Z.min())
        #print(Z.max())
        #print(Z.shape)

        x = np.linspace(-1,1, Z.shape[0])
        y = np.linspace(-1,1, Z.shape[1])
        X, Y = np.meshgrid(x,y)
        #Xc = len(x) // 2
        #Yc = len(y) //2
        self.update_graph(X = X, Y = Y, Z = Z, Xc = 0, Yc = 0)
        #if _DEBUG_:
        #    np.savetxt("./myNumpyArray.txt", 3 + 10*np.sin(np.sqrt(X**2 + Y**2)))
        return