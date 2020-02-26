#!/bin/env python

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from PyQt5.QtWidgets import QVBoxLayout, QDoubleSpinBox, QPushButton, QFormLayout, QMessageBox
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog, QLabel
from PyQt5.QtCore import Qt, pyqtSlot

import numpy as np
import sys
import os
import traceback
import gzip




class FrozenClass(object):
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "%r is a frozen class" % self )
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True


class Canvas(FigureCanvas):
    def __init__(self, parent = None, width = 28, height = 28, dpi = 100):
        fig, self.axes = plt.subplots(1, 1, figsize=(width, height), dpi=dpi, sharex=False, sharey=False)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)


    def scatter(self, xarr, yarr, Iarr, vmin = None, vmax = None):
        if vmin is None:
            vmin = min(Iarr)

        if vmax is None:
            vmax = max(Iarr)

        fig = self.figure
        axes = self.axes
        im = axes.scatter(xarr, yarr, edgecolors='none', c=Iarr, 
                            norm = LogNorm(vmin=vmin, vmax=vmax),
                            marker = 's')
        axes.set_title("GISANS map")
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label("Intensity", fontsize=20)
        plt.ylabel(r'$Q_{y}(\AA^{-1})$')
        plt.xlabel(r'$Q_{z}(\AA^{-1})$') 


class Experiment(FrozenClass):
    def __init__(self):
        self.selector_lambda = None
        self.angle_of_incidence = None

        # Define the position of the direct beam on the detector
        # and the sensitivity map file
        self.yc = 528 #530
        self.zc = 211 #220 #zc=512
        self.sens = None
        self.meansens = None
        self.monitor_counts = None
        self.y=[]
        self.z=[]
        self.I=[]
        self._freeze()
        return


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


class App(QWidget,FrozenClass):

    def __init__(self):
        super().__init__()
        self.title = 'Alexandros GISANS Viewer'
        self.left = 3000
        self.top = 500
        self.width = 640
        self.height = 480
        self.layout = QVBoxLayout()
        self.minSpinBox = QDoubleSpinBox()
        self.maxSpinBox = QDoubleSpinBox()
        self.settings = Settings()
        self.experiment = Experiment()
        self.canvas = Canvas()

        self._freeze()
        self.initUI()


    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.setLayout(self.layout)
        self.layout.setAlignment(Qt.AlignCenter)
        self.addWelcomeMessage()
        self.addCanvas()
        self.addMinMaxSpinBoxes()
        self.addOpenFileButton()
        self.show()
        #self.openFileNamesDialog()
        #self.saveFileDialog()


    def addOpenFileButton(self):
        buttonOpenDialog = QPushButton("Press here")
        buttonOpenDialog.clicked.connect(self.on_click)
        self.layout.addWidget(buttonOpenDialog)
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
        formLayout.addRow(self.tr("&Give Min Value"), self.minSpinBox)
        formLayout.addRow(self.tr("&Give Max Value"), self.maxSpinBox)
        formLayout.setFormAlignment(Qt.AlignCenter)
        self.layout.addLayout(formLayout)
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
        self.layout.addWidget(message)
        return

    def addCanvas(self):
        self.layout.addWidget(self.canvas)
        return


    @pyqtSlot()
    def on_click(self):
        print('PyQt5 button click')
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
        if fileName:
            return fileName 
        return None


    def doStuff(self):
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


    def handle_exception(self, e):
        msg = (f"Exception: {e}\n")
        msg += ("-"*60+"\n")
        msg += traceback.format_exc()
        msg += ("-"*60+"\n")
        print(msg)
        #traceback.#.print_exc(file=sys.stdout)
        pop_up = QMessageBox()
        pop_up.setWindowTitle(f"Exception: {e}\n")   
        pop_up.setText(msg)
        pop_up.setIcon(QMessageBox.Critical)
        x = pop_up.exec_()
        return


    def read_dat_file(self):
        # Open and read the dat file
        datFilePath = self.openFileNameDialog()
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
            return self.safe_parse_numpy(func, fpath, dtype='i', delimiter=' ')
        return False


    def read_intensity_file(self):
        if self.settings.gzFileName:
            fpath = self.settings.gzFilePath()
            func = self.parse_intensity_map
            return self.safe_parse_numpy(func, fpath, dtype='i', delimiter=' ')
        return False


    def parse_dat(self, file):
        for line in file:
            if line.find('omega_value')>0:
                omega_line_list = line.split()
                omega=omega_line_list[3]
                self.experiment.angle_of_incidence = omega
                print('Angle of incidence (degrees): '+str(omega))
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
        for i in range(1,1024):
             for j in range(1,1024):
                 if i>=162 and i<=862:
                     if j>=162 and j<=862:
                         meansens=meansens+float(sens[i,j])
                         num=num+1.0
        meansens=meansens/num
        self.experiment.sens = sens 
        self.experiment.meansens = meansens
        return True


    def parse_intensity_map(self, inputd):
        zc = self.experiment.zc 
        yc = self.experiment.zc 
        meansens = self.experiment.meansens
        sens = self.experiment.sens
        omega = self.experiment.angle_of_incidence
        selector_lambda = self.experiment.selector_lambda
        monitor = self.experiment.monitor_counts

        print(f"Original zc = {zc}")
        zc=zc+int((1990.0*np.tan(np.pi*float(omega)/180.0))/0.5755)
        print(f"Corrected zc = {zc}")

        y, z, I = [], [], [] 
        with open(self.settings.gisans_map_filepath(), "w") as fp:
            for i in range(1,1024):
                 for j in range(1,1024):
                     if i>=162 and i<=862:
                         if j>=162 and j<=862:
                             y.append((2*np.pi/float(selector_lambda))*
                                np.cos(np.pi*float(omega)/180.0)*
                                (0.5755*(yc-i)/np.sqrt(1990*1990+(0.5755*(yc-i))*
                                (0.5755*(yc-i)))))
                             z.append((2*np.pi/float(selector_lambda))*
                                (np.sin(np.pi*float(omega)/180.0)+
                                (0.5755*(j-zc)/np.sqrt(1990*1990+(0.5755*(zc-j))*(0.5755*(zc-j))))))
                             if float(sens[i,j]) > 0.0:
                                 I.append((meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor))
                                 fp.write(str(y[-1])+' '+str(z[-1])+' '+str(I[-1])+'\n')
        self.experiment.y = np.asarray(y)
        self.experiment.z = np.asarray(z)
        self.experiment.I = np.asarray(I)
        return True


    def show_gisans_map(self):
        if self.settings.cbar_min < 0:
            min_value = min(self.experiment.I)+float(1)/float(self.experiment.monitor_counts)
        else:
            min_value = self.settings.cbar_min

        if self.settings.cbar_max < 0:
            max_value=max(self.experiment.I)
        else:
            max_value = self.settings.cbar_max

        self.canvas.scatter(self.experiment.y,
                            self.experiment.z,
                            self.experiment.I,
                            vmin=min_value,
                            vmax=max_value
                            )

        return True

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())


exit(0)


minmax = raw_input("enter min/max (yes/no): ")
#minmax='no'
if minmax=='yes':
    gmin=input('Give min value:');
    gmax=input('Give max value:');

#qzvalue=0.0
#dqzvalue=0.05
#qzvalue=input("Give qz: ");
#dqzvalue=input("Give delta qz: ");

# Open file dialog
file_path_string = openFileNameDialog()

   
basename = file_path_string[:-4]

# Open and read the dat file
file = open(file_path_string, 'r') 

for line in file:
    if line.find('omega_value')>0:
        omega_line_list = line.split()
        omega=omega_line_list[3]
        print('Angle of incidence (degrees): '+str(omega))
    if line.find('selector_lambda_value')>0:
        lambda_line_list = line.split()
        selector_lambda=lambda_line_list[3]
        print('Neutron wavelength (Angtrom): '+str(selector_lambda))
    if line.find('.gz')>0:
        line_list = line.split()
        for word in line_list:
            if word.find('.gz')>0:
                file_path_string = word
    if line.find('.yaml')>0:
        line_list = line.split()
        for word in line_list:
            if word.find('.yaml')>0:
                file_path_stringB = word

file.close()

# Open and read the yaml file
line1=''
line2=''
line3=''
line4=''
file = open(file_path_stringB, 'r') 
for line in file:
    if line4.find('name: mon1')>0:
        line_list=line1.split()
        monitor=line_list[1]
        print('Monitor conuts: '+str(monitor)+'\n')
    line4=line3
    line3=line2
    line2=line1
    line1=line

file.close()

# Define the position of the direct beam on the detector
yc=528 #530
zc=211 #220
#zc=512

y=[]
z=[]
I=[]

print('Loading data and applying corrections... please wait..'+'\n')

# Unzip the file
os.system('cp '+file_path_string+' temp.gz')
os.system('gzip -d -q temp.gz')


sens = np.loadtxt("sensitivity_map", dtype='i', delimiter=' ')
#sens = np.loadtxt("sens.map", dtype='i', delimiter=' ')

num=0.0
meansens=0.0
for i in xrange(1,1024):
     for j in xrange(1,1024):
         if i>=162 and i<=862:
             if j>=162 and j<=862:
                 meansens=meansens+float(sens[i,j])
                 num=num+1.0

meansens=meansens/num

inputd = np.loadtxt("temp", dtype='i', delimiter=' ')
os.system('rm temp')
#os.system('rm temp.gz')

#file = open('output.dat','w') 

# Addition important... qz=0 at the horizon of the sample

print(zc)
zc=zc+int((1990.0*np.tan(np.pi*float(omega)/180.0))/0.5755)
print(zc)

file = open(basename+"_GISANS.map", "w")

for i in xrange(1,1024):
     for j in xrange(1,1024):
         if i>=162 and i<=862:
             if j>=162 and j<=862:
                 y.append((2*np.pi/float(selector_lambda))*np.cos(np.pi*float(omega)/180.0)*(0.5755*(yc-i)/np.sqrt(1990*1990+(0.5755*(yc-i))*(0.5755*(yc-i)))))
                 z.append((2*np.pi/float(selector_lambda))*(np.sin(np.pi*float(omega)/180.0)+(0.5755*(j-zc)/np.sqrt(1990*1990+(0.5755*(zc-j))*(0.5755*(zc-j))))))
                 if float(sens[i,j]) > 0.0:
                     I.append((meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor))
                     file.write(str((2*np.pi/float(selector_lambda))*np.cos(np.pi*float(omega)/180.0)*(0.5755*(yc-i)/np.sqrt(1990*1990+(0.5755*(yc-i))*(0.5755*(yc-i)))))+'   '+str((2*np.pi/float(selector_lambda))*(np.sin(np.pi*float(omega)/180.0)+(0.5755*(j-zc)/np.sqrt(1990*1990+(0.5755*(zc-j))*(0.5755*(zc-j))))))+'   '+str((meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor))+'\n')
                    #file.write(str((2*np.pi/float(selector_lambda))*np.cos(np.pi*float(omega)/162.0)*((yc-i)/np.sqrt(1990*1990+(0.6*(yc-i))*(0.6*(zc-i)))))+' '+str((-2*np.pi/float(selector_lambda))*(np.sin(np.pi*float(omega)/162.0)+((zc-j)/np.sqrt(1990*1990+(0.6*(zc-j))*(0.6*(zc-j))))))+' '+str(float(input[i,j])/float(monitor))+'\n')

file.close()

min_value=min(I)+float(1)/float(monitor)


max_value=max(I)

if minmax=='yes':
    min_value=gmin
    max_value=gmax


plt.scatter(z,y,edgecolors='none',c=I, norm = LogNorm(vmin=min_value, vmax=max_value), marker = 's')
cbar = plt.colorbar()
plt.ylabel(r'$Q_{y}(\AA^{-1})$')
plt.xlabel(r'$Q_{z}(\AA^{-1})$') 
plt.show()


linecut = raw_input("Line cut (yes/no): ");

qdirection=' '
if linecut=='yes':
    qdirection=raw_input('qy or qz? ');


if qdirection=='qy' and linecut=='yes':
    qzvalue=input("Give qz: ");
    dqzvalue=input("Give delta qz: ");

    cut_qy=[0] * len(z)
    cut_qy_xaxis=[0] * len(z)
    for i in xrange(1,1024):
        for j in xrange(1,1024):
            if i>=162 and i<=862:
                if j>=162 and j<=862:
                    cut_qy_xaxis[i]=((2*np.pi/float(selector_lambda))*np.cos(np.pi*float(omega)/180.0)*(0.5755*(yc-i)/np.sqrt(1990*1990+(0.5755*(yc-i))*(0.5755*(yc-i)))))
                    cqz=((2*np.pi/float(selector_lambda))*(np.sin(np.pi*float(omega)/180.0)+(0.5755*(j-zc)/np.sqrt(1990*1990+(0.5755*(zc-j))*(0.5755*(zc-j))))))
                    if cqz>=qzvalue-dqzvalue and cqz<=qzvalue+dqzvalue:
                        cut_qy[i]=cut_qy[i]+(meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor)
    plt.loglog(cut_qy_xaxis,cut_qy,label="cut at "+r'$Q_{z}$= '+'{:06.4f}'.format(qzvalue)+'$ \AA^{-1})$'+'\n'+r'$\Delta Q_{z}$= '+'{:06.4f}'.format(dqzvalue)+'$ \AA^{-1})$')
    plt.legend()
    plt.ylabel(r'$I(A.U.)$')
    plt.xlabel(r'$Q_{y}(\AA^{-1})$')
    plt.show()


    file = open(basename+"_line_cut_qy.out", "w")
    for i in xrange(1,1024):
        if i>=162 and i<=862:
            file.write(str(cut_qy_xaxis[i])+'   '+str(cut_qy[i])+"\n")
    file.close()

if qdirection=='qz' and linecut=='yes':
    qyvalue=input("Give qy: ");
    dqyvalue=input("Give delta qy: ");

    cut_qz=[0] * len(y)
    cut_qz_xaxis=[0] * len(y)
    for i in xrange(1,1024):
        for j in xrange(1,1024):
            if i>=162 and i<=862:
                if j>=162 and j<=862:
                    cut_qz_xaxis[j]=((2*np.pi/float(selector_lambda))*(np.sin(np.pi*float(omega)/180.0)+(0.5755*(j-zc)/np.sqrt(1990*1990+(0.5755*(zc-j))*(0.5755*(zc-j))))))
                    cqy=((2*np.pi/float(selector_lambda))*np.cos(np.pi*float(omega)/180.0)*(0.5755*(yc-i)/np.sqrt(1990*1990+(0.5755*(yc-i))*(0.5755*(yc-i)))))
                    if cqy>=qyvalue-dqyvalue and cqy<=qyvalue+dqyvalue:
                        cut_qz[j]=cut_qz[j]+(meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor)
    plt.loglog(cut_qz_xaxis,cut_qz,label="cut at "+r'$Q_{y}$= '+'{:06.4f}'.format(qyvalue)+'$ \AA^{-1})$'+'\n'+r'$\Delta Q_{y}$= '+'{:06.4f}'.format(dqyvalue)+'$ \AA^{-1})$')
    plt.legend()
    plt.ylabel(r'$I(A.U.)$')
    plt.xlabel(r'$Q_{z}(\AA^{-1})$')
    plt.show()


    file = open(basename+"_line_cut_qz.out", "w")
    for i in xrange(1,1024):
        if i>=162 and i<=862:
            omega_value=cut_qz_xaxis[i]*float(selector_lambda)*180.0/(4.0*np.pi*np.pi)
            file.write(str(omega_value)+'   '+str(cut_qz_xaxis[i])+'   '+str(cut_qz[i])+"\n")
    file.close()
