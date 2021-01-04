#!/bin/env python

import numpy as np
import os
from PyQt5.QtCore import pyqtSignal, QThread

from gisansexplorer.utils import Frozen, handle_exception


class Instrument(Frozen):
    def __init__(self, name, pixel_size_mm=None, sample_detector_distance_mm=None, data_root_dir="/"):

        if name.casefold() == "maria":
            self.name = "MARIA"
            self.pixel_size_mm = 0.5755
            self.sample_detector_distance_mm = 1990
        elif name.casefold() == "other":
            self.name = "OTHER"
            self.pixel_size_mm = 0.5755
            self.sample_detector_distance_mm = 1990
        else:
            self.name = name
            self.pixel_size_mm = pixel_size_mm
            self.sample_detector_distance_mm = sample_detector_distance_mm

        self._freeze()

        return

__instrument__ = Instrument("MARIA")

class Experiment(Frozen):
    def __init__(self, instrument=__instrument__):
        if instrument is None:
            self.instrument_name = "Default"
            self.pixel_size_mm = 0.5755
            self.sample_detector_distance_mm = 1990
        else:
            self.instrument_name = instrument.name
            self.pixel_size_mm = instrument.pixel_size_mm
            self.sample_detector_distance_mm = instrument.sample_detector_distance_mm

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
        px_size = self.pixel_size_mm
        sdd = self.sample_detector_distance_mm
        return px_size*(pixel_i-self.qyc)/np.sqrt(sdd*sdd+(px_size*(self.qyc-pixel_i))*(px_size*(self.qyc-pixel_i)))


    def sin_alpha_f(self, pixel_j):
        px_size = self.pixel_size_mm
        sdd = self.sample_detector_distance_mm
        return px_size*(pixel_j-self.qzc_corr)/np.sqrt(sdd*sdd+(px_size*(self.qzc_corr-pixel_j))*(px_size*(self.qzc_corr-pixel_j)))


    def cos_alpha_f(self, pixel_j):
        return np.sqrt(1 - self.sin_alpha_f(pixel_j)**2)



class Settings(Frozen):
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


class FileReadingThread(QThread):
    progress_signal = pyqtSignal(int)
    def __init__(self, myframe):
        super().__init__()
        self.frame = myframe
        self.retval = None
        self.datFilePath  = None

    def run(self):
        myframe = self.frame
        myframe.settings = Settings()
        myframe.experiment = Experiment()
        
        try:
            self.progress_signal.emit(0)
            if not myframe.read_dat_file(self.datFilePath):
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
        except Exception as e:
            self.retval = False, "An exception occurred"
            return