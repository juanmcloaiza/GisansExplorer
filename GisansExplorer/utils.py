#!/bin/env python

#Qt stuff:
import PyQt5.QtWidgets as qtw
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QValidator


import numpy as np
from pathlib import Path
import re
import traceback

# Modules for profiling:
import cProfile, pstats, io


#### sciSpinBox

# Regular expression to find floats. Match groups are the whole string, the
# whole coefficient, the decimal part of the coefficient, and the exponent
# part.
_float_re = re.compile(r'(([+-]?\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)')
__DEBUG__=False
__DEPLOYED__=True
__DATA_ROOT_DIR__ = "./notToVersion/Archive"

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
            return asdf
        except Exception as e:
            handle_exception(e)
        return

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

#### Frozen:
class Frozen(object):
    __isfrozen = False

    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("%r is a frozen class" % self)
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True
#### Frozen end

#### functions for profiling
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
#### end of functions for profiling


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

def is_file(path):
    my_file = Path(path)
    return my_file.is_file()