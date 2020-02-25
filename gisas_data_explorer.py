#!/bin/env python

from PyQt5.QtWidgets import QMainWindow, QApplication, QPushButton
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
 
 
 
class Window(QMainWindow):
    def __init__(self):
        super().__init__()
 
        title = "Matplotlib Embeding In PyQt5"
        top = 400
        left = 400
        width = 1000
        height = 500
 
        self.setWindowTitle(title)
        self.setGeometry(top, left, width, height)
 
        self.MyUI()
 
 
    def MyUI(self):
 
        canvas = Canvas(self, width=8, height=4)
        canvas.move(0,0)
 
        button = QPushButton("Click Me", self)
        button.move(100, 450)
 
        button2 = QPushButton("Click Me Two", self)
        button2.move(250, 450)
 
 
class Canvas(FigureCanvas):
    def __init__(self, parent = None, width = 28, height = 14, dpi = 100, I = 2, J = 3):
        fig, self.axes = plt.subplots(I, J, figsize=(width, height), dpi=dpi, sharex=False, sharey=False)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        self.plot()

    def plot(self):
        fig = self.figure
        axes = self.axes
        for i in range(axes.shape[0]):
            for j in range(axes.shape[1]):
                axes_limits = [0,10,0,10]
                axes_labels = ["axis a", "axis b"]
                image = np.random.normal(0,1,(100,100))
                title = "Title"
                im = axes[i,j].imshow(image, norm=matplotlib.colors.LogNorm(), extent=axes_limits, 
                        cmap="jet", vmin=1e-3, vmax=1e5, aspect="equal")
                axes[i,j].set_title(title)
    
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label("Intensity", fontsize=20)

        fig.text(0.5, 0.04, axes_labels[0], ha='center', fontsize=20)
        fig.text(0.04, 0.5, axes_labels[1], va='center', fontsize=20, rotation='vertical')
 
app = QApplication(sys.argv)
window = Window()
window.show()
app.exec()