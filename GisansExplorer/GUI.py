#!/bin/env python

#Qt stuff:
import PyQt5.QtWidgets as qtw
from PyQt5.QtCore import Qt, pyqtSlot
from PyQt5.QtGui import QColor

import numpy as np
import copy
import os

from gisansexplorer.utils import Frozen, handle_exception, is_file, __DEBUG__, __DATA_ROOT_DIR__
from gisansexplorer.experiment_data_handling import Settings, Experiment, FileReadingThread
from gisansexplorer.plotting import MyGraphView


class MyTabs(qtw.QTabWidget,Frozen):
    """
    Collection of tabs hosting gisans data frames

    Attributes
    ----------
    tabButton_add: QToolButton
        Adds a new tab
    tabButton_rmv: QToolButton
        Removes the current tab
    frameList: List
        List of gisans frames
    last_num: int
        index of last tab created

    Methods
    -------
    initCornerButton()
        Adds and creates connections for the add/rmv tab buttons
    addTab():
        Adds a new tab
    removeTab():
        Removes current tab
    """
    
    def __init__(self):
        """[summary]
        """
        super().__init__()
        self.tabButton_add = qtw.QToolButton()
        self.tabButton_rmv = qtw.QToolButton()
        self.frameList =[]
        self.last_num = 0
        self.addTab()
        self.initCornerButton()
        self._freeze()

    def initCornerButton(self):
        """[summary]
        """
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
        """[summary]
        """
        idx = self.currentIndex()
        del self.frameList[idx]
        super().removeTab(idx)
        if len(self.frameList) == 0:
            self.last_num = 0
            self.addTab()

        return


class MyFrame(qtw.QFrame,Frozen):
    """[summary]

    Parameters
    ----------
    qtw : [type]
        [description]
    Frozen : [type]
        [description]
    """
    def __init__(self):
        """[summary]
        """
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
        self.thread = FileReadingThread(self)

        self.initFrame()
        self._freeze()


    def initFrame(self):
        """[summary]
        """
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
        """[summary]

        Parameters
        ----------
        botSplitter : [type]
            [description]
        """
        model = qtw.QFileSystemModel()
        model.setRootPath('')
        filters = ["*.dat"]
        model.setNameFilters(filters)
        model.setNameFilterDisables(False)

        self.dirtree.setModel(model)
        self.dirtree.setRootIndex(model.index(__DATA_ROOT_DIR__))

        self.dirtree.setAnimated(True)
        self.dirtree.setIndentation(20)
        self.dirtree.setSortingEnabled(True)
        self.dirtree.doubleClicked.connect(self.on_file_double_clicked)
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
        """[summary]
        """

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
        """[summary]
        """
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
        """[summary]
        """
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
        """[summary]

        Parameters
        ----------
        value : [type]
            [description]
        """
        self.progress_bar.setValue(value)



    @pyqtSlot()
    def on_click_open_file(self):
        try:
            if __DEBUG__:
                datFilePath = os.path.join(".","notToVersion","Archive","p15496_00000992.dat")
            else:
                datFilePath = self.openFileNameDialog()

            if is_file(datFilePath):
                self.doStuff(datFilePath)
        except Exception as e:
            handle_exception(e)

    @pyqtSlot()
    def on_file_double_clicked(self):
        """[summary]
        """
        try:
            if len(self.dirtree.selectedIndexes()) < 1:
                return
            datFilePath = self.dirtree.model().filePath(self.dirtree.currentIndex())
            if is_file(datFilePath):
                self.doStuff(datFilePath)
        except Exception as e:
            handle_exception(e)


    @pyqtSlot()
    def on_subtract_checkbox_changed(self):
        """[summary]
        """
        self.update_from_selection_list()

    @pyqtSlot()
    def on_file_selection_changed(self):
        """[summary]
        """
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
        """[summary]
        """
        self.update_from_info_table()

    @pyqtSlot()
    def on_cell_changed(self):
        """[summary]
        """
        self.color_outdated()


    @pyqtSlot()
    def on_click_loglinear(self):
        """[summary]
        """
        try:
            self.graphView.update_graph(log_scale = not self.graphView.data.log_scale, reset_limits_required=False)
        except Exception as e:
            handle_exception(e)
        return


    @pyqtSlot()
    def on_click_save_png(self):
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
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
            handle_exception(e)
        return #on_click_save_png


    def build_ascii_header(self):
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
        gz_filenames = [s.text() for s in self.fileList.selectedItems()]
        dat_filenames= [self.settings_dict[s].datFileName for s in gz_filenames]
        op_str = ""
        if self.subtractCheckBox.isChecked():
            op_str = "Subtraction of:\n"
        else:
            op_str = "Addition of:\n"
        header = op_str + "\n".join([d + ": " + gz for (d,gz) in zip(dat_filenames, gz_filenames)])
        header += "\n - "+self.experiment.instrument_name+" - "
        return header


    @pyqtSlot()
    def on_click_save_ascii(self):
        """[summary]
        """
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
            handle_exception(e)


#     @pyqtSlot()
#     def on_spinbox_edit(self):
#         try:
#             self.graphView.update_graph(zmax=self.maxSpinBox.value(), zmin=self.minSpinBox.value())
#         except Exception as e:
#             handle_exception(e)
#         return


    @pyqtSlot()
    def on_graph_updated(self):
        """[summary]
        """
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
            handle_exception(e)
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
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
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
            handle_exception(e)
        return None

    def openFileNameDialog(self):
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
        try:
            options = qtw.QFileDialog.Options()
            options |= qtw.QFileDialog.DontUseNativeDialog
            fileName, _ = qtw.QFileDialog.getOpenFileName(self,"Open File", "",
            "Measurement dat file (*.dat);;All Files (*)",
            options=options)
            # self.openFileNamesDialog()
            if fileName:
                return fileName
        except Exception as e:
            handle_exception(e)
        return None


    @staticmethod
    def color_validate(table_item, value_a, value_b):
        """[summary]

        Parameters
        ----------
        table_item : [type]
            [description]
        value_a : [type]
            [description]
        value_b : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """
        green = QColor(000, 255, 0, 127)
        if value_a != value_b:
            table_item.setBackground(green)
            table_item.setSelected(False)
        else:
            table_item.setBackground(Qt.white)
        return True

    def update_single_experiment_values(self, experiment):
        """[summary]

        Parameters
        ----------
        experiment : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """
        expdict = experiment.__dict__
        for i in range(self.infoTable.rowCount()):
            key = self.infoTable.verticalHeaderItem(i).text()
            current_item = self.infoTable.item(i,0)
            if key in ['qyc', 'qzc', 'x0', 'y0', 'xf', 'yf', 'sample_detector_distance_mm']:
                value = int(current_item.text())
            elif key in ['min_intensity', "max_intensity", 'pixel_size_mm']:
                value = float(current_item.text())
            elif key in ['instrument_name']:
                value = str(current_item.text())
            else:
                continue
            expdict[key] = value
        return True

    def update_multi_experiment_values(self):
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
        selectedListEntries = self.fileList.selectedItems()
        for currentListItem in selectedListEntries:
            currentListEntry = currentListItem.text()
            self.settings = self.settings_dict[currentListEntry]
            self.experiment = self.experiment_dict[currentListEntry]
            self.update_single_experiment_values(self.experiment)
        return True


    def color_outdated(self):
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
        expdict = self.experiment.__dict__

        try:
            for i in range(self.infoTable.rowCount()):
                key = self.infoTable.verticalHeaderItem(i).text()
                current_item = self.infoTable.item(i,0)
                try:
                    if key in ['qyc', 'qzc', 'x0', 'y0', 'xf', 'yf','sample_detector_distance_mm']:
                        value = int(current_item.text())
                    elif key in ["min_intensity", "max_intensity", "pixel_size_mm"]:
                        value = float(current_item.text())
                    elif key in ['instrument_name']:
                        value = str(current_item.text())
                    else:
                        continue

                    self.color_validate(current_item, value, expdict[key])
                except Exception as e:
                    red = QColor(255, 000, 0, 127)
                    current_item.setBackground(red)
                    current_item.setSelected(False)
                    return False

        except Exception as e:
            handle_exception(e)
            return False
        return True


    def update_table(self):
        """[summary]

        Returns
        -------
        [type]
            [description]
        """

        expdict = self.experiment.__dict__

        self.infoTable.setRowCount(0)
        for k in expdict.keys():
            if k[0] == "_":
                continue

            if k in ["instrument_name", "pixel_size_mm", "sample_detector_distance_mm",
                     "selector_lambda", "qyc", "qzc", "x0", "y0", "xf", "yf",
                     "min_intensity", "max_intensity",
                     "meansens", "monitor_counts",
                     "angle_of_incidence", "qzc_corr", "qzc_spec"]:
                i = self.infoTable.rowCount()
                item_k = qtw.QTableWidgetItem(str(k))
                item_v = qtw.QTableWidgetItem(str(expdict[k]))
                self.infoTable.insertRow(i)
                self.infoTable.setVerticalHeaderItem(i,item_k)
                self.infoTable.setItem(i,0,item_v)

            if k in [#"instrument_name", "pixel_size_mm", "sample_detector_distance_mm",
                     "selector_lambda", "meansens", "monitor_counts", "angle_of_incidence",
                     "qzc_corr", "qzc_spec"]:
                current_item = self.infoTable.item(i,0)
                current_item.setBackground(QColor("lightGrey"))
                current_item.setFlags(Qt.ItemIsEnabled)

        return True


    def doStuff(self, datFilePath):
        """[summary]

        Parameters
        ----------
        datFilePath : [type]
            [description]
        """
        try:
            self.thread.datFilePath = datFilePath
            self.thread.start()
        except Exception as e:
            handle_exception(e)
            return
        return


    def safe_parse(self, parse_func, file_path):
        """[summary]

        Parameters
        ----------
        parse_func : [type]
            [description]
        file_path : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """
        try:
            with open(file_path, 'r') as fp:
                tf = parse_func(fp)
            print(f" Parsed {file_path}")
            return tf
        except Exception as e:
            handle_exception(e)
            return False


    def safe_parse_numpy(self, parse_func, file_path, dtype='i', delimiter=' '):
        """[summary]

        Parameters
        ----------
        parse_func : [type]
            [description]
        file_path : [type]
            [description]
        dtype : str, optional
            [description], by default 'i'
        delimiter : str, optional
            [description], by default ' '

        Returns
        -------
        [type]
            [description]
        """
        try:
            nparray = np.loadtxt(file_path, dtype=dtype, delimiter=delimiter)
            print(f"{file_path}:")
            print(f"Loaded array with shape: {nparray.shape}")
            tf = parse_func(nparray)
            print(f" Parsed {file_path}")
            return tf
        except Exception as e:
            handle_exception(e)
            return False



    def read_dat_file(self, datFilePath=None):
        """[summary]

        Parameters
        ----------
        datFilePath : [type], optional
            [description], by default None

        Returns
        -------
        [type]
            [description]
        """

        if datFilePath:
            path, filename = os.path.split(datFilePath)
            self.settings.datFileName = filename
            self.settings.dataDirPath = path
            return self.safe_parse(self.parse_dat, self.settings.datFilePath())
        return False


    def read_sensitivity_file(self):
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
        func = self.parse_sensitivity_map
        fpath = self.settings.sensFilePath()
        for key in self.settings.gzFileNames:
            self.experiment = self.experiment_dict[key]
            tf = self.safe_parse_numpy(func, fpath, dtype=float, delimiter=' ')
            if not tf:
                return False
        return True


    def read_yaml_file(self):
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
        # Open and read the yaml file
        for yamlFileName, key in zip(self.settings.yamlFilePaths(), self.settings.gzFileNames):
            self.experiment = self.experiment_dict[key]
            tf = self.safe_parse(self.parse_yaml, yamlFileName)
            if not tf:
                return False
        return True


    def read_intensity_file(self):
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
        func = self.parse_intensity_map
        for fpath, key in zip(self.settings.gzFilePaths(),self.settings.gzFileNames):
            self.experiment = self.experiment_dict[key]
            tf = self.safe_parse_numpy(func, fpath, dtype=float, delimiter=' ')
            if not tf:
                return False
        return True


    def parse_dat(self, file):
        """[summary]

        Parameters
        ----------
        file : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """
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
        """[summary]

        Parameters
        ----------
        fp : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """
        # Open and read the yaml file
        line1=''
        line2=''
        line3=''
        line4=''
        monitor = None
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
        """[summary]

        Parameters
        ----------
        sens : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """
        meansens = sens.astype(float).mean()
        self.experiment.sens = sens
        self.experiment.meansens = meansens
        return True


    def parse_intensity_map(self, inputd):
        """[summary]

        Parameters
        ----------
        inputd : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """
        sens = self.experiment.sens
        meansens = self.experiment.meansens
        monitor = self.experiment.monitor_counts
        quotient = np.divide(inputd, sens, out=np.zeros_like(inputd), where=sens!=0)
        self.experiment.Imatrix = float(meansens) * quotient / float(monitor)
        self.experiment.cut_Iz = self.experiment.Imatrix.sum(axis=1)
        self.experiment.cut_Iy = self.experiment.Imatrix.sum(axis=0)
        return True

    def compute_Q(self):
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
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
        """[summary]

        Parameters
        ----------
        inputd : [type]
            [description]

        Raises
        ------
        NotImplementedError
            [description]
        """
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
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
        try:
            self.update_multi_experiment_values()
            self.color_outdated()
            self.update_from_selection_list()
        except Exception as e:
            handle_exception(e)
            return False
        return True


    def update_from_selection_list(self):
        """[summary]
        """
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
            handle_exception(e)


    def sum_intensities_from_selected_files(self):
        """[summary]

        Returns
        -------
        [type]
            [description]
        """
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
        """[summary]

        Returns
        -------
        [type]
            [description]

        Raises
        ------
        ValueError
            [description]
        """
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
        """[summary]

        Returns
        -------
        [type]
            [description]

        Raises
        ------
        Exception
            [description]
        """
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
            handle_exception(e)
            return False
        finally:
            self.progress_bar.setValue(0)

        return True
