from PySide2.QtWidgets import *
from PySide2.QtCore import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backend_bases import PickEvent
from matplotlib.collections import PathCollection
import numpy as np
from nilearn import image
import os
import json
from pprint import pprint
import re


# noinspection PyCallingNonCallable
class xASL_GUI_MRIViewArtist(QWidget):
    """
    Main Widget class to handle drawing operations around seaborn's scatterplot class as well as MRI slices
    """

    def __init__(self, parent):
        super().__init__(parent=parent)
        self.mainlay = QHBoxLayout(self)
        self.parent_cw = parent
        self.manager = parent.fig_manager

        # Key Variables
        self.analysis_dir = self.parent_cw.le_analysis_dir.text()
        with open(os.path.join(self.analysis_dir, "DataPar.json")) as parms_reader:
            parms = json.load(parms_reader)
            self.subject_regex = parms["subject_regexp"]
            self.subject_regex_stripped = self.subject_regex.strip("$^")
            self.regex = re.compile(self.subject_regex)

        self.subjects_runs_dict = {}
        self.current_subject = None
        self.current_cbf_img = np.zeros((121, 145, 121))
        self.current_t1_img = np.zeros((121, 145, 121))
        self.get_filenames()

        self.plotting_fig, self.plotting_axes = plt.subplots(nrows=1, ncols=1)
        self.plotting_canvas = FigureCanvas(self.plotting_fig)
        self.plotting_canvas.mpl_connect("pick_event", self.on_pick)
        self.axial_fig, (self.axial_cbf, self.axial_t1) = plt.subplots(nrows=1, ncols=2)
        self.axial_cbf.set_facecolor("black")
        self.axial_t1.set_facecolor("black")
        self.axial_canvas = FigureCanvas(self.axial_fig)
        self.coronal_fig, (self.coronal_cbf, self.coronal_t1) = plt.subplots(nrows=1, ncols=2)
        self.coronal_cbf.set_facecolor("black")
        self.coronal_t1.set_facecolor("black")
        self.coronal_canvas = FigureCanvas(self.coronal_fig)
        self.sagittal_fig, (self.sagittal_cbf, self.sagittal_t1) = plt.subplots(nrows=1, ncols=2)
        self.sagittal_canvas = FigureCanvas(self.sagittal_fig)
        self.sagittal_cbf.set_facecolor("black")
        self.sagittal_t1.set_facecolor("black")

        self.vlay_mrifigs = QVBoxLayout()
        self.vlay_mrifigs.addWidget(self.axial_canvas)
        self.vlay_mrifigs.addWidget(self.coronal_canvas)
        self.vlay_mrifigs.addWidget(self.sagittal_canvas)
        self.mainlay.addWidget(self.plotting_canvas)
        self.mainlay.addLayout(self.vlay_mrifigs)

        self.plotupdate_axialslice(self.manager.slider_axialslice.value())
        self.plotupdate_coronalslice(self.manager.slider_coronalslice.value())
        self.plotupdate_sagittalslice(self.manager.slider_sagittalslice.value())

    def get_filenames(self):
        for subject in os.listdir(self.analysis_dir):
            match = self.regex.search(subject)
            if match:
                for run in os.listdir(os.path.join(self.analysis_dir, subject)):
                    if not os.path.isdir(os.path.join(os.path.join(self.analysis_dir, subject, run))):
                        continue
                    if os.path.exists(os.path.join(self.analysis_dir, subject, run, "CBF.nii")):
                        name = subject + "_" + run
                        t1_file = os.path.join(self.analysis_dir, "Population", f"rT1_{subject}.nii")
                        cbf_file = os.path.join(self.analysis_dir, "Population", f"qCBF_{name}.nii")
                        self.subjects_runs_dict.setdefault(name, {})
                        self.subjects_runs_dict[name]["T1"] = t1_file
                        self.subjects_runs_dict[name]["CBF"] = cbf_file

        pprint(self.subjects_runs_dict)

    def clear_axes(self):
        for ax in self.plotting_fig.axes:
            ax.clear()

    def switch_subject(self, name):
        if name == "Select a subject and run":
            self.current_subject = None
            self.current_cbf_img = np.zeros((121, 145, 121))
            self.current_t1_img = np.zeros((121, 145, 121))
        else:
            try:
                self.current_cbf_img = image.load_img(self.subjects_runs_dict[name]["CBF"]).get_fdata()
                self.current_t1_img = image.load_img(self.subjects_runs_dict[name]["T1"]).get_fdata()
                self.current_subject = name
            except KeyError:
                QMessageBox().warning(self,
                                      "Discrepancy between data and images",
                                      f"An error occurred while attempting to fetch image data for subject_run:\n"
                                      f"{name}. Please check whether see if this subject needs to be re-run",
                                      QMessageBox.Ok)
                self.current_subject = None
                self.current_cbf_img = np.zeros((121, 145, 121))
                self.current_t1_img = np.zeros((121, 145, 121))

        self.plotupdate_axialslice(self.manager.slider_axialslice.value())
        self.plotupdate_coronalslice(self.manager.slider_coronalslice.value())
        self.plotupdate_sagittalslice(self.manager.slider_sagittalslice.value())

    def on_pick(self, event: PickEvent):
        artist: PathCollection = event.artist
        coords = artist.get_offsets()
        idx = event.ind
        coord2label = dict(zip(self.plotting_axes.get_xticks(),
                               [text.get_text() for text in list(self.plotting_axes.get_xticklabels())]
                               ))
        # Index into the 2D matrix of X-coords and Y-coords, then flatten into a list
        coordinates = coords[idx].flatten()
        # print(coordinates)
        # In the case of a stripplot, you must convert the numeric x axis position into the nearest label
        if self.manager.cmb_axestype.currentText() == "Strip Plot":
            parameters = [coord2label[round(coordinates[0])], coordinates[1]]
        else:
            parameters = coordinates.tolist()

        x_var = self.manager.axes_arg_x()
        y_var = self.manager.axes_arg_y()
        # print(f"X Variable: {x_var}; Y Variable: {y_var}")

        df = self.parent_cw.loader.long_data
        subject = df.loc[(df[x_var] == parameters[0]) & (df[y_var] == parameters[1]), "SUBJECT"].tolist()
        if len(subject) == 1:
            cmb_idx = self.manager.cmb_selectsubject.findText(subject[0])
            if cmb_idx == -1:
                return
            self.manager.cmb_selectsubject.setCurrentIndex(cmb_idx)

    ##############################
    # Plotting Functions and Slots
    ##############################
    def plotupdate_axes(self):
        # If x or y axes arguments are still in their non-callable form (i.e string type), don't proceed
        if any([isinstance(self.manager.axes_arg_x, str), isinstance(self.manager.axes_arg_y, str)]):
            return

        # If x and y axes arguments are blank, don't proceed
        if any([self.manager.axes_arg_x() == '', self.manager.axes_arg_y() == '']):
            return

        # Account for a user typing the palette name, accidentally triggering this function, don't proceed
        if any([not self.manager.cmb_palette,
                self.manager.cmb_palette.currentText() not in self.manager.palettenames]):
            return

        print("PLOTUPDATE_AXES TRIGGERED")

        # Clear all axes
        self.clear_axes()

        # Otherwise, proceed
        func = self.manager.plotting_func
        x, y, hue = self.manager.axes_arg_x(), self.manager.axes_arg_y(), self.manager.axes_arg_hue()

        axes_constructor = {}
        for kwarg, call in self.manager.axes_kwargs.items():
            if call() == "":
                axes_constructor[kwarg] = None
            else:
                axes_constructor[kwarg] = call()

        # Account for a user not selecting a palette
        if axes_constructor['palette'] in ['', "None", "Default Blue", "No Palette"]:
            axes_constructor['palette'] = None

        if hue == '':
            self.plotting_axes = func(x=x, y=y, data=self.parent_cw.loader.long_data, ax=self.plotting_axes,
                                      **axes_constructor, picker=4)
        else:
            self.plotting_axes = func(x=x, y=y, hue=hue, data=self.parent_cw.loader.long_data, ax=self.plotting_axes,
                                      **axes_constructor, picker=4)
        self.plotting_canvas.draw()

    @Slot()
    def plotupdate_axescall(self):
        print("PLOTUPDATE_AXESCALL TRIGGERED")
        self.plotupdate_axes()

    def plotupdate_axialslice(self, value):

        self.axial_t1.clear()
        self.axial_t1.imshow(np.rot90(self.current_t1_img[:, :, value]), cmap="gray")
        self.axial_cbf.clear()
        self.axial_cbf.imshow(np.rot90(self.current_cbf_img[:, :, value]), cmap="gray")
        self.axial_canvas.draw()

    def plotupdate_coronalslice(self, value):

        self.coronal_t1.clear()
        self.coronal_t1.imshow(np.rot90(self.current_t1_img[:, value, :]), cmap="gray")
        self.coronal_cbf.clear()
        self.coronal_cbf.imshow(np.rot90(self.current_cbf_img[:, value, :]), cmap="gray")
        self.coronal_canvas.draw()

    def plotupdate_sagittalslice(self, value):

        self.sagittal_t1.clear()
        self.sagittal_t1.imshow(np.rot90(self.current_t1_img[value, :, :]), cmap="gray")
        self.sagittal_cbf.clear()
        self.sagittal_cbf.imshow(np.rot90(self.current_cbf_img[value, :, :]), cmap="gray")
        self.sagittal_canvas.draw()
