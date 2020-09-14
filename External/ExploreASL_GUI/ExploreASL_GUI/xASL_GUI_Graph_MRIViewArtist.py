from PySide2.QtWidgets import *
from PySide2.QtCore import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backend_bases import PickEvent
from matplotlib.collections import PathCollection
import numpy as np
from nilearn import image
from glob import glob
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
        self.current_cbf_max = np.nanmax(self.current_cbf_img)
        self.current_t1_img = np.zeros((121, 145, 121))
        self.current_t1_max = np.nanmax(self.current_t1_img)
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
                    run_path = os.path.join(self.analysis_dir, subject, run)
                    if not os.path.isdir(run_path):
                        print(f"{run_path} is not a directory")
                        continue
                    if len(glob(os.path.join(self.analysis_dir, subject, run, "CBF.nii*"))) > 0:
                        name = subject + "_" + run
                        t1_file = glob(os.path.join(self.analysis_dir, "Population", f"rT1_{subject}.nii*"))
                        cbf_file = glob(os.path.join(self.analysis_dir, "Population", f"qCBF_{name}.nii*"))
                        self.subjects_runs_dict.setdefault(name, {})
                        self.subjects_runs_dict[name]["T1"] = t1_file
                        self.subjects_runs_dict[name]["CBF"] = cbf_file
        print("\nAttempted to locate all subject_run images based on the structure of the analysis directory.\n"
              "This is what was found:")
        pprint(self.subjects_runs_dict)

    def clear_axes(self):
        for ax in self.plotting_fig.axes:
            ax.clear()

    def reset_images_to_default_black(self):
        """
        Resets the current subject back to none and the underlying images as black volumes
        """
        self.current_subject = None
        self.current_cbf_img = np.zeros((121, 145, 121))
        self.current_cbf_max = np.nanmax(self.current_cbf_img)
        self.current_t1_img = np.zeros((121, 145, 121))
        self.current_t1_max = np.nanmax(self.current_t1_img)

    def check_if_files_exist(self, name):
        """
        Wrapper for error detection around not finding a subject
        @param name: The subject_run name extracted from the dataframe (i.e sub_014_ASL_1)
        """
        subject = re.search(pattern=self.subject_regex_stripped, string=name).group()
        t1_path = os.path.join(self.analysis_dir, "Population", f"rT1_{subject}.nii")
        cbf_path = os.path.join(self.analysis_dir, "Population", f"qCBF_{name}.nii")
        check_t1 = glob(t1_path)
        check_cbf = glob(cbf_path)
        if len(check_t1) == 0:
            QMessageBox().warning(self.parent_cw,
                                  "Discrepancy between data and images",
                                  f"An error occurred while attempting to fetch the T1 data for subject {subject}\n"
                                  f"Please check if the following file exists:\n{t1_path}\n"
                                  f"If the file does not exist, you must rerun the analysis for this subject.",
                                  QMessageBox.Ok)
            return
        if len(check_cbf) == 0:
            QMessageBox().warning(self.parent_cw,
                                  "Discrepancy between data and images",
                                  f"An error occurred while attempting to fetch the CBF data for subject_run {name}\n"
                                  f"Please check if the following file exists:\n{cbf_path}\n"
                                  f"If the file does not exist, you must rerun the analysis for this subject & run",
                                  QMessageBox.Ok)
            return

    def switch_subject(self, name):
        """
        Detects that the combobox corresponding to subject_run selection has changed and updates the current_image
        variables appropriately before wrap-calling the plotupdate functions
        @param name: The subject_run name extracted from the dataframe (i.e sub_014_ASL_1)
        """
        if name == "Select a subject and run":
            self.reset_images_to_default_black()
        else:
            self.current_subject = name

            # Get the images and if an error occurs, check whether the files actually exist. If they do not exist,
            # inform the user and reset the images to be black volumes
            try:
                self.current_cbf_img = image.load_img(self.subjects_runs_dict[name]["CBF"]).get_fdata()
                self.current_cbf_max = np.nanmax(self.current_cbf_img)
                self.current_t1_img = image.load_img(self.subjects_runs_dict[name]["T1"]).get_fdata()
                self.current_t1_max = np.nanmax(self.current_t1_img)
            except ValueError:
                self.check_if_files_exist(name=name)
                self.reset_images_to_default_black()
            except KeyError:
                # Try to update the subject runs dict just in case the user has recovered the files
                self.get_filenames()
                try:
                    self.current_cbf_img = image.load_img(self.subjects_runs_dict[name]["CBF"]).get_fdata()
                    self.current_cbf_max = np.nanmax(self.current_cbf_img)
                    self.current_t1_img = image.load_img(self.subjects_runs_dict[name]["T1"]).get_fdata()
                    self.current_t1_max = np.nanmax(self.current_t1_img)
                except KeyError:
                    QMessageBox().warning(self.parent_cw,
                                          "Discrepancy between data and images",
                                          f"An error occurred while attempting to fetch image data for subject_run "
                                          f"{name}\n"
                                          f"It is possible that the images exist in the Population directory, but the "
                                          f"program could not fetch them due to a possible discrepancy between the "
                                          f"Population folder and the rest of the analysis folder.",
                                          QMessageBox.Ok)
                    self.reset_images_to_default_black()

        # Regardless of outcome, tell the canvases to update
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
        # In the case of a stripplot, you must convert the numeric x axis position into the nearest label
        if self.manager.cmb_axestype.currentText() == "Strip Plot":
            parameters = [coord2label[round(coordinates[0])], coordinates[1]]
        else:
            parameters = coordinates.tolist()

        x_var = self.manager.axes_arg_x()
        y_var = self.manager.axes_arg_y()

        df = self.parent_cw.loader.long_data
        subject = df.loc[(df[x_var] == parameters[0]) & (df[y_var] == parameters[1]), "SUBJECT"].tolist()

        # Sometimes if a dtype change has occurred and this is a stipplot, the data may still be numerical, but the
        # graph interpreted the new data correctly as categorical, resulting in strings
        if len(subject) == 0 and self.manager.cmb_axestype.currentText() == "Strip Plot":
            df[x_var] = df[x_var].values.astype(np.str)
            subject = df.loc[(df[x_var] == parameters[0]) & (df[y_var] == parameters[1]), "SUBJECT"].tolist()

        if self.parent_cw.config["DeveloperMode"]:
            print(f"Inside on_pick function in the MRIViewArtist.\n"
                  f"Based on the click made, the following subjects correspond to datapoints nearest to the click:\n"
                  f"{subject}\nSelecting the closest one.\n")

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
        self.axial_t1.imshow(np.rot90(self.current_t1_img[:, :, value]).squeeze(),
                             cmap="gray",
                             vmin=0,
                             vmax=self.current_t1_max
                             )
        self.axial_cbf.clear()
        self.axial_cbf.imshow(np.rot90(self.current_cbf_img[:, :, value]).squeeze(),
                              cmap="gray",
                              vmin=0,
                              vmax=self.current_cbf_max,
                              )
        self.axial_canvas.draw()

    def plotupdate_coronalslice(self, value):

        self.coronal_t1.clear()
        self.coronal_t1.imshow(np.rot90(self.current_t1_img[:, value, :]).squeeze(),
                               cmap="gray",
                               vmin=0,
                               vmax=self.current_t1_max
                               )
        self.coronal_cbf.clear()
        self.coronal_cbf.imshow(np.rot90(self.current_cbf_img[:, value, :]).squeeze(),
                                cmap="gray",
                                vmin=0,
                                vmax=self.current_cbf_max,
                                )
        self.coronal_canvas.draw()

    def plotupdate_sagittalslice(self, value):

        self.sagittal_t1.clear()
        self.sagittal_t1.imshow(np.rot90(self.current_t1_img[value, :, :]).squeeze(),
                                cmap="gray",
                                vmin=0,
                                vmax=self.current_t1_max
                                )
        self.sagittal_cbf.clear()
        self.sagittal_cbf.imshow(np.rot90(self.current_cbf_img[value, :, :]).squeeze(),
                                 cmap="gray",
                                 vmin=0,
                                 vmax=self.current_cbf_max,
                                 )
        self.sagittal_canvas.draw()
