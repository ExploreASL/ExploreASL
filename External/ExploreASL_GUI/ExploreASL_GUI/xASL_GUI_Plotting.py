from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit
from ExploreASL_GUI.xASL_GUI_Graph_Subsetter import xASL_GUI_Subsetter, xASL_GUI_Datatype_Indicator
from ExploreASL_GUI.xASL_GUI_Graph_Loader import xASL_GUI_Data_Loader
from ExploreASL_GUI.xASL_GUI_Graph_FacetManager import xASL_GUI_FacetManager
from ExploreASL_GUI.xASL_GUI_Graph_FacetArtist import xASL_GUI_FacetArtist
from ExploreASL_GUI.xASL_GUI_Graph_MRIViewManager import xASL_GUI_MRIViewManager
from ExploreASL_GUI.xASL_GUI_Graph_MRIViewArtist import xASL_GUI_MRIViewArtist
from ExploreASL_GUI.xASL_GUI_HelperFuncs_StringOps import set_os_dependent_text
import os


# noinspection PyCallingNonCallable
class xASL_Plotting(QMainWindow):
    """
    Central Graphing Widget which will act as a "central hub" for other widgets to be placed in
    """

    def __init__(self, parent_win=None):
        # Parent window is fed into the constructor to allow for communication with parent window devices
        super().__init__(parent=parent_win)
        self.config = self.parent().config

        # Window Size and initial visual setup
        self.setMinimumSize(1920, 1000)
        self.cw = QWidget(self)
        self.setCentralWidget(self.cw)
        self.mainlay = QVBoxLayout(self.cw)
        self.setWindowTitle("Explore ASL - Post Processing Visualization")

        # Central Classes and their connections
        self.dtype_indicator = xASL_GUI_Datatype_Indicator(self)
        self.subsetter = xASL_GUI_Subsetter(self)
        self.loader = xASL_GUI_Data_Loader(self)
        self.loader.signal_dtype_was_changed.connect(self.subsetter.add_or_remove_subsetable_field)

        # Initialize blank givens
        self.fig_manager = None
        self.fig_artist = None
        self.spacer = QSpacerItem(0, 1, QSizePolicy.Preferred, QSizePolicy.Expanding)

        # Main Widgets setup
        self.UI_Setup_Docker()

    def UI_Setup_Docker(self):
        self.dock = QDockWidget(windowTitle="Data Visualization Settings", parent=self.cw)
        self.dock.setMinimumWidth(480)
        self.dock.setFeatures(QDockWidget.AllDockWidgetFeatures)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dock)
        self.cont_maindock = QWidget(self.dock)
        self.vlay_maindock = QVBoxLayout(self.cont_maindock)
        self.dock.setWidget(self.cont_maindock)

        # Set up the directory settings - what analysis folder
        self.grp_directories = QGroupBox(title="Directory settings", parent=self.cont_maindock)
        self.formlay_directories = QFormLayout(self.grp_directories)
        self.hlay_analysis_dir, self.le_analysis_dir, self.btn_analysis_dir = self.make_droppable_clearable_le(
            btn_connect_to=self.set_analysis_dir,
            default="",
            acceptable_path_type="Directory")
        self.le_analysis_dir.setPlaceholderText("Drag & Drop an analysis directory here")
        self.le_analysis_dir.setToolTip("Indicate the analysis directory in this field\n"
                                        "(i.e C:Users\\JohnSmith\\MyStudy\\analysis)")
        self.hlay_metadata, self.le_metadata, self.btn_metadata = self.make_droppable_clearable_le(
            btn_connect_to=self.set_metadata_file,
            default="",
            acceptable_path_type="File",
            supported_extensions=[".tsv", ".csv", ".xlsx"])
        self.le_metadata.setPlaceholderText("Drag & Drop a covariates .tsv/.csv/.xlsx file")
        self.le_metadata.setToolTip("Indicate the filepath to a file containing study metadata.\n"
                                    "This should be either a .csv, .tsv, or .xlsx file")
        self.cmb_atlas_selection = QComboBox()
        self.cmb_atlas_selection.setToolTip("Select which atlas' the CBF values should be based on")
        self.cmb_atlas_selection.addItems(["MNI", "Hammers"])
        self.cmb_pvc_selection = QComboBox()
        self.cmb_pvc_selection.setToolTip("Select whether the CBF values should be based on\n"
                                          "partial-volume correction or uncorrected")
        self.cmb_pvc_selection.addItems(["Without Partial Volume Correction", "With Partial Volume Correction"])
        self.cmb_stats_selection = QComboBox()
        self.cmb_stats_selection.setToolTip("Select which statistic the CBF values should represent")
        self.cmb_stats_selection.addItems(["Mean", "Median", "Coefficient of Variation"])
        self.btn_subset_data = QPushButton("Subset Data", clicked=self.subsetter.show)
        self.btn_subset_data.setToolTip("This will open a window to allow you to subset the loaded data on based on\n"
                                        "levels in a categorical variable")
        self.btn_subset_data.setEnabled(False)
        self.btn_indicate_dtype = QPushButton("Clarify Covariate Datatype", clicked=self.dtype_indicator.show)
        self.btn_indicate_dtype.setToolTip("This will open a window to allow you to correct the program's\n"
                                           "interpretation of the datatype of a particular variable\n(i.e Sex encoded "
                                           "by 1s and 0s may be misinterpreted as a numerical type)")
        self.btn_indicate_dtype.setEnabled(False)
        self.btn_load_in_data = QPushButton("Load Data", self.grp_directories)
        self.btn_load_in_data.setToolTip("This will load in the data from the indicated analysis directory")
        self.btn_load_in_data.clicked.connect(self.loader.load_exploreasl_data)
        self.btn_load_in_data.clicked.connect(self.full_reset)

        self.formlay_directories.addRow("Analysis Directory", self.hlay_analysis_dir)
        self.formlay_directories.addRow("Metadata/Covariates file", self.hlay_metadata)
        self.formlay_directories.addRow("Which Atlas to Utilize", self.cmb_atlas_selection)
        self.formlay_directories.addRow("Which Partial-Volume Stats to View", self.cmb_pvc_selection)
        self.formlay_directories.addRow("Which Statistic to View", self.cmb_stats_selection)
        self.formlay_directories.addRow(self.btn_subset_data)
        self.formlay_directories.addRow(self.btn_indicate_dtype)
        self.formlay_directories.addRow(self.btn_load_in_data)

        # Connect the appropriate lineedits to the subsetter and dtype indicator classes
        self.le_analysis_dir.textChanged.connect(self.subsetter.clear_contents)
        self.le_analysis_dir.textChanged.connect(self.dtype_indicator.clear_contents)
        self.le_metadata.textChanged.connect(self.subsetter.clear_contents)
        self.le_metadata.textChanged.connect(self.dtype_indicator.clear_contents)

        # Setup the main Variable Viewer
        self.grp_varview = QGroupBox(title="Variables", parent=self.cont_maindock)
        self.vlay_varview = QVBoxLayout(self.grp_varview)
        self.lst_varview = QListWidget(self.grp_varview)
        self.lst_varview.setToolTip("This area houses the variables extracted from column headers found within the\n"
                                    "study's dataframes as well as any metadata loaded in by the user")
        self.lst_varview.setFixedHeight(225)
        self.lst_varview.setDragEnabled(True)
        self.vlay_varview.addWidget(self.lst_varview)

        # Setup the start of Plotting Settings
        self.grp_pltsettings = QGroupBox(title="Plotting Settings", parent=self.cont_maindock)
        self.vlay_pltsettings = QVBoxLayout(self.grp_pltsettings)
        self.cmb_figuretypeselection = QComboBox(self.grp_pltsettings)
        self.cmb_figuretypeselection.setToolTip("This will select the type of plotting submodule to use.\n"
                                                "Options include:\n"
                                                "- Facet Grid: for standard plotting (i.e Box Plots)\n"
                                                "- Plot & MRI View: for selecting ASL and T1w scans from datapoints")
        self.cmb_figuretypeselection.addItems(["Select an option", "Facet Grid", "Plot & MRI Viewer"])
        self.cmb_figuretypeselection.setEnabled(False)
        self.cmb_figuretypeselection.currentTextChanged.connect(self.set_manager_and_artist)
        self.vlay_pltsettings.addWidget(self.cmb_figuretypeselection)
        self.vlay_pltsettings.addSpacerItem(self.spacer)

        # Add the groups to the main vertical layout
        self.vlay_maindock.addWidget(self.grp_directories, 1)
        self.vlay_maindock.addWidget(self.grp_varview, 1)
        self.vlay_maindock.addWidget(self.grp_pltsettings, 2)

    def set_manager_and_artist(self, figure_type):
        self.full_reset()

        if figure_type == "Select an option":
            return

        if figure_type == "Facet Grid":
            self.fig_manager = xASL_GUI_FacetManager(self)
            self.fig_artist = xASL_GUI_FacetArtist(self)
        elif figure_type == "Plot & MRI Viewer":
            self.fig_manager = xASL_GUI_MRIViewManager(self)
            # Do NOT proceed if the manager could not be initialized appropriately
            if self.fig_manager.error_init:
                self.cmb_figuretypeselection.setCurrentIndex(0)  # This should reset the fig_manager to None by proxy
                return
            self.fig_artist = xASL_GUI_MRIViewArtist(self)
        else:
            print("set_manager_and_artist: THIS OPTION SHOULD NEVER PRINT")
            return

        self.mainlay.addWidget(self.fig_artist)
        self.vlay_pltsettings.removeItem(self.spacer)
        self.vlay_pltsettings.addWidget(self.fig_manager)
        # Create connections
        self.fig_manager.UI_Setup_ConnectManager2Artist()

    def full_reset(self):
        """
        Necessary full reset in the event that the user loads new data
        """
        if self.fig_manager is not None and self.fig_artist is not None:
            print("FULL RESET ACTIVATED")
            self.mainlay.removeWidget(self.fig_artist)
            self.vlay_pltsettings.removeWidget(self.fig_manager)
            self.vlay_pltsettings.addSpacerItem(self.spacer)
            self.fig_artist.setParent(None)
            self.fig_manager.setParent(None)
            self.fig_artist = None
            self.fig_manager = None

    def set_analysis_dir(self):
        # noinspection PyCallByClass
        analysisdir_filepath = QFileDialog.getExistingDirectory(self,
                                                                "Select the study's analysis directory",
                                                                self.config["DefaultRootDir"],
                                                                QFileDialog.ShowDirsOnly)
        # Return if the user has cancelled the operation
        if analysisdir_filepath == '':
            return

        # Quality control
        if os.path.exists(analysisdir_filepath):
            if any([not os.path.isdir(analysisdir_filepath),  # Path must be a directory
                    os.path.basename(analysisdir_filepath) != "analysis"  # The basename must be the analysis directory
                    ]):
                QMessageBox().warning(self,
                                      "Invalid Directory Selected",
                                      "Either the path you have specified is not a directory or it does not match"
                                      "ExploreASL specification (i.e select the actual 'analysis' directory of your "
                                      "study",
                                      QMessageBox.Ok)
                return
            else:
                self.le_analysis_dir.setText(analysisdir_filepath)
                set_os_dependent_text(linedit=self.le_analysis_dir,
                                      config_ossystem=self.config["Platform"],
                                      text_to_set=analysisdir_filepath)
        else:
            QMessageBox().warning(self,
                                  "The filepath you specified does not exist",
                                  "Please select an existent ExploreASL directory",
                                  QMessageBox.Ok)
            return

    def set_metadata_file(self):
        # noinspection PyCallByClass
        metadata_filepath, _ = QFileDialog.getOpenFileName(self,
                                                           "Select the filepath to your study's metadata",
                                                           self.config["DefaultRootDir"])
        if metadata_filepath == '':
            return

        print(os.path.splitext(metadata_filepath))
        print(not os.path.isfile(metadata_filepath))
        print(os.path.splitext(metadata_filepath)[1] not in [".tsv", ".csv", ".xlsx"])

        # Quality control
        if os.path.exists(metadata_filepath):
            if any([not os.path.isfile(metadata_filepath),
                    os.path.splitext(metadata_filepath)[1] not in [".tsv", ".csv", ".xlsx"]]):
                QMessageBox().warning(self,
                                      "Invalid File Selected",
                                      "The path you have specified is not a file. Supported file types include:\n"
                                      "- Comma Separated Values (.csv)\n"
                                      "- Tab Separated Values (.tsv)\n"
                                      "- Excel Spreadsheets (.xlsx)\n",
                                      QMessageBox.Ok)
                return
            else:
                self.le_metadata.setText(metadata_filepath)
                set_os_dependent_text(linedit=self.le_metadata,
                                      config_ossystem=self.config["Platform"],
                                      text_to_set=metadata_filepath)
        else:
            QMessageBox().warning(self,
                                  "The filepath you specified does not exist",
                                  "Please select an existent filepath",
                                  QMessageBox.Ok)
            return

    ############################################
    # Convenience methods for generating widgets
    ############################################

    @staticmethod
    def make_droppable_clearable_le(le_connect_to=None, btn_connect_to=None, default='', **kwargs):
        hlay = QHBoxLayout()
        le = DandD_FileExplorer2LineEdit(**kwargs)
        le.setText(default)
        le.setClearButtonEnabled(True)
        if le_connect_to is not None:
            le.textChanged.connect(le_connect_to)
        btn = QPushButton("...", )
        if btn_connect_to is not None:
            btn.clicked.connect(btn_connect_to)
        hlay.addWidget(le)
        hlay.addWidget(btn)
        return hlay, le, btn
