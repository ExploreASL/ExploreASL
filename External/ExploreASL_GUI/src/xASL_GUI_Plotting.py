from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
from src.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit
from src.xASL_GUI_Graph_Subsetter import xASL_GUI_Subsetter, xASL_GUI_Datatype_Indicator
from src.xASL_GUI_Graph_Loader import xASL_GUI_Data_Loader
from src.xASL_GUI_Graph_FacetManager import xASL_GUI_FacetManager
from src.xASL_GUI_Graph_FacetArtist import xASL_GUI_FacetArtist
from src.xASL_GUI_Graph_MRIViewManager import xASL_GUI_MRIViewManager
from src.xASL_GUI_Graph_MRIViewArtist import xASL_GUI_MRIViewArtist
from src.xASL_GUI_HelperFuncs_WidgetFuncs import make_scrollbar_area, set_formlay_options
from json import load
from pathlib import Path
from platform import system


# noinspection PyCallingNonCallable
class xASL_Plotting(QMainWindow):
    """
    Central Graphing Widget which will act as a "central hub" for other widgets to be placed in
    """

    def __init__(self, parent_win=None):
        # Parent window is fed into the constructor to allow for communication with parent window devices
        super().__init__(parent=parent_win)
        self.config = self.parent().config

        with open(Path(self.config["ProjectDir"]) / "JSON_LOGIC" / "ToolTips.json") as plot_tips_reader:
            self.plot_tips = load(plot_tips_reader)["Plotting"]

        with open(Path(self.config["ProjectDir"]) / "JSON_LOGIC" / "ErrorsListing.json") as plot_errs_reader:
            self.plot_errs = load(plot_errs_reader)

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

        # MacOS addtional actions
        if system() == "Darwin":
            set_formlay_options(self.formlay_directories, vertical_spacing=0)
            self.hlay_analysis_dir.setContentsMargins(5, 0, 5, 0)
            self.hlay_metadata.setContentsMargins(5, 0, 5, 0)
            self.vlay_directories.setSpacing(0)
            self.vlay_directories.addStretch(1)
            for widget in [self.cmb_figuretypeselection, self.cmb_stats_selection, self.cmb_atlas_selection,
                           self.cmb_pvc_selection]:
                widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)


    def resizeEvent(self, event):
        self.dock.setMaximumHeight(self.height())
        super().resizeEvent(event)

    def UI_Setup_Docker(self):
        self.dock = QDockWidget(windowTitle="Data Visualization Settings", parent=self.cw)
        self.dock.setMinimumWidth(480)
        self.dock.setMaximumHeight(self.height())
        self.dock.setFeatures(QDockWidget.AllDockWidgetFeatures)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dock)
        self.cont_maindock = QWidget(self.dock)
        self.vlay_maindock = QVBoxLayout(self.cont_maindock)
        self.dock.setWidget(self.cont_maindock)

        # Set up the directory settings - what analysis folder
        self.grp_directories = QGroupBox(title="Directory settings", parent=self.cont_maindock)
        (self.vlay_directories, self.scroll_directories,
         self.cont_directories) = make_scrollbar_area(self.grp_directories)
        self.formlay_directories = QFormLayout(self.cont_directories)
        self.formlay_directories.setContentsMargins(2, 2, 2, 2)
        self.hlay_analysis_dir, self.le_analysis_dir, self.btn_analysis_dir = self.make_droppable_clearable_le(
            btn_connect_to=self.set_analysis_dir,
            default="",
            acceptable_path_type="Directory")
        self.le_analysis_dir.setPlaceholderText("Drag & Drop an analysis directory here")
        self.le_analysis_dir.setToolTip(self.plot_tips["le_analysis_dir"])
        self.hlay_metadata, self.le_metadata, self.btn_metadata = self.make_droppable_clearable_le(
            btn_connect_to=self.set_metadata_file,
            default="",
            acceptable_path_type="File",
            supported_extensions=[".tsv", ".csv", ".xlsx"])
        self.le_metadata.setPlaceholderText("Drag & Drop a covariates .tsv/.csv/.xlsx file")
        self.le_metadata.setToolTip(self.plot_tips["le_metadata"])
        self.cmb_atlas_selection = QComboBox()
        self.cmb_atlas_selection.setToolTip(self.plot_tips["cmb_atlas_selection"])
        self.cmb_atlas_selection.addItems(["MNI", "Hammers"])
        self.cmb_pvc_selection = QComboBox()
        self.cmb_pvc_selection.setToolTip(self.plot_tips["cmb_pvc_selection"])
        self.cmb_pvc_selection.addItems(["Without PVC", "With PVC"])
        self.cmb_stats_selection = QComboBox()
        self.cmb_stats_selection.setToolTip("Select which statistic the CBF values should represent")
        self.cmb_stats_selection.addItems(["Mean", "Median", "Coefficient of Variation"])
        self.btn_subset_data = QPushButton("Subset Data", clicked=self.subsetter.show)
        self.btn_subset_data.setToolTip(self.plot_tips["btn_subset_data"])
        self.btn_subset_data.setEnabled(False)
        self.btn_indicate_dtype = QPushButton("Clarify Covariate Datatype", clicked=self.dtype_indicator.show)
        self.btn_indicate_dtype.setToolTip(self.plot_tips["btn_indicate_dtype"])
        self.btn_indicate_dtype.setEnabled(False)
        self.btn_load_in_data = QPushButton("Load Data", self.grp_directories)
        self.btn_load_in_data.setToolTip(self.plot_tips["btn_load_in_data"])
        self.btn_load_in_data.clicked.connect(self.loader.load_exploreasl_data)
        self.btn_load_in_data.clicked.connect(self.full_reset)

        self.formlay_directories.addRow("Study Directory", self.hlay_analysis_dir)
        self.formlay_directories.addRow("Metadata/Covariates file", self.hlay_metadata)
        self.formlay_directories.addRow("Which Atlas", self.cmb_atlas_selection)
        self.formlay_directories.addRow("Which Partial-Volume Statistic", self.cmb_pvc_selection)
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
        self.vlay_varview.setContentsMargins(0, 0, 0, 0)
        self.lst_varview = QListWidget(self.grp_varview)
        self.lst_varview.setToolTip(self.plot_tips["lst_varview"])
        self.lst_varview.setFixedHeight(200)
        self.lst_varview.setDragEnabled(True)
        self.vlay_varview.addWidget(self.lst_varview)

        # Setup the start of Plotting Settings
        self.grp_pltsettings = QGroupBox(title="Plotting Settings", parent=self.cont_maindock)
        self.vlay_pltsettings = QVBoxLayout(self.grp_pltsettings)
        self.vlay_pltsettings.setContentsMargins(0, 0, 0, 0)
        self.cmb_figuretypeselection = QComboBox(self.grp_pltsettings)
        self.cmb_figuretypeselection.setToolTip(self.plot_tips["cmb_figuretypeselection"])
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
        analysisdir_filepath = QFileDialog.getExistingDirectory(self, "Select the study's analysis directory",
                                                                self.config["DefaultRootDir"], QFileDialog.ShowDirsOnly)
        # Return if the user has cancelled the operation
        if analysisdir_filepath == '':
            return

        analysisdir_filepath = Path(analysisdir_filepath)
        if any([not analysisdir_filepath.exists(), not analysisdir_filepath.is_dir(),
                not (analysisdir_filepath / "Population" / "Stats").exists()]):
            QMessageBox().warning(self, self.plot_errs["BadStudyDir"][0],
                                  self.plot_errs["BadStudyDir"][1], QMessageBox.Ok)
            return
        self.le_analysis_dir.setText(str(analysisdir_filepath))

    def set_metadata_file(self):
        # noinspection PyCallByClass
        meta_path, _ = QFileDialog.getOpenFileName(self, "Select the filepath to your study's metadata",
                                                   self.config["DefaultRootDir"])
        if meta_path == '':
            return
        meta_path = Path(meta_path)
        if any([not meta_path.exists(), not meta_path.is_file(), meta_path.suffix in [".csv", ".tsv", ".xlsx"]]):
            QMessageBox().warning(self, self.plot_errs["BadMetaDataFile"][0],
                                  self.plot_errs["BadMetaDataFile"][1], QMessageBox.Ok)
            return
        self.le_metadata.setText(str(meta_path))

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
