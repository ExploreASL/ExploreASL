from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
from ExploreASL_GUI.xASL_GUI_Parms import xASL_Parms
from ExploreASL_GUI.xASL_GUI_Executor import xASL_Executor
from ExploreASL_GUI.xASL_GUI_Plotting import xASL_Plotting
from ExploreASL_GUI.xASL_GUI_Importer import xASL_GUI_Importer
from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit
from ExploreASL_GUI.xASL_GUI_HelperFuncs_StringOps import set_os_dependent_text
from ExploreASL_GUI.xASL_GUI_FileExplorer import xASL_FileExplorer
import os
import sys
import json
import platform
import subprocess
import re


# Explore ASL Main Window
# by Maurice Pasternak @ 2020
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# Quick Preface of function nomenclature for future developers
# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# UI_   --> the function is associated with building the visual interface of the GUI or connecting widgets with each
#           other; there should be minimal setting of functional attributes, processing, etc. being done here
# get_  --> the function is associated with retrieving a value
# set_  --> the function is associated with setting a value; fields and attributes are manipulated in the process
# load_ --> the function is associated with retrieving data from a file; fields and attributes may or may not be set
#           in the process
# save_ --> the function is associated with saving parameters and variable to a file
# show_ --> the function is associated with making certain widgets or widget elements visible
# on_off_ --> the function is associated purely with enabling/disabling widgets or their elements
#
# Some helpful prefix shortcuts to know:
#  act = Action; btn = PushButton; cmb = ComboBox; cont = Widget/container; dock = DockWidget; formlay = Form Layout
#  hlay = horizontal Layout; le = LineEdit; lst = ListWidget;
#  spinbox = DoubleSpinBox; vlay = vertical Layout


# noinspection PyAttributeOutsideInit
class xASL_MainWin(QMainWindow):
    def __init__(self, config):
        super().__init__()
        # Load in the master config file if it exists; otherwise, make it
        self.config = config
        self.save_config()

        self.load_tooltips()
        # Window Size and initial visual setup
        # self.setWindowFlags(Qt.WindowStaysOnTopHint)
        self.resize(self.config["ScreenSize"][0] // 2, self.config["ScreenSize"][1] // 2.5)
        self.cw = QWidget(self)
        self.setCentralWidget(self.cw)
        # Main Icon setup
        self.icon_main = QIcon(os.path.join(self.config["ProjectDir"], "media", "ExploreASL_logo.png"))
        self.setWindowIcon(self.icon_main)
        # Main Layout Setup
        self.mainlay = QGridLayout(self.cw)
        self.setWindowTitle("Explore ASL GUI")

        # Misc Players
        self.file_explorer = xASL_FileExplorer(self)

        # Set up each of the subcomponents of the main window program
        self.UI_Setup_Navigator()
        self.UI_Setup_MenuBar()
        self.UI_Setup_MainSelections()
        self.UI_Setup_Connections()

        # Pre-initialize the main players
        self.parmsmaker = xASL_Parms(self)
        self.executor = xASL_Executor(self)
        self.plotter = xASL_Plotting(self)
        self.importer = xASL_GUI_Importer(self)

    # This dockable navigator will contain the most essential parameters and will be repeatedly accessed by other
    # subwidgets within the program; should also be dockable within any of them.
    def UI_Setup_Navigator(self):
        # The main container and the main layout of the dock
        self.cont_navigator = QWidget()
        self.vlay_navigator = QVBoxLayout(self.cont_navigator)

        # Essential Widgets
        # First, the lineedit for the ExploreASL directory
        self.le_exploreasl_dir = DandD_FileExplorer2LineEdit(acceptable_path_type="Directory")
        set_os_dependent_text(linedit=self.le_exploreasl_dir,
                              config_ossystem=self.config["Platform"],
                              text_to_set=self.config["ExploreASLRoot"])
        self.le_exploreasl_dir.textChanged.connect(self.set_exploreasl_dir)
        self.le_exploreasl_dir.setReadOnly(True)
        self.le_exploreasl_dir.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.btn_exploreasl_dir = QPushButton("...", clicked=self.set_exploreasl_dir_frombtn)
        self.hlay_exploreasl_dir = QHBoxLayout()
        self.hlay_exploreasl_dir.addWidget(self.le_exploreasl_dir)
        self.hlay_exploreasl_dir.addWidget(self.btn_exploreasl_dir)

        # Second, the lineedit for the study analysis directory
        self.le_currentanalysis_dir = DandD_FileExplorer2LineEdit(acceptable_path_type="Directory")
        set_os_dependent_text(linedit=self.le_currentanalysis_dir,
                              config_ossystem=self.config["Platform"],
                              text_to_set=self.config["DefaultRootDir"])
        self.le_currentanalysis_dir.textChanged.connect(self.set_analysis_dir)
        self.le_currentanalysis_dir.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.btn_currentanalysis_dir = QPushButton("...", self.cont_navigator, clicked=self.set_analysis_dir_frombtn)
        self.hlay_currentanalysis_dir = QHBoxLayout()
        self.hlay_currentanalysis_dir.addWidget(self.le_currentanalysis_dir)
        self.hlay_currentanalysis_dir.addWidget(self.btn_currentanalysis_dir)

        # Third, the checkbox of whether to make the analysis directory list as the default display in other modules
        self.chk_makedefault_analysisdir = QCheckBox(self.cont_navigator)
        self.chk_makedefault_analysisdir.setChecked(True)
        # These aforementioned widgets will be packaged into a form layout
        self.formlay_navigator = QFormLayout()
        self.formlay_navigator.addRow("Explore ASL Directory", self.hlay_exploreasl_dir)
        self.formlay_navigator.addRow("Current Analysis Directory", self.hlay_currentanalysis_dir)
        self.formlay_navigator.addRow("Set selected directory as default?", self.chk_makedefault_analysisdir)

        self.vlay_navigator.addLayout(self.formlay_navigator)
        self.vlay_navigator.addWidget(self.file_explorer)

        # Finally add this whole widget group to the main layout
        self.mainlay.addWidget(self.cont_navigator, 0, 0, 2, 2)
        self.cont_navigator.setMinimumWidth(600)
        self.cont_navigator.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

    # Setup the standard Menu
    def UI_Setup_MenuBar(self):
        # Main menubar and main submenus setup
        self.menubar_main = QMenuBar(self)
        menu_names = ["File", "Edit", "Settings", "About"]
        self.menu_file, self.menu_edit, self.menu_settings, self.menu_about = [self.menubar_main.addMenu(name) for name
                                                                               in menu_names]
        # Setup the actions of the File menu
        # self.menu_file.addAction("Show Navigator", self.dock_navigator.show)
        self.menu_file.addAction("Save Master Config", self.save_config)
        self.menu_file.addAction("Select Analysis Directory", self.set_analysis_dir_frombtn)
        self.menu_file.addAction("About ExploreASL", self.show_AboutExploreASL)
        # Setup the actions of the Edit menu
        # Setup the actions of the Settings menu
        self.setMenuBar(self.menubar_main)

    # Setup the main selection
    def UI_Setup_MainSelections(self):

        expanding_policy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        font_policy = QFont()
        font_policy.setPointSize(16)
        font_policy.setBold(True)

        # Set up the main buttons
        self.btn_open_runexploreasl = QToolButton(self.cw)
        self.btn_open_parmsmaker = QToolButton(self.cw)
        self.btn_open_importer = QToolButton(self.cw)
        self.btn_open_postanalysis = QToolButton(self.cw)
        self.main_btns = [self.btn_open_importer,
                          self.btn_open_runexploreasl,
                          self.btn_open_parmsmaker,
                          self.btn_open_postanalysis]

        imgs = [os.path.join(self.config["ProjectDir"], "media", "importer_icon.svg"),
                os.path.join(self.config["ProjectDir"], "media", "run_exploreasl_icon.svg"),
                os.path.join(self.config["ProjectDir"], "media", "parmsmaker_icon.svg"),
                os.path.join(self.config["ProjectDir"], "media", "postrun_analysis_icon.svg")]
        txt_labels = ["Import DCM to NIFTI",
                      "Process Data",
                      "Create Parameters File",
                      "Post-Run Analysis"]

        for btn, img, text in zip(self.main_btns, imgs, txt_labels):
            btn.setSizePolicy(expanding_policy)
            btn.setIcon(QIcon(img))
            btn.setIconSize(QSize(175, 175))
            btn.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
            btn.setText(text)
            btn.setFont(font_policy)

        # Add the buttons to the layout
        self.mainlay.addWidget(self.btn_open_importer, 0, 3, 1, 1)
        self.mainlay.addWidget(self.btn_open_parmsmaker, 0, 4, 1, 1)
        self.mainlay.addWidget(self.btn_open_runexploreasl, 1, 3, 1, 1)
        self.mainlay.addWidget(self.btn_open_postanalysis, 1, 4, 1, 1)

    def resizeEvent(self, event):
        super(xASL_MainWin, self).resizeEvent(event)
        for btn in self.main_btns:
            adjusted = min([self.width(), self.height()])
            btn.setIconSize(QSize(adjusted//2.8, adjusted//2.8))

    def UI_Setup_Connections(self):
        actions = [self.show_Importer,
                   self.show_Executor,
                   self.show_ParmsMaker,
                   self.show_PostAnalysis]
        for btn, action in zip(self.main_btns, actions):
            btn.clicked.connect(action)

    # Lauch the Executor window
    def show_Executor(self):
        self.executor.show()

    # Launch the ParmsMaker window
    def show_ParmsMaker(self):
        self.parmsmaker.show()

    # Launch the Importer window
    def show_Importer(self):
        self.importer.show()

    # Launch the Post-Analysis window
    def show_PostAnalysis(self):
        self.plotter.show()

    # Launch the "About Explore ASL" window
    def show_AboutExploreASL(self):
        pass

    # Convenience function for saving to the master configuration file
    def save_config(self):
        dst = os.path.join(self.config["ProjectDir"], "JSON_LOGIC", "ExploreASL_GUI_masterconfig.json")
        with open(dst, "w") as master_config_writer:
            json.dump(self.config, master_config_writer, indent=1)
        del master_config_writer, dst

    def load_tooltips(self):
        tooltips_src = os.path.join(self.config["ProjectDir"], "JSON_LOGIC", "xASL_GUI_Tooltips.json")
        if os.path.exists(tooltips_src):
            with open(tooltips_src) as tooltips_reader:
                self.tooltips = json.load(tooltips_reader)
        else:
            self.tooltips = None

    # Sets the analysis directory the user is interested in from the push button
    def set_analysis_dir_frombtn(self):
        directory: str = QFileDialog.getExistingDirectory(QFileDialog(),
                                                          "Select analysis directory to view",  # Window title
                                                          self.config["DefaultRootDir"],  # Default dir
                                                          QFileDialog.ShowDirsOnly)  # Display options
        self.set_analysis_dir(directory)

    # The actual function that sets the analysis directory
    def set_analysis_dir(self, directory):
        # Abort if the provided directory is not a filepath
        if not os.path.exists(directory):
            return
        # Abort if the provided directory is actually not a directory
        if not os.path.isdir(directory):
            return
        # Change the lineedit
        set_os_dependent_text(linedit=self.le_currentanalysis_dir,
                              config_ossystem=self.config["Platform"],
                              text_to_set=directory)
        # Update user config to have a new default analysis dir to refer to in on future startups
        if self.chk_makedefault_analysisdir.isChecked():
            self.config["DefaultRootDir"] = directory
            self.save_config()

    # Sets the ExploreASL directory of the user from the push button
    def set_exploreasl_dir_frombtn(self):
        directory: str = QFileDialog.getExistingDirectory(self,
                                                          "Select the path to ExploreASL",  # Window title
                                                          os.getcwd(),  # Default dir
                                                          QFileDialog.ShowDirsOnly)  # Display options
        self.set_exploreasl_dir(directory)

    def set_exploreasl_dir(self, directory):
        # Abort if the provided directory is not a filepath
        if not os.path.exists(directory):
            return
        # Abort if the provided directory is actually not a directory
        if not os.path.isdir(directory):
            return
        # Change the lineedit
        set_os_dependent_text(linedit=self.le_exploreasl_dir,
                              config_ossystem=self.config["Platform"],
                              text_to_set=directory)
        # Update user config to have a new ExploreASL directory specified
        self.config["ExploreASLRoot"] = directory
        self.save_config()
