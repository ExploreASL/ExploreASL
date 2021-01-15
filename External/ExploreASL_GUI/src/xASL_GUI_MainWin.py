from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
from src.xASL_GUI_Parms import xASL_Parms
from src.xASL_GUI_Executor import xASL_Executor
from src.xASL_GUI_Plotting import xASL_Plotting
from src.xASL_GUI_Importer import xASL_GUI_Importer
from src.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit
from src.xASL_GUI_FileExplorer import xASL_FileExplorer
import os
from pathlib import Path
import json


# Explore ASL Main Window
# by Maurice Pasternak @ 2020, 2021
class xASL_MainWin(QMainWindow):
    def __init__(self, config):
        super().__init__()
        # Load in the master config file if it exists; otherwise, make it
        self.config = config

        # Window Size and initial visual setup
        # self.setWindowFlags(Qt.WindowStaysOnTopHint)
        self.resize(self.config["ScreenSize"][0] // 2.5, self.config["ScreenSize"][1] // 2.25)
        self.cw = QWidget(self)
        self.setCentralWidget(self.cw)
        # Main Icon setup
        self.icon_main = QIcon(str(Path(self.config["ProjectDir"]) / "media" / "ExploreASL_logo.png"))
        self.setWindowIcon(self.icon_main)
        # Main Layout Setup
        self.mainlay = QVBoxLayout(self.cw)
        self.mainlay.setContentsMargins(0, 0, 0, 0)
        self.setWindowTitle("Explore ASL GUI")

        # Misc Players
        self.file_explorer = xASL_FileExplorer(self)

        # Set up each of the subcomponents of the main window program
        self.UI_Setup_Navigator()
        self.le_defaultdir.setCompleter(self.file_explorer.completer_current_dir)

        # Pre-initialize the main players
        self.parmsmaker = xASL_Parms(self)
        self.executor = xASL_Executor(self)
        self.plotter = xASL_Plotting(self)
        self.importer = xASL_GUI_Importer(self)

        # Menubar has to come after, as it references the main widgets
        self.UI_Setup_MenuBar()
        self.UI_Setup_ToolBar()

        # Save the configuration
        self.save_config()

    # This dockable navigator will contain the most essential parameters and will be repeatedly accessed by other
    # subwidgets within the program; should also be dockable within any of them.
    def UI_Setup_Navigator(self):
        # The main container and the main layout of the dock
        self.cont_navigator = QWidget()
        self.vlay_navigator = QVBoxLayout(self.cont_navigator)

        # Essential Widgets
        # The lineedit for the study analysis directory
        self.le_defaultdir = DandD_FileExplorer2LineEdit(acceptable_path_type="Directory")
        self.le_defaultdir.setText(str(Path(self.config["DefaultRootDir"])))
        self.le_defaultdir.setClearButtonEnabled(True)
        self.le_defaultdir.textChanged.connect(self.set_analysis_dir)
        self.le_defaultdir.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.le_defaultdir.setToolTip("The default directory that will be suggested when searching for a file or "
                                      "folder")
        self.btn_defaultdir = QPushButton("...", self.cont_navigator, clicked=self.set_analysis_dir_frombtn,
                                          toolTip="Set the default directory")
        self.hlay_defaultdir = QHBoxLayout()
        self.hlay_defaultdir.addWidget(self.le_defaultdir)
        self.hlay_defaultdir.addWidget(self.btn_defaultdir)

        # These aforementioned widgets will be packaged into a form layout
        self.formlay_defaultdir = QFormLayout()
        self.formlay_defaultdir.addRow("Default Directory", self.hlay_defaultdir)

        # Add the form layout to the starft of the file explorer after it has been added. It isn't defined within the
        # file explorer, as it interfaces with setting the config in this
        self.vlay_navigator.addWidget(self.file_explorer)
        self.file_explorer.mainlay.insertLayout(0, self.formlay_defaultdir)
        self.vlay_navigator.setContentsMargins(0, 0, 0, 0)

        # Finally add this whole widget group to the main layout
        self.mainlay.addWidget(self.cont_navigator)
        self.cont_navigator.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

    def UI_Setup_MenuBar(self):
        # Main menubar and main submenus setup
        self.menubar_main = QMenuBar(self)
        self.setMenuBar(self.menubar_main)
        menu_names = ["File", "Modules", "About"]
        self.menu_file, self.menu_modules, self.menu_about = [self.menubar_main.addMenu(name) for name in menu_names]

        # Setup the actions of the File menu
        self.menu_file.addAction("Save Master Config", self.save_config)
        self.menu_file.addAction("Select Default Study Directory", self.set_analysis_dir_frombtn)
        self.menu_file.addSeparator()
        self.menu_file.addAction("Exit", self.close)

        # Setup the actions of the Modules menu
        self.menu_modules.addAction("Import Module", self.importer.show)
        self.menu_modules.addAction("Parameters Module", self.parmsmaker.show)
        self.menu_modules.addAction("ExploreASL Module", self.executor.show)
        self.menu_modules.addAction("Post-Processing Module", self.plotter.show)

    def UI_Setup_ToolBar(self):
        self.toolbar = QToolBar(self, orientation=Qt.Vertical, iconSize=QSize(75, 75), allowedAreas=Qt.AllToolBarAreas)
        self.addToolBar(Qt.LeftToolBarArea, self.toolbar)
        self.toolbar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        self.toolactions = {}  # Hack; keeps the toolbuttons in memory so that they can be added and remain visible
        media_path = Path(self.config["ProjectDir"]) / "media"
        paths = [media_path / "importer_icon.svg", media_path / "parmsmaker_icon.svg",
                 media_path / "run_exploreasl_icon.svg", media_path / "postrun_analysis_icon.svg"]
        funcs = [self.importer.show, self.parmsmaker.show, self.executor.show, self.plotter.show]
        descs = ["Import DCM to NIFTI", "Define Study Parameters", "Run ExploreASL", "Visualize Results"]
        tips = ["Convert a source directory containing DICOM images into NIFTI format for ExploreASL processing",
                "Define the parameters relevant to the study being analyzed and create the required DataPar.m file",
                "Analyze the data and/or modify analysis directories for specific re-runs",
                "Visualize and create figures for an analyzed dataset"]
        for path, func, desc, tip in zip(paths, funcs, descs, tips):
            self.toolactions[path] = QAction(icon=QIcon(str(path)), triggered=func, iconText=desc, toolTip=tip)
            self.toolbar.addAction(self.toolactions[path])

    # Convenience function for saving to the master configuration file
    def save_config(self):
        dst = Path(self.config["ProjectDir"]) / "JSON_LOGIC" / "ExploreASL_GUI_masterconfig.json"
        with open(dst, "w") as master_config_writer:
            json.dump(self.config, master_config_writer, indent=1)
        self.communicate_config_change()
        del master_config_writer, dst

    # Sets the analysis directory the user is interested in from the push button
    def set_analysis_dir_frombtn(self):
        directory: str = QFileDialog.getExistingDirectory(QFileDialog(),
                                                          "Select analysis directory to view",  # Window title
                                                          self.config["DefaultRootDir"],  # Default dir
                                                          QFileDialog.ShowDirsOnly)  # Display options
        if directory == "":
            return
        directory = Path(directory)
        self.le_defaultdir.setText(str(directory))
        # Update user config to have a new default analysis dir to refer to in on future startups
        self.config["DefaultRootDir"] = str(directory)
        self.save_config()

    # The actual function that sets the analysis directory
    def set_analysis_dir(self, directory: str):
        if any([directory == "", not Path(directory).exists(), not Path(directory).is_dir()]):
            return
        else:
            print(f"Updating self.le_deafultdir to {directory}")
            # Update user config to have a new default analysis dir to refer to in on future startups
            self.config["DefaultRootDir"] = str(directory)
            self.save_config()

    def communicate_config_change(self):
        self.parmsmaker.config = self.config
        self.file_explorer.config = self.config
        self.plotter.config = self.config
        self.executor.config = self.config
        self.importer.config = self.config
        self.importer.dehybridizer.config = self.config

    # This will be modified in the future to perform certain actions on end
    def closeEvent(self, event):
        super(xASL_MainWin, self).closeEvent(event)
