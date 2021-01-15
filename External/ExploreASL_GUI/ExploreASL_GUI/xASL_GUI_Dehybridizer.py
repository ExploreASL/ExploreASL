from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *

from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit
from ExploreASL_GUI.xASL_GUI_HelperFuncs_StringOps import set_os_dependent_text
from ExploreASL_GUI.xASL_GUI_AnimationClasses import xASL_ImagePlayer
import os
import shutil
import re
import sys
from glob import iglob, glob
from pprint import pprint
from platform import system
from more_itertools import divide


# noinspection PyAttributeOutsideInit
class xASL_GUI_Dehybridizer(QWidget):
    signal_startexpanding = Signal()  # If user opts not to make a backup, this will just manually force it to start

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowFlag(Qt.Window)
        self.parent_cw = parent
        self.mainlay = QVBoxLayout(self)
        if parent is not None:
            self.config = parent.config

        # Misc Variables
        self.leftside_name2le = {}
        self.rightside_name2le = {}
        self.basenames = {}
        self.documented_sourcedir = ""
        self.documented_exampledir = ""
        self.glob_path = ""
        self.threadpool = QThreadPool()
        self.paths = []
        self.translator = {}
        self.global_debt = 0

        # The setup
        self.SetupUI_UserSpecifyHybridDir()
        self.SetupUI_LeftandRightSections()
        self.expanddirs_movieplayer = xASL_ImagePlayer(os.path.join(self.config["ProjectDir"],
                                                                    "media",
                                                                    "processing.gif"))

        self.chk_makebackup = QCheckBox(text="(STRONGLY RECOMMENDED) Make a backup of the source directory "
                                             "before folder restructuring?", checked=True)

        self.btn_run_dehybridizer = QPushButton("Expand Directory Level", clicked=self.run_dehybridizer)
        self.btn_run_dehybridizer.setFixedHeight(50)
        self.btn_run_dehybridizer.setEnabled(False)
        self.mainlay.addWidget(self.chk_makebackup)
        self.mainlay.addWidget(self.btn_run_dehybridizer)

        # Connect signals
        self.signal_startexpanding.connect(self.expand_directories)

    def SetupUI_UserSpecifyHybridDir(self):
        # First, specify the root directory
        self.cont_basedirs = QWidget()
        self.formlay_basedirs = QFormLayout(self.cont_basedirs)
        self.hlay_rootdir = QHBoxLayout()
        self.le_rootdir = DandD_FileExplorer2LineEdit(acceptable_path_type="Directory")
        self.le_rootdir.setPlaceholderText("Drag and drop your study's raw directory here")
        self.le_rootdir.setToolTip("The path to the root directory that will have a backup made prior to an expansion\n"
                                   "and which tells the program where to begin looking. Usually this is the /raw/\n"
                                   "folder of your study")
        self.le_rootdir.setReadOnly(True)
        self.le_rootdir.textChanged.connect(self.can_specify_components)
        self.le_rootdir.textChanged.connect(self.can_dehybridize)
        self.btn_setrootdir = QPushButton("...", clicked=self.set_source_root_directory)
        self.btn_setrootdir.setToolTip("Specify the path to the root directory through a file dialogue rather than by\n"
                                       "drag and drop")
        self.hlay_rootdir.addWidget(self.le_rootdir)
        self.hlay_rootdir.addWidget(self.btn_setrootdir)

        # Second, specify the example directory
        self.hlay_exampledir = QHBoxLayout()
        self.le_exampledir = DandD_FileExplorer2LineEdit(acceptable_path_type="Directory")
        self.le_exampledir.setToolTip("Specify any directory that is an example of what must be expanded. This tells\n"
                                      "the program the relative folder depth at which sibling folders can be found")
        self.le_exampledir.setPlaceholderText("Drag and drop any directory existent at target depth")
        self.le_exampledir.setReadOnly(True)
        self.le_exampledir.textChanged.connect(self.can_specify_components)
        self.le_exampledir.textChanged.connect(self.can_dehybridize)
        self.btn_exampledir = QPushButton("...", clicked=self.set_example_directory)
        self.btn_exampledir.setToolTip("Specify the path to the example directory through a file dialogue rather than\n"
                                       "by drag and drop")
        self.hlay_exampledir.addWidget(self.le_exampledir)
        self.hlay_exampledir.addWidget(self.btn_exampledir)

        # Third, specify the delimiter
        self.le_delimiter = QLineEdit(placeholderText="(OPTIONAL) Specify the delimiter, if any")
        self.le_delimiter.setToolTip("Specify whether there is a character that separates two pieces of information\n"
                                     "found at the indicated directory. For example subject_session is separated by\n"
                                     "a '_' delimiter")
        self.le_delimiter.setClearButtonEnabled(True)

        # Fourth, specify the button to search fill the delimiter spaces
        self.btn_fill_spaces = QPushButton("Specify components at indicated directory level",
                                           clicked=self.specify_components)
        self.btn_fill_spaces.setToolTip("Click to extract all relevant folder names at the indicated level,\n"
                                        "separating names as appropriate if a delimiter was specified")

        # Add to layout
        self.formlay_basedirs.addRow("Source Directory", self.hlay_rootdir)
        self.formlay_basedirs.addRow("Example Directory", self.hlay_exampledir)
        self.formlay_basedirs.addRow("Delimiter", self.le_delimiter)
        self.formlay_basedirs.addRow(self.btn_fill_spaces)

        self.mainlay.addWidget(self.cont_basedirs)

    def SetupUI_LeftandRightSections(self):
        self.grp_delimiter_sections = QGroupBox(title="Split Components")
        self.hlay_delimiter_section = QHBoxLayout(self.grp_delimiter_sections)

        # Left Side
        self.grp_leftside = QGroupBox(title="Keep Left Component (default if no delimiter)", checkable=True)
        self.vlay_leftside = QVBoxLayout(self.grp_leftside)
        self.vlay_leftside.setContentsMargins(0, 0, 0, 0)
        self.hlay_leftautofill = QHBoxLayout()
        self.le_leftautofill = QLineEdit(placeholderText="Fill all rows with this text", clearButtonEnabled=True)
        self.le_leftautofill.setToolTip("Enter text in this field that you wish to be present in all rows.")
        self.btn_leftautofill = QPushButton("Fill", clicked=self.autofill_leftside)
        self.hlay_leftautofill.addWidget(self.le_leftautofill)
        self.hlay_leftautofill.addWidget(self.btn_leftautofill)
        self.scroll_leftside = QScrollArea()
        self.cont_leftside = QWidget()
        self.scroll_leftside.setWidget(self.cont_leftside)
        self.scroll_leftside.setWidgetResizable(True)
        self.formlay_leftside = QFormLayout(self.cont_leftside)
        self.vlay_leftside.addLayout(self.hlay_leftautofill)
        self.vlay_leftside.addWidget(self.scroll_leftside)

        # Right Side
        self.grp_rightside = QGroupBox(title="Keep Right Component", checkable=True)
        self.vlay_rightside = QVBoxLayout(self.grp_rightside)
        self.vlay_rightside.setContentsMargins(0, 0, 0, 0)
        self.hlay_rightautofill = QHBoxLayout()
        self.le_rightautofill = QLineEdit(placeholderText="Fill all rows with this text", clearButtonEnabled=True)
        self.le_rightautofill.setToolTip("Enter text in this field that you wish to be present in all rows.")
        self.btn_rightautofill = QPushButton("Fill", clicked=self.autofill_rightside)
        self.hlay_rightautofill.addWidget(self.le_rightautofill)
        self.hlay_rightautofill.addWidget(self.btn_rightautofill)
        self.scroll_rightside = QScrollArea()
        self.cont_rightside = QWidget()
        self.scroll_rightside.setWidget(self.cont_rightside)
        self.scroll_rightside.setWidgetResizable(True)
        self.formlay_rightside = QFormLayout(self.cont_rightside)
        self.vlay_rightside.addLayout(self.hlay_rightautofill)
        self.vlay_rightside.addWidget(self.scroll_rightside)

        # Add the sides to the hlay
        self.hlay_delimiter_section.addWidget(self.grp_leftside)
        self.hlay_delimiter_section.addWidget(self.grp_rightside)

        self.mainlay.addWidget(self.grp_delimiter_sections)

    def set_source_root_directory(self):
        dir_path = QFileDialog.getExistingDirectory(QFileDialog(),
                                                    "Select the raw directory of your study",
                                                    self.config["DefaultRootDir"],
                                                    QFileDialog.ShowDirsOnly)
        if os.path.exists(dir_path):
            set_os_dependent_text(linedit=self.le_rootdir,
                                  config_ossystem=self.config["Platform"],
                                  text_to_set=dir_path)

    def set_example_directory(self):
        dir_path = QFileDialog.getExistingDirectory(QFileDialog(),
                                                    "Select the raw directory of your study",
                                                    self.config["DefaultRootDir"],
                                                    QFileDialog.ShowDirsOnly)
        if os.path.exists(dir_path) and self.le_rootdir.text().replace('\\', '/') in dir_path.replace('\\', '/'):
            set_os_dependent_text(linedit=self.le_exampledir,
                                  config_ossystem=self.config["Platform"],
                                  text_to_set=dir_path)

    def can_specify_components(self):
        # First check, lineedits for the source and the example must not be blank
        if any([self.le_rootdir.text() == "",
                not os.path.exists(self.le_rootdir.text()),
                self.le_exampledir.text() == "",
                not os.path.exists(self.le_exampledir.text())
                ]):
            self.btn_fill_spaces.setEnabled(False)
            return

        # Make sure that they are also real directories
        if any([not os.path.isdir(self.le_rootdir.text()),
                not os.path.isdir(self.le_exampledir.text())
                ]):
            self.btn_fill_spaces.setEnabled(False)
            return

        # Make sure that the root directory is contained within the example directory
        if self.le_rootdir.text().replace('\\', '/') not in self.le_exampledir.text().replace('\\', '/'):
            self.btn_fill_spaces.setEnabled(False)
            return

        self.btn_fill_spaces.setEnabled(True)

    def can_dehybridize(self):
        # First check, filepaths should not have changed ever since the components were specified
        if any([self.le_rootdir.text() != self.documented_sourcedir,
                self.le_exampledir.text() != self.documented_exampledir
                ]):
            self.btn_run_dehybridizer.setEnabled(False)
            return

        # Second check, all lineedits in the checked sections must be filled
        if self.grp_leftside.isChecked():
            for le in self.leftside_name2le.values():
                if le.text() == "":
                    self.btn_run_dehybridizer.setEnabled(False)
                    return
        if self.grp_rightside.isChecked():
            for le in self.rightside_name2le.values():
                if le.text() == "":
                    self.btn_run_dehybridizer.setEnabled(False)
                    return

        self.btn_run_dehybridizer.setEnabled(True)

    def add_components_to_layout(self, component: str, side: str):
        if side not in ["left", "right"]:
            raise ValueError("Could not determine which side this component should be added to")
        if side == "left":
            self.leftside_name2le[component] = QLineEdit(placeholderText="Specify parent directory name",
                                                         clearButtonEnabled=True)
            self.leftside_name2le[component].setToolTip("Specify the name of the parent folder that the indicated\n "
                                                        "child folder should be placed under")
            self.leftside_name2le[component].textChanged.connect(self.can_dehybridize)
            self.formlay_leftside.addRow(component, self.leftside_name2le[component])
        else:
            self.rightside_name2le[component] = QLineEdit(placeholderText="Specify parent directory name",
                                                          clearButtonEnabled=True)
            self.rightside_name2le[component].setToolTip("Specify the name of the parent folder that the indicated\n "
                                                         "child folder should be placed under")
            self.rightside_name2le[component].textChanged.connect(self.can_dehybridize)
            self.formlay_rightside.addRow(component, self.rightside_name2le[component])

    def clear_components_from_layout(self, side: str):
        if side not in ["left", "right"]:
            raise ValueError("Could not determine which side this component should be added to")
        if side == "left":
            self.leftside_name2le.clear()
            for _ in range(self.formlay_leftside.rowCount()):
                self.formlay_leftside.removeRow(0)
        else:
            self.rightside_name2le.clear()
            for _ in range(self.formlay_rightside.rowCount()):
                self.formlay_rightside.removeRow(0)

    def autofill_leftside(self):
        for le in self.leftside_name2le.values():
            le.setText(self.le_leftautofill.text())

    def autofill_rightside(self):
        for le in self.rightside_name2le.values():
            le.setText(self.le_rightautofill.text())

    def specify_components(self):
        """
        Extracts the basenames at the indicated level, splits on the delimiter if indicated, and populates the left
        and/or right form layouts with the split basenames, allowing the user to start specifying the parent directory
        each component of a basename should be placed under
        """
        self.glob_path = ""

        if not self.le_rootdir.text() in self.le_exampledir.text():
            return

        # First step, we must determine the depth that this example directory exists at
        else:
            from_root = self.le_exampledir.text().replace(self.le_rootdir.text(), "")
            if system() == "Windows":
                from_root = from_root.strip("\\")
                nlevels = len(from_root.split("\\"))
            else:
                from_root = from_root.strip("/")
                nlevels = len(from_root.split("/"))

        if nlevels < 1:
            return

        # Second step, we must use the deduced nlevels to create a wildcard path for glob to find filepaths under
        dir_tuple = ["*"] * nlevels
        self.glob_path = os.path.join(self.le_rootdir.text(), *dir_tuple)

        # Get the basenames
        self.basenames.clear()
        for directory in iglob(self.glob_path):
            if os.path.isdir(directory):
                basename = os.path.basename(directory)
                if basename not in self.basenames:
                    self.basenames[basename] = {"left": "", "right": ""}
        if len(self.basenames) < 1:
            return

        delimiter, left_bases, right_bases = self.le_delimiter.text(), [], []
        try:
            if delimiter != "":
                for basename in self.basenames:
                    left, _, right = basename.partition(delimiter)
                    if left != "":
                        self.basenames[basename]["left"] = left
                        left_bases.append(left)
                    if right != "":
                        self.basenames[basename]["right"] = right
                        right_bases.append(right)
            else:
                for basename in self.basenames:
                    self.basenames[basename]["left"] = basename
                    left_bases.append(basename)
        except ValueError:
            QMessageBox().warning(self,
                                  "Impossible directory depth",
                                  "Could not parse through directories at that level"
                                  " Cancelling operation.",
                                  QMessageBox.Ok)
            del delimiter, dir_tuple, path, nlevels, from_root
            return

        left_bases = set(left_bases)
        right_bases = set(right_bases)
        pprint(self.basenames)
        print(f"Left bases: {left_bases}\nRight Bases: {right_bases}")

        # Exit early if there was no luck getting components
        if len(left_bases) == 0 and len(right_bases) == 0:
            return

        # Clear first
        self.clear_components_from_layout("left")
        self.clear_components_from_layout("right")

        # Fill left bases first
        if self.grp_leftside.isChecked():
            for left_component in sorted(left_bases):
                self.add_components_to_layout(component=left_component, side="left")

        # Then fill right bases
        if self.grp_rightside.isChecked():
            for right_component in sorted(right_bases):
                self.add_components_to_layout(component=right_component, side="right")

        if self.config["DeveloperMode"]:
            print(f'Left components: {self.leftside_name2le}')
            print(f'Right components: {self.rightside_name2le}')

        # These will be important for accounting for the user switching directories between determining delimiters and
        # the main function
        self.documented_sourcedir = self.le_rootdir.text()
        self.documented_exampledir = self.le_exampledir.text()

    def get_basename2dehybridized(self):
        """
        Gets a mapping of current basename to the expanded/dehybridized version of the filepath indicated
        @return: status, basename2dehybridized; status indicates whether this completed successfully and
        basename2dehybridized is the mapping
        """
        basename2dehybridized: dict = {}

        # Scenario 1: Both sides are specified
        if self.grp_leftside.isChecked() and self.grp_rightside.isChecked():
            # First, get the component2parent from each side
            leftside_component2parent = {component: self.leftside_name2le[component].text()
                                         for component in self.leftside_name2le.keys()}
            rightside_component2parent = {component: self.rightside_name2le[component].text()
                                          for component in self.rightside_name2le.keys()}
            # Then iterate over the basenames dict
            for basename in self.basenames:
                leftchild = self.basenames[basename]["left"]
                rightchild = self.basenames[basename]["right"]
                if all([leftchild in leftside_component2parent,
                        rightchild in rightside_component2parent
                        ]):
                    leftparent = leftside_component2parent[leftchild]
                    rightparent = rightside_component2parent[rightchild]
                    basename2dehybridized[basename] = os.path.join(leftparent, leftchild, rightparent, rightchild)
            return True, basename2dehybridized

        # Scenario 2: Left side only
        if self.grp_leftside.isChecked() and not self.grp_rightside.isChecked():
            # First, get the component2parent from the left side
            leftside_component2parent = {component: self.leftside_name2le[component].text()
                                         for component in self.leftside_name2le.keys()}
            # Then iterate over the basenames dict
            for basename in self.basenames:
                leftchild = self.basenames[basename]["left"]
                if leftchild in leftside_component2parent:
                    leftparent = leftside_component2parent[leftchild]
                    basename2dehybridized[basename] = os.path.join(leftparent, leftchild)
            return True, basename2dehybridized

        # Scenario 3: Right side only
        if not self.grp_leftside.isChecked() and self.grp_rightside.isChecked():
            # First, get the component2parent from the right side
            rightside_component2parent = {component: self.rightside_name2le[component].text()
                                          for component in self.rightside_name2le.keys()}
            # Then iterate over the basenames dict
            for basename in self.basenames:
                rightchild = self.basenames[basename]["right"]
                if rightchild in rightside_component2parent:
                    rightparent = rightside_component2parent[rightchild]
                    basename2dehybridized[basename] = os.path.join(rightparent, rightchild)
            return True, basename2dehybridized

        # Scenario 4: Neither are checked
        if not self.grp_leftside.isChecked() and not self.grp_rightside.isChecked():
            return False, basename2dehybridized

    def run_dehybridizer(self):
        """
        The main dehybridizer function
        """
        # The glob path is of the type raw/*/* etc. which defines how far down the tree it should located the folders
        # that need expanding
        if self.glob_path == "":
            return

        # Retrieve the mapping from the basename of a target path to the parent folder it should be placed under
        status, basename2dehybridized = self.get_basename2dehybridized()
        if not status:
            return

        # Standard warning message first to prevent accidental runs
        choice = QMessageBox().information(self,
                                           "This is a folder restructuring operation",
                                           f"The structure within the following source directory is about to be "
                                           f"altered:\n{self.le_rootdir.text()}\n\n"
                                           f"It is strongly recommended that you only proceed if you have made a "
                                           f"backup copy.\nProceed?",
                                           QMessageBox.Yes | QMessageBox.No)
        if choice == QMessageBox.No:
            return

        # Ascertain any directories that could not be processed when a delimiter is specified and warn the user with
        # a final chance to back out
        in_basename_but_not_dehybrid = set(self.basenames.keys()).difference(set(basename2dehybridized.keys()))
        if len(in_basename_but_not_dehybrid) > 0:
            names = '\n'.join(in_basename_but_not_dehybrid)
            choice = QMessageBox().warning(self,
                                           "Certain directories will be removed",
                                           f"The following folders will be removed in the process of "
                                           f"reorganizing folder structure:\n\n"
                                           f"{names}\n\nProceed?",
                                           QMessageBox.Yes | QMessageBox.No)
            if choice == QMessageBox.No:
                return

        # Disable the run button to present multi-inputs
        self.btn_run_dehybridizer.setEnabled(False)

        # Get the targets again and split them up; make as many threads as there are paths until a certain threshold
        # by which time the harddrive or SSD is fully saturated and no speedup will be gained through more threads
        paths = glob(self.glob_path)
        if len(paths) < 10:
            n_threads = len(paths)
        else:
            n_threads = 10

        # Store the divided paths and the translator under global variables for accessing by other functions
        self.paths = divide(n_threads, paths)
        self.translator = basename2dehybridized

        # Start the gif
        self.begin_movie()

        # Make the backup. The finished signal of the backup worker will indicate for the main function to begin
        if self.chk_makebackup.isChecked():
            source_dirname, source_basename = os.path.split(self.le_rootdir.text())
            backup_basename = source_basename + "_backup"
            backup_path = os.path.join(source_dirname, backup_basename)
            if os.path.exists(backup_path):
                pass
            else:
                copy_worker = Dehybridizer_CopyWorker(source_directory=self.le_rootdir.text(),
                                                      backup_directory=backup_path)
                copy_worker.signals.signal_done_copying.connect(self.expand_directories)
                self.threadpool.start(copy_worker)
                return

        # Otherwise, if no backup is specified, go straight into the folder expansion function
        else:
            self.signal_startexpanding.emit()

    @Slot()
    def expand_directories(self):
        """
        Slot; receives a signal EITHER from a Dehybridizer_CopyWorker that has completed the backup copying step OR
        from the main widget in the event that the backup step was omitted
        """
        workers = []
        for path_grp in self.paths:
            expander_worker = Dehybridizer_ExpandWorker(paths=path_grp, translator=self.translator)
            expander_worker.signals.signal_done_expanding.connect(self.re_enable_dehybridizer)
            self.global_debt -= 1
            workers.append(expander_worker)

        for worker in workers:
            self.threadpool.start(worker)

    @Slot()
    def re_enable_dehybridizer(self):
        """
        Slot; receives a signal from a Dehybridizer_ExpandWorker that has completed its job. Once all processing debt
        is paid off, the run button is re-enabled and the gif ends
        """
        self.global_debt += 1
        if self.global_debt == 0:
            self.btn_run_dehybridizer.setEnabled(True)
            self.end_movie()

    # Convenience function; adds a gif indicating processing below the textoutput and starts the gif
    def begin_movie(self):
        """
        Convenience function; inserts a gif indicating the start of folder copying and/or folder expansion
        """
        self.expanddirs_movieplayer.movie.start()
        self.expanddirs_movieplayer.setVisible(True)
        self.mainlay.insertWidget(2, self.expanddirs_movieplayer)

    # Convenience function; removes the gif below the textoutput and stops the gif
    def end_movie(self):
        """
        Convenience function; removes the gif below the mainlayout and stops it
        """
        self.expanddirs_movieplayer.movie.stop()
        self.expanddirs_movieplayer.setVisible(False)
        self.mainlay.removeWidget(self.expanddirs_movieplayer)


class Dehybridizer_CopyWorkerSignals(QObject):
    signal_done_copying = Signal()  # Signal sent by a copy worker to proceed to directory expansion


class Dehybridizer_CopyWorker(QRunnable):

    def __init__(self, source_directory, backup_directory):
        super().__init__()
        self.source_directory = source_directory
        self.backup_directory = backup_directory
        self.signals = Dehybridizer_CopyWorkerSignals()

    def run(self):
        try:
            shutil.copytree(src=self.source_directory,
                            dst=self.backup_directory,
                            dirs_exist_ok=True)
            self.signals.signal_done_copying.emit()
        except (OSError, RuntimeError) as copytree_error:
            QMessageBox().critical(self.parent(),
                                   "Failure to create backup directory",
                                   f"A critical error occurred while attempting to create a copy of:\n"
                                   f"{self.source_directory}\n\n"
                                   f"The following error was caught:\n{copytree_error}",
                                   QMessageBox.Ok)


class Dehybridizer_ExpandWorkerSignals(QObject):
    signal_done_expanding = Signal()  # Signal sent by a copy worker to reactivate the run button


class Dehybridizer_ExpandWorker(QRunnable):

    def __init__(self, paths, translator):
        super().__init__()
        self.src_paths = paths
        self.translator = translator
        self.signals = Dehybridizer_ExpandWorkerSignals()

    def run(self):
        for src in self.src_paths:
            # Skip files
            if not os.path.isdir(src):
                continue
            dirname, basename = os.path.split(src)
            try:
                dehybrid_basename = self.translator[basename]
            except KeyError:
                shutil.rmtree(path=src,
                              ignore_errors=True)
                continue
            dst = os.path.join(dirname, dehybrid_basename)
            print(f"Source: {src}\nDst: {dst}")

            # Do the copy and removal
            shutil.copytree(src=src,
                            dst=dst,
                            dirs_exist_ok=True)
            shutil.rmtree(path=src,
                          ignore_errors=True)

        self.signals.signal_done_expanding.emit()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = xASL_GUI_Dehybridizer()
    window.show()
    sys.exit(app.exec_())
