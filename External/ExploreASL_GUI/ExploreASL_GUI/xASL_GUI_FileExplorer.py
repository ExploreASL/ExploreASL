from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
from ExploreASL_GUI.xASL_GUI_HelperFuncs import set_widget_icon
from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit
from ExploreASL_GUI.xASL_GUI_HelperFuncs_StringOps import set_os_dependent_text
import os


class xASL_FileExplorer(QWidget):
    def __init__(self, parent):
        super().__init__(parent=parent)
        self.config = parent.config
        self.path_history = []
        self.path_index = 0

        # Define the lineedit of the current directory
        self.le_current_dir = DandD_FileExplorer2LineEdit(acceptable_path_type="Directory")
        self.le_current_dir.editingFinished.connect(self.go_from_text)

        # Define the buttons that will be used to navigate through the directories
        self.hlay_btns = QHBoxLayout()
        self.btn_back = QPushButton(clicked=self.go_back)
        self.btn_up = QPushButton(clicked=self.go_up)
        self.btn_forward = QPushButton(clicked=self.go_forward)
        for btn, icon_name in zip([self.btn_back,
                                   self.btn_up,
                                   self.btn_forward],
                                  ["arrow_left_encircled.svg",
                                   "arrow_up_encircled.svg",
                                   "arrow_right_encircled.svg"]):
            set_widget_icon(widget=btn, config=self.config, icon_name=icon_name, size=(25, 25))
            btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
            self.hlay_btns.addWidget(btn)

        # Define the file system model and its display container
        self.treev_file = QTreeView()
        self.treev_file.setContextMenuPolicy(Qt.CustomContextMenu)
        self.treev_file.customContextMenuRequested.connect(self.menuContextTree)
        self.model_file = QFileSystemModel()

        self.model_file.setRootPath(self.model_file.myComputer(Qt.DisplayRole))
        self.treev_file.setModel(self.model_file)
        self.treev_file.setRootIndex(self.model_file.index(self.config["DefaultRootDir"]))
        self.treev_file.header().resizeSection(0, 250)
        self.treev_file.setDragEnabled(True)
        self.treev_file.setSortingEnabled(True)
        self.treev_file.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.treev_file.sortByColumn(1, Qt.AscendingOrder)
        self.treev_file.setExpandsOnDoubleClick(False)
        self.treev_file.setAnimated(True)
        self.treev_file.doubleClicked.connect(self.go_down)
        self.treev_file.setMinimumWidth(500)
        self.treev_file.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.path_history.append(self.config["DefaultRootDir"].replace("\\", "/"))
        set_os_dependent_text(linedit=self.le_current_dir,
                              config_ossystem=self.config["Platform"],
                              text_to_set=self.config["DefaultRootDir"])

        # With the model defined, define the auto-completer class
        self.completer_current_dir = QCompleter(completionMode=QCompleter.InlineCompletion)
        self.completer_current_dir.setModel(self.model_file)
        self.le_current_dir.setCompleter(self.completer_current_dir)

        # Define main layout and add components to it
        self.mainlay = QVBoxLayout(self)
        self.mainlay.addWidget(self.le_current_dir)
        self.mainlay.addLayout(self.hlay_btns)
        self.mainlay.addWidget(self.treev_file)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

    def menuContextTree(self, point):
        # Infos about the node selected.
        index = self.treev_file.indexAt(point)

        if not index.isValid():
            return

        # We build the menu.
        menu = QMenu()
        action = menu.addAction("To be implemented to a future update")
        menu.exec_(self.treev_file.mapToGlobal(point))

    def path_change(self, newpath, current_index):
        try:
            # If the newpath is not the same as the path ahead (i.e if returning forward or up or down), clear the
            # path history ahead
            if self.path_history[current_index + 1] != newpath:
                for idx in range(current_index, len(self.path_history)):
                    del self.path_history[current_index + 1]
                self.treev_file.setRootIndex(self.model_file.index(newpath.replace('\\', '/')))
                set_os_dependent_text(linedit=self.le_current_dir,
                                      config_ossystem=self.config["Platform"],
                                      text_to_set=newpath)
                self.model_file.setRootPath(newpath)
            # Otherwise, proceed as normal
            else:
                self.treev_file.setRootIndex(self.model_file.index(newpath))
                set_os_dependent_text(linedit=self.le_current_dir,
                                      config_ossystem=self.config["Platform"],
                                      text_to_set=newpath)
                self.model_file.setRootPath(newpath)
        # If an index error was encountered, we must be at the head of the path history and there is no need to worry
        # about looking ahead
        except IndexError:
            self.treev_file.setRootIndex(self.model_file.index(newpath.replace('\\', '/')))
            set_os_dependent_text(linedit=self.le_current_dir,
                                  config_ossystem=self.config["Platform"],
                                  text_to_set=newpath)
            self.model_file.setRootPath(newpath)

    # Enter a path in the lineedit and press Enter
    def go_from_text(self):
        newpath = self.le_current_dir.text().replace('\\', '/')

        # Avoid writing into history if the Enter event is the current directory
        if newpath == self.path_history[self.path_index]:
            return

        if os.path.exists(newpath):
            if os.path.isdir(newpath):
                self.path_change(newpath=newpath, current_index=self.path_index)
                try:
                    if self.path_history[self.path_index + 1] == newpath:
                        self.path_index += 1
                    else:
                        print("go_from_text; this should never print")
                except IndexError:
                    self.path_history.append(newpath.replace('\\', '/'))
                    self.path_index += 1

                self.dev_path_print("Pressed Enter into a new directory")
            else:
                QMessageBox().warning(self,
                                      "Cannot Enter into a file",
                                      f"The File Explorer detected that the path you entered:\n{newpath}\nis a file. "
                                      f"This program is not able to open files from this location at the current time")
                # Reset the lineedit back to the text prior to the Enter event
                set_os_dependent_text(linedit=self.le_current_dir,
                                      config_ossystem=self.config["Platform"],
                                      text_to_set=self.path_history[self.path_index])
                return
        else:
            QMessageBox().warning(self,
                                  "Path does not exist",
                                  f"The File Explorer could not find the path:\n{newpath}",
                                  QMessageBox.Ok)
            # Reset the lineedit back to the text prior to the Enter event
            set_os_dependent_text(linedit=self.le_current_dir,
                                  config_ossystem=self.config["Platform"],
                                  text_to_set=self.path_history[self.path_index])
            return

    # Go up a directory
    def go_up(self):
        current_root = self.model_file.filePath(self.treev_file.rootIndex())
        dirname = os.path.dirname(current_root)
        if os.path.exists(dirname):
            if os.path.isdir(dirname) and any([os.path.dirname(dirname) != dirname, dirname != self.path_history[-1]]):
                self.path_change(newpath=dirname, current_index=self.path_index)
                try:
                    if self.path_history[self.path_index + 1] == dirname:
                        self.path_index += 1
                    else:
                        print("go_up; this should never print")
                except IndexError:
                    self.path_history.append(dirname.replace('\\', '/'))
                    self.path_index += 1

        self.dev_path_print("Pressed go_up")

    def go_down(self, filepath_modelindex):
        filepath = self.model_file.filePath(filepath_modelindex)
        if os.path.isdir(filepath):
            self.path_change(newpath=filepath, current_index=self.path_index)
            try:
                if self.path_history[self.path_index + 1] == filepath:
                    self.path_index += 1
                else:
                    print("go_down; this should never print")
            except IndexError:
                self.path_history.append(filepath.replace('\\', '/'))
                self.path_index += 1

            self.dev_path_print("Double-clicked to go down into a directory")
        else:
            print(f"The filepath: {filepath} is not a directory that can be entered into")

    def go_back(self):
        if self.path_index != 0:  # cannot be at the beginning of the history
            previous_path = self.path_history[self.path_index - 1]
            if os.path.exists(previous_path):
                if os.path.isdir(previous_path):
                    self.treev_file.setRootIndex(self.model_file.index(previous_path))
                    self.path_index -= 1
                    set_os_dependent_text(linedit=self.le_current_dir,
                                          config_ossystem=self.config["Platform"],
                                          text_to_set=previous_path)
                    self.model_file.setRootPath(previous_path)

        self.dev_path_print("Pressed go_back")

    def go_forward(self):
        if self.path_index != (len(self.path_history) - 1):
            forward_path = self.path_history[self.path_index + 1]
            if os.path.exists(forward_path):
                if os.path.isdir(forward_path):
                    self.treev_file.setRootIndex(self.model_file.index(forward_path))
                    self.path_index += 1
                    set_os_dependent_text(linedit=self.le_current_dir,
                                          config_ossystem=self.config["Platform"],
                                          text_to_set=forward_path)
                    self.model_file.setRootPath(forward_path)

        self.dev_path_print("Pressed go_forward")

    ############################################################################
    # DEVELOPER FUNCTIONS FOR WHEN RUNNING IN DEVELOPER MODE FOR TROUBLESHOOTING
    ############################################################################

    def dev_path_print(self, func_descrption=''):
        # Dev troubleshooting
        if self.config["DeveloperMode"]:
            print(func_descrption)
            print(f"self.path_history={self.path_history}")
            print(f"self.path_index={self.path_index}")
            print(f"current directory post-change: {self.model_file.filePath(self.treev_file.rootIndex())}")
            print("----------------------------\n")
