from PySide2.QtWidgets import *
from PySide2.QtGui import Qt, QFont
from PySide2.QtCore import Signal
from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit
import pandas as pd
import os
from functools import reduce
from pprint import pprint


class xASL_GUI_RerunPrep(QWidget):
    """
    Class designated to delete STATUS and other files such that a re-run of a particular module may be possible
    """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.parent = parent
        self.target_dir = self.parent.le_modjob.text()
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Explore ASL - Re-run setup")
        self.setMinimumSize(400, 720)
        self.mainlay = QVBoxLayout(self)
        self.dir_struct = self.get_directory_structure(os.path.join(self.target_dir, "lock"))

        self.lock_tree = QTreeWidget(self)
        self.lock_tree.setToolTip("Indicate which parts of the pipeline should be re-run for which \n"
                                  "modules/subjects/runs/etc.\n"
                                  "This window will delete all created .status files with the lock\n"
                                  "directory for the selected folders & files. When ExploreASL is\n"
                                  "re-run, it will detect these missing .status files and interpret\n"
                                  "that as a signal to re-run that particular section of the study's\n"
                                  "pipeline.")
        self.fill_tree(self.lock_tree.invisibleRootItem(), self.dir_struct)
        self.lock_tree.expandToDepth(2)
        self.lock_tree.itemChanged.connect(self.change_check_state)
        self.lock_tree.setHeaderLabel("Select the directories that should be redone")

        self.btn = QPushButton("Remove selected .status files \nto allow for re-run in the main widget", self,
                               clicked=self.remove_status_files)
        self.btn.setMinimumHeight(50)
        font = QFont()
        font.setPointSize(12)
        self.btn.setFont(font)

        self.mainlay.addWidget(self.lock_tree)
        self.mainlay.addWidget(self.btn)

    @staticmethod
    def get_directory_structure(rootdir):
        """
        Creates a nested dictionary that represents the folder structure of rootdir
        """
        directory = {}
        rootdir = rootdir.rstrip(os.sep)
        start = rootdir.rfind(os.sep) + 1
        for path, dirs, files in os.walk(rootdir):
            folders = path[start:].split(os.sep)
            subdir = dict.fromkeys(files)
            parent = reduce(dict.get, folders[:-1], directory)
            parent[folders[-1]] = subdir
        return directory

    def fill_tree(self, parent, d):
        if isinstance(d, dict):
            for key, value in d.items():
                it = QTreeWidgetItem()
                it.setText(0, key)
                if isinstance(value, dict):
                    parent.addChild(it)
                    it.setCheckState(0, Qt.Unchecked)
                    self.fill_tree(it, value)
                else:
                    parent.addChild(it)
                    it.setCheckState(0, Qt.Unchecked)

    @staticmethod
    def change_check_state(item: QTreeWidgetItem, col: int):
        if item.checkState(col):
            for idx in range(item.childCount()):
                item_child = item.child(idx)
                item_child.setCheckState(0, Qt.Checked)
        else:
            for idx in range(item.childCount()):
                item_child = item.child(idx)
                item_child.setCheckState(0, Qt.Unchecked)

    def return_filepaths(self):
        all_status = self.lock_tree.findItems(".status", (Qt.MatchEndsWith | Qt.MatchRecursive), 0)
        selected_status = [status for status in all_status if status.checkState(0)]
        filepaths = []

        for status in selected_status:
            filepath = [status.text(0)]
            parent = status.parent()
            while parent.text(0) != "lock":
                filepath.insert(0, parent.text(0))
                parent = parent.parent()
            filepath.insert(0, "lock")
            filepaths.append(os.path.join(self.target_dir, *filepath))

        return filepaths, selected_status

    def remove_status_files(self):
        filepaths, treewidgetitems = self.return_filepaths()
        if self.parent.config["DeveloperMode"]:
            print(f"REMOVING THE FOLLOWING STATUS FILES:")
            pprint(filepaths)

        for filepath, widget in zip(filepaths, treewidgetitems):
            os.remove(filepath)

        # Clear the tree
        self.lock_tree.clear()
        # Refresh the file structure
        self.dir_struct = self.get_directory_structure(os.path.join(self.target_dir, "lock"))
        # Refresh the tree
        self.fill_tree(self.lock_tree.invisibleRootItem(), self.dir_struct)
        self.lock_tree.expandToDepth(2)
        self.lock_tree.itemChanged.connect(self.change_check_state)

        QMessageBox.information(self.parent,
                                f"Re-run setup complete",
                                f"Successfully deleted the indicated .status files for the study:\n"
                                f"{self.target_dir}",
                                QMessageBox.Ok)


# noinspection PyAttributeOutsideInit
class xASL_GUI_TSValter(QWidget):
    """
    Class designated to alter the contents of participants.tsv
    """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.parent = parent
        self.config = parent.config
        self.target_dir = self.parent.le_modjob.text()
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Alter participants.tsv contents")
        self.mainlay = QVBoxLayout(self)

        # Misc variables
        self.default_tsvcols = []
        self.default_covcols = []

        # First section, the loading of the covariate file
        self.formlay_covfile = QFormLayout()
        self.hlay_covfile = QHBoxLayout()
        self.le_covfile = DandD_FileExplorer2LineEdit(acceptable_path_type="File",
                                                      supported_extensions=[".tsv", ".xlsx", ".csv"])
        self.le_covfile.textChanged.connect(self.load_covariates_data)
        self.btn_covfile = QPushButton("...", self)
        self.hlay_covfile.addWidget(self.le_covfile)
        self.hlay_covfile.addWidget(self.btn_covfile)
        self.formlay_covfile.addRow("Covariates File", self.hlay_covfile)

        # Second section, the Drag and Drop listviews
        self.hlay_dandcols = QHBoxLayout()
        self.vlay_covcols = QVBoxLayout()
        self.lab_covcols = QLabel(text="Covariates")
        self.list_covcols = QListWidget(self)
        self.list_covcols.setDragEnabled(True)  # This will have covariate names dragged from it
        self.list_covcols.setDefaultDropAction(Qt.MoveAction)
        self.vlay_covcols.addWidget(self.lab_covcols)
        self.vlay_covcols.addWidget(self.list_covcols)

        self.vlay_tsvcols = QVBoxLayout()
        self.lab_tsvcols = QLabel(text="TSV Variables")
        self.list_tsvcols = ColnamesDragDrop_ListWidget(self)
        self.list_tsvcols.setAcceptDrops(True)  # This will have covariate names dropped into it
        self.vlay_tsvcols.addWidget(self.lab_tsvcols)
        self.vlay_tsvcols.addWidget(self.list_tsvcols)

        self.hlay_dandcols.addLayout(self.vlay_covcols)
        self.hlay_dandcols.addLayout(self.vlay_tsvcols)

        # Third section, just the main buttons
        self.btn_altertsv = QPushButton("Commit changes", clicked=self.commit_changes)
        self.btn_reset = QPushButton("Reset to original participants.tsv", clicked=self.reset_changes)

        self.mainlay.addLayout(self.formlay_covfile)
        self.mainlay.addLayout(self.hlay_dandcols)
        self.mainlay.addWidget(self.btn_altertsv)
        self.mainlay.addWidget(self.btn_reset)

        # Once the widget UI is set up, load in the tsv data
        self.load_tsv_data()

    def reset_changes(self):
        # If changes have already been committed, remove the current participants_old.tsv and rename the _orig backup
        # back into participants.tsv, thereby restoring things to normal
        if os.path.exists(os.path.join(self.target_dir, "participants_orig.tsv")):
            os.remove(os.path.join(self.target_dir, "participants.tsv"))
            os.rename(src=os.path.join(self.target_dir, "participants_orig.tsv"),
                      dst=os.path.join(self.target_dir, "participants.tsv"))

        # Restore the lists
        self.list_covcols.clear()
        self.list_covcols.addItems(self.default_covcols)
        self.load_tsv_data()

    def load_tsv_data(self):
        tsv_file = os.path.join(self.target_dir, "participants.tsv")
        tsv_df = pd.read_csv(tsv_file, sep='\t')
        tsv_df.drop(labels=[col for col in tsv_df.columns if "Unnamed" in col], axis=1, inplace=True)
        self.tsv_df = tsv_df

        # Set the listwidget to contain the column names here
        tsv_colnames = self.tsv_df.columns.tolist()

        # Re-define the loaded data as the new defaults
        self.default_tsvcols = tsv_colnames

        # In case this is a re-load, clear the previous items in this list an try again
        self.list_tsvcols.clear()
        self.list_tsvcols.addItems(tsv_colnames)

    def load_covariates_data(self, cov_filepath: str):
        if cov_filepath == "":
            return

        # Prevent loading of the participants.tsv file itself as covariates
        if os.path.basename(cov_filepath) in ["participants.tsv", "participants_orig.tsv"]:
            self.le_covfile.clear()
            return

        # Note to self: no further QC checks needed because the le_covfile is a DandDFileExplorer2LineEdit;
        # this will already have the "File" and supported_extensions addressed

        # Load in the data
        if cov_filepath.endswith('.tsv'):
            df = pd.read_csv(cov_filepath, sep='\t')
        elif cov_filepath.endswith('.csv'):
            df = pd.read_csv(cov_filepath)
        else:
            df = pd.read_excel(cov_filepath)

        # Second set of checks: Dataframe must be the same length for proper concatenation
        if len(df) != len(self.tsv_df):
            QMessageBox().warning(self,
                                  "Incompatible Dataframes",
                                  f"The number of entires of the covariates dataframe: {len(df)}\n"
                                  f"did not equal the number of entries in participants.tsv: {len(self.tsv_df)}",
                                  QMessageBox.Ok)
            self.le_covfile.clear()
            return

        # Drop any columns already found in the participants.tsv df
        for col in df.columns:
            if col in self.tsv_df.columns.tolist():
                df.drop(labels=col, axis=1, inplace=True)

        self.cov_df = df

        # Set the listwidget to contain the column names here
        cov_colnames = self.cov_df.columns.tolist()
        self.list_covcols.clear()
        self.list_covcols.addItems(cov_colnames)
        self.default_covcols = cov_colnames

        # Also reset the tsv columns
        self.list_tsvcols.clear()
        print(f"default tsv cols: {self.default_tsvcols}")
        self.list_tsvcols.addItems(self.default_tsvcols)

    def commit_changes(self):
        commit_df = pd.DataFrame()

        # Iterate over the labels in the tsv listwidget.
        for idx in range(self.list_tsvcols.count()):
            colname = self.list_tsvcols.item(idx).text()

            # Check which dataframe it is in, then add that dataframe's column contents to the commit_df
            if colname not in self.tsv_df.columns:
                commit_df[colname] = self.cov_df[colname].values
            else:
                commit_df[colname] = self.tsv_df[colname].values

        # If the _orig backup hasn't been made from the original, make it at this point in time.
        if not os.path.exists(os.path.join(self.target_dir, "participants_orig.tsv")):
            try:
                os.replace(src=os.path.join(self.target_dir, "participants.tsv"),
                           dst=os.path.join(self.target_dir, "participants_orig.tsv"))
            except OSError:
                QMessageBox().warning(self,
                                      "Could not create a backup participants_orig.tsv from the original",
                                      "It is likely that participants.tsv is opened in another program. "
                                      "Please close that program and try committing again.",
                                      QMessageBox.Ok)
                return

        # Overwrite participants.tsv
        try:
            commit_df.to_csv(path_or_buf=os.path.join(self.target_dir, "participants.tsv"),
                             sep='\t',
                             na_rep='n/a',
                             index=False)
        except OSError:
            QMessageBox().warning(self,
                                  "Could not commit changes to participants.tsv",
                                  "It is likely that the file is opened in another program. "
                                  "Please close that program and try committing again.",
                                  QMessageBox.Ok)


class ColnamesDragDrop_ListWidget(QListWidget):
    """
    Class meant to drag and drop items between themselves
    """
    alert_regex = Signal()  # This will be called in the ParmsMaker superclass to update its regex

    def __init__(self, parent=None):
        super().__init__(parent)
        # self.setDragDropMode(QAbstractItemView.DragDrop)
        # self.setDefaultDropAction(Qt.MoveAction)  # this was the magic line
        # self.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.colnames = []

    def dragEnterEvent(self, event) -> None:
        if event.mimeData().hasFormat('application/x-qabstractitemmodeldatalist'):
            print(event.mimeData().text())
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event) -> None:
        if event.mimeData().hasFormat('application/x-qabstractitemmodeldatalist'):
            event.setDropAction(Qt.MoveAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event) -> None:
        if event.mimeData().hasFormat('application/x-qabstractitemmodeldatalist'):
            event.accept()
            for item in event.source().selectedItems():
                self.addItem(item.text())
        else:
            event.ignore()
