from PySide2.QtWidgets import *
from PySide2.QtGui import Qt, QFont, QIcon
from PySide2.QtCore import Signal, QSize, Slot
from src.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit, DandD_FileExplorer2ListWidget
from src.xASL_GUI_HelperFuncs_DirOps import *
import pandas as pd
from functools import partial
from pathlib import Path
from pprint import pprint
from json import load, dump


class xASL_GUI_MergeDirs(QWidget):
    """
    Class designated to merge processed study directories together
    """

    def __init__(self, parent=None, default_mergedir=""):
        super().__init__(parent=parent)
        self.parent = parent
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Explore ASL - Merge Study Directories")
        self.setMinimumSize(512, 512)
        self.mainlay = QVBoxLayout(self)

        # Misc Vars
        self.list_le_srcdirs: List[QLineEdit] = []
        self.list_hlay_srcdirs: List[QHBoxLayout] = []
        self.list_btn_srcdirs: List[QPushButton] = []

        # Group 1 - Overall settings
        self.grp_settings = QGroupBox(title="Merge Settings")
        self.formlay_settings = QFormLayout(self.grp_settings)
        (self.hlay_mergedst, self.le_mergedst,
         self.btn_mergedst) = self.make_le_btn_pair(le_placetxt="Specify the path to the merge directory",
                                                    btnfunc=self.select_dir, lefunc=self.can_merge)
        self.le_mergedst.setText(default_mergedir)
        self.spin_nsrcdirs = QSpinBox(minimum=2, singleStep=1)
        self.spin_nsrcdirs.valueChanged.connect(self.add_remove_srcdirs)
        self.chk_symlinks = QCheckBox(checked=True)
        self.chk_overwrite = QCheckBox(checked=False)
        self.formlay_settings.addRow(self.hlay_mergedst)
        for widget, desc, tipkey in zip([self.spin_nsrcdirs, self.chk_symlinks, self.chk_overwrite],
                                        ["Number of Studies to Merge", "Create Symlinks?", "Overwrite Existing?"],
                                        ["spin_nsrcdirs", "chk_symlinks", "chk_overwrite"]):
            self.formlay_settings.addRow(desc, widget)
            widget.setToolTip(self.parent.exec_tips["Modjob_RerunPrep"][tipkey])

        # Group 2 - Source Directories
        self.grp_srcdirs = QGroupBox(title="Directories to Merge")
        self.vlay_srcdirs_lvl1 = QVBoxLayout(self.grp_srcdirs)
        self.vlay_srcdirs_lvl1.setContentsMargins(0, 0, 0, 0)
        self.scroll_srcdirs = QScrollArea()
        self.cont_srcdirs = QWidget()
        self.vlay_srcdirs_lvl1.addWidget(self.scroll_srcdirs)
        self.scroll_srcdirs.setWidget(self.cont_srcdirs)
        self.scroll_srcdirs.setWidgetResizable(True)
        self.vlay_srcdirs_lvl2 = QVBoxLayout(self.cont_srcdirs)
        self.vlay_srcdirs_lvl2.addStretch(1)

        self.add_remove_srcdirs(2)

        # Finishing Widgets
        self.btn_mergedirs = QPushButton("Merge Directories", clicked=self.merge_dirs)
        self.btn_mergedirs.setEnabled(False)
        btn_font = QFont()
        btn_font.setPointSize(16)
        self.btn_mergedirs.setFont(btn_font)
        self.btn_mergedirs.setMinimumHeight(50)
        merge_icon = QIcon(str(Path(self.parent.config["ProjectDir"]) / "media" / "merge_ios_100ax100.png"))
        self.btn_mergedirs.setIcon(merge_icon)
        self.btn_mergedirs.setIconSize(QSize(50, 50))

        # Adding groups to main layout and adding tooltips
        self.mainlay.addWidget(self.grp_settings)
        self.mainlay.addWidget(self.grp_srcdirs)
        self.mainlay.addWidget(self.btn_mergedirs)
        for widget, tipkey in zip([self.le_mergedst, self.btn_mergedirs],
                                  ["le_mergedst", "btn_mergedirs"]):
            widget.setToolTip(self.parent.exec_tips["Modjob_RerunPrep"][tipkey])

    def add_remove_srcdirs(self, n_dirs: int):
        # Add widgets
        if n_dirs > len(self.list_le_srcdirs):
            for _ in range(abs(n_dirs - len(self.list_le_srcdirs))):
                hlay, le, btn = self.make_le_btn_pair("Specify the path to a study directory", self.select_dir,
                                                      self.can_merge)
                le.setToolTip(self.parent.exec_tips["Modjob_RerunPrep"]["le_studypath"])
                self.list_le_srcdirs.append(le)
                self.list_hlay_srcdirs.append(hlay)
                self.list_btn_srcdirs.append(btn)
                self.vlay_srcdirs_lvl2.insertLayout(len(self.list_le_srcdirs) - 1, hlay)
        # Remove wigets
        else:
            for _ in range(abs(n_dirs - len(self.list_le_srcdirs))):
                self.vlay_srcdirs_lvl2.removeItem(self.vlay_srcdirs_lvl2.itemAt(len(self.list_hlay_srcdirs) - 1))
                self.list_hlay_srcdirs[-1].setParent(None)
                self.list_le_srcdirs[-1].setParent(None)
                self.list_btn_srcdirs[-1].setParent(None)
                del self.list_hlay_srcdirs[-1]
                del self.list_le_srcdirs[-1]
                del self.list_btn_srcdirs[-1]

    @staticmethod
    def make_le_btn_pair(le_placetxt, btnfunc, lefunc, **kwargs):
        horizontal_layout = QHBoxLayout()
        lineedit = DandD_FileExplorer2LineEdit(acceptable_path_type="Directory",
                                               placeholderText=le_placetxt, clearButtonEnabled=True)
        lineedit.textChanged.connect(lefunc)
        button = QPushButton("...", clicked=partial(btnfunc, calling_le=lineedit, **kwargs))
        horizontal_layout.addWidget(lineedit)
        horizontal_layout.addWidget(button)
        return horizontal_layout, lineedit, button

    def select_dir(self, calling_le):
        directory = QFileDialog.getExistingDirectory(self, "Select Study Root",
                                                     self.parent.config["DefaultRootDir"],
                                                     QFileDialog.ShowDirsOnly)
        if directory == "":
            return
        calling_le.setText(str(Path(directory)))

    def can_merge(self):
        if self.le_mergedst.text() in {"", "/", "."}:
            self.btn_mergedirs.setEnabled(False)
            return

        merge_path = Path(self.le_mergedst.text())
        if not merge_path.is_dir():
            print("Exiting Early due to nonexistent merge destination")
            self.btn_mergedirs.setEnabled(False)
            return

        try:
            src_dirs = []
            for le in self.list_le_srcdirs:
                if le.text() in {"", "/", "."}:
                    continue
                path = Path(le.text())
                if "DataPar.json" not in {p.name for p in path.iterdir()}:
                    continue
                src_dirs.append(path)
        except FileNotFoundError:
            self.btn_mergedirs.setEnabled(False)
            return

        if len(src_dirs) <= 1:
            self.btn_mergedirs.setEnabled(False)
        else:
            self.btn_mergedirs.setEnabled(True)

    def merge_dirs(self):
        src_dirs = list({Path(le.text()) for le in self.list_le_srcdirs if le.text() not in {"", "/", "."}})
        dst_dir = Path(self.le_mergedst.text())
        merge_directories(roots=src_dirs, merge_root=dst_dir,
                          symbolic=self.chk_symlinks.isChecked(),
                          overwrite=self.chk_overwrite.isChecked())


class xASL_GUI_ModSidecars(QWidget):
    """
    Class designated to altering JSON sidecars
    """

    def __init__(self, parent=None):
        # Main Standard Setup
        super(xASL_GUI_ModSidecars, self).__init__(parent=parent)
        self.parent = parent
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Explore ASL - Modify JSON sidecars")
        self.setMinimumSize(400, 720)
        self.root_dir = Path(self.parent.le_modjob.text())
        self.mainlay = QVBoxLayout(self)

        # Grp 1 - Editing Jsons from a csv config file
        self.grp_fromfile = QGroupBox(title="Specify config from a CSV file", checkable=True, checked=True)
        self.grp_fromfile.clicked.connect(partial(self.ctrl_which_option, widget=self.grp_fromfile))
        self.grp_fromfile.clicked.connect(self.is_ready)
        self.formlay_fromfile = QFormLayout(self.grp_fromfile)
        self.hlay_fromfile = QHBoxLayout()
        self.le_fromfile = DandD_FileExplorer2LineEdit(acceptable_path_type="File",
                                                       supported_extensions=[".csv", ".tsv"])
        self.le_fromfile.setClearButtonEnabled(True)
        self.le_fromfile.textChanged.connect(self.is_ready)
        self.btn_fromfile = QPushButton("...", clicked=self.select_file)
        self.hlay_fromfile.addWidget(self.le_fromfile)
        self.hlay_fromfile.addWidget(self.btn_fromfile)
        self.formlay_fromfile.addRow(self.hlay_fromfile)

        # Grp 2 - Editing Jsons from a list of subjects and the indicated key + value
        self.grp_fromlist = QGroupBox(title="Specify subject list", checkable=True, checked=False)
        self.grp_fromlist.clicked.connect(partial(self.ctrl_which_option, widget=self.grp_fromlist))
        self.grp_fromlist.clicked.connect(self.is_ready)
        self.formlay_fromlist = QFormLayout(self.grp_fromlist)
        self.lab_subs = QLabel(text="Drag and drop the directories of the subjects\n"
                               "whose json sidecars should be altered")
        self.lst_subs = DandD_FileExplorer2ListWidget()
        self.lst_subs.itemsAdded.connect(self.le_fromfile.clear)
        self.lst_subs.itemsAdded.connect(self.is_ready)
        self.btn_clearsubjects = QPushButton("Clear the above list", clicked=self.lst_subs.clear)
        self.btn_clearsubjects.clicked.connect(self.is_ready)
        self.le_key = QLineEdit(placeholderText="Specify the name of the field to be changed", clearButtonEnabled=True)
        self.le_key.textChanged.connect(self.is_ready)
        self.le_value = QLineEdit(placeholderText="Specify the value, if applicable", clearButtonEnabled=True)
        for widget in [self.lab_subs, self.lst_subs, self.btn_clearsubjects, self.le_key, self.le_value]:
            self.formlay_fromlist.addRow(widget)

        # Grp 3 - Other Settings
        self.grp_runsettings = QGroupBox(title="Run Settings")
        self.formlay_runsettings = QFormLayout(self.grp_runsettings)
        self.cmb_actiontype = QComboBox()
        self.cmb_actiontype.addItems(["Add/Edit a field", "Remove a field"])
        self.cmb_actiontype.currentIndexChanged.connect(self.is_ready)
        self.chk_asl = QCheckBox(checked=True)
        self.chk_asl.stateChanged.connect(self.is_ready)
        self.chk_m0 = QCheckBox(checked=False)
        self.chk_m0.stateChanged.connect(self.is_ready)
        self.chk_t1 = QCheckBox(checked=False)
        self.chk_t1.stateChanged.connect(self.is_ready)
        self.btn_run = QPushButton("Alter JSON sidecars", clicked=self.alter_json_sidecars)
        self.btn_run.setEnabled(False)
        for widget, desc in zip([self.cmb_actiontype,
                                 self.chk_asl, self.chk_m0, self.chk_t1, self.btn_run],
                                ["Which action to perform",
                                 "Do this for ASL JSONs", "Do this for M0 JSONs", "Do this for T1 JSONs", ""]):
            self.formlay_runsettings.addRow(widget) if desc == "" else self.formlay_runsettings.addRow(desc, widget)

        # Put it all together
        for grp in [self.grp_fromfile, self.grp_fromlist, self.grp_runsettings]:
            self.mainlay.addWidget(grp)

        # Add tooltips
        for tipkey, tiptext in self.parent.exec_tips["Modjob_ModSidecars"].items():
            getattr(self, tipkey).setToolTip(tiptext)

    def ctrl_which_option(self, widget):
        if widget == self.grp_fromlist:
            self.grp_fromfile.setChecked(False)
        elif widget == self.grp_fromfile:
            self.grp_fromlist.setChecked(False)

    def select_file(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select the CSV Config File", self.parent.config["DefaultRootDir"],
                                              "Comma/Tab-Separated Values (*.csv | *.tsv)")
        if file == "":
            return
        self.le_fromfile.setText(str(Path(file)))

    def is_ready(self):
        if any([all([not self.chk_asl.isChecked(), not self.chk_t1.isChecked(), not self.chk_m0.isChecked()]),
                self.grp_fromfile.isChecked() and self.le_fromfile.text() in ["~", "/", ".", ""],
                self.grp_fromfile.isChecked() and not Path(self.le_fromfile.text()).exists(),
                self.grp_fromlist.isChecked() and self.lst_subs.count() == 0,
                self.grp_fromlist.isChecked() and self.le_key.text() == ""]):
            self.btn_run.setEnabled(False)
        else:
            self.btn_run.setEnabled(True)

    def alter_json_sidecars(self):
        QApplication.setOverrideCursor(Qt.WaitCursor)

        # Ascertain which types of scans to search for regardless of method
        which_scans = []
        for scan_type, widget in zip(["asl", "t1", "m0"], [self.chk_asl, self.chk_t1, self.chk_m0]):
            if widget.isChecked():
                which_scans.append(scan_type)
        del scan_type, widget

        # User wishes to alter sidecars using a pre-configured csv file
        if self.grp_fromfile.isChecked():
            for scan_type in which_scans:
                alter_sidecars(root_dir=self.root_dir, subjects=self.le_fromfile.text(), which_scan=scan_type,
                               action="remove" if self.cmb_actiontype.currentText() == "Remove a field" else "alter")
        # User wishes to alter sidecars using the drag & drop list
        else:
            # Get the list of subjects
            subs = []
            for idx in range(self.lst_subs.count()):
                sub_name = self.lst_subs.item(idx).text()
                if (self.root_dir / sub_name).exists():
                    subs.append(sub_name)
            if len(subs) == 0:
                QApplication.restoreOverrideCursor()
                QMessageBox.warning(self, self.parent.exec_errs["SubjectsNotFound"][0],
                                    self.parent.exec_errs["SubjectsNotFound"][1], QMessageBox.Ok)
                return

            for scan_type in which_scans:
                alter_sidecars(root_dir=self.root_dir, subjects=subs, which_scan=scan_type,
                               action="remove" if self.cmb_actiontype.currentText() == "Remove a field" else "alter",
                               key=self.le_key.text(), value=interpret_value(self.le_value.text()))
        # Afterwards...
        QMessageBox.information(self, "Finished json sidecar operation",
                                "Completed the requested json sidecar operation on the indicated subjects",
                                QMessageBox.Ok)
        QApplication.restoreOverrideCursor()
        return


class xASL_GUI_RerunPrep(QWidget):
    """
    Class designated to delete STATUS and other files such that a re-run of a particular module may be possible
    """

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.parent = parent
        self.root_dir = Path(self.parent.le_modjob.text())
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Explore ASL - Re-run setup")
        self.setMinimumSize(400, 720)
        self.mainlay = QVBoxLayout(self)
        self.directory_struct = dict()
        self.directory_struct["lock"] = self.get_path_directory_structure(self.root_dir / "lock")

        self.lock_tree = QTreeWidget(self)
        self.lock_tree.setToolTip(self.parent.exec_tips["Modjob_RerunPrep"]["lock_tree"])
        self.fill_tree(self.lock_tree.invisibleRootItem(), self.directory_struct)
        self.lock_tree.expandToDepth(2)
        self.lock_tree.itemChanged.connect(self.change_check_state)
        self.lock_tree.setHeaderLabel("Select the directories that should be redone")

        self.btn = QPushButton("Remove selected .status files", self, clicked=self.remove_status_files)
        self.btn.setMinimumHeight(50)
        font = QFont()
        font.setPointSize(14)
        self.btn.setFont(font)

        self.mainlay.addWidget(self.lock_tree)
        self.mainlay.addWidget(self.btn)

    def get_path_directory_structure(self, rootdir: Path):
        directory = {}
        for path in sorted(rootdir.iterdir()):
            if path.is_dir():
                directory[path.name] = self.get_path_directory_structure(path)
            else:
                directory[path.name] = None
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
        selected_status: List[QTreeWidgetItem] = [status for status in all_status if status.checkState(0)]
        filepaths: List[Path] = []

        for status in selected_status:
            filepath: List[str] = [status.text(0)]
            parent = status.parent()
            while parent.text(0) != "lock":
                filepath.insert(0, parent.text(0))
                parent = parent.parent()
            filepath.insert(0, "lock")
            filepaths.append(self.root_dir.joinpath(*filepath))

        return filepaths, selected_status

    def remove_status_files(self):
        filepaths, treewidgetitems = self.return_filepaths()
        if self.parent.config["DeveloperMode"]:
            print(f"REMOVING THE FOLLOWING STATUS FILES:")
            pprint(filepaths)

        for filepath in filepaths:
            filepath.unlink(missing_ok=True)

        # Clear the tree
        self.lock_tree.clear()
        # Refresh the file structure
        self.directory_struct.clear()
        self.directory_struct["lock"] = self.get_path_directory_structure(self.root_dir / "lock")
        # Refresh the tree
        self.fill_tree(self.lock_tree.invisibleRootItem(), self.directory_struct)
        self.lock_tree.expandToDepth(2)
        self.lock_tree.itemChanged.connect(self.change_check_state)

        QMessageBox.information(self.parent,
                                f"Re-run setup complete",
                                f"Successfully deleted the indicated .status files for the study:\n"
                                f"{self.root_dir}",
                                QMessageBox.Ok)


class xASL_GUI_TSValter(QWidget):
    """
    Class designated to alter the contents of participants.tsv
    """

    def __init__(self, parent=None):
        super().__init__(parent=parent)

        # Main Widget attributes
        self.parent = parent
        self.config = parent.config
        self.target_dir: str = self.parent.le_modjob.text()
        self.setWindowFlag(Qt.Window)
        self.setWindowTitle("Alter participants.tsv contents")
        self.resize(600, 600)

        # Misc attributes
        self.default_tsvcols = []
        self.default_metacols = []
        self.df_metadata = pd.DataFrame()
        self.df_parttsv = pd.DataFrame()
        self.df_backupmetadata = pd.DataFrame()
        self.df_backupparttsv = pd.DataFrame()
        self.metadatafile_duringalter: str = ""

        # Create the Widget
        self.Setup_UI_TSValter()

        # Followup
        self.btn_reset.setEnabled(False)
        self.btn_altertsv.setEnabled(False)
        if not (Path(self.target_dir) / "participants_orig.tsv").exists():
            self.btn_revert.setEnabled(False)
        self.load_parttsvfile()

    def Setup_UI_TSValter(self):
        self.mainlay = QVBoxLayout(self)
        self.hlay_listgrps = QHBoxLayout()  # Container holding the

        # Group for loading in metadata
        self.grp_metadatafile = QGroupBox(title="Load Metadata")
        self.formlay_metadatafile = QFormLayout(self.grp_metadatafile)
        self.hlay_metadatafile = QHBoxLayout()
        self.le_metadatafile = DandD_FileExplorer2LineEdit(acceptable_path_type="File",
                                                           supported_extensions=[".csv", ".tsv", ".xlsx"])
        self.btn_metadatafile = QPushButton("...", clicked=self.get_metadatafile)
        self.hlay_metadatafile.addWidget(self.le_metadatafile)
        self.hlay_metadatafile.addWidget(self.btn_metadatafile)
        self.btn_loadmetadata = QPushButton("Load Metadata", clicked=self.load_metadatafile)
        self.formlay_metadatafile.addRow("Metadata Filepath", self.hlay_metadatafile)
        self.formlay_metadatafile.addRow(self.btn_loadmetadata)

        # Group for holding the widgets associated with the metadata columns
        self.grp_metadatacols = QGroupBox(title="Metadata Columns")
        self.vlay_metadatacols = QVBoxLayout(self.grp_metadatacols)
        self.lst_metadata = QListWidget()
        self.lst_metadata.setDragEnabled(True)
        self.lst_metadata.setDefaultDropAction(Qt.MoveAction)
        self.le_metadata_subjectcol = QLineEdit(placeholderText="Indicate the metadata subject column",
                                                clearButtonEnabled=True)
        self.le_metadata_subjectcol.textChanged.connect(self.is_ready_merge)
        self.vlay_metadatacols.addWidget(self.lst_metadata)
        self.vlay_metadatacols.addWidget(self.le_metadata_subjectcol)

        # Group for holding the widgets associated with the participants.tsv columns
        self.grp_parttsvcols = QGroupBox(title="Participants.tsv Columns")
        self.vlay_parttsvcols = QVBoxLayout(self.grp_parttsvcols)
        self.lst_parttsv = ColnamesDragDrop_ListWidget()
        self.lst_parttsv.setAcceptDrops(True)
        self.lst_parttsv.signal_received_item.connect(self.set_reset_enabled)
        self.lst_parttsv.signal_received_item.connect(self.is_ready_merge)
        self.le_parttsv_subjectcol = QLineEdit(text="participant_id",
                                               placeholderText="Indicate the participants.tsv subject column",
                                               clearButtonEnabled=True)
        self.le_parttsv_subjectcol.textChanged.connect(self.is_ready_merge)
        self.vlay_parttsvcols.addWidget(self.lst_parttsv)
        self.vlay_parttsvcols.addWidget(self.le_parttsv_subjectcol)

        # Remainder for holding buttons to perform certain functions
        self.btn_reset = QPushButton("Reset", clicked=self.reset_colnames)
        self.btn_altertsv = QPushButton("Alter participants.tsv", clicked=self.alter_tsv)
        self.btn_revert = QPushButton("Emergency Revert", clicked=self.revert)

        # Add widgets to their respective layouts
        self.mainlay.addWidget(self.grp_metadatafile)
        self.hlay_listgrps.addWidget(self.grp_metadatacols)
        self.hlay_listgrps.addWidget(self.grp_parttsvcols)
        self.mainlay.addLayout(self.hlay_listgrps)
        self.mainlay.addWidget(self.btn_reset)
        self.mainlay.addWidget(self.btn_altertsv)
        self.mainlay.addWidget(self.btn_revert)

        # Add tooltips
        for tipkey, tiptext in self.parent.exec_tips["Modjob_TSValter"].items():
            getattr(self, tipkey).setToolTip(tiptext)

    def load_parttsvfile(self, from_reset: bool = False):
        if from_reset:
            path_parttsvfile = Path(self.target_dir, "participants_orig.tsv")
        else:
            path_parttsvfile = Path(self.target_dir, "participants.tsv")
        self.df_parttsv: pd.DataFrame = robust_read_csv(path_parttsvfile)
        self.lst_parttsv.clear()
        self.lst_parttsv.addItems(self.df_parttsv.columns)
        self.default_tsvcols = self.df_parttsv.columns.tolist()
        self.is_ready_merge()

    def get_metadatafile(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select the Metadata File",
                                              self.target_dir, "Data Files (*.tsv *.csv *.xlsx)")
        if path == "":
            return
        path = Path(path).resolve()
        self.le_metadatafile.setText(str(path))

    def load_metadatafile(self, from_reset: bool = False):
        if from_reset:
            meta_path = Path(self.metadatafile_duringalter).resolve()
        else:
            meta_path = Path(self.le_metadatafile.text()).resolve()
        if any([not meta_path.exists(), meta_path.is_dir()]):
            QMessageBox.warning(self, "Could not load metadata",
                                "You have indicated either a non-existent path or a directory", QMessageBox.Ok)
            return
        self.df_metadata: pd.DataFrame = robust_read_csv(meta_path)
        self.lst_metadata.clear()
        self.lst_metadata.addItems(self.df_metadata.columns)
        self.default_metacols = self.df_metadata.columns.tolist()
        self.is_ready_merge()

    def is_ready_merge(self):
        meta_path = Path(self.le_metadatafile.text()).resolve()

        if any([not meta_path.exists(), meta_path.is_dir(),
                self.le_metadata_subjectcol.text() not in self.df_metadata.columns,
                self.le_parttsv_subjectcol.text() not in self.df_parttsv.columns,
                len(self.get_current_parttsvcolnames() - set(self.default_tsvcols)) == 0]):
            self.btn_altertsv.setEnabled(False)
            return
        self.btn_altertsv.setEnabled(True)

    def get_current_metacolnames(self):
        return {self.lst_metadata.item(idx).text() for idx in range(self.lst_metadata.count())}

    def get_current_parttsvcolnames(self):
        return {self.lst_parttsv.item(idx).text() for idx in range(self.lst_parttsv.count())}

    @Slot()
    def set_reset_enabled(self):
        if not self.btn_reset.isEnabled():
            self.btn_reset.setEnabled(True)

    def reset_colnames(self):
        print(f"Resetting back to established default metadata columns and default participants.tsv columns")
        self.lst_parttsv.clear()
        self.lst_parttsv.addItems(self.default_tsvcols)
        self.lst_metadata.clear()
        self.lst_metadata.addItems(self.default_metacols)
        self.le_metadata_subjectcol.clear()
        self.btn_reset.setEnabled(False)  # Prevent False resetting

    def alter_tsv(self):

        # Last Sanity Check - Does the choice of SUBJECT columns make sense?
        subs_from_meta = set(self.df_metadata[self.le_metadata_subjectcol.text()].unique())
        subs_from_parttsv = set(self.df_parttsv[self.le_parttsv_subjectcol.text()].unique())
        if len(subs_from_meta.intersection(subs_from_parttsv)) == 0:
            QMessageBox.warning(self, "Imminent Bad Merge",
                                "No subjects from the indicated metadata column could be found in the participants.tsv "
                                "subject column. Cancelling update.\n\nPlease ensure you have listed the correct "
                                "subject columns or that the subject labelling is consistent between your metadata "
                                "and participants.tsv.", QMessageBox.Ok)
            return

        # Prepare certain variables in the event of a revert
        self.metadatafile_duringalter = self.le_metadatafile.text()
        self.df_backupmetadata: pd.DataFrame = self.df_metadata.copy()
        self.df_backupparttsv: pd.DataFrame = self.df_parttsv.copy()

        # Get the columns that are needed for merging, subset the metadata, then merge
        which_transfer_cols = self.get_current_parttsvcolnames() - set(self.default_tsvcols)
        meta_subset = self.df_metadata.loc[:, list(which_transfer_cols) + [self.le_metadata_subjectcol.text()]]
        merge_df = pd.merge(left=self.df_parttsv, right=meta_subset,
                            left_on=self.le_parttsv_subjectcol.text(),
                            right_on=self.le_metadata_subjectcol.text(),
                            how="left", suffixes=("_from_meta", "_from_participants"))
        merge_df.drop(self.le_metadata_subjectcol.text(), axis=1, inplace=True)

        # Save backup and new merge to their respective files
        self.df_backupparttsv.to_csv(Path(self.target_dir) / "participants_orig.tsv", sep="\t", index=False)
        merge_df.to_csv(Path(self.target_dir) / "participants.tsv", sep="\t", index=False)

        QMessageBox.information(self, "Successfully updated participants.tsv",
                                "participants.tsv has been updated with the indicated metadata columns. "
                                "A backup of the previous content has been saved as participants_orig.tsv",
                                QMessageBox.Ok)

        self.btn_altertsv.setEnabled(False)
        self.btn_revert.setEnabled(True)

    def revert(self):
        if not (Path(self.target_dir) / "participants_orig.tsv").exists():
            QMessageBox.warning(self, "No participants_orig.tsv file found",
                                "Could not perform a revert without the expected backup file", QMessageBox.Ok)
            return
        # Restore the old file
        (Path(self.target_dir) / "participants.tsv").unlink()

        # Load from the restored file
        self.load_parttsvfile(from_reset=True)
        self.load_metadatafile(from_reset=True)
        self.le_metadata_subjectcol.clear()
        (Path(self.target_dir) / "participants_orig.tsv").rename(target=Path(self.target_dir) / "participants.tsv")

        # Disabled buttons as necessary
        self.metadatafile_duringalter = ""
        self.btn_altertsv.setEnabled(False)
        self.btn_revert.setEnabled(False)
        self.btn_reset.setEnabled(False)
        self.is_ready_merge()

        QMessageBox.information(self, "participants.tsv has been reset",
                                "participants_orig.tsv was used to restore the previous iteration of participants.tsv",
                                QMessageBox.Ok)


class ColnamesDragDrop_ListWidget(QListWidget):
    """
    Class meant to drag and drop items between themselves
    """
    signal_received_item = Signal()  # Signal sent whenever the widget accepts

    def __init__(self, parent=None, **kwargs):
        super().__init__(parent, **kwargs)
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
            self.signal_received_item.emit()
        else:
            event.ignore()
