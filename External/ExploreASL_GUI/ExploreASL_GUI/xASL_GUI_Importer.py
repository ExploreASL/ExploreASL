from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit
from ExploreASL_GUI.xASL_GUI_DCM2BIDS import get_dicom_directories, asldcm2bids_onedir, create_import_summary, \
    bids_m0_followup
from ExploreASL_GUI.xASL_GUI_HelperFuncs_StringOps import set_os_dependent_text
from glob import iglob
from tdda import rexpy
from pprint import pprint
from collections import OrderedDict
from more_itertools import divide
import json
import os
import platform


class Importer_WorkerSignals(QObject):
    """
    Class for handling the signals sent by an ExploreASL worker
    """
    signal_send_summaries = Signal(list)  # Signal sent by worker to process the summaries of imported files
    signal_send_errors = Signal(list)  # Signal sent by worker to indicate the file where something has failed


class Importer_Worker(QRunnable):
    """
    Worker thread for running the import for a particular group.
    """

    def __init__(self, dcm_dirs, config, use_legacy_mode):
        self.dcm_dirs = dcm_dirs
        self.import_config = config
        self.use_legacy_mode = use_legacy_mode
        super().__init__()
        self.signals = Importer_WorkerSignals()
        self.import_summaries = []
        self.failed_runs = []
        print("Initialized Worker with args:\n")
        pprint(self.import_config)

    # This is called by the threadpool during threadpool.start(worker)
    def run(self):
        for dicom_dir in self.dcm_dirs:
            result, last_job, import_summary = asldcm2bids_onedir(dcm_dir=dicom_dir,
                                                                  config=self.import_config,
                                                                  legacy_mode=self.use_legacy_mode)
            if result:
                self.import_summaries.append(import_summary)
            else:
                self.failed_runs.append((dicom_dir, last_job))

        self.signals.signal_send_summaries.emit(self.import_summaries)
        if len(self.failed_runs) > 0:
            self.signals.signal_send_errors.emit(self.failed_runs)


# noinspection PyCallingNonCallable
class xASL_GUI_Importer(QMainWindow):
    def __init__(self, parent_win=None):
        # Parent window is fed into the constructor to allow for communication with parent window devices
        super().__init__(parent=parent_win)
        if parent_win is not None:
            self.config = self.parent().config
        else:
            with open("ExploreASL_GUI_masterconfig.json") as config_reader:
                self.config = json.load(config_reader)

        # Misc and Default Attributes
        self.labfont = QFont()
        self.labfont.setPointSize(16)
        self.rawdir = ''
        self.subject_regex = None
        self.visit_regex = None
        self.scan_regex = None
        self.visit_aliases = OrderedDict()
        self.scan_aliases = dict.fromkeys(["ASL4D", "T1", "M0", "FLAIR"])
        self.cmb_visitaliases_dict = {}
        self.threadpool = QThreadPool()
        self.import_summaries = []
        self.failed_runs = []

        # Window Size and initial visual setup
        self.setMinimumSize(540, 720)
        self.resize(600, 720)
        self.setWindowTitle("ExploreASL ASL2BIDS Importer")
        self.cw = QWidget(self)
        self.setCentralWidget(self.cw)
        self.mainlay = QVBoxLayout(self.cw)
        self.mainsplit = QSplitter(Qt.Vertical, self.cw)
        handle_path = os.path.join(self.config["ProjectDir"], "media", "3_dots_horizontal.svg").replace('\\', '/')
        handle_style = 'QSplitter::handle {image: url(' + handle_path + ');}'
        self.mainsplit.setStyleSheet(handle_style)
        self.mainsplit.setHandleWidth(20)
        self.cw.setLayout(self.mainlay)
        self.mainlay.addWidget(self.mainsplit)

        self.Setup_UI_UserSpecifyDirStuct()
        self.Setup_UI_UserSpecifyScanAliases()
        self.Setup_UI_UserSpecifyVisitAliases()

        self.btn_run_importer = QPushButton("Run ASL2BIDS", self.cw, clicked=self.run_importer)
        self.btn_run_importer.setFont(self.labfont)
        self.btn_run_importer.setFixedHeight(50)
        self.btn_run_importer.setEnabled(False)
        self.mainsplit.addWidget(self.btn_run_importer)

        self.mainsplit.setSizes([150, 225, 300, 50])

    def Setup_UI_UserSpecifyDirStuct(self):
        self.grp_dirstruct = QGroupBox(title="Specify Directory Structure")
        self.grp_dirstruct.setMaximumHeight(225)
        self.vlay_dirstruct = QVBoxLayout(self.grp_dirstruct)

        # First specify the root directory
        self.formlay_rootdir = QFormLayout()
        self.hlay_rootdir = QHBoxLayout()
        self.le_rootdir = DandD_FileExplorer2LineEdit(acceptable_path_type="Directory")
        self.le_rootdir.setPlaceholderText("Drag and drop your study's raw directory here")
        self.le_rootdir.setReadOnly(True)
        self.le_rootdir.textChanged.connect(self.set_rootdir_variable)
        self.le_rootdir.textChanged.connect(self.clear_widgets)
        self.le_rootdir.textChanged.connect(self.is_ready_import)
        self.btn_setrootdir = QPushButton("...", self.grp_dirstruct, clicked=self.set_import_root_directory)
        self.hlay_rootdir.addWidget(self.le_rootdir)
        self.hlay_rootdir.addWidget(self.btn_setrootdir)
        self.chk_uselegacy = QCheckBox(checked=True)
        self.formlay_rootdir.addRow("Raw Root Directory", self.hlay_rootdir)
        self.formlay_rootdir.addRow("Use Legacy Import", self.chk_uselegacy)

        # Next specify the QLabels that can be dragged to have their text copied elsewhere
        self.hlay_placeholders = QHBoxLayout()
        self.lab_holdersub = DraggableLabel("Subject", self.grp_dirstruct)
        self.lab_holdersess = DraggableLabel("Visit", self.grp_dirstruct)
        self.lab_holderscan = DraggableLabel("Scan", self.grp_dirstruct)
        self.lab_holderdummy = DraggableLabel("Dummy", self.grp_dirstruct)
        for lab in [self.lab_holdersub, self.lab_holdersess, self.lab_holderscan, self.lab_holderdummy]:
            self.hlay_placeholders.addWidget(lab)

        # Next specify the QLineEdits that will be receiving the dragged text
        self.hlay_receivers = QHBoxLayout()
        self.lab_rootlabel = QLabel(text="raw")
        self.lab_rootlabel.setFont(self.labfont)
        self.levels = {}
        for idx, (level, func) in enumerate(zip(["Level1", "Level2", "Level3", "Level4", "Level5"],
                                                [self.get_nth_level_dirs] * 5)):
            le = DandD_Label2LineEdit(self, self.grp_dirstruct, idx)
            le.modified_text.connect(self.get_nth_level_dirs)
            le.textChanged.connect(self.update_sibling_awareness)
            le.textChanged.connect(self.is_ready_import)
            self.levels[level] = le

        self.hlay_receivers.addWidget(self.lab_rootlabel)
        if platform.system() == "Windows":
            separator = '\\'
        else:
            separator = '/'
        lab_sep = QLabel(text=separator)
        lab_sep.setFont(self.labfont)
        self.hlay_receivers.addWidget(lab_sep)
        for ii, level in enumerate(self.levels.values()):
            level.setFont(self.labfont)
            self.hlay_receivers.addWidget(level)
            if ii < 4:
                lab_sep = QLabel(text=separator)
                lab_sep.setFont(self.labfont)
                self.hlay_receivers.addWidget(lab_sep)

        # Include the button that will clear the current structure for convenience
        self.btn_clear_receivers = QPushButton("Clear the fields", self.grp_dirstruct, clicked=self.clear_receivers)

        # Organize layouts
        self.vlay_dirstruct.addLayout(self.formlay_rootdir, 1)
        self.vlay_dirstruct.addLayout(self.hlay_placeholders, 2)
        self.vlay_dirstruct.addLayout(self.hlay_receivers, 2)
        self.vlay_dirstruct.addWidget(self.btn_clear_receivers, 2)

        self.mainsplit.addWidget(self.grp_dirstruct)

    def Setup_UI_UserSpecifyScanAliases(self):
        # Next specify the scan aliases
        self.grp_scanaliases = QGroupBox(title="Specify Scan Aliases")
        self.cmb_scanaliases_dict = dict.fromkeys(["ASL4D", "T1", "M0", "FLAIR"])
        self.formlay_scanaliases = QFormLayout(self.grp_scanaliases)
        for description, scantype in zip(["ASL scan alias:\n(Mandatory)",
                                          "T1 scan alias:\n(Mandatory)",
                                          "M0 scan alias:\n(Optional)",
                                          "FLAIR scan alias:\n(Optional)"],
                                         self.cmb_scanaliases_dict.keys()):
            cmb = QComboBox(self.grp_scanaliases)
            cmb.addItems(["Select an alias"])
            cmb.currentTextChanged.connect(self.update_scan_aliases)
            cmb.currentTextChanged.connect(self.is_ready_import)
            self.cmb_scanaliases_dict[scantype] = cmb
            self.formlay_scanaliases.addRow(description, cmb)

        self.mainsplit.addWidget(self.grp_scanaliases)

    def Setup_UI_UserSpecifyVisitAliases(self):
        # Define the groupbox and its main layout
        self.grp_visitaliases = QGroupBox(title="Specify Visit Aliases and Ordering")
        self.vlay_visitaliases = QVBoxLayout(self.grp_visitaliases)
        self.vlay_visitaliases.setContentsMargins(0, 0, 0, 0)
        self.scroll_visitaliases = QScrollArea(self.grp_visitaliases)
        self.cont_visitaliases = QWidget()
        self.scroll_visitaliases.setWidget(self.cont_visitaliases)
        self.scroll_visitaliases.setWidgetResizable(True)

        # Arrange widgets and layouts
        self.le_visitaliases_dict = dict()
        self.formlay_visitaliases = QFormLayout(self.cont_visitaliases)
        self.vlay_visitaliases.addWidget(self.scroll_visitaliases)
        self.mainsplit.addWidget(self.grp_visitaliases)

    # Purpose of this function is to set the directory of the root path lineedit based on the adjacent pushbutton
    @Slot()
    def set_import_root_directory(self):
        dir_path = QFileDialog.getExistingDirectory(QFileDialog(),
                                                    "Select the raw directory of your study",
                                                    self.parent().config["DefaultRootDir"],
                                                    QFileDialog.ShowDirsOnly)
        if os.path.exists(dir_path):
            set_os_dependent_text(linedit=self.le_rootdir,
                                  config_ossystem=self.parent().config["Platform"],
                                  text_to_set=dir_path)

    # Purpose of this function is to change the value of the rawdir attribute based on the current text
    @Slot()
    def set_rootdir_variable(self, path):
        if path == '' or not os.path.exists(path):
            return

        if all([os.path.isdir(path),
                os.path.basename(path) == "raw"]
               ):
            self.rawdir = self.le_rootdir.text()

    def get_nth_level_dirs(self, dir_type: str, level: int):
        """
        :param dir_type: whether this is a subject, visit, or scan
        :param level: which lineedit, in python index terms, emitted this signal
        """
        # Requirements to proceed
        if any([self.rawdir == '',  # Raw dir must be specified
                not os.path.exists(self.rawdir),  # Raw dir must exist
                os.path.basename(self.rawdir) != 'raw',  # Raw dir's basename must be raw
                ]):
            return

        # Check if a reset is needed
        self.check_if_reset_needed()

        # If this was a clearing, the dir_type will be an empty string and the function should exit after any resetting
        # has been performed
        if dir_type == '':
            return

        # Get the directories at the depth according to which lineedit's text was changed
        dir_tuple = ["*"] * (level + 1)
        path = os.path.join(self.rawdir, *dir_tuple)
        try:
            directories, basenames = zip(*[(directory, os.path.basename(directory)) for directory in iglob(path)
                                           if os.path.isdir(directory)])
        except ValueError:
            QMessageBox().warning(self,
                                  "Impossible directory depth",
                                  "The directory depth you've indicated does not have "
                                  "directories present at that level."
                                  " Cancelling operation.",
                                  QMessageBox.Ok)
            # Clear the appropriate lineedit that called this function after the error message
            list(self.levels.values())[level].clear()
            return

        # Do not proceed if no directories were found and clear the linedit that emitted the textChanged signal
        if len(directories) == 0:
            idx = list(self.levels.keys())[level]
            print(f"idx: {idx}")
            self.levels[idx].clear()
            return

        # Otherwise, make the appropriate adjustment depending on which label was dropped in
        if dir_type == "Subject":
            self.subject_regex = self.infer_regex(list(basenames))
            print(f"Subject regex: {self.subject_regex}")
            del directories, basenames

        elif dir_type == "Visit":
            self.visit_regex = self.infer_regex(list(set(basenames)))
            print(f"Visit regex: {self.visit_regex}")
            self.reset_visit_aliases(basenames=list(set(basenames)))
            del directories, basenames

        elif dir_type == "Scan":
            self.scan_regex = self.infer_regex(list(set(basenames)))
            print(f"Scan regex: {self.scan_regex}")
            self.reset_scan_alias_cmbs(basenames=list(set(basenames)))
            del directories, basenames

        elif dir_type == "Dummy":
            del directories, basenames
            return

        else:
            del directories, basenames
            print("Error. This should never print")
            return

    #####################################
    # SECTION - RESET AND CLEAR FUNCTIONS
    #####################################

    def clear_widgets(self):
        """
        Raw reset. Resets all important variables upon a change in the indicated raw directory text.
        """
        # Resets everything back to normal
        self.subject_regex = None
        self.visit_regex = None
        self.scan_regex = None
        self.clear_receivers()
        self.clear_visit_alias_cmbs_and_les()
        self.reset_scan_alias_cmbs(basenames=[])
        self.visit_aliases = OrderedDict()
        self.scan_aliases = dict.fromkeys(["ASL4D", "T1", "M0", "FLAIR"])

        if self.config["DeveloperMode"]:
            print("clear_widgets engaged due to a change in the indicated Raw directory")

    def check_if_reset_needed(self):
        """
        More specialized reset function. If any of the drop-enabled lineedits has their field change,
        this function will accomodate that change by resetting the variable that may have been removed
        during the drop
        """
        used_directories = [le.text() for le in self.levels.values()]
        # If subjects is not in the currently-specified structure and the regex has been already set
        if "Subject" not in used_directories and self.subject_regex is not None:
            self.subject_regex = None

        # If visits is not in the currently-specified structure and the regex has been already set
        if "Visit" not in used_directories and self.visit_regex is not None:
            self.visit_regex = None
            self.visit_aliases.clear()
            self.clear_visit_alias_cmbs_and_les()  # This clears the visitaliases dict and the widgets

        if "Scan" not in used_directories and self.scan_regex is not None:
            self.scan_regex = None
            self.scan_aliases = dict.fromkeys(["ASL4D", "T1", "M0", "FLAIR"])
            self.reset_scan_alias_cmbs(basenames=[])

    def clear_receivers(self):
        """
        Convenience function for resetting the drop-enabled lineedits
        """
        for le in self.levels.values():
            le.clear()

    def reset_scan_alias_cmbs(self, basenames=None):
        """
        Resets all comboboxes in the scans section and repopulates them with new options
        :param basenames: filepath basenames to populate the comboboxes with
        """
        if basenames is None:
            basenames = []

        # Must first disconnect the combobox or else update_scan_aliases goes berserk because the index
        # will be reset for each combobox in the process. Reconnect after changes.
        for key, cmb in self.cmb_scanaliases_dict.items():
            cmb.currentTextChanged.disconnect(self.update_scan_aliases)
            cmb.clear()
            cmb.addItems(["Select an alias"] + basenames)
            cmb.currentTextChanged.connect(self.update_scan_aliases)
            cmb.currentTextChanged.connect(self.is_ready_import)

    def update_scan_aliases(self):
        """
        Updates the scan aliases global variable as comboboxes in the scans section are selected
        """
        for key, value in self.cmb_scanaliases_dict.items():
            if value.currentText() != "Select an alias":
                self.scan_aliases[key] = value.currentText()
            else:
                self.scan_aliases[key] = None

    def clear_visit_alias_cmbs_and_les(self):
        """
        Removes all row widgets from the visits section. Clears the lineedits dict linking directory names to user-
        preferred aliases. Clears the comboboxes dictionary specifying order.
        """
        for idx in range(self.formlay_visitaliases.rowCount()):
            self.formlay_visitaliases.removeRow(0)
        self.le_visitaliases_dict.clear()
        self.cmb_visitaliases_dict.clear()

    def reset_visit_aliases(self, basenames=None):
        """
        Resets the entire visit section. Clears previous rows if necessary. Resets the global variables for the
        lineedits and comboboxes containing mappings of the basename to the row widgets.
        :param basenames: filepath basenames to populate the row labels with and establish alias mappings with
        """
        if basenames is None:
            basenames = []

        # If this is an update, remove the previous widgets and clear the dict
        if len(self.le_visitaliases_dict) > 0:
            self.clear_visit_alias_cmbs_and_les()

        # Generate the new dict mappings of directory basename to preferred alias name and mapping
        self.le_visitaliases_dict = dict.fromkeys(basenames)
        self.cmb_visitaliases_dict = dict.fromkeys(basenames)

        # Repopulate the format layout, and establish mappings for the lineedits and the comboboxes
        for ii, key in enumerate(self.le_visitaliases_dict):
            hlay = QHBoxLayout()
            cmb = QComboBox()
            nums_to_add = [str(num) for num in range(1, len(self.le_visitaliases_dict) + 1)]
            cmb.addItems(nums_to_add)
            cmb.setCurrentIndex(ii)
            cmb.currentIndexChanged.connect(self.is_ready_import)
            le = QLineEdit()
            le.setPlaceholderText("(Optional) Specify the alias for this visit")
            hlay.addWidget(le)
            hlay.addWidget(cmb)
            self.formlay_visitaliases.addRow(key, hlay)
            # This is where the mappings are re-established
            self.le_visitaliases_dict[key] = le
            self.cmb_visitaliases_dict[key] = cmb

    ##########################
    # SECTION - MISC FUNCTIONS
    ##########################

    @staticmethod
    def infer_regex(list_of_strings):
        """
        Self-explanatory: deduces a regex string to match a provided list of strings
        :param list_of_strings: the list of string to be matched
        :return: The inferred regex string matching the all the items in the list of strings
        """
        extractor = rexpy.Extractor(list_of_strings)
        extractor.extract()
        regex = extractor.results.rex[0]
        return regex

    @Slot()
    def update_sibling_awareness(self):
        """
        Updates the awareness of what each drop-enabled lineedits contain such that certain variables cannot be dropped
        in for multiple lineedits
        """
        current_texts = [le.text() for le in self.levels.values()]
        for le in self.levels.values():
            le.sibling_awareness = current_texts

    @Slot()
    def is_ready_import(self):
        """
        Quality controls several conditions required in order to be able to run the Importer.
        """
        current_texts = [le.text() for le in self.levels.values()]
        # First requirement; raw directory must be an existent directory
        if os.path.exists(self.le_rootdir.text()):
            if not os.path.isdir(self.le_rootdir.text()):
                self.btn_run_importer.setEnabled(False)
                return
        else:
            return

        # Next requirement; a minimum of "Subject" and "Scan" must be present in the lineedits
        if not all(["Subject" in current_texts, "Scan" in current_texts]):
            self.btn_run_importer.setEnabled(False)
            return

        # Next requirement; a minimum of "ASL4D" and "T1" must have their aliases specified
        if any([self.scan_aliases["ASL4D"] is None, self.scan_aliases["T1"] is None]):
            self.btn_run_importer.setEnabled(False)
            return

        # Next requirement; if Visit is indicated, the aliases and ordering must both be unique
        if "Visit" in current_texts and len(self.cmb_visitaliases_dict) > 0:
            current_visit_aliases = [le.text() for le in self.le_visitaliases_dict.values() if le.text() != '']
            current_visit_ordering = [cmb.currentText() for cmb in self.cmb_visitaliases_dict.values()]
            if any([
                len(set(current_visit_aliases)) != len(current_visit_aliases),  # unique aliases requires
                len(set(current_visit_ordering)) != len(current_visit_ordering)  # unique ordering required
            ]):
                self.btn_run_importer.setEnabled(False)
                return

        self.btn_run_importer.setEnabled(True)

    def set_widgets_on_or_off(self, state: bool):
        """
        Convenience function for turning off widgets during an import run and then re-enabling them afterwards
        :param state: the boolean state of whether the widgets should be enabled or not
        """
        self.btn_run_importer.setEnabled(state)
        self.btn_setrootdir.setEnabled(state)
        for le in self.levels.values():
            le.setEnabled(state)

    ##################################
    # SECTION - RETRIEVAL OF VARIABLES
    ##################################

    # Returns the directory structure in preparation of running the import
    def get_directory_structure(self):
        dirnames = [le.text() for le in self.levels.values()]
        valid_dirs = []
        encountered_nonblank = False
        # Iterate backwards to remove false
        for name in reversed(dirnames):
            # Cannot have blank lines existing between the important directories
            if name == '' and encountered_nonblank:
                QMessageBox().warning(self,
                                      "Invalid directory structure entered",
                                      "You must indicate filler directories occuring between"
                                      "\nSubject/Visit/Scan directories using the Dummy label provided",
                                      QMessageBox.Ok)
                return False, []
            elif name == '' and not encountered_nonblank:
                continue
            else:
                encountered_nonblank = True
                valid_dirs.append(name)

        # Sanity check for false user input
        if any(["Subject" not in valid_dirs,
                "Scan" not in valid_dirs]):
            QMessageBox().warning(self,
                                  "Invalid directory structure entered",
                                  "A minimum of Visit and Scan directories must be present in your study for"
                                  "ExploreASL to import data correctly.")
            return False, []

        valid_dirs = list(reversed(valid_dirs))
        # print(valid_dirs)
        return True, valid_dirs

    def get_scan_aliases(self):
        """
        Retrieves a mapping of the standard scan name for ExploreASL (i.e ASL4D) and the user-specifed corresponding
        scan directory
        @return: status, whether the operation was a success; scan_aliases, the mapping
        """
        try:
            if any([self.scan_aliases["ASL4D"] is None,
                    self.scan_aliases["T1"] is None]):
                QMessageBox().warning(self,
                                      "Invalid scan aliases entered",
                                      "At minimum, the aliases corresponding to the ASL and T1-weighted scans "
                                      "should be specified",
                                      QMessageBox.Ok)
                return False, None
        except KeyError as e:
            print(f'ENCOUNTERED KEYERROR: {e}')
            return False, None

        # Filter out scans that have None to avoid problems down the line
        scan_aliases = {key: value for key, value in self.scan_aliases.items() if value is not None}

        return True, scan_aliases

    def get_visit_aliases(self):
        """
        Retrieves a mapping of the visit alias names and the user-specified preferred name
        @return: status, whether the operation was a success; visit_aliases, the mapping
        """

        visit_aliases = OrderedDict()

        # If the visit aliases dict is empty, simply return the empty dict, as visits are not mandatory to outline
        if len(self.cmb_visitaliases_dict) == 0:
            return True, visit_aliases

        # First, make sure that every number is unique:
        current_orderset = [cmb.currentText() for cmb in self.cmb_visitaliases_dict.values()]
        if len(current_orderset) != len(set(current_orderset)):
            QMessageBox().warning(self,
                                  "Invalid visits alias ordering entered",
                                  "Please check for accidental doublings",
                                  QMessageBox.Ok)
            return False, visit_aliases

        basename_keys = list(self.le_visitaliases_dict.keys())
        aliases = list(le.text() for le in self.le_visitaliases_dict.values())
        orders = list(cmb.currentText() for cmb in self.cmb_visitaliases_dict.values())

        if self.config["DeveloperMode"]:
            print(f"Inside get_visit_aliases, the following variable values were in play prior to generating the "
                  f"visit aliases dict:\n"
                  f"basename_keys: {basename_keys}\n"
                  f"aliases: {aliases}\n"
                  f"orders: {orders}")

        for num in range(1, len(orders) + 1):
            idx = orders.index(str(num))
            current_alias = aliases[idx]
            current_basename = basename_keys[idx]
            if current_alias == '':
                visit_aliases[current_basename] = f"ASL_{num}"
            else:
                visit_aliases[current_basename] = current_alias

        return True, visit_aliases

    # Utilizes the other get_ functions above to create the import parameters file
    def get_import_parms(self):
        import_parms = {}.fromkeys(["Regex", "Directory Structure", "Scan Aliases", "Ordered Visit Aliases"])
        # Get the directory structure, the scan aliases, and the visit aliases
        directory_status, valid_directories = self.get_directory_structure()
        scanalias_status, scan_aliases = self.get_scan_aliases()
        visitalias_status, visit_aliases = self.get_visit_aliases()
        if any([self.subject_regex == '',  # Subject regex must be established
                self.scan_regex == '',  # Scan regex must be established
                not directory_status,  # Getting the directory structure must have been successful
                not scanalias_status,  # Getting the scan aliases must have been successful
                not visitalias_status  # Getting the visit aliases must have been successful
                ]):
            return None

        # Otherwise, green light to create the import parameters
        import_parms["RawDir"] = self.le_rootdir.text()
        import_parms["Regex"] = [self.subject_regex, self.visit_regex, self.scan_regex]
        import_parms["Directory Structure"] = valid_directories
        import_parms["Scan Aliases"] = scan_aliases
        import_parms["Ordered Visit Aliases"] = visit_aliases

        # Save a copy of the import parms to the raw directory in question
        with open(os.path.join(self.le_rootdir.text(), "ImportConfig.json"), 'w') as w:
            json.dump(import_parms, w, indent=3)

        return import_parms

    #############################################
    # SECTION - CONCURRENT AND POST-RUN FUNCTIONS
    #############################################

    @Slot(list)
    def create_import_summary_file(self, signalled_summaries: list):
        """
        Creates the summary file. Increments the "debt" due to launching workers back towards zero. Resets widgets
        once importer workers are done.
        :param signalled_summaries: A list of dicts, each dict being all the relevant DICOM and NIFTI parameters of
        a converted directory
        """
        # Stockpile the completed summaries and increment the "debt" back towards zero
        self.import_summaries.extend(signalled_summaries)
        self.n_import_workers -= 1

        # Don't proceed until all importer workers are finished
        if self.n_import_workers > 0 or self.import_parms is None:
            return

        # Otherwise, proceed to post-import processing
        self.import_postprocessing()

    @Slot(list)
    def update_failed_runs_log(self, signalled_failed_runs: list):
        """
        Updates the attribute failed_runs in order to write the json file summarizing failed runs once everything is
        complete
        :param signalled_failed_runs: A list of dicts, each dict being the name of the DICOM directory attempted for
        conversion and the value being a description of the step in DCM2BIDS that it failed on.
        """
        self.failed_runs.extend(signalled_failed_runs)

    def import_postprocessing(self):
        """
        Performs the bulk of the post-import work, especially if the import type was specified to be BIDS
        """
        analysis_dir = os.path.join(os.path.dirname(self.import_parms["RawDir"]), "analysis")

        # Re-enable widgets and change the cwd back to the scripts directory
        self.set_widgets_on_or_off(state=True)
        os.chdir(self.config["ScriptsDir"])

        # Create the import summary
        create_import_summary(import_summaries=self.import_summaries, config=self.import_parms)

        # If there were any failures, write them to disk now
        if len(self.failed_runs) > 0:
            with open(os.path.join(analysis_dir, "import_summary_failed.json"), 'w') as failed_writer:
                json.dump(dict(self.failed_runs), failed_writer, indent=3)

        # If the settings is BIDS...
        if not self.chk_uselegacy.isChecked():
            # Ensure all M0 jsons have the appropriate "IntendedFor" field if this is in BIDS
            bids_m0_followup(analysis_dir=analysis_dir)

            # Create the template for the dataset description
            self.create_dataset_description_template(analysis_dir)

            # Create the "bidsignore" file
            with open(os.path.join(analysis_dir, ".bidsignore"), 'w') as ignore_writer:
                ignore_writer.writelines(["import_summary.tsv\n", "DataPar.json\n"])

    @staticmethod
    def create_dataset_description_template(analysis_dir):
        """
        Creates a template for the dataset description file for the user to complete at a later point in time
        :param analysis_dir: The analysis directory where the dataset description will be saved to.
        """
        template = {
            "BIDSVersion": "0.1.0",
            "License": "CC0",
            "Name": "A multi-subject, multi-modal human neuroimaging dataset",
            "Authors": [],
            "Acknowledgements": "",
            "HowToAcknowledge": "This data was obtained from [owner]. "
                                "Its accession number is [id number]'",
            "ReferencesAndLinks": ["https://www.ncbi.nlm.nih.gov/pubmed/25977808",
                                   "https://openfmri.org/dataset/ds000117/"],
            "Funding": ["UK Medical Research Council (MC_A060_5PR10)"]
        }
        with open(os.path.join(analysis_dir, "dataset_description.json"), 'w') as dataset_writer:
            json.dump(template, dataset_writer, indent=3)

    ########################
    # SECTION - RUN FUNCTION
    ########################
    def run_importer(self):
        """
        First confirms that all import parameters are set, then runs ASL2BIDS using multi-threading
        """
        # Set (or reset if this is another run) the essential variables
        self.n_import_workers = 0
        self.import_parms = None
        self.import_summaries = []
        workers = []

        # Disable the run button to prevent accidental re-runs
        self.set_widgets_on_or_off(state=False)

        # Ensure the dcm2niix path is visible
        os.chdir(os.path.join(self.config["ProjectDir"], "External", "DCM2NIIX", f"DCM2NIIX_{platform.system()}"))

        # Get the import parameters
        self.import_parms = self.get_import_parms()
        if self.import_parms is None:
            # Reset widgets back to normal and change the directory back
            self.set_widgets_on_or_off(state=True)
            os.chdir(self.config["ScriptsDir"])
            return

        # Get the dicom directories
        dicom_dirs = get_dicom_directories(config=self.import_parms)

        if self.config["DeveloperMode"]:
            print("Detected the following dicom directories:")
            pprint(dicom_dirs)
            print('\n')

        # Create workers
        dicom_dirs = list(divide(4, dicom_dirs))
        for ddirs in dicom_dirs:
            worker = Importer_Worker(ddirs,  # The list of dicom directories
                                     self.import_parms,  # The import parameters
                                     self.chk_uselegacy.isChecked())  # Whether to use legacy mode or not
            worker.signals.signal_send_summaries.connect(self.create_import_summary_file)
            worker.signals.signal_send_errors.connect(self.update_failed_runs_log)
            workers.append(worker)
            self.n_import_workers += 1

        # Launch them
        for worker in workers:
            self.threadpool.start(worker)

        print(r"""
  ______            _                          _____ _         _____ _    _ _____ 
 |  ____|          | |                  /\    / ____| |       / ____| |  | |_   _|
 | |__  __  ___ __ | | ___  _ __ ___   /  \  | (___ | |      | |  __| |  | | | |  
 |  __| \ \/ / '_ \| |/ _ \| '__/ _ \ / /\ \  \___ \| |      | | |_ | |  | | | |  
 | |____ >  <| |_) | | (_) | | |  __// ____ \ ____) | |____  | |__| | |__| |_| |_ 
 |______/_/\_\ .__/|_|\___/|_|  \___/_/    \_\_____/|______|  \_____|\____/|_____|
             | |                                                                  
             |_|                                                                  
  _____                            _                                              
 |_   _|                          | |                                             
   | |  _ __ ___  _ __   ___  _ __| |_                                            
   | | | '_ ` _ \| '_ \ / _ \| '__| __|                                           
  _| |_| | | | | | |_) | (_) | |  | |_                                            
 |_____|_| |_| |_| .__/ \___/|_|   \__|                                           
                 | |                                                              
                 |_|                                                             """)


class DraggableLabel(QLabel):
    """
    Modified QLabel to support dragging out the text content
    """

    def __init__(self, text='', parent=None):
        super(DraggableLabel, self).__init__(parent)
        self.setText(text)
        self.setStyleSheet("border-style: solid;"
                           "border-width: 2px;"
                           "border-color: black;"
                           "border-radius: 10px;"
                           "background-color: white;"
                           "margin-bottom: 5px;")
        font = QFont()
        font.setPointSize(16)
        self.setFont(font)
        self.setMinimumHeight(75)
        self.setMaximumHeight(100)
        self.setAlignment(Qt.AlignCenter)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.drag_start_position = event.pos()

    def mouseMoveEvent(self, event):
        if not (event.buttons() & Qt.LeftButton):
            return
        if (event.pos() - self.drag_start_position).manhattanLength() < QApplication.startDragDistance():
            return
        drag = QDrag(self)
        mimedata = QMimeData()
        mimedata.setText(self.text())
        drag.setMimeData(mimedata)
        drag.setHotSpot(event.pos())
        drag.exec_(Qt.CopyAction | Qt.MoveAction)


class DandD_Label2LineEdit(QLineEdit):
    """
    Modified QLineEdit to support accepting text drops from a QLabel with Drag enabled
    """

    modified_text = Signal(str, int)

    def __init__(self, superparent, parent=None, identification: int = None):
        super().__init__(parent)

        self.setAcceptDrops(True)
        self.setReadOnly(True)
        self.superparent = superparent  # This is the Importer Widget itself
        self.sibling_awareness = ['', '', '', '', '']
        self.id = identification  # This is the python index of which level after ..\\raw does this lineedit represent
        self.textChanged.connect(self.modifiedtextChanged)

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:
        if event.mimeData().hasText():
            if all([event.mimeData().text() not in self.sibling_awareness,
                    self.superparent.le_rootdir.text() != '']) or event.mimeData().text() == "Dummy":
                event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event: QDragMoveEvent) -> None:
        if event.mimeData().hasText():
            if all([event.mimeData().text() not in self.sibling_awareness,
                    self.superparent.le_rootdir.text() != '']) or event.mimeData().text() == "Dummy":
                event.accept()
                event.setDropAction(Qt.CopyAction)
        else:
            event.ignore()

    def dropEvent(self, event: QDropEvent) -> None:
        if event.mimeData().hasText():
            if all([event.mimeData().text() not in self.sibling_awareness,
                    self.superparent.le_rootdir.text() != '']) or event.mimeData().text() == "Dummy":
                event.accept()
                self.setText(event.mimeData().text())
        else:
            event.ignore()

    def modifiedtextChanged(self):
        self.modified_text.emit(self.text(), self.id)
