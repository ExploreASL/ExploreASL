from PySide2.QtWidgets import *
from PySide2.QtGui import *
from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit, DandD_FileExplorer2ListWidget
from ExploreASL_GUI.xASL_GUI_HelperFuncs_StringOps import set_os_dependent_text
import json
import os
from glob import iglob, glob
from tdda import rexpy
from more_itertools import peekable


class xASL_Parms(QMainWindow):
    def __init__(self, parent_win=None):
        # Parent window is fed into the constructor to allow for communication with parent window devices
        super().__init__(parent=parent_win)
        if parent_win is not None:
            self.config = self.parent().config
        else:
            with open("ExploreASL_GUI_masterconfig.json") as f:
                self.config = json.load(f)

        # Window Size and initial visual setup
        self.setMinimumSize(512, 960)
        self.cw = QWidget(self)
        self.setCentralWidget(self.cw)
        self.mainlay = QVBoxLayout(self.cw)
        self.setLayout(self.mainlay)
        self.setWindowTitle("Explore ASL - Parameter File Maker")

        # Buttons for executing the fundamental functions
        btn_font = QFont()
        btn_font.setPointSize(16)
        self.btn_make_parms = QPushButton("Run", self.cw, clicked=self.saveparms2json)
        self.btn_load_parms = QPushButton("Load from existing Json", self.cw, clicked=self.loadjson2parms)
        for btn in [self.btn_make_parms, self.btn_load_parms]:
            btn.setFont(btn_font)
            btn.setMinimumHeight(50)

        # TabWidget Setup and containers
        self.tab_main = QTabWidget(self.cw)
        self.mainlay.addWidget(self.tab_main)
        self.cont_basic = QWidget()
        self.cont_advanced = QWidget()
        self.tab_main.addTab(self.cont_basic, "Basic Settings")
        self.tab_main.addTab(self.cont_advanced, "Advanced Settings")
        self.mainlay.addWidget(self.btn_make_parms)
        self.mainlay.addWidget(self.btn_load_parms)

        # Misc Players
        self.import_error_logger = []
        self.asl_json_sidecar_data = {}
        self.can_update_slicereadouttime = False

        self.UI_Setup_Basic()
        self.UI_Setup_Advanced()

        # After all UI is set up, make certain connections
        self.le_study_dir.textChanged.connect(self.update_asl_json_sidecar_data)

    def UI_Setup_Basic(self):
        self.formlay_basic = QFormLayout(self.cont_basic)
        self.hlay_easl_dir, self.le_easl_dir, self.btn_easl_dir = self.make_droppable_clearable_le(
            btn_connect_to=self.set_exploreasl_dir,
            default=''
        )
        self.le_studyname = QLineEdit(text="My Study")
        self.hlay_study_dir, self.le_study_dir, self.btn_study_dir = self.make_droppable_clearable_le(
            btn_connect_to=self.set_study_dir,
            default=''
        )
        self.le_subjectregex = QLineEdit(text='\\d+')
        self.le_subjectregex.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.lst_included_subjects = DandD_FileExplorer2ListWidget()
        self.lst_included_subjects.alert_regex.connect(self.update_regex)
        self.btn_included_subjects = QPushButton("Clear Subjects", clicked=self.clear_included)
        self.lst_excluded_subjects = DandD_FileExplorer2ListWidget()
        self.btn_excluded_subjects = QPushButton("Clear Excluded", clicked=self.clear_excluded)
        self.le_run_names = QLineEdit(text="ASL_1",
                                      placeholderText="Indicate run names, each separated by a comma and space")
        self.le_run_options = QLineEdit(placeholderText="Indicate option names, each separated by a comma and space")
        self.cmb_vendor = self.make_cmb_and_items(["Siemens", "Philips", "GE", "GE_WIP"])
        self.cmb_sequencetype = self.make_cmb_and_items(["3D GRaSE", "2D EPI", "3D Spiral"])
        self.cmb_labelingtype = self.make_cmb_and_items(["Q2 TIPS PASL", "pCASL/CASL"])
        self.cmb_labelingtype.currentTextChanged.connect(self.autocalc_slicereadouttime)
        self.cmb_m0_isseparate = self.make_cmb_and_items(["Proton density scan (M0) was acquired",
                                                          "Use mean control ASL as proton density mimic"])
        self.cmb_m0_posinasl = self.make_cmb_and_items(
            ["M0 exists as a separate scan", "M0 is the first ASL control-label pair",
             "M0 is the first ASL scan volume", "M0 is the second ASL scan volume"])
        self.cmb_quality = self.make_cmb_and_items(["Low", "High"])

        for desc, widget in zip(["ExploreASL Directory", "Name of Study", "Study Directory", "Subject Regex",
                                 "Subjects to Assess\n(Drag and Drop Directories)",
                                 "Subjects to Exclude\n(Drag and Drop Directories)",
                                 "Run Names", "Run Options", "Vendor", "Sequence Type", "Labelling Type",
                                 "M0 was acquired?", "M0 Position in ASL", "Quality"],
                                [self.hlay_easl_dir, self.le_studyname, self.hlay_study_dir,
                                 self.le_subjectregex, self.lst_included_subjects, self.lst_excluded_subjects,
                                 self.le_run_names, self.le_run_options, self.cmb_vendor, self.cmb_sequencetype,
                                 self.cmb_labelingtype, self.cmb_m0_isseparate, self.cmb_m0_posinasl,
                                 self.cmb_quality]):
            self.formlay_basic.addRow(desc, widget)
        self.formlay_basic.insertRow(5, "", self.btn_included_subjects)
        self.formlay_basic.insertRow(7, "", self.btn_excluded_subjects)

    def UI_Setup_Advanced(self):
        # First, set up the groupboxes and add them to the advanced tab layout
        self.vlay_advanced = QVBoxLayout(self.cont_advanced)
        self.grp_sequenceparms = QGroupBox(title="Sequence Parameters")
        self.grp_quantparms = QGroupBox(title="Quantification Parameters")
        self.grp_m0parms = QGroupBox(title="M0 Parameters")
        self.grp_procparms = QGroupBox(title="Processing Parameters")
        self.grp_envparms = QGroupBox(title="Environment Parameters")
        for grp in [self.grp_sequenceparms, self.grp_quantparms, self.grp_m0parms,
                    self.grp_procparms, self.grp_envparms]:
            self.vlay_advanced.addWidget(grp)

        # Set up the Sequence Parameters
        self.formlay_sequenceparms = QFormLayout(self.grp_sequenceparms)
        self.cmb_nsup_pulses = self.make_cmb_and_items(["0", "2", "4", "5"], 1)
        self.cmb_readout_dim = self.make_cmb_and_items(["3D", "2D"])
        self.spinbox_initialpld = QDoubleSpinBox(maximum=2500, minimum=0, value=1800)
        self.spinbox_initialpld.valueChanged.connect(self.autocalc_slicereadouttime)
        self.spinbox_labdur = QDoubleSpinBox(maximum=2000, minimum=0, value=800)
        self.spinbox_labdur.valueChanged.connect(self.autocalc_slicereadouttime)
        self.hlay_slice_readout = QHBoxLayout()
        self.cmb_slice_readout = self.make_cmb_and_items(["Use Indicated Value", "Use Shortest TR"])
        self.spinbox_slice_readout = QDoubleSpinBox(maximum=1000, minimum=0, value=37)
        self.hlay_slice_readout.addWidget(self.cmb_slice_readout)
        self.hlay_slice_readout.addWidget(self.spinbox_slice_readout)
        for description, widget in zip(["Number of Suppression Pulses", "Readout Dimension",
                                        "Initial Post-Labeling Delay (ms)", "Labeling Duration (ms)",
                                        "Slice Readout Time (ms)"],
                                       [self.cmb_nsup_pulses, self.cmb_readout_dim, self.spinbox_initialpld,
                                        self.spinbox_labdur, self.hlay_slice_readout]):
            self.formlay_sequenceparms.addRow(description, widget)

        # Set up the Quantification Parameters
        self.formlay_quantparms = QFormLayout(self.grp_quantparms)
        self.spinbox_lambda = QDoubleSpinBox(maximum=1, minimum=0, value=0.9, singleStep=0.01)
        self.spinbox_artt2 = QDoubleSpinBox(maximum=100, minimum=0, value=50, singleStep=0.1)
        self.spinbox_bloodt1 = QDoubleSpinBox(maximum=2000, minimum=0, value=1650, singleStep=0.1)
        self.spinbox_tissuet1 = QDoubleSpinBox(maximum=2000, minimum=0, value=1240, singleStep=0.1)
        self.cmb_ncomparts = self.make_cmb_and_items(["1", "2"], 0)
        self.le_quantset = QLineEdit(text="1 1 1 1 1")
        for description, widget in zip(["Lambda", "Arterial T2*", "Blood T1",
                                        "Tissue T1", "Number of Compartments", "Quantification Settings"],
                                       [self.spinbox_lambda, self.spinbox_artt2, self.spinbox_bloodt1,
                                        self.spinbox_tissuet1, self.cmb_ncomparts, self.le_quantset]):
            self.formlay_quantparms.addRow(description, widget)

        # Set up the remaining M0 Parameters
        self.formlay_m0parms = QFormLayout(self.grp_m0parms)
        self.cmb_m0_algorithm = self.make_cmb_and_items(["New Image Processing", "Standard Processing"], 0)
        self.spinbox_gmscale = QDoubleSpinBox(maximum=100, minimum=0.01, value=1, singleStep=0.01)
        for description, widget in zip(["M0 Processing Algorithm", "GM Scale Factor"],
                                       [self.cmb_m0_algorithm, self.spinbox_gmscale]):
            self.formlay_m0parms.addRow(description, widget)

        # Set up the Processing Parameters
        self.vlay_procparms = QVBoxLayout(self.grp_procparms)
        self.scroll_procparms = QScrollArea()
        self.scroll_procparms.setWidgetResizable(True)
        self.vlay_procparms.addWidget(self.scroll_procparms)
        self.cont_procparms = QWidget()
        self.scroll_procparms.setWidget(self.cont_procparms)
        self.formlay_procparms = QFormLayout(self.cont_procparms)
        self.chk_removespikes = QCheckBox(checked=True)
        self.spinbox_spikethres = QDoubleSpinBox(maximum=1, minimum=0, value=0.01, singleStep=0.01)
        self.chk_motioncorrect = QCheckBox(checked=True)
        self.chk_deltempfiles = QCheckBox(checked=True)
        self.chk_skipnoflair = QCheckBox(checked=False)
        self.chk_skipnoasl = QCheckBox(checked=True)
        self.chk_skipnom0 = QCheckBox(checked=False)
        self.chk_uset1dartel = QCheckBox(checked=True)
        self.chk_usepwidartel = QCheckBox(checked=False)
        self.chk_usebilatfilter = QCheckBox(checked=False)
        self.cmb_segmethod = self.make_cmb_and_items(["CAT12", "SPM12"], 0)
        self.cmb_imgcontrast = self.make_cmb_and_items(["Automatic", "Control --> T1w", "CBF --> pseudoCBF",
                                                        "Force CBF --> pseudoCBF"], 0)
        self.cmb_affineregbase = self.make_cmb_and_items(["Based on PWI CoV", "Enabled", "Disabled"])
        self.chk_regm0toasl = QCheckBox(checked=True)
        self.chk_usemniasdummy = QCheckBox(checked=False)
        for description, widget in zip(["Remove Spikes", "Spike Removal Threshold", "Correct for Motion",
                                        "Delete Temporary Files", "Skip Subjects without FLAIR",
                                        "Skip Subjects without ASL", "Skip subjects without M0", "Use T1 DARTEL",
                                        "Use PWI DARTEL", "Use Bilateral Filter", "Segmentation Method",
                                        "Image Contrast used for", "Use Affine Registration", "Register M0 to ASL",
                                        "Use MNI as Dummy Template"],
                                       [self.chk_removespikes, self.spinbox_spikethres, self.chk_motioncorrect,
                                        self.chk_deltempfiles, self.chk_skipnoflair, self.chk_skipnoasl,
                                        self.chk_skipnom0, self.chk_uset1dartel, self.chk_usepwidartel,
                                        self.chk_usebilatfilter, self.cmb_segmethod, self.cmb_imgcontrast,
                                        self.cmb_affineregbase, self.chk_regm0toasl, self.chk_usemniasdummy]):
            self.formlay_procparms.addRow(description, widget)

        # Set up the Environment Parameters
        self.formlay_envparms = QFormLayout(self.grp_envparms)
        self.chk_detectfsl = QCheckBox(checked=True)
        self.chk_outputcbfmaps = QCheckBox(checked=False)
        for desc, widget in zip(["Detect FSL Automatically", "Output CBF native space maps?"],
                                [self.chk_detectfsl, self.chk_outputcbfmaps]):
            self.formlay_envparms.addRow(desc, widget)

    ################################
    # Json Sidecar Related Functions
    ################################
    def update_asl_json_sidecar_data(self, analysis_dir_text):
        """
        Receives a signal from the le_study_dir lineedit and will accordingly update several fields
        @param analysis_dir_text: the text updated from the analysis directory
        """

        # First set of checks
        if any([analysis_dir_text == '',  # Do not react to blank lines
                not os.path.exists(analysis_dir_text),  # Do not react to nonexistent paths
                ]):
            return

        # Second set of checks
        if any([not os.path.isdir(analysis_dir_text),  # Must be a directory
                os.path.basename(analysis_dir_text) != "analysis"
                ]):
            return

        if self.config["DeveloperMode"]:
            print(f"Detected an update to the specified analysis directory. Attempting to find asl json sidecars and "
                  f"infer appropriate field values from within.\n")

        # Retrieve any asl json sidecar
        asl_sides_legacy = iglob(os.path.join(analysis_dir_text, "**", "ASL4D.json"), recursive=True)
        asl_sides_legacy = peekable(asl_sides_legacy)
        asl_sides_bids = iglob(os.path.join(analysis_dir_text, "**", "*_asl.json"), recursive=True)
        asl_sides_bids = peekable(asl_sides_bids)

        # Disengage if there is no luck finding any sidecar
        if not asl_sides_legacy:
            if not asl_sides_bids:
                return
            else:
                asl_sidecar = next(asl_sides_bids)
        else:
            asl_sidecar = next(asl_sides_legacy)
        try:
            with open(asl_sidecar) as sidecar_reader:
                self.asl_json_sidecar_data = json.load(sidecar_reader)
        except json.decoder.JSONDecodeError as json_e:
            QMessageBox.warning(self.parent(),
                                "Json sidecars not in proper json format",
                                f"ExploreASL GUI has detected that the json sidecars present for this dataset "
                                f"are not in appropriate format. The following inconsistency was found:\n{json_e}",
                                QMessageBox.Ok)
            return

        # Now we can update a couple of fields
        # First, the vendor
        try:
            idx = self.cmb_vendor.findText(self.asl_json_sidecar_data["Manufacturer"])
            if idx != -1:
                self.cmb_vendor.setCurrentIndex(idx)
        except KeyError:
            if self.config["DeveloperMode"]:
                print(f"Warning in update_asl_json_sidecar_data. The field: Manufacturer was not present in the "
                      f"detected asl json sidecar.\n")

        # Next, the readout dimension
        try:
            idx = self.cmb_readout_dim.findText(self.asl_json_sidecar_data["MRAcquisitionType"])
            if idx != -1:
                self.cmb_readout_dim.setCurrentIndex(idx)
        except KeyError:
            if self.config["DeveloperMode"]:
                print(f"Warning in update_asl_json_sidecar_data. The field: MRAcquisitionType was not present in the "
                      f"detected asl json sidecar.\n")

        # Next the inversion time (i.e Post-Label Duration)
        try:
            value = self.asl_json_sidecar_data["InversionTime"]
            self.spinbox_labdur.setValue(value * 1000)
        except KeyError:
            if self.config["DeveloperMode"]:
                print(f"Warning in update_asl_json_sidecar_data. The field: InversionTime was not present in the "
                      f"detected asl json sidecar.\n")

        # Next get a few essentials for auto-calculating the SliceReadoutTime
        try:
            has_tr = "RepetitionTime" in self.asl_json_sidecar_data
            has_nslices = "NumberOfSlices" in self.asl_json_sidecar_data
            print(has_tr, has_nslices)
            if has_tr and has_nslices:
                self.can_update_slicereadouttime = True
            else:
                self.can_update_slicereadouttime = False
        except KeyError:
            pass

        # Retrieve any M0 json sidecar
        m0_sides_legacy = iglob(os.path.join(analysis_dir_text, "**", "M0.json"), recursive=True)
        m0_sides_legacy = peekable(m0_sides_legacy)
        m0_sides_bids = iglob(os.path.join(analysis_dir_text, "**", "*_m0scan.json"), recursive=True)
        m0_sides_bids = peekable(m0_sides_bids)

        # Disengage if there is no luck finding any m0 sidecar
        if not m0_sides_legacy:
            if not m0_sides_bids:
                return
            else:
                m0_sidecar = next(m0_sides_bids)
        else:
            m0_sidecar = next(m0_sides_legacy)

        if m0_sidecar:
            idx = self.cmb_m0_isseparate.findText("Proton density scan (M0) was acquired")
            if idx != -1:
                self.cmb_m0_isseparate.setCurrentIndex(idx)

            # If there was an M0 sidecar, it stands to reason there was also background suppression and that field
            # should be set appropriately
            try:
                manufac = self.asl_json_sidecar_data["Manufacturer"]
                acq_type = self.asl_json_sidecar_data["MRAcquisitionType"]
                if manufac == "GE":
                    idx = self.cmb_nsup_pulses.findText('5')
                elif manufac == "Philips" and acq_type == "2D":
                    idx = self.cmb_nsup_pulses.findText("2")
                elif manufac == "Philips" and acq_type == "3D":
                    idx = self.cmb_nsup_pulses.findText("4")
                elif manufac == "Siemens":
                    idx = self.cmb_nsup_pulses.findText("2")
                else:
                    if self.config["DeveloperMode"]:
                        print(f"Warning in update_asl_json_sidecar_data. An M0 json sidecar was detected, but its "
                              f"Manufacturer field was not present, preventing the appropriate setting for the number "
                              f"of background suppression pulses to be established.\n")
                    return
                self.cmb_nsup_pulses.setCurrentIndex(idx)

            except KeyError:
                pass

    def autocalc_slicereadouttime(self):
        if not self.can_update_slicereadouttime or self.cmb_labelingtype.currentText() != "pCASL/CASL":
            return

        tr = self.asl_json_sidecar_data["RepetitionTime"]*1000
        labdur = self.spinbox_labdur.value()
        ini_pld = self.spinbox_initialpld.value()
        nslices = self.asl_json_sidecar_data["NumberOfSlices"]
        readouttime = round((tr - labdur - ini_pld) / nslices, 2)

        self.spinbox_slice_readout.setValue(readouttime)

    ################
    # Misc Functions
    ################
    # Clears the currently-included subjects list and resets the regex
    def clear_included(self):
        self.lst_included_subjects.clear()
        self.le_subjectregex.clear()

    # Clears the currently-excluded subjects list
    def clear_excluded(self):
        self.lst_study_parms_exclusions.clear()

    # Updates the current recognized regex
    def update_regex(self):
        n_subjects = self.lst_included_subjects.count()
        if n_subjects == 0:
            return
        subject_list = [self.lst_included_subjects.item(idx).text() for idx in range(n_subjects)]
        extractor = rexpy.Extractor(subject_list)
        extractor.extract()
        inferred_regex = extractor.results.rex[0]
        del extractor
        if inferred_regex:
            self.le_subjectregex.setText(inferred_regex)

    def set_exploreasl_dir(self):
        exploreasl_filepath = QFileDialog.getExistingDirectory(self,
                                                               "Select ExploreASL directory",
                                                               self.config["DefaultRootDir"],
                                                               QFileDialog.ShowDirsOnly)
        # Return if the user has cancelled the operation
        if exploreasl_filepath == '':
            return

        # Quality control
        if os.path.exists(exploreasl_filepath):
            if any([not os.path.isdir(exploreasl_filepath),  # Path must be a directory
                    "ExploreASL_Master.m" not in os.listdir(exploreasl_filepath)  # Must contain the master script
                    ]):
                QMessageBox().warning(self,
                                      "Invalid Directory Selected",
                                      "Either the path you have specified is not a directory or if it is an ExploreASL"
                                      "directory, it does not contain the required scripts for processing a study",
                                      QMessageBox.Ok)
                return
            else:
                set_os_dependent_text(linedit=self.le_easl_dir,
                                      config_ossystem=self.config["Platform"],
                                      text_to_set=exploreasl_filepath)
        else:
            QMessageBox().warning(self,
                                  "The filepath you specified does not exist",
                                  "Please select an existent ExploreASL directory",
                                  QMessageBox.Ok)
            return

    def set_study_dir(self):
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
                set_os_dependent_text(linedit=self.le_study_dir,
                                      config_ossystem=self.config["Platform"],
                                      text_to_set=analysisdir_filepath)
        else:
            QMessageBox().warning(self,
                                  "The filepath you specified does not exist",
                                  "Please select an existent ExploreASL directory",
                                  QMessageBox.Ok)
            return

    #######################################################################
    # Main Functions for this module - saving to json and loading from json
    #######################################################################
    def saveparms2json(self):
        json_parms = {
            "MyPath": self.le_easl_dir.text(),
            "name": self.le_studyname.text(),
            "D": {"ROOT": self.le_study_dir.text()},
            "subject_regexp": self.le_subjectregex.text(),
            "SESSIONS": self.prep_run_names(),
            "session": {"options": self.prep_run_options()},
            "exclusion": [self.lst_excluded_subjects.item(row).text() for row in
                          range(self.lst_excluded_subjects.count())],
            "M0_conventionalProcessing":
                {"New Image Processing": 0, "Standard Processing": 1}[self.cmb_m0_algorithm.currentText()],
            "M0": {"Proton density scan (M0) was acquired": "separate_scan",
                   "Use mean control ASL as proton density mimic": "UseControlAsM0"}
            [self.cmb_m0_isseparate.currentText()],
            "M0_GMScaleFactor": float(self.spinbox_gmscale.value()),
            "readout_dim": self.cmb_readout_dim.currentText(),
            "Vendor": self.cmb_vendor.currentText(),
            "Sequence": {"3D Spiral": "3D_spiral", "3D GRaSE": "3D_GRASE", "2D EPI": "2D_EPI"}
            [self.cmb_sequencetype.currentText()],
            "Quality": {"High": 1, "Low": 0}[self.cmb_quality.currentText()],
            "DELETETEMP": int(self.chk_deltempfiles.isChecked()),
            "SPIKE_REMOVAL": int(self.chk_removespikes.isChecked()),
            "SpikeRemovalThreshold": float(self.spinbox_spikethres.value()),
            "SkipIfNoFlair": int(self.chk_skipnoflair.isChecked()),
            "SkipIfNoM0": int(self.chk_skipnom0.isChecked()),
            "SkipIfNoASL": int(self.chk_skipnoasl.isChecked()),
            "T1_DARTEL": int(self.chk_uset1dartel.isChecked()),
            "PWI_DARTEL": int(self.chk_usepwidartel.isChecked()),
            "BILAT_FILTER": int(self.chk_usebilatfilter.isChecked()),
            "motion_correction": int(self.chk_motioncorrect.isChecked()),
            "SegmentSPM12": {"SPM12": 1, "CAT12": 0}[self.cmb_segmethod.currentText()],
            "bRegistrationContrast": {"Automatic": 2, "Control --> T1w": 0, "CBF --> pseudoCBF": 1,
                                      "Force CBF --> pseudoCBF": 3}[self.cmb_imgcontrast.currentText()],
            "bAffineRegistration": {"Based on PWI CoV": 2, "Enabled": 1, "Disabled": 0}
            [self.cmb_affineregbase.currentText()],
            "bRegisterM02ASL": int(self.chk_removespikes.isChecked()),
            "bUseMNIasDummyStructural": int(self.chk_usemniasdummy.isChecked()),
            "ApplyQuantification": self.prep_quantparms(),
            "Q": {
                "BackGrSupprPulses": int(self.cmb_nsup_pulses.currentText()),
                "LabelingType": {"Q2 TIPS PASL": "PASL", "pCASL/CASL": "CASL"}[self.cmb_labelingtype.currentText()],
                "Initial_PLD": float(self.spinbox_initialpld.value()),
                "LabelingDuration": float(self.spinbox_labdur.value()),
                "SliceReadoutTime": float(self.spinbox_slice_readout.value()),
                "Lambda": float(self.spinbox_lambda.value()),
                "T2art": float(self.spinbox_artt2.value()),
                "TissueT1": float(self.spinbox_tissuet1.value()),
                "nCompartments": int(self.cmb_ncomparts.currentText())
            }
        }

        if self.cmb_m0_posinasl.currentText() != "M0 exists as a separate scan":
            parms_m0_pos_translate = {"M0 is the first ASL control-label pair": "[1 2]",
                                      "M0 is the first ASL scan volume": 1, "M0 is the second ASL scan volume": 2}
            json_parms["M0PositionInASL4D"] = parms_m0_pos_translate[self.cmb_m0_posinasl.currentText()]
        try:
            with open(os.path.join(self.le_study_dir.text(), "DataPar.json"), 'w') as w:
                json.dump(json_parms, w, indent=3)

            # Also, if this is BIDS, write to the root level asl.json
            asl_json = os.path.join(self.le_study_dir.text(), "asl.json")
            if os.path.exists(asl_json):
                asl_parms = {
                    "LabelingType": json_parms["Q"]["LabelingType"],
                    "PostLabelingDelay": json_parms["Q"]["Initial_PLD"],
                    "BackgroundSuppression": json_parms["Q"]["BackGrSupprPulses"] == 0}
                with open(asl_json, 'w') as asl_json_writer:
                    json.dump(asl_parms, asl_json_writer, indent=3)

        except FileNotFoundError:
            QMessageBox().warning(self,
                                  "Could not save parameters to json",
                                  f"Check whether the following path exists in your computer:\n"
                                  f"{self.le_study_dir.text()}",
                                  QMessageBox.Ok)
            return

        QMessageBox().information(self,
                                  "DataPar.json successfully saved",
                                  f"The parameter file was successfully saved to:\n"
                                  f"{self.le_study_dir.text()}",
                                  QMessageBox.Ok)

    def loadjson2parms(self):
        self.import_error_logger.clear()
        json_filepath, _ = QFileDialog.getOpenFileName(QFileDialog(),
                                                       "Select the JSON parameters file",
                                                       self.config["DefaultRootDir"],
                                                       "Json files (*.json)")
        if json_filepath == '':
            return
        if not os.path.exists(json_filepath) and json_filepath != '':
            QMessageBox().warning(self,
                                  "Incorrect File Selected",
                                  "The file selected was either not a json file "
                                  "or did not contain the essential parameters",
                                  QMessageBox.Ok)
            return

        with open(json_filepath, 'r') as reader:
            parms: dict = json.load(reader)

        json_setter = {
            "MyPath": self.le_easl_dir.setText,
            "name": self.le_studyname.setText,
            "D": {"ROOT": self.le_study_dir.setText},
            "subject_regexp": self.le_subjectregex.setText,
            "SESSIONS": self.get_run_names,
            "session": {"options": self.get_run_options},
            "exclusion": self.lst_excluded_subjects.addItems,
            "M0_conventionalProcessing": self.get_m0_algo,
            "M0": self.get_m0_isseparate,
            "M0_GMScaleFactor": self.spinbox_gmscale.setValue,
            "M0PositionInASL4D": self.get_m0_posinasl,
            "readout_dim": self.get_readout_dim,
            "Vendor": self.get_vendor,
            "Sequence": self.get_sequence,
            "Quality": self.get_quality,
            "DELETETEMP": self.chk_deltempfiles.setChecked,
            "SPIKE_REMOVAL": self.chk_removespikes.setChecked,
            "SpikeRemovalThreshold": self.spinbox_spikethres.setValue,
            "SkipIfNoFlair": self.chk_skipnoflair.setChecked,
            "SkipIfNoM0": self.chk_skipnom0.setChecked,
            "SkipIfNoASL": self.chk_skipnoasl.setChecked,
            "T1_DARTEL": self.chk_uset1dartel.setChecked,
            "PWI_DARTEL": self.chk_usepwidartel.setChecked,
            "BILAT_FILTER": self.chk_usebilatfilter.setChecked,
            "motion_correction": self.chk_motioncorrect.setChecked,
            "SegmentSPM12": self.get_segmethod,
            "bRegistrationContrast": self.get_imgconstrast,
            "bAffineRegistration": self.get_affinereg,
            "bRegisterM02ASL": self.chk_regm0toasl.setChecked,
            "bUseMNIasDummyStructural": self.chk_usemniasdummy.setChecked,
            "ApplyQuantification": self.get_quantparms,
            "Q": {
                "BackGrSupprPulses": self.get_backgroundpulses,
                "LabelingType": self.get_labelingtype,
                "Initial_PLD": self.spinbox_initialpld.setValue,
                "LabelingDuration": self.spinbox_labdur.setValue,
                "SliceReadoutTime": self.spinbox_slice_readout.setValue,
                "Lambda": self.spinbox_lambda.setValue,
                "T2art": self.spinbox_artt2.setValue,
                "TissueT1": self.spinbox_tissuet1.setValue,
                "nCompartments": self.get_ncompartments
            }
        }

        to_bool = ["DELETETEMP", "SPIKE_REMOVAL", "SkipIfNoFlair", "SkipIfNoM0", "SkipIfNoASL", "T1_DARTEL",
                   "PWI_DARTEL", "BILAT_FILTER", "motion_correction", "bRegisterM02ASL", "bUseMNIasDummyStructural"]

        for key, call in parms.items():
            try:
                if isinstance(call, dict):
                    for subkey, subcall in call.items():
                        try:
                            json_setter[key][subkey](subcall)
                        except TypeError as te:
                            print(f"TypeError encountered with key {key} and subkey {subkey}")
                            print(te)
                        except KeyError as ke:
                            print(f"TypeError encountered with key {key} and subkey {subkey}")
                            print(ke)
                else:
                    if key in to_bool:
                        json_setter[key](bool(call))
                    else:
                        try:
                            json_setter[key](call)
                        except KeyError as ke:
                            print(f"TypeError encountered with key {key}")
                            print(ke)

            except TypeError as te:
                print(f"TypeError encountered with key: {key}")
                print(te)

        if len(self.import_error_logger) > 0:
            errors = "\n -".join(self.import_error_logger)
            QMessageBox().warning(self,
                                  "Errors were encountered importing certain values",
                                  f"The following fields could not be properly imported:\n -{errors}",
                                  QMessageBox.Ok)

    #############################################
    # Convenience methods for translation to json
    #############################################
    def prep_run_names(self):
        if "," not in self.le_run_names.text():
            return [self.le_run_names.text()]
        else:
            return self.le_run_names.text().split(", ")

    def prep_run_options(self):
        if "," not in self.le_run_options.text():
            return [self.le_run_options.text()]
        else:
            return self.le_run_options.text().split(", ")

    def prep_quantparms(self):
        parms_logvec = self.le_quantset.text().split(" ")
        if all([len(parms_logvec) == 5,  # Must be 5 1s or 0s
                all([x in ['1', '0'] for x in parms_logvec])]):  # Check that all are 1s or 0s
            return [int(option) for option in parms_logvec]
        else:
            QMessageBox().warning(self,
                                  "Incorrect Input for Quantification Settings",
                                  "Must be a series of five 1s or 0s separated by single spaces",
                                  QMessageBox.Ok)

    ###############################################
    # Convenience methods for translation from json
    ###############################################
    def get_run_names(self, value):
        if isinstance(value, list) and len(value) > 1:
            self.le_run_names.setText(", ".join(value))
        else:
            self.le_run_names.setText(value[0])

    def get_run_options(self, value):
        if isinstance(value, list) and len(value) > 1:
            self.le_run_options.setText(", ".join(value))
        else:
            self.le_run_options.setText(value[0])

    def get_m0_algo(self, value):
        translator = {0: "New Image Processing", 1: "Standard Processing"}
        text = translator[value]
        index = self.cmb_m0_algorithm.findText(text)
        if index == -1:
            self.import_error_logger.append("M0 Processing Algorithm")
            return
        self.cmb_m0_algorithm.setCurrentIndex(index)

    def get_m0_isseparate(self, value):
        translator = {"separate_scan": "Proton density scan (M0) was acquired",
                      "UseControlAsM0": "Use mean control ASL as proton density mimic"}
        text = translator[value]
        index = self.cmb_m0_isseparate.findText(text)
        if index == -1:
            self.import_error_logger.append("M0 was acquired?")
            return
        self.cmb_m0_isseparate.setCurrentIndex(index)

    def get_m0_posinasl(self, value):
        translator = {"[1 2]": "M0 is the first ASL control-label pair", 1: "M0 is the first ASL scan volume",
                      2: "M0 is the second ASL scan volume"}
        text = translator[value]
        index = self.cmb_m0_posinasl.findText(text)
        if index == -1:
            self.import_error_logger.append("M0 Position in ASL")
            return
        self.cmb_m0_posinasl.setCurrentIndex(index)

    def get_readout_dim(self, value):
        index = self.cmb_readout_dim.findText(value)
        if index == -1:
            self.import_error_logger.append("Readout Dimension")
            return
        self.cmb_readout_dim.setCurrentIndex(index)

    def get_vendor(self, value):
        index = self.cmb_vendor.findText(value)
        if index == -1:
            self.import_error_logger.append("Vendor")
            return
        self.cmb_vendor.setCurrentIndex(index)

    def get_sequence(self, value):
        translator = {"3D_spiral": "3D Spiral", "3D_GRASE": "3D GRaSE", "2D_EPI": "2D EPI"}
        text = translator[value]
        index = self.cmb_sequencetype.findText(text)
        if index == -1:
            self.import_error_logger.append("Sequence Type")
            return
        self.cmb_sequencetype.setCurrentIndex(index)

    def get_quality(self, value):
        translator = {0: "Low", 1: "High"}
        text = translator[value]
        index = self.cmb_quality.findText(text)
        if index == -1:
            self.import_error_logger.append("Quality")
            return
        self.cmb_quality.setCurrentIndex(index)

    def get_segmethod(self, value):
        translator = {1: "SPM12", 0: "CAT12"}
        text = translator[value]
        index = self.cmb_segmethod.findText(text)
        if index == -1:
            self.import_error_logger.append("Segmentation Method")
            return
        self.cmb_segmethod.setCurrentIndex(index)

    def get_imgconstrast(self, value):
        translator = {2: "Automatic", 0: "Control --> T1w", 1: "CBF --> pseudoCBF", 3: "Force CBF --> pseudoCBF"}
        text = translator[value]
        index = self.cmb_imgcontrast.findText(text)
        if index == -1:
            self.import_error_logger.append("Image Contrast used for")
            return
        self.cmb_imgcontrast.setCurrentIndex(index)

    def get_affinereg(self, value):
        translator = {0: "Disabled", 1: "Enabled", 2: "Based on PWI CoV"}
        text = translator[value]
        index = self.cmb_affineregbase.findText(text)
        if index == -1:
            self.import_error_logger.append("Use Affine Registration")
            return
        self.cmb_affineregbase.setCurrentIndex(index)

    def get_quantparms(self, value):
        if any([len(value) != 5, sum([str(val).isdigit() for val in value]) != 5]):
            self.import_error_logger.append("Quantification Settings")
            return
        str_value = " ".join([str(val) for val in value])
        self.le_quantset.setText(str_value)

    def get_backgroundpulses(self, value):
        index = self.cmb_nsup_pulses.findText(str(value))
        if index == -1:
            self.import_error_logger.append("Number of Suppression Pulses")
            return
        self.cmb_nsup_pulses.setCurrentIndex(index)

    def get_labelingtype(self, value):
        translator = {"PASL": "Q2 TIPS PASL", "CASL": "pCASL/CASL"}
        text = translator[value]
        index = self.cmb_labelingtype.findText(text)
        if index == -1:
            self.import_error_logger.append("Labeling Type")
            return
        self.cmb_labelingtype.setCurrentIndex(index)

    def get_ncompartments(self, value):
        index = self.cmb_ncomparts.findText(str(value))
        if index == -1:
            self.import_error_logger.append("Number of Compartments")
            return
        self.cmb_ncomparts.setCurrentIndex(index)

    ############################################
    # Convenience methods for generating widgets
    ############################################

    @staticmethod
    def make_droppable_clearable_le(le_connect_to=None, btn_connect_to=None, default=''):
        """
        Convenience function for creating a lineedit-button pair within a HBoxlayout
        @param le_connect_to: Optional; the function that the linedit's textChanged signal should connect to
        @param btn_connect_to: Optional; the function that the button's clicked signal should connect to
        @param default: the default text to specify
        @return: HBoxlayout, DandD_FileExplorer2LineEdit, QPushButton, in that order
        """
        hlay = QHBoxLayout()
        le = DandD_FileExplorer2LineEdit()
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

    @staticmethod
    def make_cmb_and_items(items, default=None):
        cmb = QComboBox()
        cmb.addItems(items)
        if default is not None:
            cmb.setCurrentIndex(default)
        return cmb
