from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtGui import *
from xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit
from xASL_GUI_HelperFuncs import set_widget_icon
from pprint import pprint
from glob import glob
import sys
import json
import os
from tdda import rexpy


# noinspection PyAttributeOutsideInit
class xASL_ParmsMaker(QMainWindow):
    def __init__(self, parent_win=None):
        # Parent window is fed into the constructor to allow for communication with parent window devices
        super().__init__(parent=parent_win)
        if parent_win is not None:
            self.config = self.parent().config
        else:
            with open("ExploreASL_GUI_masterconfig.json") as f:
                self.config = json.load(f)

        # Window Size and initial visual setup
        self.setMinimumSize(1080, 480)
        self.cw = QWidget(self)
        self.setCentralWidget(self.cw)
        self.mainlay = QHBoxLayout(self.cw)
        self.setLayout(self.mainlay)
        self.setWindowTitle("Explore ASL - Parameter File Maker")
        self.setWindowIcon(QIcon(os.path.join(os.getcwd(), "media", "ExploreASL_logo.png")))
        # Buttons for executing the fundamental functions
        self.btn_make_parms = QPushButton("Run", self.cw, clicked=self.save_parms2json)
        self.btn_load_parms = QPushButton("Load from existing Json", self.cw, clicked=self.load_json2parms)
        font = QFont()
        font.setPointSize(20)
        self.btn_make_parms.setFont(font)
        self.btn_load_parms.setFont(font)
        # Set up each of the components of the widget
        self.UI_Setup_GroupBoxes()
        self.UI_Setup_Layouts()
        self.UI_Setup_EnvParmsSection()
        self.UI_Setup_StudyParmsSection()
        self.UI_Setup_M0ParmsSection()
        self.UI_Setup_SequenceParmsSection()
        self.UI_Setup_QuantParmsSection()
        self.UI_Setup_ProcParmsSection()
        self.UI_Setup_Connections()
        self.UI_Setup_ToolTips()

    def UI_Setup_GroupBoxes(self):
        # Create groupbox widgets for each of the parameter types
        self.grp_env_parms = QGroupBox("Environment Parameters", self.cw)
        self.grp_study_parms = QGroupBox("Study Parameters", self.cw)
        self.grp_m0_parms = QGroupBox("M0 Parameters", self.cw)
        self.grp_seq_parms = QGroupBox("Sequence Parameters", self.cw)
        self.grp_quant_parms = QGroupBox("Quantification Parameters", self.cw)
        self.grp_proc_parms = QGroupBox("Processing Parameters", self.cw)

    def UI_Setup_Layouts(self):
        # Divide the main layout into 3 sections; left, middle, and right
        self.left_vlay = QVBoxLayout()
        self.mid_vlay = QVBoxLayout()
        self.right_vlay = QVBoxLayout()
        self.mainlay.addLayout(self.left_vlay)
        self.mainlay.addLayout(self.mid_vlay)
        self.mainlay.addLayout(self.right_vlay)
        # Add the appropriate groups to their sections
        self.left_vlay.addWidget(self.grp_env_parms)
        self.left_vlay.addWidget(self.grp_study_parms)
        self.mid_vlay.addWidget(self.grp_seq_parms)
        self.mid_vlay.addWidget(self.grp_quant_parms)
        self.mid_vlay.addWidget(self.grp_m0_parms)
        self.right_vlay.addWidget(self.grp_proc_parms)
        self.right_vlay.addWidget(self.btn_make_parms)
        self.right_vlay.addWidget(self.btn_load_parms)
        # Initialize the form layouts and bind them to their respective groupbox
        self.frmlay_env_parms = QFormLayout(self.grp_env_parms)
        self.frmlay_study_parms = QFormLayout(self.grp_study_parms)
        self.frmlay_m0_parms = QFormLayout(self.grp_m0_parms)
        self.frmlay_seq_parms = QFormLayout(self.grp_seq_parms)
        self.frmlay_quant_parms = QFormLayout(self.grp_quant_parms)
        self.frmlay_proc_parms = QFormLayout(self.grp_proc_parms)

    def UI_Setup_EnvParmsSection(self):
        # Set up the environment parameters group
        self.chk_env_parms_fsl = QCheckBox(self.grp_env_parms, checkable=True, checked=True)
        self.chk_env_parms_mk4dcm = QCheckBox(self.grp_env_parms, checkable=True, checked=False)
        self.frmlay_env_parms.addRow("Detect FSL Automatically?", self.chk_env_parms_fsl)
        self.frmlay_env_parms.addRow("Output CBF native space maps?", self.chk_env_parms_mk4dcm)
        self.collection_env_parms = [self.chk_env_parms_fsl, self.chk_env_parms_mk4dcm]

    def UI_Setup_StudyParmsSection(self):
        # Set up the study parameters group
        self.hlay_study_parms_exploreasl = QHBoxLayout(self.grp_study_parms)
        self.le_study_parms_exploreasl = DandD_FileExplorer2LineEdit(self.grp_study_parms)
        self.le_study_parms_exploreasl.setText(self.config['ExploreASLRoot'])
        self.btn_study_parms_exploreasl = QPushButton("...", self.grp_study_parms,
                                                      clicked=self.set_exploreasl_dir)
        self.hlay_study_parms_exploreasl.addWidget(self.le_study_parms_exploreasl)
        self.hlay_study_parms_exploreasl.addWidget(self.btn_study_parms_exploreasl)
        self.le_study_parms_name = QLineEdit("My Study", self.grp_study_parms)
        self.hlay_study_parms_analysisdir = QHBoxLayout(self.grp_study_parms)
        self.le_study_parms_analysisdir = DandD_FileExplorer2LineEdit(self.grp_study_parms)
        self.le_study_parms_analysisdir.setText(self.parent().le_currentanalysis_dir.text())
        self.btn_study_parms_analysisdir = QPushButton("...", self.grp_study_parms,
                                                       clicked=self.set_currentanalysis_dir)
        self.hlay_study_parms_analysisdir.addWidget(self.le_study_parms_analysisdir, 6)
        self.hlay_study_parms_analysisdir.addWidget(self.btn_study_parms_analysisdir, 1)
        self.le_study_parms_subregexp = QLineEdit(r"\d+", self.grp_study_parms)
        self.lst_study_parms_includedsubjects = DirectoryDragDrop_ListWidget(self.grp_study_parms)
        self.btn_study_parms_clearincluded = QPushButton("Clear Subjects",
                                                         self.grp_study_parms,
                                                         clicked=self.clear_included)
        self.lst_study_parms_exclusions = DirectoryDragDrop_ListWidget(self.grp_study_parms)
        self.btn_study_parms_clearexclusions = QPushButton("Clear Exclusions",
                                                           self.grp_study_parms,
                                                           clicked=self.clear_exclusions)
        self.le_study_parms_sessions = QLineEdit("ASL_1", self.grp_study_parms)
        self.le_study_parms_session_opts = QLineEdit(self.grp_study_parms)
        self.le_study_parms_session_opts.setPlaceholderText("Indicate option names separated by a commas")
        self.frmlay_study_parms.addRow("ExploreASL Directory", self.hlay_study_parms_exploreasl)
        self.frmlay_study_parms.addRow("Name of Study", self.le_study_parms_name)
        self.frmlay_study_parms.addRow("Root Directory", self.hlay_study_parms_analysisdir)
        self.frmlay_study_parms.addRow("Subject Regexp", self.le_study_parms_subregexp)
        self.frmlay_study_parms.addRow("Subjects to Assess\n(Drag and Drop\nDirectories)",
                                       self.lst_study_parms_includedsubjects)
        self.frmlay_study_parms.addRow("", self.btn_study_parms_clearincluded)
        self.frmlay_study_parms.addRow("Exclusions\n(Drag and Drop\nDirectories)", self.lst_study_parms_exclusions)
        self.frmlay_study_parms.addRow("", self.btn_study_parms_clearexclusions)
        self.frmlay_study_parms.addRow("Session Names", self.le_study_parms_sessions)
        self.frmlay_study_parms.addRow("Session Options", self.le_study_parms_session_opts)
        self.collection_study_parms = self.grp_study_parms.findChildren(QLineEdit) + \
                                      self.grp_study_parms.findChildren(QListWidget)

    def UI_Setup_M0ParmsSection(self):
        # Set up the M0 parameters group
        self.cmb_m0_parms_proctype = QComboBox(self.grp_m0_parms)
        self.cmb_m0_parms_proctype.addItems(["New Image Processing", "Standard Processing"])
        self.cmb_m0_parms_proctype.setCurrentIndex(0)
        self.cmb_m0_parms_source = QComboBox(self.grp_m0_parms)
        self.cmb_m0_parms_source.addItems(["M0 exists as separate scan", "Use Mean ASL Control as M0"])
        self.cmb_m0_parms_source.setCurrentIndex(0)
        self.spinbox_m0_parms_scale = QDoubleSpinBox(self.grp_m0_parms)
        self.spinbox_m0_parms_scale.setValue(1)
        self.cmb_m0_parms_pos = QComboBox(self.grp_m0_parms)
        self.cmb_m0_parms_pos.addItems(["M0 is independent of ASL", "First Control-Label Pair",
                                        "First Image", "Second Image"])
        self.cmb_m0_parms_pos.setCurrentIndex(0)
        self.frmlay_m0_parms.addRow("M0 Processing Type", self.cmb_m0_parms_proctype)
        self.frmlay_m0_parms.addRow("M0 Source", self.cmb_m0_parms_source)
        self.frmlay_m0_parms.addRow("GM Scale Factor", self.spinbox_m0_parms_scale)
        self.frmlay_m0_parms.addRow("M0 Position in ASL", self.cmb_m0_parms_pos)
        self.collection_m0_parms = self.grp_m0_parms.findChildren(QComboBox) + \
                                   self.grp_m0_parms.findChildren(QDoubleSpinBox)

    def UI_Setup_SequenceParmsSection(self):
        # Set up the sequence parameters group
        self.cmb_seq_parms_npulses = QComboBox(self.grp_seq_parms)
        self.cmb_seq_parms_npulses.addItems(["0", "2", "4", "5"])
        self.cmb_seq_parms_npulses.setCurrentIndex(1)
        self.cmb_seq_parms_readdim = QComboBox(self.grp_seq_parms)
        self.cmb_seq_parms_readdim.addItems(["2D", "3D"])
        self.cmb_seq_parms_readdim.setCurrentIndex(1)
        self.cmb_seq_parms_vendor = QComboBox(self.grp_seq_parms)
        self.cmb_seq_parms_vendor.addItems(["GE Product", "GE WIP", "Philips", "Siemens"])
        self.cmb_seq_parms_vendor.setCurrentIndex(3)
        self.cmb_seq_parms_seqtype = QComboBox(self.grp_seq_parms)
        self.cmb_seq_parms_seqtype.addItems(["3D Spiral", "3D GRaSE", "2D EPI"])
        self.cmb_seq_parms_seqtype.setCurrentIndex(1)
        self.cmb_seq_parms_labeltype = QComboBox(self.grp_seq_parms)
        self.cmb_seq_parms_labeltype.addItems(["Q2 TIPS PASL", "pCASL or CASL"])
        self.cmb_seq_parms_labeltype.setCurrentIndex(0)
        self.spinbox_seq_parms_iniPLD = QDoubleSpinBox(self.grp_seq_parms)
        self.spinbox_seq_parms_iniPLD.setRange(0, 5000)
        self.spinbox_seq_parms_iniPLD.setValue(1800)
        self.spinbox_seq_parms_labdur = QDoubleSpinBox(self.grp_seq_parms)
        self.spinbox_seq_parms_labdur.setRange(0, 5000)
        self.spinbox_seq_parms_labdur.setValue(800)
        self.cmb_seq_parms_slicerdt = QComboBox(self.grp_seq_parms)
        self.cmb_seq_parms_slicerdt.addItems(["Use Shortest TR", "Use Indicated Value"])
        self.cmb_seq_parms_slicerdt.setCurrentIndex(1)
        self.cmb_seq_parms_slicerdt.currentTextChanged.connect(self.on_off_slicereadtime_box)
        self.spinbox_seq_parms_slicerdtime = QDoubleSpinBox(self.grp_seq_parms)
        self.spinbox_seq_parms_slicerdtime.setValue(37)
        self.hbox_seq_parms_slicerdtime = QHBoxLayout(self.grp_seq_parms)
        self.hbox_seq_parms_slicerdtime.addWidget(self.cmb_seq_parms_slicerdt)
        self.hbox_seq_parms_slicerdtime.addWidget(self.spinbox_seq_parms_slicerdtime)
        self.frmlay_seq_parms.addRow("Number of Suppression Pulses", self.cmb_seq_parms_npulses)
        self.frmlay_seq_parms.addRow("Readout Dimension", self.cmb_seq_parms_readdim)
        self.frmlay_seq_parms.addRow("Vendor", self.cmb_seq_parms_vendor)
        self.frmlay_seq_parms.addRow("Sequence Type", self.cmb_seq_parms_seqtype)
        self.frmlay_seq_parms.addRow("Labeling Type", self.cmb_seq_parms_labeltype)
        self.frmlay_seq_parms.addRow("Initial Post-Labeling Delay", self.spinbox_seq_parms_iniPLD)
        self.frmlay_seq_parms.addRow("Labeling Duration", self.spinbox_seq_parms_labdur)
        self.frmlay_seq_parms.addRow("Slice Readout Time", self.hbox_seq_parms_slicerdtime)
        self.collection_seq_parms = [self.cmb_seq_parms_npulses, self.cmb_seq_parms_readdim, self.cmb_seq_parms_vendor,
                                     self.cmb_seq_parms_seqtype, self.cmb_seq_parms_labeltype,
                                     self.spinbox_seq_parms_iniPLD, self.spinbox_seq_parms_labdur,
                                     self.cmb_seq_parms_slicerdt, self.spinbox_seq_parms_slicerdtime]

    def UI_Setup_QuantParmsSection(self):
        # Set up Quantification Parameters
        self.spinbox_quant_parms_lambda = QDoubleSpinBox(self.grp_quant_parms)
        self.spinbox_quant_parms_lambda.setRange(0, 100)
        self.spinbox_quant_parms_lambda.setSingleStep(0.01)
        self.spinbox_quant_parms_lambda.setValue(0.9)
        self.spinbox_quant_parms_artT2 = QDoubleSpinBox(self.grp_quant_parms)
        self.spinbox_quant_parms_artT2.setRange(0, 5000)
        self.spinbox_quant_parms_artT2.setValue(50)
        self.spinbox_quant_parms_bloodT1 = QDoubleSpinBox(self.grp_quant_parms)
        self.spinbox_quant_parms_bloodT1.setRange(0, 5000)
        self.spinbox_quant_parms_bloodT1.setValue(1650)
        self.spinbox_quant_parms_tissueT1 = QDoubleSpinBox(self.grp_quant_parms)
        self.spinbox_quant_parms_tissueT1.setRange(0, 5000)
        self.spinbox_quant_parms_tissueT1.setValue(1240)
        self.cmb_quant_parms_ncomparts = QComboBox(self.grp_quant_parms)
        self.cmb_quant_parms_ncomparts.addItems(["1", "2"])
        self.cmb_quant_parms_ncomparts.setCurrentIndex(0)
        self.le_quant_parms_logvec = QLineEdit(self.grp_quant_parms)
        self.le_quant_parms_logvec.setText("1 1 1 1 1")
        self.le_quant_parms_logvec.setPlaceholderText("Specify a vector of 0s or 1s separated by single spaces")
        self.frmlay_quant_parms.addRow("Lambda", self.spinbox_quant_parms_lambda)
        self.frmlay_quant_parms.addRow("Arterial T2*", self.spinbox_quant_parms_artT2)
        self.frmlay_quant_parms.addRow("Blood T1", self.spinbox_quant_parms_bloodT1)
        self.frmlay_quant_parms.addRow("Tissue T1", self.spinbox_quant_parms_tissueT1)
        self.frmlay_quant_parms.addRow("Number of Compartments", self.cmb_quant_parms_ncomparts)
        self.frmlay_quant_parms.addRow("Quantification Settings", self.le_quant_parms_logvec)
        self.collection_quant_parms = [self.spinbox_quant_parms_lambda, self.spinbox_quant_parms_artT2,
                                       self.spinbox_quant_parms_bloodT1, self.spinbox_quant_parms_tissueT1,
                                       self.cmb_quant_parms_ncomparts, self.le_quant_parms_logvec]

    def UI_Setup_ProcParmsSection(self):
        # Set up Processing Parameters
        self.chk_proc_parms_spikes = QCheckBox(self.grp_proc_parms, checkable=True, checked=True)
        self.spinbox_proc_parms_spikethres = QDoubleSpinBox(self.grp_proc_parms)
        self.spinbox_proc_parms_spikethres.setValue(0.01)
        self.spinbox_proc_parms_spikethres.setSingleStep(0.001)
        self.chk_proc_parms_motion = QCheckBox(self.grp_proc_parms, checkable=True, checked=True)
        self.cmb_proc_parms_quality = QComboBox(self.grp_proc_parms)
        self.cmb_proc_parms_quality.addItems(["Low", "High"])
        self.cmb_proc_parms_quality.setCurrentIndex(1)
        self.chk_proc_parms_deltemp = QCheckBox(self.grp_proc_parms, checkable=True, checked=True)
        self.chk_proc_parms_skipnoflair = QCheckBox(self.grp_proc_parms, checkable=True, checked=False)
        self.chk_proc_parms_skipnoasl = QCheckBox(self.grp_proc_parms, checkable=True, checked=True)
        self.chk_proc_parms_skipnom0 = QCheckBox(self.grp_proc_parms, checkable=True, checked=False)
        self.chk_proc_parms_T1dartel = QCheckBox(self.grp_proc_parms, checkable=True, checked=True)
        self.chk_proc_parms_PWIdartel = QCheckBox(self.grp_proc_parms, checkable=True, checked=False)
        self.chk_proc_parms_bilatfilter = QCheckBox(self.grp_proc_parms, checkable=True, checked=False)
        self.cmb_proc_parms_segtype = QComboBox(self.grp_proc_parms)
        self.cmb_proc_parms_segtype.addItems(["CAT12", "SPM12"])
        self.cmb_proc_parms_segtype.setCurrentIndex(0)
        self.cmb_proc_parms_imgcontrastsrc = QComboBox(self.grp_proc_parms)
        self.cmb_proc_parms_imgcontrastsrc.addItems(["Control --> T1w",
                                                     "CBF --> pseudoCBF",
                                                     "Automatic",
                                                     "Forced CBF --> pseudoCBF"])
        self.cmb_proc_parms_imgcontrastsrc.setCurrentIndex(2)
        self.cmb_proc_parms_affinereg = QComboBox(self.grp_proc_parms)
        self.cmb_proc_parms_affinereg.addItems(["Disabled", "Enabled", "Based on PWI CoV"])
        self.cmb_proc_parms_affinereg.setCurrentIndex(2)
        self.chk_proc_parms_m02aslreg = QCheckBox(self.grp_proc_parms, checkable=True, checked=True)
        self.chk_proc_parms_useMNIasdummy = QCheckBox(self.grp_proc_parms, checkable=True, checked=False)
        self.frmlay_proc_parms.addRow("Remove Spikes", self.chk_proc_parms_spikes)
        self.frmlay_proc_parms.addRow("Spike Removal Threshold", self.spinbox_proc_parms_spikethres)
        self.frmlay_proc_parms.addRow("Correct for Motion", self.chk_proc_parms_motion)
        self.frmlay_proc_parms.addRow("Quality", self.cmb_proc_parms_quality)
        self.frmlay_proc_parms.addRow("Delete Temporary Files", self.chk_proc_parms_deltemp)
        self.frmlay_proc_parms.addRow("Skip Subjects without FLAIR", self.chk_proc_parms_skipnoflair)
        self.frmlay_proc_parms.addRow("Skip Subjects without ASL", self.chk_proc_parms_skipnoasl)
        self.frmlay_proc_parms.addRow("Skip Subjects without M0", self.chk_proc_parms_skipnom0)
        self.frmlay_proc_parms.addRow("Use T1 DARTEL", self.chk_proc_parms_T1dartel)
        self.frmlay_proc_parms.addRow("Use PWI DARTEL", self.chk_proc_parms_PWIdartel)
        self.frmlay_proc_parms.addRow("Use Bilateral Filter", self.chk_proc_parms_bilatfilter)
        self.frmlay_proc_parms.addRow("Segmentation Method", self.cmb_proc_parms_segtype)
        self.frmlay_proc_parms.addRow("Image Contrast used for", self.cmb_proc_parms_imgcontrastsrc)
        self.frmlay_proc_parms.addRow("Use Affine Registration", self.cmb_proc_parms_affinereg)
        self.frmlay_proc_parms.addRow("Register M0 to ASL", self.chk_proc_parms_m02aslreg)
        self.frmlay_proc_parms.addRow("Use MNI as a Dummy Template", self.chk_proc_parms_useMNIasdummy)
        self.collection_proc_parms = self.grp_proc_parms.findChildren(QCheckBox) + \
                                     self.grp_proc_parms.findChildren(QComboBox) + \
                                     self.grp_proc_parms.findChildren(QDoubleSpinBox)

    # Setup of all inner signals here that are not included within constructors (ex. QPushButton 'clicked' arguments)
    def UI_Setup_Connections(self):
        self.lst_study_parms_includedsubjects.alert_regex.connect(self.get_regex_from_subjectlist)

    # Setup hover tooltips over all relevant widgets
    def UI_Setup_ToolTips(self):
        for widget_group, section in zip([self.collection_env_parms, self.collection_study_parms,
                                          self.collection_m0_parms, self.collection_seq_parms,
                                          self.collection_quant_parms, self.collection_proc_parms],
                                         ["EnvParms", "StudyParms", "M0Parms", "SeqParms", "QuantParms", "ProcParms"]):
            for widget, tip in zip(widget_group, self.parent().tooltips["ParmsMaker"][section]):
                widget.setToolTip(tip)

    # Clears the currently-included subjects list and resets the regex
    def clear_included(self):
        self.lst_study_parms_includedsubjects.clear()
        self.le_study_parms_subregexp.clear()

    # Clears the currently-excluded subjects list
    def clear_exclusions(self):
        self.lst_study_parms_exclusions.clear()

    # Retrieve the directory ExploreASL with user interaction; set the appropriate field and save to config file
    def set_exploreasl_dir(self):
        result = QFileDialog.getExistingDirectory(QFileDialog(),
                                                  "Select the ExploreASL root directory",
                                                  os.getcwd(),
                                                  QFileDialog.ShowDirsOnly)
        if '/' in result and "windows" in self.parent().config["Platform"].lower():
            result = result.replace('/', '\\')
        if result:
            self.le_study_parms_exploreasl.setText(result)

    # Retrieve the analysis directory with user interaction; set the appropriate field
    def set_currentanalysis_dir(self):
        result = QFileDialog.getExistingDirectory(QFileDialog(),
                                                  "Select the study analysis directory",
                                                  os.getcwd(),
                                                  QFileDialog.ShowDirsOnly)
        if '/' in result and "windows" in self.parent().config["Platform"].lower():
            result = result.replace('/', '\\')
        if result:
            self.le_study_parms_analysisdir.setText(result)

    # Updates the regex based on the text patterns placed each time the content of includedsubjects changes
    def get_regex_from_subjectlist(self):
        n_subjects = self.lst_study_parms_includedsubjects.count()
        if n_subjects == 0:
            return ""
        subject_list = [self.lst_study_parms_includedsubjects.item(idx).text() for idx in range(n_subjects)]
        extractor = rexpy.Extractor(subject_list)
        extractor.batch_extract(subject_list)
        inferred_regex = extractor.results.rex[0]
        del extractor
        if inferred_regex:
            self.le_study_parms_subregexp.setText(inferred_regex)

    def on_off_slicereadtime_box(self, text):
        if text == "Use Shortest TR":
            self.spinbox_seq_parms_slicerdtime.setEnabled(False)
        elif text == "Use Indicated Value":
            self.spinbox_seq_parms_slicerdtime.setEnabled(True)

    def save_parms2json(self):
        json_parms = {}
        # Only export if the analysis study directory is known
        if not os.path.exists(self.le_study_parms_analysisdir.text()):
            QMessageBox.warning(QMessageBox(),
                                "The Study Root does not exist",
                                "Please select the correct directory",
                                QMessageBox.Ok)
            return

        json_parms["MyPath"] = self.le_study_parms_exploreasl.text()
        json_parms["name"] = self.le_study_parms_name.text()
        json_parms["D"] = {}
        json_parms["D"]["ROOT"] = self.le_study_parms_analysisdir.text()
        json_parms["subject_regexp"] = self.le_study_parms_subregexp.text()
        parms_sessions = self.le_study_parms_sessions.text()
        if "," in parms_sessions:
            parms_sessions = [sess.strip() for sess in parms_sessions.split(",")]
        else:
            parms_sessions = [parms_sessions]
        json_parms["SESSIONS"] = parms_sessions
        json_parms["session"] = {}
        parms_sessions_opts = self.le_study_parms_session_opts.text()
        if "," in parms_sessions_opts:
            parms_sessions_opts = [opt.strip() for opt in parms_sessions_opts.split(",")]
        else:
            parms_sessions_opts = [parms_sessions_opts]
        json_parms["session"]["options"] = parms_sessions_opts
        exclusions = [self.lst_study_parms_exclusions.item(row).text() for row in
                      range(self.lst_study_parms_exclusions.count())]
        json_parms["exclusion"] = exclusions
        json_parms["M0_conventionalProcessing"] = self.cmb_m0_parms_proctype.currentIndex()
        parms_source_translate = {"M0 exists as separate scan": "separate_scan",
                                  "Use Mean ASL Control as M0": "UseControlAsM0"}
        json_parms["M0"] = parms_source_translate[self.cmb_m0_parms_source.currentText()]
        json_parms["M0_GMScaleFactor"] = float(self.spinbox_m0_parms_scale.value())
        parms_m0_pos_translate = {"First Control-Label Pair": "[1 2]", "First Image": 1, "Second Image": 2}
        if self.cmb_m0_parms_pos.currentText() != "M0 is independent of ASL":
            json_parms["M0PositionInASL4D"] = parms_m0_pos_translate[self.cmb_m0_parms_pos.currentText()]
        json_parms["readout_dim"] = self.cmb_seq_parms_readdim.currentText()
        parms_vendor_translate = {"GE Product": "GE_Product", "GE WIP": "GE_WIP", "Philips": "Philips",
                                  "Siemens": "Siemens"}
        json_parms["Vendor"] = parms_vendor_translate[self.cmb_seq_parms_vendor.currentText()]
        parms_sequence_translate = {"3D Spiral": "3D_Spiral", "3D GRaSE": "3D_GRASE", "2D EPI": "2D_EPI"}
        json_parms["Sequence"] = parms_sequence_translate[self.cmb_seq_parms_seqtype.currentText()]
        parms_quality_translate = {"Low": 0, "High": 1}
        json_parms["Quality"] = parms_quality_translate[self.cmb_proc_parms_quality.currentText()]
        json_parms["DELETETEMP"] = int(self.chk_proc_parms_deltemp.isChecked())
        json_parms["SPIKE_REMOVAL"] = int(self.chk_proc_parms_spikes.isChecked())
        json_parms["SpikeRemovalThreshold"] = float(self.spinbox_proc_parms_spikethres.value())
        json_parms["SkipIfNoFlair"] = int(self.chk_proc_parms_skipnoflair.isChecked())
        json_parms["SkipIfNoM0"] = int(self.chk_proc_parms_skipnom0.isChecked())
        json_parms["SkipIfNoASL"] = int(self.chk_proc_parms_skipnoasl.isChecked())
        json_parms["T1_DARTEL"] = int(self.chk_proc_parms_T1dartel.isChecked())
        json_parms["PWI_DARTEL"] = int(self.chk_proc_parms_PWIdartel.isChecked())
        json_parms["BILAT_FILTER"] = int(self.chk_proc_parms_bilatfilter.isChecked())
        json_parms["motion_correction"] = int(self.chk_proc_parms_motion.isChecked())
        parms_segtype_translate = {"SPM12": 1, "CAT12": 0}
        json_parms["SegmentSPM12"] = parms_segtype_translate[self.cmb_proc_parms_segtype.currentText()]
        json_parms["bRegistrationContrast"] = int(self.cmb_proc_parms_imgcontrastsrc.currentIndex())
        json_parms["bAffineRegistration"] = int(self.cmb_proc_parms_affinereg.currentIndex())
        json_parms["bRegisterM02ASL"] = int(self.chk_proc_parms_m02aslreg.isChecked())
        json_parms["bUseMNIasDummyStructural"] = int(self.chk_proc_parms_useMNIasdummy.isChecked())
        # Handling the Quantification logical vector
        parms_logvec = self.le_quant_parms_logvec.text().split(" ")
        if all([len(parms_logvec) == 5,  # Must be 5 1s or 0s
                all([x in ['1', '0'] for x in parms_logvec])]):  # Check that all are 1s or 0s
            json_parms["ApplyQuantification"] = [int(option) for option in parms_logvec]
        else:
            QMessageBox.warning(QMessageBox(),
                                "Incorrect Input for Quantification Settings",
                                "Must be a series of five 1s or 0s separated by single spaces",
                                QMessageBox.Ok)
            return
        json_parms["Q"] = {}
        json_parms["Q"]["BackGrSupprPulses"] = int(self.cmb_seq_parms_npulses.currentText())
        parms_labeltype_translate = {"Q2 TIPS PASL": "PASL", "pCASL or CASL": "CASL"}
        json_parms["Q"]["LabelingType"] = parms_labeltype_translate[self.cmb_seq_parms_labeltype.currentText()]
        json_parms["Q"]["Initial_PLD"] = float(self.spinbox_seq_parms_iniPLD.value())
        json_parms["Q"]["LabelingDuration"] = float(self.spinbox_seq_parms_labdur.value())
        parms_slicereadtime_translate = {"Use Shortest TR": "shortestTR",
                                         "Use Indicated Value": float(self.spinbox_seq_parms_slicerdtime.value())}
        json_parms["Q"]["SliceReadoutTime"] = parms_slicereadtime_translate[self.cmb_seq_parms_slicerdt.currentText()]
        json_parms["Q"]["Lambda"] = float(self.spinbox_quant_parms_lambda.value())
        json_parms["Q"]["T2art"] = float(self.spinbox_quant_parms_artT2.value())
        json_parms["Q"]["TissueT1"] = float(self.spinbox_quant_parms_tissueT1.value())
        json_parms["Q"]["nCompartments"] = int(self.cmb_quant_parms_ncomparts.currentText())

        with open(os.path.join(self.le_study_parms_analysisdir.text(), "DataPar.json"), 'w') as w:
            json.dump(json_parms, w, indent=1)
        QMessageBox.information(QMessageBox(),
                                "DataPar.json successfully saved",
                                f"The parameter file was successfully saved to:\n"
                                f"{self.le_study_parms_analysisdir.text()}",
                                QMessageBox.Ok)

    def load_json2parms(self):
        errors_log = {}
        json_filepath, _ = QFileDialog.getOpenFileName(QFileDialog(),
                                                       "Select the JSON parameters file",
                                                       os.getcwd(),
                                                       "Json files (*.json)")
        if json_filepath == '': return
        if not os.path.exists(json_filepath) and json_filepath != '':
            QMessageBox.warning(QMessageBox(),
                                "Incorrect File Selected",
                                "The file selected was either not a json file "
                                "or did not contain the essential parameters",
                                QMessageBox.Ok)
            return

        with open(json_filepath, 'r') as reader:
            parms: dict = json.load(reader)

        for key, value in parms.items():
            try:
                if key == "MyPath":
                    self.le_study_parms_exploreasl.setText(value)
                elif key == "name":
                    self.le_study_parms_name.setText(value)
                elif key == "D":
                    self.le_study_parms_analysisdir.setText(value["ROOT"])
                elif key == "subject_regexp":
                    self.le_study_parms_subregexp.setText(value)
                elif key == "SESSIONS":
                    self.le_study_parms_sessions.setText(", ".join([item for item in value]))
                elif key == "session":
                    if value["options"] == [""]:
                        self.le_study_parms_session_opts.setText("")
                    else:
                        self.le_study_parms_session_opts.setText(", ".join([item for item in value["options"]]))
                elif key == "exclusion":
                    self.lst_study_parms_exclusions.clear()
                    self.lst_study_parms_exclusions.addItems([item for item in value])
                elif key == "M0_conventionalProcessing":
                    self.cmb_m0_parms_proctype.setCurrentIndex(value)
                elif key == "M0":
                    parms_source_translate = {"separate_scan": "M0 exists as separate scan",
                                              "UseControlAsM0": "Use Mean ASL Control as M0"}
                    self.cmb_m0_parms_proctype.setCurrentIndex(
                        self.cmb_m0_parms_proctype.findText(parms_source_translate[value]))
                elif key == "M0_GMScaleFactor":
                    self.spinbox_m0_parms_scale.setValue(value)
                elif key == "M0PositionInASL4D":
                    parms_m0_pos_translate = {"[1 2]": "First Control-Label Pair", 1: "First Image", 2: "Second Image"}
                    self.cmb_m0_parms_pos.setCurrentIndex(self.cmb_m0_parms_pos.findText(parms_m0_pos_translate[value]))
                elif key == "readout_dim":
                    self.cmb_seq_parms_readdim.setCurrentIndex(self.cmb_seq_parms_readdim.findText(value))
                elif key == "Vendor":
                    parms_vendor_translate = {"GE_Product": "GE Product", "GE_WIP": "GE WIP", "Philips": "Philips",
                                              "Siemens": "Siemens"}
                    self.cmb_seq_parms_vendor.setCurrentIndex(
                        self.cmb_seq_parms_vendor.findText(parms_vendor_translate[value]))
                elif key == "Sequence":
                    parms_sequence_translate = {"3D_Spiral": "3D Spiral", "3D_GRaSE": "3D GRASE", "2D_EPI": "2D EPI"}
                    self.cmb_seq_parms_seqtype.setCurrentIndex(
                        self.cmb_seq_parms_seqtype.findText(parms_sequence_translate[value]))
                elif key == "Quality":
                    parms_quality_translate = {0: "Low", 1: "High"}
                    self.cmb_proc_parms_quality.setCurrentIndex(
                        self.cmb_proc_parms_quality.findText(parms_quality_translate[value]))
                elif key == "DELETETEMP":
                    self.chk_proc_parms_deltemp.setChecked(bool(value))
                elif key == "SPIKE_REMOVAL":
                    self.chk_proc_parms_spikes.setChecked(bool(value))
                elif key == "SpikeRemovalThreshold":
                    self.spinbox_proc_parms_spikethres.setValue(value)
                elif key == "SkipIfNoFlair":
                    self.chk_proc_parms_skipnoflair.setChecked(bool(value))
                elif key == "SkipIfNoM0":
                    self.chk_proc_parms_skipnom0.setChecked(bool(value))
                elif key == "SkipIfNoASL":
                    self.chk_proc_parms_skipnoasl.setChecked(bool(value))
                elif key == "T1_DARTEL":
                    self.chk_proc_parms_T1dartel.setChecked(bool(value))
                elif key == "PWI_DARTEL":
                    self.chk_proc_parms_PWIdartel.setChecked(bool(value))
                elif key == "BILAT_FILTER":
                    self.chk_proc_parms_bilatfilter.setChecked(bool(value))
                elif key == "motion_correction":
                    self.chk_proc_parms_motion.setChecked(bool(value))
                elif key == "SegmentSPM12":
                    parms_segtype_translate = {1: "SPM12", 0: "CAT12"}
                    self.cmb_proc_parms_segtype.setCurrentIndex(
                        self.cmb_proc_parms_segtype.findText(parms_segtype_translate[value]))
                elif key == "bRegistrationContrast":
                    self.cmb_proc_parms_imgcontrastsrc.setCurrentIndex(value)
                elif key == "bAffineRegistration":
                    self.cmb_proc_parms_affinereg.setCurrentIndex(value)
                elif key == "bRegisterM02ASL":
                    self.chk_proc_parms_m02aslreg.setChecked(bool(value))
                elif key == "bUseMNIasDummyStructural":
                    self.chk_proc_parms_useMNIasdummy.setChecked(bool(value))
                elif key == "ApplyQuantification":
                    self.le_quant_parms_logvec.setText(" ".join(str(num) for num in value))
                elif key == "Q":
                    for sub_key, sub_val in value.items():
                        try:
                            if sub_key == "BackGrSupprPulses":
                                self.cmb_seq_parms_npulses.setCurrentIndex(
                                    self.cmb_seq_parms_npulses.findText(str(sub_val)))
                            elif sub_key == "LabelingType":
                                parms_labeltype_translate = {"PASL": "Q2 TIPS PASL", "CASL": "pCASL or CASL"}
                                self.cmb_seq_parms_labeltype.setCurrentIndex(
                                    self.cmb_seq_parms_labeltype.findText(parms_labeltype_translate[sub_val]))
                            elif sub_key == "Initial_PLD":
                                self.spinbox_seq_parms_iniPLD.setValue(sub_val)
                            elif sub_key == "LabelingDuration":
                                self.spinbox_seq_parms_labdur.setValue(sub_val)
                            elif sub_key == "SliceReadoutTime":
                                self.spinbox_seq_parms_slicerdtime.setValue(sub_val)
                            elif sub_key == "Lambda":
                                self.spinbox_quant_parms_lambda.setValue(sub_val)
                            elif sub_key == "T2art":
                                self.spinbox_quant_parms_artT2.setValue(sub_val)
                            elif sub_key == "TissueT1":
                                self.spinbox_quant_parms_tissueT1.setValue(sub_val)
                            elif sub_key == "nCompartments":
                                self.cmb_quant_parms_ncomparts.setCurrentIndex(sub_val - 1)
                            else:
                                print(f"{sub_key} is not a valid quantification field")
                        except:
                            errors_log[sub_key] = sub_val
                else:
                    print(f"{key} is not a valid field")
            except:
                errors_log[key] = value

        if len(errors_log) > 0:
            warning_string = "\n".join([f"{key}: {value}" for key, value in errors_log.items()])
            QMessageBox.warning(QMessageBox(),
                                "Errors Occured at certain fields",
                                f"The following parameters did not load correctly:\n{warning_string}\n"
                                f"Please ascertain the validity of the values and/or labels for these parameters",
                                QMessageBox.Ok)

    # This will be utilized in a future setting to control the navigator movement and prevent its deletion if it is
    # docked within the
    def closeEvent(self, a0) -> None:
        super(xASL_ParmsMaker, self).closeEvent(a0)


class DirectoryDragDrop_ListWidget(QListWidget):
    """
    Class meant to accept MULTIPLE directory inputs and add them to the underlying QListWidget
    """
    alert_regex = Signal()  # This will be called in the ParmsMaker superclass to update its regex

    def __init__(self, parent=None):
        super(DirectoryDragDrop_ListWidget, self).__init__(parent)
        self.setAcceptDrops(True)

    def dragEnterEvent(self, event) -> None:
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event: QDragMoveEvent) -> None:
        if event.mimeData().hasUrls():
            event.setDropAction(Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event: QDropEvent) -> None:
        if event.mimeData().hasUrls():
            event.accept()
            subject_directories = []
            for url in event.mimeData().urls():
                if url.isLocalFile():
                    subject_directories.append(str(url.toLocalFile()))
            # Only filter for the basenames of directories; also, avoid bringing in unnecessary directories like lock
            basenames = [os.path.basename(directory) for directory in subject_directories if os.path.isdir(directory)
                         and directory not in ['lock', "Population"]]
            # Avoid double-dipping the names
            current_names = [self.item(idx).text() for idx in range(self.count())]
            filtered_basenames = [name for name in basenames if name not in current_names]
            self.addItems(filtered_basenames)
            self.alert_regex.emit()
        else:
            event.ignore()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    maker = xASL_ParmsMaker()
    maker.show()
    sys.exit(app.exec())
