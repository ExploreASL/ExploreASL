import json
import os
import re
from glob import glob, iglob
from itertools import chain
from time import sleep
from datetime import date
from more_itertools import peekable
from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit
from ExploreASL_GUI.xASL_GUI_Executor_ancillary import (initialize_all_lock_dirs, calculate_anticipated_workload,
                                                        calculate_missing_STATUS, interpret_statusfile_errors)
from ExploreASL_GUI.xASL_GUI_AnimationClasses import xASL_ImagePlayer
from ExploreASL_GUI.xASL_GUI_Executor_Modjobs import xASL_GUI_RerunPrep, xASL_GUI_TSValter
from ExploreASL_GUI.xASL_GUI_HelperFuncs import (set_widget_icon, make_droppable_clearable_le)
from ExploreASL_GUI.xASL_GUI_HelperFuncs_StringOps import set_os_dependent_text
from pprint import pprint
from collections import defaultdict
import subprocess
from shutil import rmtree


class ExploreASL_WorkerSignals(QObject):
    """
    Class for handling the signals sent by an ExploreASL worker
    """
    started_processing = Signal()
    finished_processing = Signal()  # Signal sent by worker to watcher to help indicate it when to stop watching
    stdout_processing_error = Signal(dict)  # Signal sent by a worker to update the main stdout errors dict
    stderr_processing_error = Signal(dict)  # Signal sent by a worker to update the main stderr errors dict
    encountered_fatal_error = Signal()  # Signal sent by worker to raise a dialogue informing the user of an error


class ExploreASL_Worker(QRunnable):
    """
    Worker thread for running lauching an ExploreASL MATLAB session with the given arguments
    """

    def __init__(self, *args):
        self.args = args
        super().__init__()
        self.signals = ExploreASL_WorkerSignals()
        self.is_running = False
        print(f"Initialized Worker with args {self.args}")

    # This is called by the threadpool during threadpool.start(worker)
    # noinspection RegExpRedundantEscape
    def run(self):
        ##################################################
        # PREPARE ARGUMENTS AND RUN THE UNDERLYING PROGRAM
        ##################################################

        exploreasl_path, par_path, process_data, skip_pause, iworker, nworkers, imodules = self.args[0:-1]
        subject_regex_str = self.args[-1]

        # Generate the string that the command line will feed into the complied MATLAB session
        func_line = f"('{par_path}', " \
                    f"{process_data}, " \
                    f"{skip_pause}, " \
                    f"{iworker}, " \
                    f"{nworkers}, " \
                    f"[{' '.join([str(item) for item in imodules])}])"

        cmds = ["matlab",
                "-nodesktop",
                "-nosplash",
                "-batch",
                f"cd('{exploreasl_path}'); ExploreASL_Master{func_line}"]

        # # Compile and run MATLAB session from command line
        # result = subprocess.run(cmds, capture_output=True, text=True)
        # stdout = result.stdout
        # stderr = result.stderr
        # print(f"RESULTs for IWorker {iworker} of {nworkers}:"
        #       f"\n\tReturn code = {result.returncode}")

        self.proc = subprocess.Popen(cmds, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.is_running = True

        # Must keep it in a loop to prevent code from proceeding until the process is definitely done for purposely
        # terminated
        # while self.proc.poll() is None:
        #     print(f"Worker {iworker} of {nworkers} for study {par_path} is still running.")
        #     print(f"{self.proc.poll()=}")
        #     sleep(10)

        stdout, stderr = self.proc.communicate()
        print(f"RESULTs for IWorker {iworker} of {nworkers}:"
              f"\n\tReturn code = {self.proc.returncode}")
        # print(f"{stdout=}")
        # print(f"{stderr=}")
        self.is_running = False

        worker_analysis_dir = re.search(r'.*analysis', par_path).group()
        stdout_error_dict = {}
        stderr_error_dict = {}

        ####################################################
        # ATTEMPT TO CATCH STDOUT ERRORS FOR TROUBLESHOOTING
        ####################################################
        if imodules == [3]:  # Population module
            stdout_error_regex = re.compile(r"ERROR: Job iteration terminated![\n\s]+ans ="
                                            r"([\n\s\w'\/()|;,><=%^&*$#!@~`:.\[\]]+)"
                                            r"CONT: but continue with next iteration")
            captured_stdout_errors = stdout_error_regex.findall(stdout)
            if len(captured_stdout_errors) == 1:
                stdout_error_dict[worker_analysis_dir] = {"Population": captured_stdout_errors[0]}
                self.signals.stdout_processing_error.emit(stdout_error_dict)

        else:  # Anything other than the Population module
            stdout_error_regex = re.compile(
                f"({subject_regex_str})" +
                r"(?:.*)"
                r"ERROR: Job iteration terminated![\n\s]+ans ="
                r"([\n\s\w'\/()|;,><=%^&*$#!@~`:.\[\]]+)"
                r"CONT: but continue with next iteration", re.DOTALL)

            captured_stdout_errors = stdout_error_regex.findall(stdout)
            if len(captured_stdout_errors) > 0:
                captured_stdout_errors = dict(captured_stdout_errors)
                stdout_error_dict[worker_analysis_dir] = captured_stdout_errors
                self.signals.stdout_processing_error.emit(stdout_error_dict)

        #################################
        # SEND THE APPROPRIATE END SIGNAL
        #################################
        # if result.returncode == 0:
        if self.proc.returncode == 0:
            print(f"WORKER {iworker} IS EMITTING FINISHED PROCESSING SIGNAL")
            self.signals.finished_processing.emit()
        else:
            print(f"WORKER {iworker} IS EMMITING ENCOUNTERED ERROR SIGNAL")
            stderr_error_dict[re.search(r'.*analysis', par_path).group()] = stderr
            self.signals.stderr_processing_error.emit(stderr_error_dict)
            self.signals.encountered_fatal_error.emit()

    @Slot()
    def terminate_run(self):
        print("RECEIVED TERMINATE CLICK. CHECKING FOR WHETHER THE WORKER IS STILL RUNNING!")
        if self.is_running:
            if self.proc.poll() is None:
                print("THE WORKER WAS FOUND TO STILL BE RUNNING. ATTEMPTING TO TERMINATE THE PROCESS NOW!")
                self.proc.terminate()


# noinspection PyCallingNonCallable,PyAttributeOutsideInit,PyCallByClass
class xASL_Executor(QMainWindow):
    cont_nstudies: QWidget

    def __init__(self, parent_win):
        # Parent window is fed into the constructor to allow for communication with parent window devices
        super().__init__(parent=parent_win)
        self.config = self.parent().config

        # Window Size and initial visual setup
        self.setMinimumSize(self.config["ScreenSize"][0] // 2, self.config["ScreenSize"][1] // 2)
        self.resize(self.config["ScreenSize"][0] // 2, self.config["ScreenSize"][1] // 2)
        self.cw = QWidget(self)
        self.setCentralWidget(self.cw)
        self.mainlay = QHBoxLayout(self.cw)
        self.setLayout(self.mainlay)
        self.setWindowTitle("Explore ASL - Executor")
        self.setWindowIcon(QIcon(os.path.join(self.config["ProjectDir"], "media", "ExploreASL_logo.png")))

        # Other instance variables
        self.threadpool = QThreadPool()
        self.processing_movieplayer = xASL_ImagePlayer(os.path.join(self.config["ProjectDir"],
                                                                    "media",
                                                                    "processing.gif"))

        # MISC VARIABLES
        self.red_palette = QPalette()
        self.red_palette.setColor(QPalette.Highlight, Qt.red)
        self.green_palette = QPalette()
        self.green_palette.setColor(QPalette.Highlight, Qt.green)
        self.total_process_dbt = 0
        with open(os.path.join(self.config["ProjectDir"],
                               "JSON_LOGIC",
                               "ExecutorTranslators.json")) as translator_reader:
            self.executor_translators = json.load(translator_reader)

        self.UI_Setup_Layouts_and_Groups()
        self.UI_Setup_TaskScheduler()
        self.UI_Setup_TextFeedback_and_Executor()
        self.UI_Setup_ProcessModification()

    def UI_Setup_Layouts_and_Groups(self):
        self.splitter_leftside = QSplitter(Qt.Vertical, self.cw)
        self.splitter_rightside = QSplitter(Qt.Vertical, self.cw)

        # Group Boxes
        self.grp_taskschedule = QGroupBox(title="Task Scheduler")
        self.grp_textoutput = QGroupBox(title="Output")
        self.grp_procmod = QGroupBox(title="Process Modifier")

        # Run Button
        self.frame_runExploreASL = QFrame()
        self.frame_runExploreASL.setFrameStyle(QFrame.StyledPanel | QFrame.Raised)
        self.frame_runExploreASL.setLineWidth(0)
        self.vlay_frame_runExploreASL = QVBoxLayout(self.frame_runExploreASL)
        self.btn_runExploreASL = QPushButton("Run Explore ASL", clicked=self.run_Explore_ASL)
        self.btn_runExploreASL.setEnabled(False)
        run_icon_font = QFont()
        run_icon_font.setPointSize(24)
        self.btn_runExploreASL.setFont(run_icon_font)
        set_widget_icon(self.btn_runExploreASL, self.config, "run_icon.svg", (75, 75))
        self.btn_runExploreASL.setMinimumHeight(80)
        self.btn_runExploreASL.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.vlay_frame_runExploreASL.addWidget(self.btn_runExploreASL)

        # Add main players to the appropriate splitters
        self.splitter_leftside.addWidget(self.grp_taskschedule)
        self.splitter_leftside.addWidget(self.grp_procmod)
        self.splitter_rightside.addWidget(self.grp_textoutput)
        self.splitter_rightside.addWidget(self.frame_runExploreASL)

        # Adjust splitter spacing, handle width, and display
        self.splitter_rightside.setSizes([self.height() - 100, 100])
        self.splitter_leftside.setSizes([self.height() - 200, 200])
        self.splitter_rightside.setHandleWidth(25)
        self.splitter_leftside.setHandleWidth(25)
        handle_path = os.path.join(self.config["ProjectDir"], "media", "3_dots_horizontal.svg").replace('\\', '/')
        handle_style = 'QSplitter::handle {image: url(' + handle_path + ');}'
        self.splitter_rightside.setStyleSheet(handle_style)
        self.splitter_leftside.setStyleSheet(handle_style)

        self.mainlay.addWidget(self.splitter_leftside)
        self.mainlay.addWidget(self.splitter_rightside)

    # Left side setup; define the number of studies
    def UI_Setup_TaskScheduler(self):
        self.vlay_scrollholder = QVBoxLayout(self.grp_taskschedule)
        self.scroll_taskschedule = QScrollArea()
        self.cont_taskschedule = QWidget()
        self.scroll_taskschedule.setWidget(self.cont_taskschedule)
        self.scroll_taskschedule.setWidgetResizable(True)
        self.vlay_scrollholder.addWidget(self.scroll_taskschedule)

        self.vlay_taskschedule = QVBoxLayout(self.cont_taskschedule)
        self.lab_coresinfo = QLabel(text=f"CPU Count: A total of {os.cpu_count() // 2} "
                                         f"processors are available on this machine")
        self.ncores_left = os.cpu_count() // 2
        self.lab_coresleft = QLabel(text=f"You are permitted to set up to {self.ncores_left} more core(s)")
        self.cont_nstudies = QWidget()
        self.hlay_nstudies = QHBoxLayout(self.cont_nstudies)
        self.lab_nstudies = QLabel(text=f"Indicate the number of studies you wish to process:")
        self.cmb_nstudies = QComboBox(self.cont_nstudies)
        self.nstudies_options = ["Select"] + list(map(str, range(1, (os.cpu_count() // 2 + 1))))
        self.cmb_nstudies.addItems(self.nstudies_options)
        self.cmb_nstudies.currentTextChanged.connect(self.UI_Setup_TaskScheduler_FormUpdate)
        self.cmb_nstudies.currentTextChanged.connect(self.set_ncores_left)
        self.cmb_nstudies.currentTextChanged.connect(self.is_ready_to_run)
        self.hlay_nstudies.addWidget(self.lab_nstudies)
        self.hlay_nstudies.addWidget(self.cmb_nstudies)

        # self.cont_filler = QWidget(self.grp_taskschedule)
        # self.formlay_filler = QFormLayout(self.cont_filler)
        # self.formlay_filler.addRow("Cores to Allocate ||", QLabel(text="Filepaths to Analysis Directories"))

        self.cont_tasks = QWidget(self.grp_taskschedule)
        self.formlay_tasks = QFormLayout(self.cont_tasks)

        self.cont_progbars = QWidget(self.grp_taskschedule)
        self.formlay_progbars = QFormLayout(self.cont_progbars)

        # Need python lists to keep track of row additions/removals; findChildren's ordering is incorrect
        self.formlay_lineedits_list = []
        self.formlay_buttons_list = []
        self.formlay_cmbs_ncores_list = []
        self.formlay_cmbs_runopts_list = []
        self.formlay_nrows = 0
        self.formlay_progbars_list = []
        self.formlay_stopbtns_list = []

        self.vlay_taskschedule.addWidget(self.lab_coresinfo)
        self.vlay_taskschedule.addWidget(self.lab_coresleft)
        self.vlay_taskschedule.addWidget(self.cont_nstudies)
        # self.vlay_taskschedule.addWidget(self.cont_filler)
        self.vlay_taskschedule.addWidget(self.cont_tasks)
        self.vlay_taskschedule.addWidget(self.cont_progbars)
        self.vlay_taskschedule.addStretch(2)

        self.cmb_nstudies.setCurrentIndex(1)

    # Right side setup; this will have a text editor to display feedback coming from ExploreASL or any future watchers
    # that are installed. Also, the Run buttons will be set up here.
    def UI_Setup_TextFeedback_and_Executor(self):
        self.vlay_textoutput = QVBoxLayout(self.grp_textoutput)
        self.textedit_textoutput = QTextEdit(self.grp_textoutput)
        self.textedit_textoutput.setPlaceholderText("Processing Progress will appear within this window")
        self.vlay_textoutput.addWidget(self.textedit_textoutput)

    # Rare exception of a UI function that is also technically a setter; this will dynamically alter the number of
    # rows present in the task scheduler form layout to allow for ExploreASL analysis of multiple studies at once
    def UI_Setup_TaskScheduler_FormUpdate(self, n_studies):
        if n_studies == "Select":
            return  # Don't do anything if the user selects Select again by accident
        n_studies = int(n_studies)
        diff = n_studies - self.formlay_nrows  # The difference between the current n_rows and n_studies

        # Addition of rows
        if diff > 0:
            for ii in range(diff):
                self.formlay_nrows += 1

                # Combobox to specify the number of cores to allocate
                inner_cmb_ncores = QComboBox()
                inner_cmb_ncores.setMinimumWidth(140)
                inner_cmb_ncores.setToolTip("Specify the number of cores to allocate to this study.\n"
                                            "Important points:\n"
                                            "\t-DO NOT specify more cores than there are subjects for the study\n"
                                            "\t-DO NOT specify more than one core for a study that will have the \n"
                                            "\tPopulation Module run on it")
                inner_cmb_ncores.addItems(list(map(str, range(1, os.cpu_count() // 2 + 1))))
                inner_cmb_ncores.currentTextChanged.connect(self.set_ncores_left)
                inner_cmb_ncores.currentTextChanged.connect(self.set_ncores_selectable)
                inner_cmb_ncores.currentTextChanged.connect(self.set_nstudies_selectable)
                inner_cmb_ncores.currentTextChanged.connect(self.is_ready_to_run)

                # Lineedit to specify
                inner_le = DandD_FileExplorer2LineEdit(acceptable_path_type="Directory")
                inner_le.setToolTip("Specify the filepath to the analysis folder of your study.\nFor example:\n"
                                    "C:\\Users\\JohnSmith\\MyStudy\\analysis")
                inner_le.setClearButtonEnabled(True)
                inner_le.setPlaceholderText("Select the analysis directory to your study")
                inner_le.textChanged.connect(self.is_ready_to_run)
                inner_btn_browsedirs = RowAwareQPushButton(self.formlay_nrows, "...")
                inner_btn_browsedirs.row_idx_signal.connect(self.set_analysis_directory)

                # Combobox to specify when module to run
                inner_cmb_procopts = QComboBox()
                inner_cmb_procopts.setToolTip("Specify which ExploreASL module to run:\n"
                                              "\t-Structural: Structural Module for processing T1w and FLAIR scans\n"
                                              "\t-ASL: ASL Module for processing ASL and M0 scans\n"
                                              "\t-Both: Run both the Structural and ASL modules\n"
                                              "\t-Population: Population module for determining statistics,\n"
                                              "\tstudywide masks, etc.")
                inner_cmb_procopts.addItems(["Structural", "ASL", "Both", "Population"])
                inner_cmb_procopts.setCurrentIndex(2)
                inner_cmb_procopts.currentTextChanged.connect(self.is_ready_to_run)

                # Stop Button
                stop_icon_path = os.path.join(self.config["ProjectDir"], "media", "stop_processing.png")
                inner_btn_stop = QPushButton(QIcon(stop_icon_path), "", enabled=False)
                inner_btn_stop.setToolTip("Stop this study being run")

                # Package it all in another FormLayout
                inner_grp = QGroupBox(f"Study {inner_btn_browsedirs.row_idx}")
                inner_formlay = QFormLayout(inner_grp)
                inner_formlay.addRow("Assigned Number of Cores", inner_cmb_ncores)
                inner_hbox = QHBoxLayout()
                inner_hbox.addWidget(inner_le)
                inner_hbox.addWidget(inner_btn_browsedirs)
                inner_formlay.addRow("Path to Analysis Folder", inner_hbox)
                inner_formlay.addRow("Which Modules to Run", inner_cmb_procopts)
                inner_formlay.addRow("Terminate Processing", inner_btn_stop)

                # Package it all in an HBoxLayout
                # inner_hbox = QHBoxLayout()
                # inner_hbox.addWidget(inner_le)
                # inner_hbox.addWidget(inner_btn_browsedirs)
                # inner_hbox.addWidget(inner_cmb_procopts)
                # inner_hbox.addWidget(inner_btn_stop)

                # Add progressbars
                inner_progbar = QProgressBar(orientation=Qt.Horizontal, value=0, maximum=100, minimum=0)
                inner_progbar.setPalette(self.green_palette)

                # Update format layouts through addition of the appropriate row
                # self.formlay_tasks.addRow(inner_cmb_ncores, inner_hbox)
                self.formlay_tasks.addRow(inner_grp)
                self.formlay_progbars.addRow(f"Study {inner_btn_browsedirs.row_idx}", inner_progbar)

                # Add widgets to their respective containers
                self.formlay_cmbs_ncores_list.append(inner_cmb_ncores)
                self.formlay_lineedits_list.append(inner_le)
                self.formlay_buttons_list.append(inner_btn_browsedirs)
                self.formlay_cmbs_runopts_list.append(inner_cmb_procopts)
                self.formlay_progbars_list.append(inner_progbar)
                self.formlay_stopbtns_list.append(inner_btn_stop)

        # Removal of rows
        elif diff < 0:
            for ii in range(abs(diff)):
                row_to_remove = self.formlay_nrows - 1

                # Update format layouts through removal of the appropriate row
                self.formlay_tasks.removeRow(row_to_remove)
                self.formlay_progbars.removeRow(row_to_remove)

                # Remove widgets from their respective containers starting from the latest addition
                self.formlay_cmbs_ncores_list.pop()
                self.formlay_lineedits_list.pop()
                self.formlay_buttons_list.pop()
                self.formlay_cmbs_runopts_list.pop()
                self.formlay_progbars_list.pop()
                self.formlay_stopbtns_list.pop()
                self.formlay_nrows -= 1

        # Adjust the number of cores selectable in each of the comboboxes
        self.set_ncores_left()
        self.set_ncores_selectable()
        self.set_nstudies_selectable()
        self.is_ready_to_run()

    # Left side setup; launches additional windows for specialized jobs such as modifying participants.tsv, preparing
    # a study for re-run in certain subjects, etc.
    def UI_Setup_ProcessModification(self):
        self.vlay_procmod = QVBoxLayout(self.grp_procmod)
        self.formlay_promod = QFormLayout()
        # Set up the widgets in this section
        self.cmb_modjob = QComboBox(self.grp_procmod)
        self.cmb_modjob.addItems(["Re-run a study", "Alter participants.tsv"])
        self.cmb_modjob.setToolTip("Specify the type of re-run or pre-processing modification you'd like to perform.\n"
                                   "Currently the following options are avaliable:\n"
                                   "\t'Re-run a study': Re-run parts of a previously-run study\n"
                                   "\t'Alter participants.tsv': Add metadata to the tsv file such that biasfields for\n"
                                   "\tthat metadata may be created when running the Population module")
        (self.hlay_modjob,
         self.le_modjob,
         self.btn_modjob) = make_droppable_clearable_le(btn_connect_to=self.set_modjob_analysis_dir,
                                                        acceptable_path_type="Directory")
        self.le_modjob.setPlaceholderText("Drag & Drop analysis directory here")
        self.le_modjob.setToolTip("Specify the filepath to the analysis folder of your study.\nFor example:\n"
                                  "C:\\Users\\JohnSmith\\MyStudy\\analysis")
        self.btn_runmodjob = QPushButton("Modify for Re-run", self.grp_procmod, clicked=self.run_modjob)

        modjob_icon_font = QFont()
        modjob_icon_font.setPointSize(24)
        self.btn_runmodjob.setFont(modjob_icon_font)
        set_widget_icon(self.btn_runmodjob, self.config, "run_modjob_icon.svg", (75, 75))
        self.btn_runmodjob.setMinimumHeight(80)
        self.btn_runmodjob.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.MinimumExpanding)

        # Add the widgets to the form layout
        self.formlay_promod.addRow("Which modification job to run", self.cmb_modjob)
        self.formlay_promod.addRow("Path to analysis dir to modify", self.hlay_modjob)
        self.vlay_procmod.addLayout(self.formlay_promod)
        self.vlay_procmod.addWidget(self.btn_runmodjob)

    # Runs the selected modification widget
    # noinspection PyCallByClass
    def run_modjob(self):
        selected_job = self.cmb_modjob.currentText()
        modjob_widget = None
        if selected_job == "Alter participants.tsv":
            # Requirements to launch the widget
            try:
                if any([self.le_modjob.text() == '',
                        "participants.tsv" not in os.listdir(self.le_modjob.text())  # participants.tsv must be present
                        ]):
                    QMessageBox().warning(self,
                                          "participants.tsv not found",
                                          f"participants.tsv was not located in the study directory you provided:\n"
                                          f"{self.le_modjob.text()}\n"
                                          f"Please run the Population module for this study "
                                          f"to generate the required file",
                                          QMessageBox.Ok)
                    return
            except FileNotFoundError:
                QMessageBox().warning(self,
                                      "participants.tsv not found",
                                      f"participants.tsv was not located in the study directory you provided:\n"
                                      f"{self.le_modjob.text()}\n"
                                      f"Please run the Population module for this study to generate the required file",
                                      QMessageBox.Ok)
                return
            modjob_widget = xASL_GUI_TSValter(self)
        elif selected_job == "Re-run a study":
            # Requirements to launch the widget
            try:
                if any([self.le_modjob.text() == '',
                        'lock' not in os.listdir(self.le_modjob.text())  # lock dirs must be initialized
                        ]):
                    QMessageBox().warning(self,
                                          "lock directory system not found",
                                          f"The study you provided does not have a lock directory system in play.\n"
                                          f"This would be automatically created if a study was run. Please run your\n"
                                          f"study first before trying to modify it for a re-run.",
                                          QMessageBox.Ok)
                    return
            except FileNotFoundError:
                QMessageBox().warning(self,
                                      "lock directory system not found",
                                      f"The study you provided does not have a lock directory system in play.\n"
                                      f"This would be automatically created if a study was run. Please run your\n"
                                      f"study first before trying to modify it for a re-run.",
                                      QMessageBox.Ok)
                return
            modjob_widget = xASL_GUI_RerunPrep(self)

        if modjob_widget is not None:
            modjob_widget.show()

    ##########################
    # TASK SCHEDULER FUNCTIONS
    ##########################

    # Function responsible for adjusting the label of how many cores are still accessible
    def set_ncores_left(self):
        self.ncores_left = os.cpu_count() // 2 - sum([int(cmb.currentText()) for cmb in self.formlay_cmbs_ncores_list])
        if self.ncores_left > 0:
            self.lab_coresleft.setText(f"You are permitted to set up to {self.ncores_left} more core(s)")
        elif self.ncores_left == 0:
            self.lab_coresleft.setText(f"No more cores are avaliable for allocation")
        else:
            self.lab_coresleft.setText(f"Something went terribly wrong")

    # Function responsible for adjusting the choices avaliable within each of the comboboxes of a given task row
    def set_ncores_selectable(self):
        self.ncores_left = os.cpu_count() // 2 - sum([int(cmb.currentText()) for cmb in self.formlay_cmbs_ncores_list])
        cores_left = self.ncores_left
        for box in self.formlay_cmbs_ncores_list:
            current_selection = int(box.currentText())
            max_cores_allowed = current_selection + cores_left
            for idx in range(box.count()):
                val_at_idx = int(box.itemText(idx))
                if val_at_idx <= max_cores_allowed:
                    box.model().item(idx).setEnabled(True)
                else:
                    box.model().item(idx).setEnabled(False)

    # Function responsible for adjusting the number of studies still permitted (assuming 1 core will be initially
    # allocated to it)
    def set_nstudies_selectable(self):
        self.ncores_left = os.cpu_count() // 2 - sum([int(cmb.currentText()) for cmb in self.formlay_cmbs_ncores_list])
        current_n_studies = int(self.cmb_nstudies.currentText())
        max_studies_allowed = current_n_studies + self.ncores_left
        for idx in range(self.cmb_nstudies.count()):
            val_at_idx = self.cmb_nstudies.itemText(idx)
            if not val_at_idx.isdigit():
                continue
            val_at_idx = int(val_at_idx)
            if val_at_idx <= max_studies_allowed:
                self.cmb_nstudies.model().item(idx).setEnabled(True)
            else:
                self.cmb_nstudies.model().item(idx).setEnabled(False)

    # This slot is responsible for setting the correct analysis directory to a given task row's lineedit correcting
    @Slot(int)
    def set_analysis_directory(self, row_idx):
        """
        :param row_idx: The index of the row from the pushbutton calling this slot. It is a non-pythonic index, so it
        must be reduced by 1 to index properly
        """
        dir_path = QFileDialog.getExistingDirectory(self.cw,
                                                    "Select the analysis directory of your study",
                                                    self.parent().config["DefaultRootDir"],
                                                    QFileDialog.ShowDirsOnly)
        set_os_dependent_text(linedit=self.formlay_lineedits_list[row_idx - 1],
                              config_ossystem=self.config["Platform"],
                              text_to_set=dir_path)

    def set_modjob_analysis_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self.cw,
                                                    "Select the analysis directory of your study",
                                                    self.config["DefaultRootDir"],
                                                    QFileDialog.ShowDirsOnly)
        set_os_dependent_text(linedit=self.le_modjob,
                              config_ossystem=self.config["Platform"],
                              text_to_set=dir_path)

    # Define whether the run Explore ASL button should be enabled
    def is_ready_to_run(self):
        # First check: all lineedits must have an appropriate analysis directory with a par file
        checks = []
        for le in self.formlay_lineedits_list:
            directory = le.text()
            if os.path.exists(directory):
                if all([os.path.isdir(directory),
                        len(glob(os.path.join(directory, "*Par*.json")))]):
                    checks.append(True)
                else:
                    checks.append(False)
            else:
                checks.append(False)

        # Second check: for any study that has its module set as Population, the number of cores must be 1
        for cmb_cores, le_analysisdir, cmb_runopt in zip(self.formlay_cmbs_ncores_list,
                                                         self.formlay_lineedits_list,
                                                         self.formlay_cmbs_runopts_list):
            if cmb_runopt.currentText() == "Population" and cmb_cores.currentText() != "1":
                checks.append(False)

            if cmb_runopt.currentText() != "Population":
                try:
                    with open(os.path.join(le_analysisdir.text(), "DataPar.json")) as datapar_reader:
                        regex = json.load(datapar_reader)["subject_regexp"]
                    n_subjects = sum([True if re.search(regex, subject) else False
                                      for subject in os.listdir(le_analysisdir.text())])
                    if n_subjects >= int(cmb_cores.currentText()):
                        checks.append(True)
                    else:
                        checks.clear()
                        self.btn_runExploreASL.setEnabled(False)
                        return
                except (KeyError, FileNotFoundError):
                    checks.clear()
                    self.btn_runExploreASL.setEnabled(False)
                    return

            if all(checks):
                self.btn_runExploreASL.setEnabled(True)
            else:
                self.btn_runExploreASL.setEnabled(False)

        ###################################
        # CONCURRENT AND POST-RUN FUNCTIONS
        ###################################

    def post_run_main(self):
        """
        Main wrapper function for all post-run followup
        """
        # Re-activate all relevant widgets
        self.set_widgets_activation_states(True)

        # Stop the movie
        self.end_movie()

        # Check if all expected operations took place
        self.post_run_statusfile_and_error_assessment()

        # Check the progressbars
        for progbar in self.formlay_progbars_list:
            if progbar.value() != progbar.maximum():
                progbar.setPalette(self.red_palette)

    def post_run_statusfile_and_error_assessment(self):
        """
        1) Assesses any stdout errors that were captured during the run. That is, processes that finish successfully
        with code 0 may have still experienced errors that stopped an iteration. These must be captured.
        2) Assesses any stderr errors that were captured during the run. That is, processes that did not finish
        correctly (i.e MATLAB crashed)
        3) Assess the STATUS files that were generated against the expected files
        4) Joins all the aforementioned together in an error file for that study
        """
        # First, assess the status of the stdout_errordicts_list
        master_stdout_err_dict = defaultdict(list)
        if len(self.stdout_errordicts_list) > 0:
            # Must convert from a list of dicts to a dict of lists
            for stdout_errordict in self.stdout_errordicts_list:
                for study_path, study_errors in stdout_errordict.items():
                    master_stdout_err_dict[study_path].append(study_errors)

        # Same deal for the stderr_errordicts_list
        master_stderr_err_dict = defaultdict(list)
        if len(self.stderr_errordicts_list) > 0:
            # Must convert from a list of dicts to a dict of lists
            for stderr_errordict in self.stderr_errordicts_list:
                for study_path, study_errors in stderr_errordict.items():
                    master_stderr_err_dict[study_path].append(study_errors)

        # Next, assess the status files
        postrun_diagnosis = {}
        for study_dir, expected_files in self.expected_status_files.items():
            all_completed, incomplete_status_files = calculate_missing_STATUS(study_dir, expected_files)

            if self.config["DeveloperMode"] and len(incomplete_status_files) > 0:
                print(f"INCOMPLETE STATUS FILES FOR STUDY DIR: {study_dir}:")
                pprint(incomplete_status_files)

            if all_completed:
                continue
            else:
                postrun_diagnosis[study_dir] = incomplete_status_files

        # If at least one study had an issue, print out a warning and generate a log file in that study directory
        if len(postrun_diagnosis) > 0:
            for study_dir, incomplete_files in postrun_diagnosis.items():
                (structmod_msgs,
                 aslmod_msgs,
                 popmod_msgs) = interpret_statusfile_errors(incomplete_files, self.executor_translators)

                specific_stdout_msg = []
                # Extract the specific error messages
                forward_study_dir = study_dir.replace("\\", "/")
                if forward_study_dir in list(master_stdout_err_dict.keys()):
                    study_errs = master_stdout_err_dict[forward_study_dir]
                    for err_dict in study_errs:
                        for subject, err in err_dict.items():
                            specific_stdout_msg.append(f"{subject}:\n{err.strip()}\n\n")

                specific_stderr_msg = []
                # Extract the specific stderr messages
                if forward_study_dir in list(master_stderr_err_dict.keys()):
                    study_errs = master_stderr_err_dict[forward_study_dir]
                    for std_err_msg in study_errs:
                        specific_stderr_msg.append(std_err_msg)

                msg = ["The following could not be completed during the run:\n"] + \
                      ["\n########################", "\nIn the Structural Module:\n"] + \
                      [file + '\n' for file in structmod_msgs] + \
                      ["\n########################", "\nIn the ASL Module:\n"] + \
                      [file + '\n' for file in aslmod_msgs] + \
                      ["\n########################", "\nIn the Population Module:\n"] + \
                      [file + '\n' for file in popmod_msgs] + \
                      ["\n########################", "\nSpecific Stdout Errors:\n"] + specific_stdout_msg + \
                      ["\n########################", "\nSpecific Stderr Errors:\n"] + specific_stderr_msg

                with open(os.path.join(study_dir, f"{date.today()} ExploreASL run - errors.txt"), 'w') as writer:
                    writer.writelines(msg)

            dirs_string = "\n".join(list(postrun_diagnosis.keys()))
            QMessageBox().warning(self,
                                  f"Errors were encountered during ExploreASL run",
                                  f"The following study directories encountered at least 1 error during their run:\n"
                                  f"{dirs_string}\n"
                                  f"Please review the generated [Date of run] ExploreASL run - errors.txt \n"
                                  f"within each of the listed studies to see which subjects could not be processed",
                                  QMessageBox.Ok)
        else:
            QMessageBox().information(self,
                                      f"Successful ExploreASL run",
                                      f"All expected operations successfully took place.\n"
                                      f"Many thanks for using ExploreASL. If this has been helpful in your study "
                                      f"please don't forget to cite this program in your manuscript.",
                                      QMessageBox.Ok)

    # Based on signal outout, this receives messages from watchers and outputs text feedback to the user
    @Slot(str)
    def update_text(self, msg):
        self.textedit_textoutput.append(msg)

    @Slot(dict)
    def update_stdout_error_dicts(self, stdout_error_dict):
        self.stdout_errordicts_list.append(stdout_error_dict)
        if self.config["DeveloperMode"]:
            print("update_stdout_error_dicts got a signal")

    @Slot(dict)
    def update_stderr_error_dicts(self, stderr_error_dict):
        self.stderr_errordicts_list.append(stderr_error_dict)
        if self.config["DeveloperMode"]:
            print("update_stderr_error_dicts got a signal")

    # This slot is responsible for updating the progressbar based on signals set from the watcher
    @Slot(int, int)
    def update_progressbar(self, val_to_inc_by, study_idx):
        """
        :param val_to_inc_by: the value of the .STATUS file that was just completed
        :param study_idx: the index of the progressbar contained within formlay_progbars_list to be selected
        """
        selected_progbar: QProgressBar = self.formlay_progbars_list[study_idx]
        if self.config["DeveloperMode"]:
            print(
                f"update_progressbar received a signal from watcher {study_idx} to increase the progressbar value by "
                f"{val_to_inc_by}")
            print(f"The progressbar's value before update: {selected_progbar.value()} "
                  f"out of maximum {selected_progbar.maximum()}")
        selected_progbar.setValue(selected_progbar.value() + val_to_inc_by)
        if self.config["DeveloperMode"]:
            print(f"The progressbar's value after update: {selected_progbar.value()} "
                  f"out of maximum {selected_progbar.maximum()}")

    # This slot is responsible for re-activating the Run ExploreASL button
    @Slot()
    def update_process_debt_and_check_done(self):
        """
        Receives a signal that a process has finished, regardless of whether it was a success or crash, updates the
        total debt accordingly, and launches the post-run function if the debt is cleared
        """
        self.total_process_dbt += 1
        if self.config["DeveloperMode"]:
            print(f"update_process_debt_and_check_done received a signal and incremented the debt. "
                  f"The debt at this time is: {self.total_process_dbt}")

        # If the debt is cleared, proceed to the post-run assessments and widget processing
        if self.total_process_dbt == 0:
            self.post_run_main()

    # Convenience function; deactivates all widgets associated with running exploreASL
    def set_widgets_activation_states(self, state: bool):
        self.btn_runExploreASL.setEnabled(state)
        self.cmb_nstudies.setEnabled(state)
        for core_cmb, le, browse_btn, runopt_cmb, stop_btn in zip(self.formlay_cmbs_ncores_list,
                                                                  self.formlay_lineedits_list,
                                                                  self.formlay_buttons_list,
                                                                  self.formlay_cmbs_runopts_list,
                                                                  self.formlay_stopbtns_list):
            core_cmb.setEnabled(state)
            le.setEnabled(state)
            browse_btn.setEnabled(state)
            runopt_cmb.setEnabled(state)
            stop_btn.setEnabled(not state)

    # Convenience function; adds a gif indicating processing below the textoutput and starts the gif
    def begin_movie(self):
        """
        Convenience function; adds a gif indicating processing below the textoutput and starts the gif
        """
        self.processing_movieplayer.movie.start()
        self.processing_movieplayer.setVisible(True)
        self.vlay_textoutput.addWidget(self.processing_movieplayer)

    # Convenience function; removes the gif below the textoutput and stops the gif
    def end_movie(self):
        """
        Convenience function; removes the gif below the textoutput and stops the gif
        """
        self.processing_movieplayer.movie.stop()
        self.processing_movieplayer.setVisible(False)
        self.vlay_textoutput.removeWidget(self.processing_movieplayer)

    ###################################################################################################################
    #                                              THE MAIN RUN FUNCTION
    ###################################################################################################################
    def run_Explore_ASL(self):

        # Immediately abandon this if the MATLAB version is not newer
        version_checks = []
        for version in ["R2019a", "R2019b", "R2020a", "R2020b", "R2021a"]:
            if version in os.path.basename(self.config["MATLABROOT"]):
                version_checks.append(True)
            else:
                version_checks.append(False)

        # If none of those versions could be detected, exit out
        if not any(version_checks):
            QMessageBox().warning(self,
                                  f"Incompatible MATLAB version on your machine",
                                  f"The program has detected that you have MATLAB version:\n"
                                  f"{os.path.basename(self.config['MATLABROOT'])}\n"
                                  f"This program requires a MATLAB installation of 2019a or later.",
                                  QMessageBox.Ok)
            del version_checks, version
            return
        else:  # Otherwise, garbage collect
            del version_checks, version

        if self.config["DeveloperMode"]:
            print("%" * 60)
        translator = {"Structural": [1], "ASL": [2], "Both": [1, 2], "Population": [3]}
        self.workers = []
        self.watchers = []
        self.total_process_dbt = 0
        self.expected_status_files = {}
        self.stdout_errordicts_list = []
        self.stderr_errordicts_list = []

        # Disable widgets
        self.set_widgets_activation_states(False)

        # Clear the textoutput each time
        self.textedit_textoutput.clear()

        # Outer for loop; loops over the studies
        for study_idx, (box, path, run_opts, progressbar, stop_btn) in enumerate(
                zip(self.formlay_cmbs_ncores_list,  # Comboboxes for number of cores
                    self.formlay_lineedits_list,  # Lineedits for analysis directory path
                    self.formlay_cmbs_runopts_list,  # Comboboxes for the modules to run
                    self.formlay_progbars_list,  # Progressbars for the study
                    self.formlay_stopbtns_list)):  # Stop Buttons for that study

            #########################################
            # INNER FOR LOOP - For a particular study
            #########################################
            # Create a container to keep a block of workers in for a particular study
            # Also create a debt counter to be fed to the watcher so that it can eventually disengage
            inner_worker_block = []
            debt = 0

            # %%%%%%%%%%%%%%%%%%%%%%
            # Step 1 - For the study, load in the parameters
            try:
                parms_file = glob(os.path.join(path.text(), "*Par*.json"))[0]
            except IndexError:
                QMessageBox().warning(self,
                                      f"Problem prior to starting ExploreASL",
                                      f"No DataPar file found within study:\n{path.text()}\n"
                                      f"Please ensure that the study's parameters file is present within the directory",
                                      QMessageBox.Ok)

                # Remember to re-activate widgets
                self.set_widgets_activation_states(True)
                return

            # Get a few essentials that will be needed for the workers and the watcher for this study
            try:
                with open(parms_file) as f:
                    parms = json.load(f)
                    try:
                        str_regex: str = parms["subject_regexp"]
                        explore_asl_path = parms["MyPath"]
                        regex = re.compile(str_regex)
                        sess_names = parms["SESSIONS"]
                        if str_regex.startswith("^") or str_regex.endswith('$'):
                            str_regex = str_regex.strip('^$')
                    except KeyError:
                        QMessageBox().warning(self,
                                              f"Problem prior to starting ExploreASL",
                                              f"Essential parameters not present within the DataPar file for study:"
                                              f"\n{path.text()}\n"
                                              f"Please ensure that parameter file contains fields detailing:\n"
                                              f"1) A MyPath key detailing the path to the Explore ASL directory\n"
                                              f"2) A subject_regex key representing subjects in the study\n"
                                              f"3) A SESSIONS key detailing the sessions in the study, if any",
                                              QMessageBox.Ok)

                        # Remember to re-activate widgets
                        self.set_widgets_activation_states(True)
                        return
            except json.decoder.JSONDecodeError as json_read_error:
                QMessageBox().warning(self.parent(),
                                      f"Problem prior to starting ExploreASL - Bad Json syntax",
                                      f"The DataPar.json file in {os.path.dirname(parms_file)}\n"
                                      f"could not be read.\n"
                                      f"Specifically, the error capture is as follows:\n"
                                      f"{json_read_error}\n"
                                      f"Please either correct this error or ensure the DataPar.json file has come "
                                      f"from the GUI module for defining study-level run parameters.",
                                      QMessageBox.Ok)
                # Remember to re-activate widgets
                self.set_widgets_activation_states(True)
                return

            # With the parms now loaded, additional checks can be made to prevent false runs
            if any([not os.path.exists(parms["MyPath"]),  # ExploreASL dir must exist
                    sum([bool(regex.search(file)) for file in os.listdir(path.text())
                         if file not in ["lock", "Population"]]) == 0,  # The regex must match at least 1 subject
                    str_regex == '',  # The string used to make the regex cannot be blank
                    parms["D"]["ROOT"] != path.text()  # The study provided must match the one in the parms file
                    ]):
                QMessageBox().warning(self,
                                      f"Problem prior to starting ExploreASL",
                                      f"An error was encountered while preparing study:\n{path.text()}\n"
                                      f"Please ensure:\n"
                                      f"1) That the filepath to the study exists\n"
                                      f"2) That the filepath to the study in the parms file matches the one provided to"
                                      f" the Task Scheduler\n"
                                      f"3) That the ExploreASL filepath specified within the parms file exists\n"
                                      f"4) The subjects within the study are the same as when the parms file was "
                                      f"created",
                                      QMessageBox.Ok)

                # Remember to re-activate widgets
                self.set_widgets_activation_states(True)
                return

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Step 2 - Prepare the workers for that study
            # Inner for loop: loops over the range of the num of cores within a study. Each core will be an iWorker
            for ii in range(box.count()):
                if ii < int(box.currentText()):
                    worker = ExploreASL_Worker(
                        explore_asl_path,  # The exploreasl filepath
                        glob(os.path.join(path.text(), "*Par*.json"))[0].replace('\\', '/'),  # The filepath to parms
                        1,  # Process the data
                        1,  # Skip the pause
                        ii + 1,  # iWorker
                        int(box.currentText()),  # nWorkers
                        translator[run_opts.currentText()],  # Processing options
                        str_regex  # For catching errors for particular subjects
                    )

                    inner_worker_block.append(worker)
                    debt -= 1
                    self.total_process_dbt -= 1

            # Add the block to the main workers argument
            self.workers.append(inner_worker_block)

            # %%%%%%%%%%%%%%%%%%%%%%%%%
            # Step 3 - Create all lock directories in advance; this will give the watcher full oversight
            initialize_all_lock_dirs(analysis_dir=path.text(),
                                     regex=regex,
                                     run_options=run_opts.currentText(),
                                     session_names=sess_names)

            # Also delete any directories called "locked" in the study
            locked_dirs = iglob(os.path.join(path.text(), "**", "locked"), recursive=True)
            locked_dirs = peekable(locked_dirs)
            if locked_dirs:
                if self.config["DeveloperMode"]:
                    print(f"Detected locked direcorties in {path.text()} prior to starting ExploreASL. "
                          f"Removing them first.")
                for lock_dir in locked_dirs:
                    try:
                        os.rmdir(lock_dir)
                    except OSError as lock_err:  # Just in case a user tampers with the lock directory
                        print(f"{lock_err}...but proceeding to recursive delete")
                        rmtree(path=lock_dir, ignore_errors=True)

            # %%%%%%%%%%%%%%%%%%%%%%%%%
            # Step 4 - Calculate the anticipated workload based on missing .STATUS files; adjust the progressbar's
            # maxvalue from that
            workload, expected_status_files = calculate_anticipated_workload(parmsdict=parms,
                                                                             run_options=run_opts.currentText(),
                                                                             translators=self.executor_translators)
            if not workload or len(expected_status_files) == 0:
                # Remember to re-activate widgets
                self.set_widgets_activation_states(True)
                return
            progressbar.reset()
            progressbar.setMaximum(workload)
            progressbar.setMinimum(0)
            progressbar.setValue(0)
            progressbar.setPalette(self.green_palette)
            del workload
            # Save the expected status files to the dict container; these will be iterated over after workers are done
            self.expected_status_files[path.text()] = expected_status_files
            if self.config["DeveloperMode"]:
                print(f"EXPECTED STATUS FILES TO BE GENERATED FOR STUDY: {path.text()}")
                pprint(expected_status_files)

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Step 5 - Create a Watcher for that study
            watcher = ExploreASL_Watcher(target=path.text(),  # the analysis directory
                                         regex=str_regex,  # the regex used to recognize subjects
                                         watch_debt=debt,  # the debt used to determine when to stop watching
                                         study_idx=study_idx,
                                         # the identifier used to know which progressbar to signal
                                         translators=self.executor_translators,
                                         config=self.config
                                         )
            self.textedit_textoutput.append(f"Setting a Watcher thread on {path.text()}")

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Step 6 - Set up watcher connections
            # Connect the watcher to signal to the text output
            watcher.signals.update_text_output_signal.connect(self.update_text)
            # Connect the watcher to signal to the progressbar
            watcher.signals.update_progbar_signal.connect(self.update_progressbar)

            # Finally, add the watcher to the container
            self.watchers.append(watcher)

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Step 7 - Set up worker connections
            # Worker connections setup within a study
            for idx, worker in enumerate(inner_worker_block):
                # Connect the finished signal to the watcher debt to help it understand when it should stop watching
                worker.signals.finished_processing.connect(watcher.increment_debt)
                worker.signals.encountered_fatal_error.connect(watcher.increment_debt)
                # Connect the finished signal to the run btn reactivator so that the run button may be reactivated
                # via debt repayment counter
                worker.signals.finished_processing.connect(self.update_process_debt_and_check_done)
                worker.signals.encountered_fatal_error.connect(self.update_process_debt_and_check_done)
                # Connect the error encountered signal to the slot that updates the executor's track record of the
                # accounted for errors
                worker.signals.stdout_processing_error.connect(self.update_stdout_error_dicts)
                worker.signals.stderr_processing_error.connect(self.update_stderr_error_dicts)
                # Connect the stop button for that study to its terminator function
                stop_btn.clicked.connect(worker.terminate_run)

                self.textedit_textoutput.append(
                    f"Preparing Worker {idx + 1} of {len(inner_worker_block)} for study: "
                    f"{path.text()}")

        ######################################
        # THIS IS NOW OUTSIDE OF THE FOR LOOPS

        # self.watchers is nested at this point; we need to flatten it
        self.workers = list(chain(*self.workers))

        # Launch all threads in one go
        for runnable in self.workers + self.watchers:
            self.threadpool.start(runnable)

        self.begin_movie()
        begin_msg = r"""
##################################################################################
______            _                          _____ _         _____ _    _ _____ 
|  ____|          | |                  /\    / ____| |       / ____| |  | |_   _|
| |__  __  ___ __ | | ___  _ __ ___   /  \  | (___ | |      | |  __| |  | | | |  
|  __| \ \/ / '_ \| |/ _ \| '__/ _ \ / /\ \  \___ \| |      | | |_ | |  | | | |  
| |____ >  <| |_) | | (_) | | |  __// ____ \ ____) | |____  | |__| | |__| |_| |_ 
|______/_/\_\ .__/|_|\___/|_|  \___/_/    \_\_____/|______|  \_____|\____/|_____|
         | |                                                                  
         |_|                                                                  
______                     _                                                    
|  ____|                   | |                                                   
| |__  __  _____  ___ _   _| |_ ___  _ __                                        
|  __| \ \/ / _ \/ __| | | | __/ _ \| '__|                                       
| |____ >  <  __/ (__| |_| | || (_) | |                                          
|______/_/\_\___|\___|\__,_|\__\___/|_|                                          """
        print(begin_msg)


class ExploreASL_WatcherSignals(QObject):
    """
    Defines the signals avaliable from a running watcher thread.
    """
    update_progbar_signal = Signal(int, int)
    update_text_output_signal = Signal(str)
    update_debt_signal = Signal(int)


# noinspection PyCallingNonCallable
class ExploreASL_Watcher(QRunnable):
    """
    Modified file system watcher. Will monitor the appearance of STATUS files within the lock dirs of the analysis
    directory. If it detects a STATUS file, it will emit signals to:
    1) update the progress bars
    2) inform the text editor view of which STATUS file was made so as to give user feedback
    3)
    """

    def __init__(self, target, regex, watch_debt, study_idx, translators, config):
        super().__init__()
        self.signals = ExploreASL_WatcherSignals()
        self.dir_to_watch = os.path.join(target, "lock")
        self.subject_regex = re.compile(regex)
        self.module_regex = re.compile('module_(ASL|Structural|Population)')
        self.watch_debt = watch_debt
        self.study_idx = study_idx
        self.config = config

        self.pop_mod_started = False
        self.struct_mod_started = False
        self.asl_mod_started = False

        self.sessions_seen_struct = []
        self.sessions_seen_asl = []

        self.observer = Observer()
        self.event_handler = ExploreASL_EventHandler()
        self.event_handler.signals.inform_file_creation.connect(self.process_message)
        self.observer.schedule(event_handler=self.event_handler,
                               path=self.dir_to_watch,
                               recursive=True)
        self.stuct_status_file_translator = translators["Structural_Module_Filename2Description"]
        self.asl_status_file_translator = translators["ASL_Module_Filename2Description"]
        self.pop_status_file_translator = translators["Population_Module_Filename2Description"]
        self.workload_translator = translators["ExploreASL_Filename2Workload"]
        if self.config["DeveloperMode"]:
            print(
                f"Initialized a watcher for the directory {self.dir_to_watch} "
                f"and will communicate with the progressbar at Python idx: {self.study_idx}")

    # Processes the information sent from the event hander and emits signals to update widgets in the main Executor
    @Slot(str)
    def process_message(self, created_path):
        if self.config["DeveloperMode"]:
            print(f"Watcher process_message received message: {created_path}")

        detected_subject = self.subject_regex.search(created_path)
        detected_module = self.module_regex.search(created_path)
        msg = None
        workload_val = None

        if os.path.isdir(created_path):  # Lock dir
            if detected_module.group(
                    1) == "Structural" and not detected_subject.group() in self.sessions_seen_struct:
                self.sessions_seen_struct.append(detected_subject.group())
                msg = f"Structural Module has started for subject: {detected_subject.group()}"
            elif detected_module.group(1) == "ASL" and not detected_subject.group() in self.sessions_seen_asl:
                self.sessions_seen_asl.append(detected_subject.group())
                msg = f"ASL Module has started for subject: {detected_subject.group()}"
            elif detected_module.group(1) == "Population" and not self.pop_mod_started:
                self.pop_mod_started = True
                msg = f"Population Module has started"
            else:
                pass

        elif os.path.isfile(created_path):  # Status file
            basename = os.path.basename(created_path)
            if detected_module.group(1) == "Structural" and detected_subject:
                msg = f"Completed {self.stuct_status_file_translator[basename]} in the Structural module " \
                      f"for subject: {detected_subject.group()}"
            elif detected_module.group(1) == "ASL" and detected_subject:
                msg = f"Completed {self.asl_status_file_translator[basename]} in the ASL module " \
                      f"for subject: {detected_subject.group()}"
            elif detected_module.group(1) == "Population":
                msg = f"Completed {self.pop_status_file_translator[basename]} in the Population module"

            workload_val = self.workload_translator[basename]

        else:
            print("Neither a file nor a directory was detected")

        # Emit the message to inform the user of the most recent progress
        if msg:
            self.signals.update_text_output_signal.emit(msg)

        # Emit the workload value associated with the completion of that status file as well as the study idx so that
        # the appropriate progressbar is updated
        if workload_val:
            self.signals.update_progbar_signal.emit(workload_val, self.study_idx)

    # Must use slots system, as it is thread-safe
    @Slot()
    def increment_debt(self):
        self.watch_debt += 1

    def run(self):
        self.observer.start()
        if self.config["DeveloperMode"]:
            print(f"THE WATCHER FOR {self.dir_to_watch} HAS STARTED")
        while self.watch_debt < 0:
            sleep(10)
        self.observer.stop()
        self.observer.join()
        if self.config["DeveloperMode"]:
            print(f"THE WATCHER FOR {self.dir_to_watch} IS SHUTTING DOWN")
        return


class ExploreASL_EventHanderSignals(QObject):
    """
    Defines the signals used by the EventHandler class
    """
    inform_file_creation = Signal(str)


class ExploreASL_EventHandler(FileSystemEventHandler):
    """
    The real watcher behind the scenes
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.signals = ExploreASL_EventHanderSignals()

    def on_created(self, event):
        self.signals.inform_file_creation.emit(event.src_path)


class RowAwareQPushButton(QPushButton):
    """
    A subset of QPushButton that has awareness of which row within Task Scheduler it is located in at all times. This
    is a convenience measure for easily communicating to self.set_analysis_directory and knowing which lineedit to
    alter.
    """
    row_idx_signal = Signal(int)

    def __init__(self, row_idx, text, parent=None):
        super().__init__(text=text, parent=parent)
        self.row_idx = row_idx

    def mousePressEvent(self, e):
        self.row_idx_signal.emit(self.row_idx)
