import os
from glob import glob
import re
from ExploreASL_GUI.xASL_GUI_HelperClasses import DandD_FileExplorer2LineEdit
from PySide2.QtWidgets import QWidget, QLabel, QVBoxLayout
from PySide2.QtGui import QMovie, Qt
from PySide2.QtCore import QSize, QByteArray
from platform import system
from itertools import chain
from pprint import pprint


def initialize_all_lock_dirs(analysis_dir, regex, run_options, session_names):
    """
    Convenience function for creating the lock directories in advance of a run such that a file system watcher
    could effectively be set at the root and detect any downstream changes
    :param analysis_dir: the root analysis directory of the study ex. User\\Study_Name\\analysis
    :param regex: the regex used to identify subjects
    :param run_options: the type of run (ex. ASL, Structural, Both, Population)
    :param session_names: the expected session names that should be encountered (ex. ASL_1, ASL_2, etc.)
    """

    def dirnames_for_asl(analysis_directory, study_sessions, within_session_names):
        dirnames = [os.path.join(analysis_directory, "lock", f"xASL_module_ASL", sess, f"xASL_module_ASL_{name}")
                    for sess in study_sessions for name in within_session_names]
        return dirnames

    def dirnames_for_structural(analysis_directory, study_sessions):
        dirnames = [
            os.path.join(analysis_directory, "lock", f"xASL_module_Structural", sess, "xASL_module_Structural")
            for sess in study_sessions]
        return dirnames

    def dirnames_for_population(analysis_directory):
        return [os.path.join(analysis_directory, "lock", "xASL_module_Population", "xASL_module_Population")]

    print(f"Generating the appropriate lock dirs for {analysis_dir}")
    sessions = [session for session in os.listdir(analysis_dir) if regex.search(session)]

    # Prepare the list of names of lockdirs to expect to create based on the run options, session names and detected
    # sessions in the analysis root directory
    if run_options == "Both":
        struc_dirs = dirnames_for_structural(analysis_dir, sessions)
        asl_dirs = dirnames_for_asl(analysis_dir, sessions, session_names)
        lock_dirs = struc_dirs + asl_dirs
    elif run_options == "ASL":
        lock_dirs = dirnames_for_asl(analysis_dir, sessions, session_names)
    elif run_options == "Structural":
        lock_dirs = dirnames_for_structural(analysis_dir, sessions)
    elif run_options == "Population":
        lock_dirs = dirnames_for_population(analysis_dir)
    else:
        raise ValueError("Impossible outcome in initialize_all_lock_dirs")

    # Create empty directories where applicable
    for lock_dir in lock_dirs:
        if not os.path.exists(lock_dir):
            os.makedirs(lock_dir)


def calculate_anticipated_workload(parmsdict, run_options, translators):
    """
    Convenience function for calculating the anticipated workload
    :param parmsdict: the parameter file of the study; given parameters such as the regex are used from this
    :param run_options: "Structural", "ASL", "Both" or "Population"; which module is being run
    :param translators: The ExecutorTranslators, primarily for calculating the workload
    :return: workload; a numerical representation of the cumulative value of all status files made; these will be
    used to determine the appropriate maximum value for the progressbar
    """

    def get_structural_workload(analysis_directory, study_subjects, structuralmod_dict,
                                skip_if_no_asl, skip_if_no_m0, skip_if_no_flair, workload_translator):
        default_workload = ["010_LinearReg_T1w2MNI.status",
                            "060_Segment_T1w.status",
                            "080_Resample2StandardSpace.status",
                            "090_GetVolumetrics.status",
                            "100_VisualQC_Structural.status",
                            "999_ready.status"]
        flair_workload = ["020_LinearReg_FLAIR2T1w.status",
                          "030_FLAIR_BiasfieldCorrection.status",
                          "040_LST_Segment_FLAIR_WMH.status",
                          "050_LST_T1w_LesionFilling_WMH.status",
                          "070_CleanUpWMH_SEGM.status"]
        status_files = []
        for subject in study_subjects:

            # Ascertain if there are any FLAIR images and compare their presence/absence to the FLAIR skip flag
            path_to_flair = glob(os.path.join(analysis_directory, subject, "*FLAIR.nii*"))
            try:
                has_flair_img = os.path.exists(path_to_flair[0])
            except IndexError:
                has_flair_img = False

            path_to_m0 = glob(os.path.join(analysis_directory, subject, "*", "*M0.nii*"))
            try:
                has_m0_img = os.path.exists(path_to_m0[0])
            except IndexError:
                has_m0_img = False
            path_to_asl = glob(os.path.join(analysis_directory, subject, "*", "*ASL4D.nii*"))
            try:
                has_asl_img = os.path.exists(path_to_asl[0])
            except IndexError:
                has_asl_img = False

            # Do not proceed to if a particular subject/session is missing the required image based on parms
            if any([skip_if_no_m0 and not has_m0_img,
                    skip_if_no_asl and not has_asl_img,
                    skip_if_no_flair and not has_flair_img
                    ]):
                continue

            directory = os.path.join(analysis_directory, "lock", "xASL_module_Structural", subject,
                                     "xASL_module_Structural")
            current_status_files = os.listdir(directory)
            ###########################################
            # Keeping this here in case of switch back
            if has_flair_img:
                workload = default_workload + flair_workload
            else:
                workload = default_workload
            # workload = default_workload + flair_workload  # The full workload is assumed now every time
            ############################################
            # Filter out any anticipated status files that are already present in the lock dirs
            filtered_workload = set(workload).difference(set(current_status_files))
            # Append the filepaths; these will be used after analysis to check for incompleted STATUS workloads
            status_files.append([os.path.join(directory, status_file) for status_file in filtered_workload])
            numerical_representation = sum([workload_translator[stat_file] for stat_file in filtered_workload])
            structuralmod_dict[subject] = numerical_representation

        # Flatten the nested list that is status_files
        status_files = list(chain(*status_files))

        return structuralmod_dict, status_files

    def get_asl_workload(analysis_directory, study_subjects, session_names, aslmod_dict,
                         skip_if_no_asl, skip_if_no_m0, skip_if_no_flair, workload_translator):
        default_workload = ["020_RealignASL.status",
                            "030_RegisterASL.status",
                            "040_ResampleASL.status",
                            "050_PreparePV.status",
                            "060_ProcessM0.status",
                            "070_Quantification.status",
                            "080_CreateAnalysisMask.status",
                            "090_VisualQC_ASL.status",
                            "999_ready.status"]

        # Must iterate through both the subject level listing AND the session level (ASL_1, ASL_2, etc.) listing
        status_files = []
        for subject in study_subjects:
            for session in session_names:
                path_to_m0 = glob(os.path.join(analysis_directory, subject, session, "*M0.nii*"))
                try:
                    has_m0_img = os.path.exists(path_to_m0[0])
                except IndexError:
                    has_m0_img = False
                path_to_asl = glob(os.path.join(analysis_directory, subject, session, "*ASL4D.nii*"))
                try:
                    has_asl_img = os.path.exists(path_to_asl[0])
                except IndexError:
                    has_asl_img = False
                path_to_flair = glob(os.path.join(analysis_directory, subject, "*FLAIR.nii*"))
                try:
                    has_flair_img = os.path.exists(path_to_flair[0])
                except IndexError:
                    has_flair_img = False

                # Do not proceed to if a particular subject/session is missing the required image based on parms
                if any([skip_if_no_m0 and not has_m0_img,
                        skip_if_no_asl and not has_asl_img,
                        skip_if_no_flair and not has_flair_img
                        ]):
                    continue

                directory = os.path.join(analysis_directory, "lock", "xASL_module_ASL", subject,
                                         f"xASL_module_ASL_{session}")
                current_status_files = os.listdir(directory)
                workload = default_workload
                # Filter out any anticipated status files that are already present in the lock dirs
                filtered_workload = set(workload).difference(set(current_status_files))
                # Append the filepaths; these will be used after analysis to check for incompleted STATUS workloads
                status_files.append([os.path.join(directory, status_file) for status_file in filtered_workload])
                # Calculate the numerical representation of the STATUS files workload
                numerical_representation = sum([workload_translator[stat_file] for stat_file in filtered_workload])
                aslmod_dict[subject][session] = numerical_representation

        # Flatten the nested list that is status_files
        status_files = list(chain(*status_files))

        return aslmod_dict, status_files

    def get_population_workload(analysis_directory, workload_translator):
        default_workload = ["010_CreatePopulationTemplates.status",
                            "020_CreateAnalysisMask.status",
                            "030_CreateBiasfield.status",
                            "040_GetDICOMStatistics.status",
                            "050_GetVolumeStatistics.status",
                            "060_GetMotionStatistics.status",
                            "070_GetROIstatistics.status",
                            "080_SortBySpatialCoV.status",
                            "090_DeleteAndZip.status",
                            "999_ready.status"]

        directory = os.path.join(analysis_directory, "lock", "xASL_module_Population", "xASL_module_Population")
        current_status_files = os.listdir(directory)

        workload = default_workload
        filtered_workload = set(workload).difference(set(current_status_files))
        status_files = [os.path.join(directory, status_file) for status_file in filtered_workload]
        numerical_representation = sum([workload_translator[stat_file] for stat_file in filtered_workload])
        # No need for flattening the status_files for this one; not nested
        return numerical_representation, status_files

    # Define the individual translators
    filename2workload = translators["ExploreASL_Filename2Workload"]

    # First get all the subjects
    analysis_dir: str = parmsdict["D"]["ROOT"]
    regex = re.compile(parmsdict["subject_regexp"])
    sess_names: list = parmsdict["SESSIONS"]
    skipifnoasl: bool = parmsdict["SkipIfNoASL"]
    skipifnom0: bool = parmsdict["SkipIfNoM0"]
    skipifnoflair: bool = parmsdict["SkipIfNoFlair"]
    subjects = [subject for subject in os.listdir(analysis_dir) if
                all([regex.search(subject),  # regex must fit
                     os.path.isdir(os.path.join(analysis_dir, subject)),  # must be a directory
                     subject not in ["Population", "lock"]  # can't accidentally be the non-subject directories
                     ])]

    # Use a dict to keep track of everything
    struct_dict = {subject: 0 for subject in subjects}
    asl_dict = {subject: {} for subject in subjects}

    # Update the dicts as appropriate
    if run_options == "Both":
        struct_dict, struct_status = get_structural_workload(analysis_dir,
                                                             subjects,
                                                             struct_dict,
                                                             skipifnoasl,
                                                             skipifnom0,
                                                             skipifnoflair,
                                                             workload_translator=filename2workload)
        asl_dict, asl_status = get_asl_workload(analysis_dir,
                                                subjects,
                                                sess_names,
                                                asl_dict,
                                                skipifnoasl,
                                                skipifnom0,
                                                skipifnoflair,
                                                workload_translator=filename2workload)

        struct_totalworkload = sum(struct_dict.values())
        asl_totalworkload = {subject: sum(asl_dict[subject].values()) for subject in subjects}
        asl_totalworkload = sum(asl_totalworkload.values())
        print(f"Structural Calculated Workload: {struct_totalworkload}")
        print(f"ASL Calculated Workload: {asl_totalworkload}")
        # Return the numerical sum of the workload and the combined list of the expected status files
        return struct_totalworkload + asl_totalworkload, struct_status + asl_status

    elif run_options == "ASL":
        asl_dict, asl_status = get_asl_workload(analysis_dir,
                                                subjects,
                                                sess_names,
                                                asl_dict,
                                                skipifnoasl,
                                                skipifnom0,
                                                skipifnoflair,
                                                workload_translator=filename2workload)
        asl_totalworkload = {asl_subject: sum(asl_dict[asl_subject].values()) for asl_subject in subjects}
        asl_totalworkload = sum(asl_totalworkload.values())
        print(f"ASL Calculated Workload: {asl_totalworkload}")
        # Return the numerical sum of the workload and the list of expected status files
        return asl_totalworkload, asl_status

    elif run_options == "Structural":
        struct_dict, struct_status = get_structural_workload(analysis_dir,
                                                             subjects,
                                                             struct_dict,
                                                             skipifnoasl,
                                                             skipifnom0,
                                                             skipifnoflair,
                                                             workload_translator=filename2workload)
        # pprint(struct_dict)
        struct_totalworkload = sum(struct_dict.values())
        print(f"Structural Calculated Workload: {struct_totalworkload}")
        # Return the numerical sum of the workload and the list of expected status files
        return struct_totalworkload, struct_status

    elif run_options == "Population":
        pop_totalworkload, pop_status = get_population_workload(analysis_dir, workload_translator=filename2workload)
        print(f"Population Calculated Workload: {pop_totalworkload}")
        # Return the numerical sum of the workload and the list of expected status files
        return pop_totalworkload, pop_status

    else:
        print("THIS SHOULD NEVER PRINT AS YOU HAVE SELECTED AN IMPOSSIBLE WORKLOAD OPTION")


# Called after processing is done to compare the present status files against the files that were expected to be created
# at the time the run was initialized
def calculate_missing_STATUS(analysis_dir, expected_status_files):
    postrun_status_files = glob(os.path.join(analysis_dir, 'lock', "**", "*.status"), recursive=True)
    incomplete = [file for file in expected_status_files if file not in postrun_status_files]
    if len(incomplete) == 0:
        return True, incomplete
    else:
        return False, incomplete


def interpret_statusfile_errors(incomplete_files, translators: dict):
    """
    Interprets the step in the ExploreASL pipeline for which particular subjects/sessions/etc. must have failed and
    returns these interpretations as messages for each of the modules
    :param incomplete_files: the list of status files that were not generated in the pipeline
    :param translators: the translators (from JSON_LOGIC directory) used to convert filenames to their descriptions for
    generating the correct error message
    :return: 3 lists of interpreted error messages, one for each ExploreASL module
    """
    if system() == "Windows":
        delimiter = '\\\\'
    else:
        delimiter = '/'

    # Prepare containers and translators
    struct_dict = {}
    asl_dict = {}
    pop_list = []
    asl_msgs = []
    struct_msgs = []
    pop_msgs = []
    stuct_status_file_translator = translators["Structural_Module_Filename2Description"]
    asl_status_file_translator = translators["ASL_Module_Filename2Description"]
    population_file_translator = translators["Population_Module_Filename2Description"]

    # Prepare regex detectors
    asl_regex = re.compile(f"(?:.*){delimiter}lock{delimiter}xASL_module_(?:Structural|ASL){delimiter}(.*){delimiter}"
                           f"xASL_module_(?:Structural|ASL)_?(.*)?{delimiter}(.*\\.status)")
    pop_regex = re.compile(f"(?:.*){delimiter}lock{delimiter}xASL_module_Population{delimiter}"
                           f"xASL_module_Population{delimiter}(.*\\.status)")

    # Use the regex to extract subject, session, filename fields from the status filepaths, then organize the status
    # file basenames into the appropriate dictionary structure
    for file in incomplete_files:
        asl_match = asl_regex.search(file)
        pop_match = pop_regex.search(file)
        if asl_match:
            subject, session, file = asl_match.groups()
            # Structural
            if session == "":
                struct_dict.setdefault(subject, [])
                struct_dict[subject].append(file)
            # ASL
            else:
                asl_dict.setdefault(subject, {})
                asl_dict[subject].setdefault(session, [])
                asl_dict[subject][session].append(file)
        # Population
        elif pop_match:
            file = pop_match.group(1)
            pop_list.append(file)
        else:
            return

    # For each of the dictionaries corresponding to an ExploreASL module, sort the basenames and create the appropriate
    # error message for that subject/session
    if len(asl_dict) > 0:
        for subject, inner_1 in asl_dict.items():
            for session, files in inner_1.items():
                sorted_files = sorted(files)
                msg = f"Subject: {subject}; Session: {session}; Failed in the ASL module prior to: " \
                      f"{asl_status_file_translator[sorted_files[0]]}; " \
                      f"STATUS file failed to be created: {sorted_files[0]}"
                asl_msgs.append(msg)

    if len(struct_dict) > 0:
        for subject, files in struct_dict.items():
            sorted_files = sorted(files)
            msg = f"Subject: {subject}; Failed in the Structural module prior to: " \
                  f"{stuct_status_file_translator[sorted_files[0]]}; " \
                  f"STATUS file failed to be created: {sorted_files[0]}"
            struct_msgs.append(msg)

    if len(pop_list) > 0:
        sorted_files = sorted(pop_list)
        msg = f"Population module failed at: {population_file_translator[sorted_files[0]]}; " \
              f"STATUS file failed to be created: {sorted_files[0]}"
        pop_msgs = [msg]

    return struct_msgs, asl_msgs, pop_msgs
