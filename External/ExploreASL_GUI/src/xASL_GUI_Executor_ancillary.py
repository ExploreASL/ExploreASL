from pathlib import Path
from more_itertools import peekable
import re
from platform import system
from itertools import chain
from pprint import pprint
from typing import List, Set, Tuple


def initialize_all_lock_dirs(analysis_dir: Path, regex: re.Pattern, run_options: str, run_names: list):
    """
    Convenience function for creating the lock directories in advance of a run such that a file system watcher
    could effectively be set at the root and detect any downstream changes
    :param analysis_dir: the root analysis directory of the study ex. User\\Study_Name\\analysis
    :param regex: the regex used to identify subjects
    :param run_options: the type of run (ex. ASL, Structural, Both, Population)
    :param run_names: the expected session names that should be encountered (ex. ASL_1, ASL_2, etc.)
    """

    def dirnames_for_asl(analysis_directory: Path, subjects, runs):
        dirnames = [analysis_directory / "lock" / f"xASL_module_ASL" / subject / f"xASL_module_ASL_{run}"
                    for subject in subjects for run in runs]
        return dirnames

    def dirnames_for_structural(analysis_directory: Path, subjects):
        dirnames = [analysis_directory / "lock" / f"xASL_module_Structural" / subject / "xASL_module_Structural"
                    for subject in subjects]
        return dirnames

    def dirnames_for_population(analysis_directory):
        return [analysis_directory / "lock" / "xASL_module_Population" / "xASL_module_Population"]

    print(f"Generating the appropriate lock dirs for {analysis_dir}")
    subject_names = [subject.name for subject in analysis_dir.iterdir() if regex.search(str(subject))]

    # Prepare the list of names of lockdirs to expect to create based on the run options, session names and detected
    # sessions in the analysis root directory
    if run_options == "Both":
        struc_dirs = dirnames_for_structural(analysis_dir, subject_names)
        asl_dirs = dirnames_for_asl(analysis_dir, subject_names, run_names)
        lock_dirs = struc_dirs + asl_dirs
    elif run_options == "ASL":
        lock_dirs = dirnames_for_asl(analysis_dir, subject_names, run_names)
    elif run_options == "Structural":
        lock_dirs = dirnames_for_structural(analysis_dir, subject_names)
    elif run_options == "Population":
        lock_dirs = dirnames_for_population(analysis_dir)
    else:
        raise ValueError("Impossible outcome in initialize_all_lock_dirs")

    # Create empty directories where applicable
    for lock_dir in lock_dirs:
        if not lock_dir.exists():
            lock_dir.mkdir(parents=True)


def calculate_anticipated_workload(parmsdict, run_options, translators):
    """
    Convenience function for calculating the anticipated workload
    :param parmsdict: the parameter file of the study; given parameters such as the regex are used from this
    :param run_options: "Structural", "ASL", "Both" or "Population"; which module is being run
    :param translators: The ExecutorTranslators, primarily for calculating the workload
    :return: workload; a numerical representation of the cumulative value of all status files made; these will be
    used to determine the appropriate maximum value for the progressbar
    """

    def get_structural_workload(analysis_directory: Path, study_subjects, structuralmod_dict,
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
            # Do not proceed to if a particular subject/session is missing the required image based on parms
            try:
                has_flair_img = next((analysis_directory / subject).glob("*FLAIR.nii*")).exists()
            except StopIteration:
                has_flair_img = False
            try:
                has_m0_img = next((analysis_directory / subject).glob("*/*M0.nii*")).exists()
            except StopIteration:
                has_m0_img = False
            try:
                has_asl_img = next((analysis_directory / subject).glob("*/*ASL*.nii*")).exists()
            except StopIteration:
                has_asl_img = False

            if any([skip_if_no_m0 and not has_m0_img,
                    skip_if_no_asl and not has_asl_img,
                    skip_if_no_flair and not has_flair_img]):
                continue

            directory = analysis_directory / "lock" / "xASL_module_Structural" / subject / "xASL_module_Structural"
            workload = set(default_workload + flair_workload)  # The full workload is assumed now every time
            filtered_workload = [directory / name for name in workload if not (directory / name).exists()]
            # Filter out any anticipated status files that are already present in the lock dirs
            status_files.extend(filtered_workload.copy())
            num_repr = sum([workload_translator[stat_file.name] for stat_file in filtered_workload])
            structuralmod_dict[subject] = num_repr

        return structuralmod_dict, status_files

    def get_asl_workload(analysis_directory, study_subjects, session_names, aslmod_dict,
                         skip_if_no_asl, skip_if_no_m0, skip_if_no_flair, workload_translator,
                         conditions: List[Tuple[str, bool]] = None):
        if conditions is None:
            conditions = []
        workload = ["020_RealignASL.status",
                    "030_RegisterASL.status",
                    "040_ResampleASL.status",
                    "050_PreparePV.status",
                    "060_ProcessM0.status",
                    "070_CreateAnalysisMask.status",
                    "080_Quantification.status",
                    "090_VisualQC_ASL.status",
                    "999_ready.status"]

        # conditions is a list of tuples whose first element is a workload filename that may be impacted and whose
        # second element is a boolean that defines whether to remove it or not
        for condition in conditions:
            filename, to_remove = condition
            if to_remove:
                workload.remove(filename)

        # Must iterate through both the subject level listing AND the session level (ASL_1, ASL_2, etc.) listing
        status_files = []
        for subject in study_subjects:
            for run in session_names:
                try:
                    has_flair_img = next((analysis_directory / subject).glob("*FLAIR.nii*")).exists()
                except StopIteration:
                    has_flair_img = False
                try:
                    has_m0_img = next((analysis_directory / subject / run).glob("*M0.nii*")).exists()
                except StopIteration:
                    has_m0_img = False
                try:
                    has_asl_img = next((analysis_directory / subject / run).glob("*ASL*.nii*")).exists()
                except StopIteration:
                    has_asl_img = False

                if any([skip_if_no_m0 and not has_m0_img,
                        skip_if_no_asl and not has_asl_img,
                        skip_if_no_flair and not has_flair_img]):
                    continue

                directory = analysis_directory / "lock" / "xASL_module_ASL" / subject / f"xASL_module_ASL_{run}"
                # Filter out any anticipated status files that are already present in the lock dirs
                filtered_workload = [directory / name for name in workload if not (directory / name).exists()]
                status_files.extend(filtered_workload)
                # Calculate the numerical representation of the STATUS files workload
                num_repr = sum([workload_translator[stat_file.name] for stat_file in filtered_workload])
                aslmod_dict[subject][run] = num_repr

        return aslmod_dict, status_files

    def get_population_workload(analysis_directory, workload_translator):
        workload = ["010_CreatePopulationTemplates.status",
                    "020_CreateAnalysisMask.status",
                    "030_CreateBiasfield.status",
                    "040_GetDICOMStatistics.status",
                    "050_GetVolumeStatistics.status",
                    "060_GetMotionStatistics.status",
                    "065_GetRegistrationStatistics.status",
                    "070_GetROIstatistics.status",
                    "080_SortBySpatialCoV.status",
                    "090_DeleteAndZip.status",
                    "999_ready.status"]
        directory = analysis_directory / "lock" / "xASL_module_Population" / "xASL_module_Population"
        status_files = [directory / name for name in workload if not (directory / name).exists()]
        numerical_representation = sum([workload_translator[stat_file.name] for stat_file in status_files])
        return numerical_representation, status_files

    # Define the individual translators
    filename2workload = translators["ExploreASL_Filename2Workload"]

    # First get all the subjects
    analysis_dir = Path(parmsdict["D"]["ROOT"])
    regex: re.Pattern = re.compile(parmsdict["subject_regexp"])
    sess_names: List[str] = parmsdict["SESSIONS"]
    skipnoasl, skipnom0, skipnoflair = parmsdict["SkipIfNoASL"], parmsdict["SkipIfNoM0"], parmsdict["SkipIfNoFlair"]
    subjects: List[str] = [sub.name for sub in analysis_dir.iterdir() if
                           all([regex.search(str(sub.name)), sub.is_dir(), sub.name not in ["lock", "Population"]])]

    # Account for conditions that influence whether a .status file is to be removed from the expected workload or not
    asl_conditions = []
    # for statfile, corresponding_key, default in zip(["085_PVCorrection.status"],
    #                                                 ["bPVCorrectionNativeSpace"],
    #                                                 [True]):
    #     try:
    #         asl_conditions.append((statfile, not parmsdict[corresponding_key]))
    #     except KeyError:
    #         asl_conditions.append((statfile, default))

    # Use a dict to keep track of everything
    struct_dict = {subject: 0 for subject in subjects}
    asl_dict = {subject: {} for subject in subjects}
    # Update the dicts as appropriate
    if run_options == "Both":
        s_res = get_structural_workload(analysis_dir, subjects, struct_dict, skipnoasl, skipnom0, skipnoflair,
                                        workload_translator=filename2workload)
        struct_dict, struct_status = s_res
        a_res = get_asl_workload(analysis_dir, subjects, sess_names, asl_dict, skipnoasl, skipnom0, skipnoflair,
                                 workload_translator=filename2workload, conditions=asl_conditions)
        asl_dict, asl_status = a_res

        struct_totalworkload = sum(struct_dict.values())
        asl_totalworkload = {subject: sum(asl_dict[subject].values()) for subject in subjects}
        asl_totalworkload = sum(asl_totalworkload.values())
        print(f"Structural Calculated Workload: {struct_totalworkload}")
        print(f"ASL Calculated Workload: {asl_totalworkload}")
        # Return the numerical sum of the workload and the combined list of the expected status files
        return struct_totalworkload + asl_totalworkload, sorted(struct_status + asl_status)

    elif run_options == "ASL":
        a_res = get_asl_workload(analysis_dir, subjects, sess_names, asl_dict, skipnoasl, skipnom0, skipnoflair,
                                 workload_translator=filename2workload, conditions=asl_conditions)
        asl_dict, asl_status = a_res
        asl_totalworkload = {asl_subject: sum(asl_dict[asl_subject].values()) for asl_subject in subjects}
        asl_totalworkload = sum(asl_totalworkload.values())
        print(f"ASL Calculated Workload: {asl_totalworkload}")
        # Return the numerical sum of the workload and the list of expected status files
        return asl_totalworkload, asl_status

    elif run_options == "Structural":
        s_res = get_structural_workload(analysis_dir, subjects, struct_dict, skipnoasl, skipnom0, skipnoflair,
                                        workload_translator=filename2workload)
        struct_dict, struct_status = s_res
        struct_totalworkload = sum(struct_dict.values())
        print(f"Structural Calculated Workload: {struct_totalworkload}")
        # Return the numerical sum of the workload and the list of expected status files
        return struct_totalworkload, sorted(struct_status)

    elif run_options == "Population":
        pop_totalworkload, pop_status = get_population_workload(analysis_dir, workload_translator=filename2workload)
        print(f"Population Calculated Workload: {pop_totalworkload}")
        # Return the numerical sum of the workload and the list of expected status files
        return pop_totalworkload, pop_status

    else:
        print("THIS SHOULD NEVER PRINT AS YOU HAVE SELECTED AN IMPOSSIBLE WORKLOAD OPTION")


# Called after processing is done to compare the present status files against the files that were expected to be created
# at the time the run was initialized
def calculate_missing_STATUS(analysis_dir: Path, expected_status_files: List[Path]):
    postrun_status_files = set((analysis_dir / "lock").rglob("*.status"))
    incomplete = [file for file in expected_status_files if file not in postrun_status_files]
    if len(incomplete) == 0:
        return True, incomplete
    else:
        return False, incomplete


def interpret_statusfile_errors(incomplete_files: List[Path], translators: dict):
    """
    Interprets the step in the ExploreASL pipeline for which particular subjects/sessions/etc. must have failed and
    returns these interpretations as messages for each of the modules
    :param incomplete_files: the list of Path objects pointing to status files that were not generated in the pipeline
    :param translators: the translators (from JSON_LOGIC directory) used to convert filenames to their descriptions for
    generating the correct error message
    :return: 3 lists of interpreted error messages, one for each ExploreASL module
    """

    # Prepare containers and translators
    struct_dict = {}  # dict whose keys are subject names and whose values are Path objects to .status files
    asl_dict = {}  # dict whose keys are subject names and whose value are dicts of keys runs and values Path objects
    pop_list = []  # list whose elements are Path objects to .status files
    asl_msgs = []  # list whose elements are string messages pertaining to the ASL module
    struct_msgs = []  # list whose elements are string messages pertaining to the Structural module
    pop_msgs = []  # list whose elements are string messages pertaining to the Population
    stuct_status_file_translator = translators["Structural_Module_Filename2Description"]
    asl_status_file_translator = translators["ASL_Module_Filename2Description"]
    population_file_translator = translators["Population_Module_Filename2Description"]

    # Prepare regex detectors
    delimiter = "\\\\" if system() == "Windows" else "/"
    asl_regex = re.compile(f"(?:.*){delimiter}lock{delimiter}xASL_module_(?:Structural|ASL){delimiter}(.*){delimiter}"
                           f"xASL_module_(?:Structural|ASL)_?(.*)?{delimiter}(.*\\.status)")
    pop_regex = re.compile(f"(?:.*){delimiter}lock{delimiter}xASL_module_Population{delimiter}"
                           f"xASL_module_Population{delimiter}(.*\\.status)")

    # Use the regex to extract subject, session, filename fields from the status filepaths, then organize the status
    # file basenames into the appropriate dictionary structure
    incomplete_file: Path
    for incomplete_file in incomplete_files:
        asl_match = asl_regex.search(str(incomplete_file))
        pop_match = pop_regex.search(str(incomplete_file))
        if asl_match:
            subject, run, file_basename = asl_match.groups()
            # Structural
            if run == "":
                struct_dict.setdefault(subject, [])
                struct_dict[subject].append(file_basename)
            # ASL
            else:
                asl_dict.setdefault(subject, {})
                asl_dict[subject].setdefault(run, [])
                asl_dict[subject][run].append(file_basename)
        # Population
        elif pop_match:
            file_basename = pop_match.group(1)
            pop_list.append(file_basename)
        else:
            return

    # For each of the dictionaries corresponding to an ExploreASL module, sort the basenames and create the appropriate
    # error message for that subject/session
    if len(asl_dict) > 0:
        for subject, inner_run_dict in asl_dict.items():
            for run, file_basenames in inner_run_dict.items():
                sorted_basenames: List[str] = sorted(file_basenames)
                msg = f"Subject: {subject}; Run: {run}; Failed in the ASL module prior to: " \
                      f"{asl_status_file_translator[sorted_basenames[0]]}; " \
                      f"STATUS file failed to be created: {sorted_basenames[0]}"
                asl_msgs.append(msg)

    if len(struct_dict) > 0:
        for subject, file_basenames in struct_dict.items():
            sorted_basenames: List[str] = sorted(file_basenames)
            msg = f"Subject: {subject}; Failed in the Structural module prior to: " \
                  f"{stuct_status_file_translator[sorted_basenames[0]]}; " \
                  f"STATUS file failed to be created: {sorted_basenames[0]}"
            struct_msgs.append(msg)

    if len(pop_list) > 0:
        sorted_basenames: List[str] = sorted(pop_list)
        msg = f"Population module failed at: {population_file_translator[sorted_basenames[0]]}; " \
              f"STATUS file failed to be created: {sorted_basenames[0]}"
        pop_msgs = [msg]

    return struct_msgs, asl_msgs, pop_msgs
