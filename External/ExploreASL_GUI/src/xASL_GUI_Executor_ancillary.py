from pathlib import Path
from more_itertools import peekable
import re
from platform import system
from itertools import chain
from pprint import pprint
from typing import List, Set, Tuple, Union


def is_earlier_version(easl_dir: Union[Path, str], threshold_higher: int = 140, threshold_lower: int = 0):
    """
    Helper function to determine whether a given ExploreASL directory is between some set of integer thresholds
    representing the versions
    """
    ver_regex = re.compile(r"VERSION_(.*)")
    easl_dir = Path(easl_dir).resolve()
    try:
        ver_file = next(easl_dir.rglob("VERSION_*"))
    except StopIteration:
        return True
    ver_str = ver_regex.search(str(ver_file)).group(1).replace(".", "")
    if threshold_lower < int(ver_str) < threshold_higher:
        return True
    else:
        return False


def is_valid_for_analysis(path: Path, parms: dict, glob_dict: dict):
    """
    Helper function. Given a subject or session, the parameters from DataPar.json, and a dict of glob_patterns to use,
    determine whether this path should be skipped.
    """
    try:
        has_flair_img = next(path.glob(glob_dict["FLAIR"])).exists()
    except StopIteration:
        has_flair_img = False
    try:
        has_m0_img = next(path.glob(glob_dict["M0"])).exists()
    except StopIteration:
        has_m0_img = False
    try:
        has_asl_img = next(path.glob(glob_dict["ASL"])).exists()
    except StopIteration:
        has_asl_img = False

    if any([parms["SkipIfNoM0"] and not has_m0_img,
            parms["SkipIfNoASL"] and not has_asl_img,
            parms["SkipIfNoFlair"] and not has_flair_img]):
        return False
    return True


def calculate_anticipated_workload(parmsdict, run_options, translators):
    """
    Convenience function for calculating the anticipated workload
    :param parmsdict: the parameter file of the study; given parameters such as the regex are used from this
    :param run_options: "Structural", "ASL", "Both" or "Population"; which module is being run
    :param translators: The ExecutorTranslators, primarily for calculating the workload
    :return: workload; a numerical representation of the cumulative value of all status files made; these will be
    used to determine the appropriate maximum value for the progressbar
    """

    def get_structural_workload(analysis_directory: Path, parms: dict, incl_regex: re.Pattern,
                                workload_translator: dict):
        structuralmod_dict = {}
        status_files = []
        workload = {"010_LinearReg_T1w2MNI.status", "020_LinearReg_FLAIR2T1w.status",
                    "030_FLAIR_BiasfieldCorrection.status", "040_LST_Segment_FLAIR_WMH.status",
                    "050_LST_T1w_LesionFilling_WMH.status", "060_Segment_T1w.status", "070_CleanUpWMH_SEGM.status",
                    "080_Resample2StandardSpace.status", "090_GetVolumetrics.status", "100_VisualQC_Structural.status",
                    "999_ready.status"}
        glob_dictionary = {"ASL": "*/*ASL*.nii*", "FLAIR": "*FLAIR.nii*", "M0": "*/*M0.nii*"}

        for subject_path in analysis_directory.iterdir():
            # Disregard files, standard directories, subjects that fail regex, and subjects that are to be excluded
            if any([subject_path.is_file(), subject_path.name in ["Population", "lock"],
                    subject_path.name in parms["exclusion"], not incl_regex.search(subject_path.name)]):
                continue
            # Account for SkipIfNo flags
            if not is_valid_for_analysis(path=subject_path, parms=parms, glob_dict=glob_dictionary):
                continue

            # Account for version 1.2.1 and earlier
            has_flair = len(list(subject_path.glob(glob_dictionary["FLAIR"]))) > 0
            if is_earlier_version(parms["MyPath"], threshold_higher=130) and not has_flair:
                workload = {"010_LinearReg_T1w2MNI.status", "060_Segment_T1w.status",
                            "080_Resample2StandardSpace.status", "090_GetVolumetrics.status",
                            "100_VisualQC_Structural.status", "999_ready.status"}

            lock_dir: Path = analysis_directory / "lock" / "xASL_module_Structural" / subject_path.name / \
                             "xASL_module_Structural"
            if not lock_dir.exists():
                lock_dir.mkdir(parents=True)
            filtered_workload = [lock_dir / name for name in workload if not (lock_dir / name).exists()]
            # Filter out any anticipated status files that are already present in the lock dirs
            status_files.extend(filtered_workload.copy())
            num_repr = sum([workload_translator[stat_file.name] for stat_file in filtered_workload])
            structuralmod_dict[subject_path.name] = num_repr

        return structuralmod_dict, status_files

    def get_asl_workload(analysis_directory, parms: dict, workload_translator: dict, incl_regex: re.Pattern,
                         conditions: List[Tuple[str, bool]] = None):
        aslmod_dict = {}
        status_files = []
        glob_dictionary = {"ASL": "*ASL*.nii*", "FLAIR": "*FLAIR.nii*", "M0": "*M0.nii*"}
        # OLD EXPECTATION
        if is_earlier_version(parms["MyPath"], 140):
            workload = {"020_RealignASL.status", "030_RegisterASL.status", "040_ResampleASL.status",
                        "050_PreparePV.status", "060_ProcessM0.status", "070_Quantification.status",
                        "080_CreateAnalysisMask.status", "090_VisualQC_ASL.status", "999_ready.status"}
        # NEW EXPECTATION
        else:
            workload = {"020_RealignASL.status", "030_RegisterASL.status", "040_ResampleASL.status",
                        "050_PreparePV.status", "060_ProcessM0.status", "070_CreateAnalysisMask.status",
                        "080_Quantification.status", "090_VisualQC_ASL.status", "999_ready.status"}

        # conditions is a list of tuples whose first element is a workload filename that may be impacted and whose
        # second element is a boolean that defines whether to remove it or not
        if conditions is None:
            conditions = []
        for condition in conditions:
            filename, to_remove = condition
            if to_remove:
                workload.remove(filename)

        # Must iterate through both the subject level listing AND the session level (ASL_1, ASL_2, etc.) listing
        for subject_path in analysis_directory.iterdir():
            # Disregard files, standard directories, subjects that fail regex and subjects that are to be excluded
            if any([subject_path.is_file(), subject_path.name in ["Population", "lock"],
                    subject_path.name in parms["exclusion"], not incl_regex.search(subject_path.name)]):
                continue
            aslmod_dict[subject_path.name] = {}
            for run_path in subject_path.iterdir():
                if run_path.is_file():  # This is kept separate since many files are expected
                    continue
                if not is_valid_for_analysis(path=run_path, parms=parms, glob_dict=glob_dictionary):
                    continue

                # Deduce the lock dir path and make it if it doesn't exist
                lock_dir: Path = analysis_directory / "lock" / "xASL_module_ASL" / subject_path.name / \
                                 f"xASL_module_ASL_{run_path.name}"
                if not lock_dir.exists():
                    lock_dir.mkdir(parents=True)

                # Filter out any anticipated status files that are already present in the lock dirs
                filtered_workload = [lock_dir / name for name in workload if not (lock_dir / name).exists()]
                status_files.extend(filtered_workload)
                # Calculate the numerical representation of the STATUS files workload
                num_repr = sum([workload_translator[stat_file.name] for stat_file in filtered_workload])
                aslmod_dict[subject_path.name][run_path.name] = num_repr

        return aslmod_dict, status_files

    def get_population_workload(analysis_directory, workload_translator):
        workload = {"010_CreatePopulationTemplates.status", "020_CreateAnalysisMask.status",
                    "030_CreateBiasfield.status", "040_GetDICOMStatistics.status", "050_GetVolumeStatistics.status",
                    "060_GetMotionStatistics.status", "065_GetRegistrationStatistics.status",
                    "070_GetROIstatistics.status", "080_SortBySpatialCoV.status", "090_DeleteAndZip.status",
                    "999_ready.status"}
        directory = analysis_directory / "lock" / "xASL_module_Population" / "xASL_module_Population"
        if not directory.exists():
            directory.mkdir(parents=True)
        status_files = [directory / name for name in workload if not (directory / name).exists()]
        numerical_representation = sum([workload_translator[stat_file.name] for stat_file in status_files])
        return numerical_representation, status_files

    # Define the individual translators and analysis directory
    filename2workload = translators["ExploreASL_Filename2Workload"]
    analysis_dir = Path(parmsdict["D"]["ROOT"])
    subject_regex = re.compile(parmsdict["subject_regexp"])

    # Account for conditions that influence whether a .status file is to be removed from the expected workload or not
    asl_conditions = []
    # for statfile, corresponding_key, default in zip(["085_PVCorrection.status"],
    #                                                 ["bPVCorrectionNativeSpace"],
    #                                                 [True]):
    #     try:
    #         asl_conditions.append((statfile, not parmsdict[corresponding_key]))
    #     except KeyError:
    #         asl_conditions.append((statfile, default))

    # Update the dicts as appropriate
    if run_options == "Both":
        s_res = get_structural_workload(analysis_dir, parms=parmsdict, workload_translator=filename2workload,
                                        incl_regex=subject_regex)
        struct_dict, struct_status = s_res
        a_res = get_asl_workload(analysis_dir, parms=parmsdict, workload_translator=filename2workload,
                                 conditions=asl_conditions, incl_regex=subject_regex)
        asl_dict, asl_status = a_res

        struct_totalworkload = sum(struct_dict.values())
        asl_totalworkload = sum([sum(subject_dict.values()) for subject_dict in asl_dict.values()])

        print(f"Structural Calculated Workload: {struct_totalworkload}")
        print(f"ASL Calculated Workload: {asl_totalworkload}")
        # Return the numerical sum of the workload and the combined list of the expected status files
        return struct_totalworkload + asl_totalworkload, sorted(struct_status + asl_status)

    elif run_options == "ASL":
        a_res = get_asl_workload(analysis_dir, parms=parmsdict, workload_translator=filename2workload,
                                 conditions=asl_conditions, incl_regex=subject_regex)
        asl_dict, asl_status = a_res
        asl_totalworkload = sum([sum(subject_dict.values()) for subject_dict in asl_dict.values()])
        print(f"ASL Calculated Workload: {asl_totalworkload}")
        # Return the numerical sum of the workload and the list of expected status files
        return asl_totalworkload, asl_status

    elif run_options == "Structural":
        s_res = get_structural_workload(analysis_dir, parms=parmsdict, workload_translator=filename2workload,
                                        incl_regex=subject_regex)
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
