import os
import numpy as np
import shutil
import subprocess
from glob import glob, iglob
import nibabel as nib
import pydicom
import json
import pandas as pd
from more_itertools import peekable
import struct
from ast import literal_eval
from nilearn import image

pd.set_option("display.width", 600)
pd.set_option("display.max_columns", 15)


def get_dicom_directories(config: dict):
    """
    Convenience function for globbing the dicom directories from the config file
    :param config: the configuration file that specifies the directory structure
    :return: dcm_firs: the list of filepaths to directories containing the dicom files
    """
    raw_dir = config["RawDir"]
    n_levels = ["*"] * len(config["Directory Structure"])
    dcm_dirs = glob(os.path.join(raw_dir, *n_levels))
    return dcm_dirs


def get_manufacturer(dcm_dir: str):
    """
    Returns the string suggesting which manufacturer the dicoms in the provided dicom directory belong to.
    :param dcm_dir: the absolute path to the directory containing the dicoms to be converted.
    :return: returns the string "Siemens", "Philips", or "GE". Otherwise returns None.
    """
    dcm_files = iglob(os.path.join(dcm_dir, "*.dcm"))
    dcm_files = peekable(dcm_files)
    if not dcm_files:
        return None

    dcm_data = pydicom.read_file(dcm_files[0])
    detected_manufac = []
    manufac_tags = [(0x0008, 0x0070), (0x0019, 0x0010)]
    for tag in manufac_tags:
        try:
            detected_manufac.append(f"{dcm_data[tag].value}".upper())
        except KeyError:
            detected_manufac.append("")

    if any(["SIEMENS" in result for result in detected_manufac]):
        return "Siemens"
    elif any(["PHILIPS" in result for result in detected_manufac]):
        return "Philips"
    elif any(["GE" in result for result in detected_manufac]):
        return "GE"
    else:
        return


def get_dicom_value(data: pydicom.Dataset, tags: list, default=None):
    """
    Convenience function for retrieving the value of a dicom tag. Otherwise, returns the indicated default
    :param data: the dicom data as a Pydicom Dataset
    :param tags: a list of tuples, each tuple usually being 2 elements (0x####, 0x####). 4 element tuples are used
    to delve into nested dicom structures
    :param default: the default value to return if nothing can be found
    :return: value: the first valid value associated with the tag
    """
    detected_values = []
    for tag in tags:
        # For tags that are nested
        if len(tag) == 4:
            try:
                first = (tag[0], tag[1])
                second = (tag[2], tag[3])
                detected_values.append(f"{data[first][0][second].value}")
            except (KeyError, TypeError):
                detected_values.append(None)
        # For base level tags
        else:
            try:
                detected_values.append(f"{data[tag].value}")
            except KeyError:
                detected_values.append(None)

    # Additional for loop for types
    types = [type(value) for value in detected_values]

    if str in types and float in types:
        idx = types.index(float)
        return detected_values[idx]

    if str in types and int in types:
        idx = types.index(int)
        return detected_values[idx]

    for value in detected_values:
        if value is not None:
            if isinstance(value, str):
                if value.isnumeric():
                    return float(value)
                elif value.startswith("b'") and value.endswith("'"):
                    return struct.unpack('f', literal_eval(value))[0]
                else:
                    return value

            return value

    return default


def get_additional_dicom_parms(dcm_dir: str, manufacturer: str):
    """
    Retrieves some additional important dicom headers that dcm2niix may not capture and which may be important for
    processing
    :param manufacturer: the string detailing which MRI scanner manufacturer corresponds to this DICOM directory
    :param dcm_dir: absolute path to the dicom directory, where the first dicom will be used to determine parameters
    :return: additional_dcm_info: a dict of the additional parameters as keys and their values as the dict values.
    """
    print(f"PROCESSING {dcm_dir}!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    dcm_files = iglob(os.path.join(dcm_dir, "*.dcm"))
    dcm_files = peekable(dcm_files)
    if not dcm_files:
        return None

    dcm_data = pydicom.read_file(next(dcm_files))

    if manufacturer == 'Philips':
        tags_dict = {
            "AcquisitionMatrix": [(0x0018, 0x1310)],
            "SoftwareVersions": [(0x0018, 0x1020)],
            "NumberOfAverages": [(0x0018, 0x0083)],
            "RescaleSlope": [(0x0028, 0x1053), (0x2005, 0x110A)],
            "RescaleIntercept": [(0x0028, 0x1052)],
            "MRScaleSlope": [(0x2005, 0x120E), (0x2005, 0x110E), (0x2005, 0x100E)],
            "RealWorldValueSlope": [(0x0040, 0x9096, 0x0040, 0x9225)],
        }
        defaults = [None, None, 1, 1, 0, 1, None, None]
    else:
        tags_dict = {
            "AcquisitionMatrix": [(0x0018, 0x1310)],
            "SoftwareVersions": [(0x0018, 0x1020)],
            "NumberOfAverages": [(0x0018, 0x0083)],
            "RescaleSlope": [(0x0028, 0x1053), (0x2005, 0x110A)],
            "RescaleIntercept": [(0x0028, 0x1052)]
        }
        defaults = [None, None, 1, 1, 0]

    additional_dcm_info = {}.fromkeys(tags_dict.keys())
    for (key, value), default in zip(tags_dict.items(),
                                     defaults  # Default values
                                     ):
        result = get_dicom_value(dcm_data, value, default)

        # Additional processing for specific keys
        if key == "AcquisitionMatrix" and result is not None:
            result = result.strip('[]').split(", ")
            result = [int(number) for number in result]

        # Convert any lingering strings to float
        if key in ["NumberOfAverages", "RescaleIntercept", "RescaleSlope", "MRScaleSlope", "RealWorldValueSlope"]:
            if result is not None and not isinstance(result, list):
                try:
                    result = float(result)
                except ValueError:
                    pass

        additional_dcm_info[key] = result

        # Final corrections for Philips scans in particular
        if manufacturer == "Philips":
            # First correction - disagreeing values between RescaleSlope and RealWorldValueSlope if they ended up in the
            # same dicom. Choose the small value of the two and set it for both
            if all([additional_dcm_info["RescaleSlope"] is not None,
                    additional_dcm_info["RealWorldValueSlope"] is not None,
                    additional_dcm_info["RescaleSlope"] != 1,
                    additional_dcm_info["RealWorldValueSlope"] != 1,
                    additional_dcm_info["RescaleSlope"] != additional_dcm_info["RealWorldValueSlope"]
                    ]):
                additional_dcm_info["RescaleSlope"] = min([additional_dcm_info["RescaleSlope"],
                                                           additional_dcm_info["RealWorldValueSlope"]])
                additional_dcm_info["RealWorldValueSlope"] = min([additional_dcm_info["RescaleSlope"],
                                                                  additional_dcm_info["RealWorldValueSlope"]])

            # Second correction - just to ease things on the side of ExploreASL; if RescaleSlope could not be determined
            # while "RealWorldValueSlope" could be, copy over the latter's value for the former
            if all([additional_dcm_info["RealWorldValueSlope"] is not None,
                    additional_dcm_info["RealWorldValueSlope"] != 1,
                    additional_dcm_info["RescaleSlope"] == 1]):
                additional_dcm_info["RescaleSlope"] = additional_dcm_info["RealWorldValueSlope"]

    # remove the "RealWorldValueSlope" as it is no longer needed
    try:
        del additional_dcm_info["RealWorldValueSlope"]
    except KeyError:
        pass

    return additional_dcm_info


def get_structure_components(dcm_dir: str, config: dict):
    """
    Returns the essential subject, session, and scan names of the currently-assessed dicom directory.
    :param dcm_dir: the absolute path to the directory containing the dicoms to be converted
    :param config: a dict with essentials details the raw directory structure and mappings of aliases
    :return: returns the current subject, session, and scan names in their alias form
    """
    subject, session, scan = None, None, None
    dirname = dcm_dir
    for dir_type in reversed(config["Directory Structure"]):
        if dir_type == "Subject":
            dirname, basename = os.path.split(dirname)
            subject = basename
        elif dir_type == "Session":
            dirname, basename = os.path.split(dirname)
            session = basename
        elif dir_type == "Scan":
            dirname, basename = os.path.split(dirname)
            scan = basename
        else:
            dirname, _ = os.path.split(dirname)

    scan_translator = {value: key for key, value in config["Scan Aliases"].items()}
    scan_dst_name = scan_translator[scan]

    if session is not None:
        session_translator = config["Ordered Session Aliases"]
        session_dst_name = session_translator[session]
    else:
        session_dst_name = None

    return subject, session_dst_name, scan_dst_name


def get_dst_dirname(raw_dir: str, subject: str, session: str, scan: str, legacy_mode: bool = False):
    """
    Creates the essential destination directory for nifti and json files to be created in during the conversion process
    :param raw_dir: the absolute path to the raw folder directory
    :param subject: the string representing the current subject
    :param session: the string representing the current session
    :param scan: the string representing the scan. Is either ASL4D, T1, M0, FLAIR, or WMH_SEGM
    :param legacy_mode: whether to format things in the old ExploreASL analysis directory structure or not
    :return: a string representation of the output directory for dcm2niix to create nifti files in
    """
    # NON-BIDS FORMAT
    if legacy_mode:
        try:
            analysis_dir = os.path.join(os.path.dirname(raw_dir), "analysis")
            if session is None:
                if scan not in ["T1", "FLAIR"]:
                    dst_dir = os.path.join(analysis_dir, subject, "ASL_1", "TEMP")
                else:
                    dst_dir = os.path.join(analysis_dir, subject, "TEMP")
            else:
                if scan not in ["T1", "FLAIR"]:
                    dst_dir = os.path.join(analysis_dir, subject, session, "TEMP")
                else:
                    dst_dir = os.path.join(analysis_dir, subject, "TEMP")
        except NotADirectoryError:
            return False, None

    # BIDS FORMAT
    else:
        try:
            analysis_dir = os.path.join(os.path.dirname(raw_dir), "analysis")
            # Get rid of illegal characters for subject
            if "-" in subject:
                subject = subject.replace("-", "")
            if "_" in subject:
                subject = subject.replace("_", "")

            if session is None:
                if scan not in ["T1", "FLAIR"]:
                    dst_dir = os.path.join(analysis_dir, f"sub-{subject}", "perf", "TEMP")
                else:
                    dst_dir = os.path.join(analysis_dir, f"sub-{subject}", "anat", "TEMP")
            else:
                # Get rid of illegal characters for session
                if "-" in session:
                    session = session.replace("-", "")
                if "_" in session:
                    session = session.replace("_", "")
                # Determine output directory based on the scan type
                if scan not in ["T1", "FLAIR"]:
                    dst_dir = os.path.join(analysis_dir, f"sub-{subject}", f"ses-{session}", "perf", "TEMP")
                else:
                    dst_dir = os.path.join(analysis_dir, f"sub-{subject}", f"ses-{session}", "anat", "TEMP")
        except NotADirectoryError:
            return False, None

    try:
        os.makedirs(dst_dir, exist_ok=True)
    except OSError:
        return False, None

    return True, dst_dir


def run_dcm2niix(temp_dir: str, dcm_dir: str, subject: str, session, scan: str):
    """
    Runs the dcm2niix program as a subprocess and generates the appropriate nifti and json files
    :param temp_dir: the TEMP dst where nifti and json files will be deposited
    :param dcm_dir: the directory where the dicoms of interest are being held
    :param subject: the string representing the current subject
    :param session: the string representing the current session alias
    :param scan: the string representing the scan. It is either ASL4D, T1, M0, FLAIR, or WMH_SEGM
    :return status: whether the operation was a success or not
    """

    print(f"Inside run_dcn2niix with args:\n"
          f"temp_dir: {temp_dir}\n"
          f"dcm_dir: {dcm_dir}\n"
          f"subject: {subject}\n"
          f"session: {session}\n"
          f"scan: {scan}\n")

    if session is None:
        output_filename_format = f"{subject}_{scan}_%s"
    else:
        output_filename_format = f"{subject}_{session}_{scan}_%s"

    # Must ensure that no files currently exist within the destination
    try:
        if len(os.listdir(temp_dir)) > 0:
            for path in os.listdir(temp_dir):
                fullpath = os.path.join(temp_dir, path)
                if os.path.isfile(fullpath):
                    os.remove(fullpath)
                elif os.path.isdir(fullpath):
                    os.removedirs(fullpath)
                else:
                    continue
    except FileNotFoundError:
        return False

    # Arguments used:
    # -b : Create a BIDS json sidecar; value y indicates yes
    # -z : Whether to compress the nifti or not; value n indicates no
    # -X : Whether to crop the nifti or not; value n indicates no
    # -t : Whether to include private patient details; value n indicates no
    # -m : Whether to merge 2D slices from the same series; value n indicates no
    # -f : The format of the output file
    # -o : The directory to output to
    # -s : Whether to export as "single file" mode; value n indicates no
    # -v : Whether to be very verbose; value n indicates no (although its still somewhat verbose...)
    # final argument is the source directory where the dicoms are located

    command = f"dcm2niix -b y -z n -x n -t n -m n -s n -v n " \
              f"-f {output_filename_format} " \
              f"-o {temp_dir} " \
              f"{dcm_dir}"

    return_code = subprocess.run(command.split(" "))
    print(f"Return code: {return_code.returncode}")
    if return_code.returncode == 0:
        return True
    else:
        return False


def clean_niftis_in_temp(temp_dir: str, add_parms: dict, subject: str, session: str, scan: str,
                         legacy_mode: bool = False):
    """
    Concatenates the niftis, deletes the previous ones, and moves the concatenated one out of the temp dir
    :param temp_dir: the absolute filepath to the TEMP directory where the niftis are present
    :param add_parms: additional parameters retrieved such as Acquisition Matrix used to deal with mosaics that did
    not get converted correctly
    :param subject: the string representing the current subject
    :param session: the string representing the current session alias
    :param scan: the string representing the scan. It is either ASL4D, T1, M0, FLAIR, or WMH_SEGM
    :param legacy_mode: whether to adjust file naming to the old ExploreASL format or not.
    :return: status: whether the operation was a success or not; the import summary parameters, and the filepath to the
    new nifti created
    """
    jsons = iglob(os.path.join(temp_dir, "*.json"))
    reorganized_niftis = []

    ###################################
    # PART 1 - ORGANIZING DCM2NII FILES
    ###################################
    # Must go over the jsons first to gain an understanding of the series
    json_data = {"SeriesNumber": {}, "AcquisitionTime": {}, "AcquisitionNumber": {}}
    for json_file in jsons:
        with open(json_file) as json_reader:
            sidecar_data = json.load(json_reader)
            for parm in json_data.keys():
                json_data[parm][json_file.replace(".json", ".nii")] = sidecar_data[parm]

    # Philips case; same series and probably same acq number too; opt for acq_time as the differentiating factor
    if len(set(json_data["SeriesNumber"].values())) == 1:
        reorganized_data = {key: value for key, value in sorted(json_data["AcquisitionTime"].items(),
                                                                key=lambda x: x[1])}
        reorganized_niftis = list(reorganized_data.keys())

    # Better Siemens scenario, usually has sequential series that increments
    elif len(set(json_data["SeriesNumber"].values())) > 1:
        reorganized_data = {key: value for key, value in sorted(json_data["SeriesNumber"].items(),
                                                                key=lambda x: x[1])}
        reorganized_niftis = list(reorganized_data.keys())

    # Hail Mary backup, maybe the acquisition number is by chance different
    elif len(set(json_data["AcquisitionNumber"].values())) > 1:
        reorganized_data = {key: value for key, value in sorted(json_data["AcquisitionNumber"].items(),
                                                                key=lambda x: x[1])}
        reorganized_niftis = list(reorganized_data.keys())

    import_summary = dict.fromkeys(["subject", "session", "scan", "filename",
                                    "dx", "dy", "dz", "nx", "ny", "nz", "nt"])

    # Get the acquisition matrix
    acq_matrix = add_parms["AcquisitionMatrix"]
    if acq_matrix[0] == 0:
        acq_rows = int(acq_matrix[1])
        acq_cols = int(acq_matrix[2])
    else:
        acq_rows = int(acq_matrix[0])
        acq_cols = int(acq_matrix[3])

    ########################################
    # PART 2 PROCESSING THE ORGANIZED NIFTIS
    ########################################
    # Must process niftis differently depending on the scan and the number present after conversion
    # Scenario: ASL4D
    if len(reorganized_niftis) > 1 and scan == "ASL4D":
        nii_objs = []
        for nifti in reorganized_niftis:
            nii_obj: nib.Nifti1Image = nib.load(nifti)

            # dcm2niix error: imports a 4D NIFTI instead of a 3D one. Solution: must be split first and concatenated
            # with the others at a later step
            if len(nii_obj.shape) == 4:
                volumes = nib.funcs.four_to_three(nii_obj)
                for volume in volumes:
                    nii_objs.append(volume)

            # Otherwise, correct 3D import
            elif len(nii_obj.shape) == 3:

                # dcm2niix error: imports a 3D mosaic. Solution: reformat as a 3D stack
                if nii_obj.shape[2] == 1:
                    nii_obj = fix_mosaic(mosaic_nifti=nii_obj, acq_dims=(acq_rows, acq_cols))
                nii_objs.append(nii_obj)
            else:
                return False, import_summary, None, None

        final_nifti_obj = nib.funcs.concat_images(nii_objs)

    # Scenario: multiple M0; will take their mean as final
    elif len(reorganized_niftis) > 1 and scan == "M0":
        nii_objs = [nib.load(nifti) for nifti in reorganized_niftis]
        final_nifti_obj = image.mean_img(nii_objs)

    # Scenario: single M0
    elif len(reorganized_niftis) == 1 and scan == "M0":
        final_nifti_obj = nib.load(reorganized_niftis[0])

    # Scenario: one of the structural types
    elif len(reorganized_niftis) == 1 and scan in ["T1", "FLAIR", "WMH_SEGM"]:
        final_nifti_obj = nib.load(reorganized_niftis[0])

    # Otherwise, something went wrong and the operation should stop
    else:
        return False, import_summary, None, None

    # Take the oppurtunity to get more givens for the import summary
    zooms = final_nifti_obj.header.get_zooms()
    shape = final_nifti_obj.shape
    import_summary["subject"], import_summary["session"], import_summary["scan"] = subject, session, scan
    import_summary["filename"] = scan + ".nii"

    if len(zooms) == 4:
        import_summary["dx"], import_summary["dy"], import_summary["dz"] = zooms[0:3]
    else:
        import_summary["dx"], import_summary["dy"], import_summary["dz"] = zooms

    if len(shape) == 4:
        import_summary["nx"], import_summary["ny"], \
        import_summary["nz"], import_summary["nt"] = shape
    else:
        import_summary["nx"], import_summary["ny"], \
        import_summary["nz"], import_summary["nt"] = shape[0], shape[1], shape[2], 1

    ################################
    # PART 3 NAMING AND MOVING FILES
    ################################
    if session is None:
        session = ""
    else:
        # Get rid of illegal characters for session
        if "-" in session:
            session = session.replace("-", "")
        if "_" in session:
            session = session.replace("_", "")
        session = f"ses-{session}_"

    if legacy_mode:
        final_nifti_filename = os.path.join(os.path.dirname(temp_dir), f"{scan}.nii")
        final_json_filename = os.path.join(os.path.dirname(temp_dir), f"{scan}.json")
    else:
        scan = {"ASL4D": "asl", "M0": "m0scan", "T1": "T1w", "FLAIR": "FLAIR"}[scan]
        # Get rid of illegal characters for subject
        if "-" in subject:
            subject = subject.replace("-", "")
        if "_" in subject:
            subject = subject.replace("_", "")

        final_nifti_filename = os.path.join(os.path.dirname(temp_dir), f"sub-{subject}_{session}{scan}.nii")
        final_json_filename = os.path.join(os.path.dirname(temp_dir), f"sub-{subject}_{session}{scan}.json")

    nib.save(final_nifti_obj, final_nifti_filename)

    jsons = iglob(os.path.join(temp_dir, "*.json"))
    json_file = next(jsons)

    # Json can bring up errors; os.replace should handle the rename if os.rename fails
    try:
        os.rename(json_file, final_json_filename)
    except FileExistsError:
        os.replace(json_file, final_json_filename)

    return True, import_summary, final_nifti_filename, final_json_filename


def update_json_sidecar(json_file: str, dcm_parms: dict):
    """
    Updates the json sidecar to include additional DICOM parameters not extracted by dcm2niix
    :param json_file: the absolute filepath to the json sidecar
    :param dcm_parms: a dict of the additional extracted DICOM headers that should be added to the json
    :return: 2 items: whether the operation was a success; the updated parameters of the json file (to be used in the
    summary file generation)
    """
    try:
        with open(json_file) as json_sidecar_reader:
            json_sidecar_parms: dict = json.load(json_sidecar_reader)

        json_sidecar_parms.update(dcm_parms)

        # First, rename certain elements
        for old_name, new_name in {"EstimatedEffectiveEchoSpacing": "EffectiveEchoSpacing",
                                   "EstimatedTotalReadoutTime": "TotalReadoutTime"}.items():
            if old_name in json_sidecar_parms.keys():
                json_sidecar_parms[new_name] = json_sidecar_parms.pop(old_name)

        # Next, take care of Philips keys generated in dcm2niix
        if "PhilipsRescaleSlope" in list(json_sidecar_parms.keys()) and \
                "RescaleSlope" in list(json_sidecar_parms.keys()):
            if all([json_sidecar_parms["RescaleSlope"] == 1,
                    isinstance(json_sidecar_parms["PhilipsRescaleSlope"], (float, int)),
                    json_sidecar_parms["PhilipsRescaleSlope"] != 1
                    ]):
                json_sidecar_parms["RescaleSlope"] = json_sidecar_parms.pop("PhilipsRescaleSlope")
            else:
                del json_sidecar_parms["PhilipsRescaleSlope"]

        if "PhilipsRescaleIntercept" in json_sidecar_parms.keys() and "RescaleIntercept" in json_sidecar_parms.keys():
            if all([json_sidecar_parms["RescaleIntercept"] == 0,
                    isinstance(json_sidecar_parms["PhilipsRescaleIntercept"], (float, int)),
                    json_sidecar_parms["PhilipsRescaleIntercept"] != 1
                    ]):
                json_sidecar_parms["RescaleIntercept"] = json_sidecar_parms.pop("PhilipsRescaleIntercept")
            else:
                del json_sidecar_parms["PhilipsRescaleIntercept"]

        if "PhilipsScaleSlope" in json_sidecar_parms.keys() and "MRScaleSlope" in json_sidecar_parms.keys():
            if all([json_sidecar_parms["MRScaleSlope"] == 1,
                    isinstance(json_sidecar_parms["PhilipsScaleSlope"], (float, int)),
                    json_sidecar_parms["PhilipsScaleSlope"] != 1
                    ]):
                json_sidecar_parms["MRScaleSlope"] = json_sidecar_parms.pop("PhilipsScaleSlope")
            else:
                del json_sidecar_parms["PhilipsScaleSlope"]

        with open(json_file, 'w') as w:
            json.dump(json_sidecar_parms, w, indent=3)
    except FileNotFoundError:
        return False, None

    return True, json_sidecar_parms


def asldcm2bids_onedir(dcm_dir: str, config: dict, legacy_mode: bool = False):
    """
    The main function for most ASL-processing steps centered around processing a single dicom directory immediately
    preceding the dicom files.
    :param dcm_dir: the absolute path to the dicom directory
    :param config: the configuration file containing the essential parameters for importing
    :param legacy_mode: whether to adjust file naming to the old ExploreASL format or not.
    :returns: 3 items: whether the operation was a success; the most recent step performed, and the import summary (or
    None if the import for this directory failed)
    """
    # Get the subject, session, and scan for that directory
    subject, session, scan = get_structure_components(dcm_dir=dcm_dir, config=config)

    # Get the manufacturer tag
    manufacturer = get_manufacturer(dcm_dir=dcm_dir)

    # Retrieve the additional DICOM parameters and include the Philips rescale slope indicator
    addtional_dcm_parameters = get_additional_dicom_parms(dcm_dir=dcm_dir, manufacturer=manufacturer)
    # if manufacturer == "Philips":
    #     addtional_dcm_parameters["UsePhilipsFloatNotDisplayScaling"] = 1

    # Generate the directories for dumping dcm2niix output
    successful_run, temp_dst_dir = get_dst_dirname(raw_dir=config["RawDir"],
                                                   subject=subject,
                                                   session=session,
                                                   scan=scan,
                                                   legacy_mode=legacy_mode)
    if not successful_run:
        print(f"FAILURE ENCOUNTERED AT THE DCM2NIIX STEP")
        return False, f"SUBJECT: {subject}; " \
                      f"SCAN: {scan}; " \
                      f"ERROR: Failed at Temp folder generation", None

    # Run the main program
    successful_run = run_dcm2niix(temp_dir=temp_dst_dir, dcm_dir=dcm_dir,
                                  subject=subject, session=session, scan=scan)
    if not successful_run:
        print(f"FAILURE ENCOUNTERED AT THE DCM2NIIX STEP")
        return False, f"SUBJECT: {subject}; " \
                      f"SCAN: {scan}; " \
                      f"ERROR: Failed at DCM2NIIX conversion", None

    # Clean the niftis in the TEMP directory
    successful_run, \
    nifti_parmameters, \
    nifti_filepath, \
    json_filepath = clean_niftis_in_temp(temp_dir=temp_dst_dir,
                                         add_parms=addtional_dcm_parameters,
                                         subject=subject,
                                         session=session,
                                         scan=scan,
                                         legacy_mode=legacy_mode)
    if not successful_run:
        print(f"FAILURE ENCOUNTERED AT CLEANING THE NIFTIs IN THE TEMP FOLDER")
        return False, f"SUBJECT: {subject}; " \
                      f"SCAN: {scan}; " \
                      f"ERROR: Failed at post-conversion nifti cleanup", None

    # For ASL4D and M0 scans, the JSON sidecar from dcm2niix must include additional parameters
    successful_run, json_parameters = update_json_sidecar(json_file=json_filepath, dcm_parms=addtional_dcm_parameters)
    if not successful_run:
        print(f"FAILURE ENCOUNTERED AT CORRECTION THE JSON SIDECAR STEP")
        return False, f"SUBJECT: {subject}; " \
                      f"SCAN: {scan}; " \
                      f"ERROR: Failed at updating json sidecar with additional parms", None

    # If everything was a success, prepare the dict for creating the import summary
    import_summary = {}
    import_summary.update(addtional_dcm_parameters)
    import_summary.update(nifti_parmameters)
    import_summary.update(json_parameters)

    # Finally, delete the TEMP folder
    try:
        shutil.rmtree(temp_dst_dir, ignore_errors=True)
        return True, f"SUBJECT: {subject}; " \
                     f"SCAN: {scan}; " \
                     f"Message: Successful conversion", import_summary
    except FileNotFoundError:
        return True, f"SUBJECT: {subject}; " \
                     f"SCAN: {scan}; " \
                     f"ERROR: Failed at deleting the TEMP folder", import_summary


def create_import_summary(import_summaries: list, config: dict):
    """
    Given a list of individual summaries of each subject/session/scan, this function will bring all those givens
    together into a single dataframe for easy viewing
    :param import_summaries: a list of dicts, with each dict being the parameters of that subject-session-scan
    :param config: the import configuration file generated by the GUI to help locate the analysis directory
    """
    analysis_dir = os.path.join(os.path.dirname(config["RawDir"]), "analysis")
    df = pd.concat([pd.Series(import_summary) for import_summary in import_summaries], axis=1, sort=True).T
    df["dt"] = df["RepetitionTime"]
    appropriate_ordering = ['subject', 'session', 'scan', 'dx', 'dy', 'dz', 'dt', 'nx', 'ny', 'nz', 'nt',
                            "RepetitionTime", "EchoTime", "NumberOfAverages", "RescaleSlope", "RescaleIntercept",
                            "MRScaleSlope", "AcquisitionTime",
                            "AcquisitionMatrix", "TotalReadoutTime", "EffectiveEchoSpacing"]
    df = df.reindex(columns=appropriate_ordering)
    df = df.sort_values(by=["scan", "subject"]).reset_index(drop=True)

    print(df)

    try:
        df.to_csv(os.path.join(analysis_dir, "import_summary.tsv"), sep='\t', index=False, na_rep='n/a')
    except PermissionError:
        df.to_csv(os.path.join(analysis_dir, "import_summary_backup.tsv"), sep='\t', index=False, na_rep='n/a')


def bids_m0_followup(analysis_dir):
    """
    In a BIDS import, this function will run through the imported dataset and adjust any BIDS-standard fields that
    should be present in the m0scan.json sidecar, such as "IntendedFor"
    :param analysis_dir: the absolute path to the analysis directory
    """
    m0_jsons = iglob(os.path.join(analysis_dir, "**", "*_m0scan.json"), recursive=True)
    for m0_json in m0_jsons:
        asl_json = m0_json.replace("_m0scan.json", "_asl.json")
        asl_nifti = m0_json.replace("_m0scan.json", "_asl.nii")

        # Ensure that the asl json sidecar and nifti images actually exist adjacent to the m0scan.json
        if os.path.exists(asl_json) and os.path.exists(asl_nifti):

            # BIDS standard: the "IntendedFor" filepath must be relative to the subject and contain forward slashes
            truncated_asl_json = asl_nifti.replace(analysis_dir, "")
            if '\\' in truncated_asl_json:
                truncated_asl_json = truncated_asl_json.replace("\\", "/")

            with open(m0_json) as m0_json_reader:
                m0_parms = json.load(m0_json_reader)
            m0_parms["IntendedFor"] = truncated_asl_json
            with open(m0_json, 'w') as m0_json_writer:
                json.dump(m0_parms, m0_json_writer, indent=3)


def fix_mosaic(mosaic_nifti: nib.Nifti1Image, acq_dims: tuple):
    """
    Fixes incorrectly-processed NIFTIs by dcm2niix where they still remain mosaics due to a lack of
    NumberOfImagesInMosaic header. This function implements a hack to
    :param mosaic_nifti: the nifti image object that needs to be fixed. Should be of shape m x n x 1
    the sliding window algorithm
    :param acq_dims: the (row, col) acquisition dimensions for rows and columns from the AcquisitionMatrix DICOM field.
    Used to determine the appropriate kernel size to use for the sliding window algorithm
    :return: new_nifti; a 3D NIFTI that is no longer mosaic
    """
    acq_rows, acq_cols = acq_dims

    # Get the shape and array values of the mosaic (flatten the latter into a 2D array)
    img_shape = mosaic_nifti.shape
    # noinspection PyTypeChecker
    img_data = np.rot90(np.squeeze(mosaic_nifti.get_fdata()))

    # If this is a square, and the rows perfectly divides the mosaic
    if img_shape[0] == img_shape[1] and img_shape[0] % acq_rows == 0:
        nsplits_w, nsplits_h = img_shape[0] / acq_rows, img_shape[0] / acq_rows
        kernel_w, kernel_h = acq_rows, acq_rows
    # If this is a square, and the cols perfectly divides the mosaic
    elif img_shape[0] == img_shape[1] and img_shape[0] % acq_cols == 0:
        nsplits_w, nsplits_h = img_shape[0] / acq_cols, img_shape[0] / acq_cols
        kernel_w, kernel_h = acq_cols, acq_cols
    # If this is a rectangle
    elif all([img_shape[0] != img_shape[1],
              img_shape[0] % acq_rows == 0,
              img_shape[1] % acq_cols == 0
              ]):
        nsplits_w, nsplits_h = img_shape[0] / acq_rows, img_shape[1] / acq_cols
        kernel_w, kernel_h = acq_rows, acq_cols
    else:
        return

    # Initialize the data that will house the split mosaic into slices
    new_img_data = np.zeros(shape=(kernel_w, kernel_h, int(nsplits_w * nsplits_h)))
    slice_num = 0

    # Sliding Window algorithm
    for ii in range(int(nsplits_w)):
        for jj in range(int(nsplits_h)):
            x_start, x_end = ii * kernel_w, (ii + 1) * kernel_w
            y_start, y_end = jj * kernel_h, (jj + 1) * kernel_h
            img_slice = img_data[x_start:x_end, y_start:y_end]

            # Disregard slices that are only zeros
            if np.nanmax(img_slice) == 0:
                continue
            # Otherwise update the zeros array at the appropriate slice with the new values
            else:
                new_img_data[:, :, slice_num] = img_slice
                slice_num += 1

    # Filter off slices that had only zeros
    new_img_data = np.rot90(new_img_data[:, :, 0:slice_num], 3)
    new_nifti = image.new_img_like(mosaic_nifti, new_img_data, affine=mosaic_nifti.affine)
    return new_nifti
