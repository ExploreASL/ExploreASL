import numpy as np
import shutil
import subprocess
import nibabel as nib
import pydicom
from pydicom.errors import InvalidDicomError
import json
import pandas as pd
from more_itertools import peekable
import struct
from ast import literal_eval
from nilearn import image
from platform import system
from pathlib import Path
from typing import Union, List
from datetime import datetime

pd.set_option("display.width", 600)
pd.set_option("display.max_columns", 15)


def get_dicom_directories(config: dict) -> List[Path]:
    """
    Convenience function for globbing the dicom directories from the config file
    :param config: the configuration file that specifies the directory structure
    :return: dcm_firs: the list of filepaths to directories containing the dicom files
    """
    raw_dir = Path(config["RawDir"])
    n_levels = ["*"] * len(config["Directory Structure"])
    delimiter = "\\" if system() == "Windows" else "/"
    # dcm_dirs = [directory for directory in raw_dir.glob(delimiter.join(n_levels)) if directory.is_dir()]
    dcm_dirs = []
    directory: Path
    for directory in raw_dir.glob(delimiter.join(n_levels)):
        hit = set(config["Scan Aliases"].values()).intersection(directory.parts)
        if len(hit) > 0:
            dcm_dirs.append(directory)

    return dcm_dirs


def get_structure_components(dcm_dir: Path, config: dict):
    """
    Returns the essential subject, visit, run, and scan names of the currently-assessed dicom directory.
    :param dcm_dir: the absolute path to the directory containing the dicoms to be converted
    :param config: a dict with essentials details the raw directory structure and mappings of aliases
    :return: returns the current subject, visit, run, and scan names in their alias form
    """
    subject, visit, run, scan = None, None, None, None
    for path_partname, dir_type in zip(reversed(dcm_dir.parts), reversed(config["Directory Structure"])):
        if dir_type == "Subject":
            subject = path_partname
        elif dir_type == "Visit":
            visit = path_partname
        elif dir_type == "Run":
            run = path_partname
        elif dir_type == "Scan":
            scan = path_partname
        else:
            pass

    try:
        scan_translator = {value: key for key, value in config["Scan Aliases"].items()}
        scan_dst_name = scan_translator[scan]
    except KeyError:
        return False, subject, visit, run, scan

    if run is not None:
        run_dst_name = config["Ordered Run Aliases"][run]
    else:
        run_dst_name = None

    return True, subject, visit, run_dst_name, scan_dst_name


def get_manufacturer(dcm_dir: Path) -> (bool, str):
    """
    Returns the string suggesting which manufacturer the dicoms in the provided dicom directory belong to.
    :param dcm_dir: the absolute path to the directory containing the dicoms to be converted.
    :return: returns the string "Siemens", "Philips", or "GE". Otherwise returns None.
    """
    dcm_files = peekable(dcm_dir.glob("*"))
    if not dcm_files:
        return False, None

    dcm_data = None
    for dcm_file in dcm_files:
        if dcm_file.is_dir():
            continue
        try:
            dcm_data = pydicom.read_file(dcm_file)
            break
        except InvalidDicomError:
            print(f"DICOM dir bad bad file: {dcm_file}")
            continue
        except PermissionError:
            continue

    if dcm_data is None:
        return False, None

    detected_manufac = []
    manufac_tags = [(0x0008, 0x0070), (0x0019, 0x0010)]
    for tag in manufac_tags:
        try:
            detected_manufac.append(f"{dcm_data[tag].value}".upper())
        except KeyError:
            detected_manufac.append("")

    if any(["SIEMENS" in result.upper() for result in detected_manufac]):
        return True, "Siemens"
    elif any(["PHILIPS" in result.upper() for result in detected_manufac]):
        return True, "Philips"
    elif any(["GE" in result.upper() for result in detected_manufac]):
        return True, "GE"
    else:
        return False, None


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


def get_additional_dicom_parms(dcm_dir: Path) -> (bool, dict):
    """
    Retrieves some additional important dicom headers that dcm2niix may not capture and which may be important for
    processing
    :param dcm_dir: absolute path to the dicom directory, where the first dicom will be used to determine parameters
    :return: additional_dcm_info: a dict of the additional parameters as keys and their values as the dict values.
    """
    dcm_files = peekable(dcm_dir.glob("*"))
    if not dcm_files:
        return False, None

    dcm_data = None
    for dcm_file in dcm_files:
        try:
            dcm_data = pydicom.read_file(dcm_file)
            break
        except InvalidDicomError:
            continue
        except PermissionError:
            continue

    if dcm_data is None:
        return False, None

    manufacturer = get_dicom_value(data=dcm_data, tags=[(0x0008, 0x0070), (0x0019, 0x0010)], default=None)
    if manufacturer is None:
        return False, None
    else:
        if "SIEMENS" in manufacturer.upper():
            manufacturer = "Siemens"
        elif "PHILIPS" in manufacturer.upper():
            manufacturer = "Philips"
        elif "GE" in manufacturer.upper():
            manufacturer = "GE"
        else:
            return False, None

    if manufacturer == 'Philips':
        tags_dict = {
            "AcquisitionMatrix": [(0x0018, 0x1310)],
            "SoftwareVersions": [(0x0018, 0x1020)],
            "NumberOfAverages": [(0x0018, 0x0083)],
            "RescaleSlope": [(0x0028, 0x1053), (0x2005, 0x110A), (0x2005, 0x140A)],
            "RescaleIntercept": [(0x0028, 0x1052)],
            "MRScaleSlope": [(0x2005, 0x120E), (0x2005, 0x110E), (0x2005, 0x100E)],
            "RealWorldValueSlope": [(0x0040, 0x9096, 0x0040, 0x9225)],
            "NumberOfSlices": [(0x0054, 0x0081)],
            "AcquisitionTime": [(0x0008, 0x0032)]
        }
        defaults = [None,  # AcquisitionMatrix
                    None,  # SoftwareVersions
                    1,  # NumberOfAverages
                    1,  # RescaleSlope
                    0,  # RescaleIntercept
                    1,  # MRScaleSlope
                    None,  # RealWorldValueSlope
                    None,  # NumberOfSlices
                    0  # AcquisitionTime
                    ]
    else:
        tags_dict = {
            "AcquisitionMatrix": [(0x0018, 0x1310)],
            "SoftwareVersions": [(0x0018, 0x1020)],
            "NumberOfAverages": [(0x0018, 0x0083)],
            "RescaleSlope": [(0x0028, 0x1053), (0x2005, 0x110A)],
            "RescaleIntercept": [(0x0028, 0x1052)],
            "NumberOfSlices": [(0x0054, 0x0081)],
            "AcquisitionTime": [(0x0008, 0x0032)]
        }
        defaults = [None,  # AcquisitionMatrix
                    None,  # SoftwareVersions
                    1,  # NumberOfAverages
                    1,  # RescaleSlope
                    0,  # RescaleIntercept
                    None,  # NumberOfSlices
                    0  # AcquisitionTime
                    ]

    additional_dcm_info = {}.fromkeys(tags_dict.keys())
    additional_dcm_info["Manufacturer"] = manufacturer
    for (key, value), default in zip(tags_dict.items(), defaults):
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

    return True, additional_dcm_info


def get_tempdst_dirname(raw_dir: Path, subject: str, visit: str, run: str, scan: str, legacy_mode: bool = False):
    """
    Creates the essential TEMP destination directory for nifti and json files to be created in during the conversion
    process
    :param raw_dir: the absolute path to the raw folder directory
    :param subject: the string representing the current subject
    :param visit: the string representing the current visit
    :param run: the string representing the current run
    :param scan: the string representing the scan. Is either ASL4D, T1, M0, FLAIR, or WMH_SEGM
    :param legacy_mode: whether to format things in the old ExploreASL analysis directory structure or not
    :return: a string representation of the output directory for dcm2niix to create nifti files in
    """
    # NON-BIDS FORMAT
    # For legacy import, this will be in the form raw/subject[_visit]/run
    if legacy_mode:
        try:
            if visit is None:
                visit = ""
            else:
                visit = f"_{visit}"

            analysis_dir = raw_dir.parent / "analysis"
            if run is None:
                run = "ASL_1"
            if scan not in ["T1", "FLAIR"]:
                dst_dir = analysis_dir / f"{subject}{visit}" / run / "TEMP"
            else:
                dst_dir = analysis_dir / f"{subject}{visit}" / "TEMP"
        except NotADirectoryError:
            return False, None

    # BIDS FORMAT
    # For BIDS import, this will be in the form raw/subject/[visit]/perf or anat
    else:
        try:
            analysis_dir = raw_dir.parent / "analysis"
            # Get rid of illegal characters for subject
            if "-" in subject:
                subject = subject.replace("-", "")
            if "_" in subject:
                subject = subject.replace("_", "")

            if visit is None:
                if scan not in ["T1", "FLAIR"]:
                    dst_dir = analysis_dir / f"sub-{subject}" / "perf" / "TEMP"
                else:
                    dst_dir = analysis_dir / f"sub-{subject}" / "anat" / "TEMP"
            else:
                # Get rid of illegal characters for visit
                if "-" in visit:
                    visit = visit.replace("-", "")
                if "_" in visit:
                    visit = visit.replace("_", "")
                # Determine output directory based on the scan type
                if scan not in ["T1", "FLAIR"]:
                    dst_dir = analysis_dir / f"sub-{subject}" / f"ses-{visit}" / "perf" / "TEMP"
                else:
                    dst_dir = analysis_dir / f"sub-{subject}" / f"ses-{visit}" / "anat" / "TEMP"
        except NotADirectoryError:
            return False, None

    try:
        dst_dir.mkdir(parents=True, exist_ok=True)
    except OSError:
        return False, None

    return True, dst_dir


def run_dcm2niix(temp_dir: Path, dcm_dir: Path, subject: str, visit, run: str, scan: str):
    """
    Runs the dcm2niix program as a subprocess and generates the appropriate nifti and json files
    :param temp_dir: the TEMP dst where nifti and json files will be deposited
    :param dcm_dir: the directory where the dicoms of interest are being held
    :param subject: the string representing the current subject
    :param visit: the string representing the current visit alias
    :param run: the string representing the current run alias
    :param scan: the string representing the scan. It is either ASL4D, T1, M0, FLAIR, or WMH_SEGM
    :return status: whether the operation was a success or not
    """

    print(f"Inside run_dcn2niix with args:\n"
          f"temp_dir: {temp_dir}\n"
          f"dcm_dir: {dcm_dir}\n"
          f"subject: {subject}\n"
          f"visit: {visit}\n"
          f"run: {run}\n"
          f"scan: {scan}\n")

    if visit is None:
        visit = ""
    else:
        visit = f"_{visit}"

    if run is None:
        run = "_1"
    else:
        run = f"_{run}"

    scan = f"_{scan}"

    if visit is None:
        output_filename_format = f"{subject}{scan}{run}_%s"
    else:
        output_filename_format = f"{subject}{visit}{scan}{run}_%s"

    # Must ensure that no files currently exist within the destination
    try:
        for path in temp_dir.iterdir():
            if path.is_dir():
                path.rmdir()
            else:
                path.unlink(missing_ok=True)
    except FileNotFoundError as ferr:
        print(f"Encountered Error while processing subject: {subject}; visit: {visit}; run: {run}; scan: {scan}"
              f"\n{ferr}")
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

    if system() == "Windows":
        command = f"dcm2niix -b y -z n -x n -t n -m n -s n -v n " \
                  f"-f {output_filename_format} " \
                  f"-o {temp_dir} " \
                  f"{dcm_dir}"

        return_code = subprocess.run(command.split(" "))
        if return_code.returncode == 0:
            return True
        else:
            print(f"DCM2NIIX Did not exit gracefully!!!\nStd Err:\n{return_code.stderr}")
            return False
    else:
        command = f"./dcm2niix -b y -z n -x n -t n -m n -s n -v n " \
                  f"-f {output_filename_format} " \
                  f"-o {temp_dir} " \
                  f"{dcm_dir}".strip()
        p = subprocess.Popen(command, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                             text=True, shell=True)
        p.wait()
        stdout, stderr = p.communicate()
        return_code = p.returncode
        if return_code == 0:
            return True
        else:
            print(f"DCM2NIIX Did not exit gracefully!!!\nStd Err:\n{stderr}")
            return False


def remove_illegal_chars(item):
    if "-" in item:
        item = item.replace("-", "")
    if "_" in item:
        item = item.replace("_", "")
    return item


def clean_niftis_in_temp(temp_dir: Path, add_parms: dict, subject: str, run: str, visit: str, scan: str,
                         legacy_mode: bool = False):
    """
    Concatenates the niftis, deletes the previous ones, and moves the concatenated one out of the temp dir
    :param temp_dir: the absolute filepath to the TEMP directory where the niftis are present
    :param add_parms: additional parameters retrieved such as Acquisition Matrix used to deal with mosaics that did
    not get converted correctly
    :param subject: the string representing the current subject
    :param visit: the string representing the current visit alias
    :param run: the string representing the current run alias
    :param scan: the string representing the scan. It is either ASL4D, T1, M0, FLAIR, or WMH_SEGM
    :param legacy_mode: whether to adjust file naming to the old ExploreASL format or not.
    :return: status: whether the operation was a success or not; the import summary parameters, and the filepath to the
    new nifti created
    """
    reorganized_niftis = []
    import_summary = dict.fromkeys(["subject", "visit", "scan", "filename",
                                    "dx", "dy", "dz", "nx", "ny", "nz", "nt"])
    jsons = peekable(temp_dir.glob("*.json"))
    if not jsons:
        print(f"Did not locate any json files in {temp_dir}")
        return False, import_summary, None, None

    ###################################
    # PART 1 - ORGANIZING DCM2NII FILES
    ###################################
    # Must go over the jsons first to gain an understanding of the series
    json_data = {"SeriesNumber": {}, "AcquisitionTime": {}, "AcquisitionNumber": {}}
    for json_file in jsons:
        with open(json_file) as json_reader:
            sidecar_data = json.load(json_reader)
            for parm in json_data.keys():
                try:
                    json_data[parm][json_file.with_suffix(".nii")] = sidecar_data[parm]
                except KeyError:
                    if parm == "AcquisitionTime":
                        json_data[parm][json_file.with_suffix(".nii")] = 0
                        if "TriggerDelayTime" in list(sidecar_data.keys()):
                            json_data[parm][json_file.with_suffix(".nii")] = sidecar_data["TriggerDelayTime"]

    # Philips case; same series and probably same acq number too; opt for acq_time as the differentiating factor
    if len(set(json_data["SeriesNumber"].values())) == 1 and len(set(json_data["AcquisitionTime"].values())) > 1:
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

    else:
        if len(reorganized_niftis) == 0:
            reorganized_data = {key: value for key, value in sorted(json_data["AcquisitionNumber"].items(),
                                                                    key=lambda x: x[1])}
            reorganized_niftis = list(set(reorganized_data.keys()))
        else:
            return False, import_summary, None, None

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
        # Must correct for bad headers under BIDS specification
        if not legacy_mode and final_nifti_obj.ndim < 4:
            pixdim_copy = final_nifti_obj.header["pixdim"].copy()
            final_nifti_obj.header["dim"][0] = 4
            final_nifti_obj.header["pixdim"] = pixdim_copy
            final_nifti_obj = nib.Nifti1Image(np.expand_dims(image.get_data(final_nifti_obj), axis=-1),
                                              final_nifti_obj.affine,
                                              final_nifti_obj.header)

    # Scenario: single M0 or single ASL4D
    elif any([len(reorganized_niftis) == 1 and scan == "M0",  # Single independent M0
              len(reorganized_niftis) == 1 and scan == "ASL4D"]):  # Could be a CBF or delta image
        final_nifti_obj: nib.Nifti1Image = nib.load(reorganized_niftis[0])
        # Must correct for bad headers under BIDS specification
        if not legacy_mode and final_nifti_obj.ndim < 4:
            pixdim_copy = final_nifti_obj.header["pixdim"].copy()
            final_nifti_obj.header["dim"][0] = 4
            final_nifti_obj.header["pixdim"] = pixdim_copy
            final_nifti_obj = nib.Nifti1Image(np.expand_dims(image.get_data(final_nifti_obj), axis=-1),
                                              final_nifti_obj.affine,
                                              final_nifti_obj.header)

    # Scenario: one of the structural types
    elif len(reorganized_niftis) == 1 and scan in ["T1", "FLAIR"]:
        final_nifti_obj = nib.load(reorganized_niftis[0])

    # Scenario: multiple T1 acquisitions...take the mean
    elif len(reorganized_niftis) > 1 and scan in ["T1", "FLAIR"]:
        nii_objs = [nib.load(nifti) for nifti in reorganized_niftis]
        final_nifti_obj = image.mean_img(nii_objs)

    # Otherwise, something went wrong and the operation should stop
    else:
        print(f"Error in clean_niftis_in_temp while attempting to process subject {subject}; visit {visit};"
              f" run {run}; scan {scan}\n"
              f"Reorganized_niftis did not fit into any of the foreseen scenarios\n"
              f"Length of reorganized_niftis: {len(reorganized_niftis)}")
        return False, import_summary, None, None

    # Take the oppurtunity to get more givens for the import summary
    zooms = final_nifti_obj.header.get_zooms()
    shape = final_nifti_obj.shape
    import_summary["subject"] = subject
    import_summary["visit"] = visit
    import_summary["run"] = run
    import_summary["scan"] = scan
    import_summary["filename"] = scan + ".nii"

    if len(zooms) >= 4:
        import_summary["dx"], import_summary["dy"], import_summary["dz"] = zooms[0:3]
    else:
        import_summary["dx"], import_summary["dy"], import_summary["dz"] = zooms

    if len(shape) == 4:
        (import_summary["nx"], import_summary["ny"],
         import_summary["nz"], import_summary["nt"]) = shape
    else:
        (import_summary["nx"], import_summary["ny"],
         import_summary["nz"], import_summary["nt"]) = shape[0], shape[1], shape[2], 1

    ################################
    # PART 3 NAMING AND MOVING FILES
    ################################
    if run is None:
        run = ""
    else:
        # Get rid of illegal characters for visit
        run = remove_illegal_chars(run)
        run = f"run-{run}_"

    if visit is None:
        visit = ""
    else:
        # Get rid of illegal character for visit
        visit = remove_illegal_chars(visit)
        visit = f"ses-{visit}_"

    if legacy_mode:
        final_nifti_filename = temp_dir.parent / f"{scan}.nii"
        final_json_filename = temp_dir.parent / f"{scan}.json"

    else:
        scan = {"ASL4D": "asl", "M0": "m0scan", "T1": "T1w", "FLAIR": "FLAIR"}[scan]
        # Get rid of illegal characters for subject
        if "-" in subject:
            subject = subject.replace("-", "")
        if "_" in subject:
            subject = subject.replace("_", "")

        final_nifti_filename = temp_dir.parent / f"sub-{subject}_{visit}{run}{scan}.nii"
        final_json_filename = temp_dir.parent / f"sub-{subject}_{visit}{run}{scan}.json"

    nib.save(final_nifti_obj, final_nifti_filename)
    jsons = peekable(temp_dir.glob("*json"))
    if not jsons:
        print("Error in clean_niftis_in_temp while attempting to rename remaining json files")
        return False, import_summary, None, None
    json_file = next(jsons)
    json_file.replace(final_json_filename)

    return True, import_summary, final_nifti_filename, final_json_filename


def update_json_sidecar(json_file: Path, nifti_file: Path, dcm_parms: dict):
    """
    Updates the json sidecar to include additional DICOM parameters not extracted by dcm2niix
    :param json_file: the absolute filepath to the json sidecar
    :param nifti_file: the absolute filepath to the nifti corresponding to the json
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

        # Next, must see if Philips-related fixes post-DCM2NIIX are necessary
        manufac = json_sidecar_parms.get("Manufacturer", None)
        if manufac == "Philips":
            # One possibility: Array values are Stored Values and must be corrected to Philips Floating Point
            if all(["PhilipsRescaleSlope" in json_sidecar_parms,
                    "PhilipsRescaleIntercept" in json_sidecar_parms,
                    "PhilipsScaleSlope" in json_sidecar_parms,
                    "PhilipsRWVSlope" not in json_sidecar_parms,
                    json_sidecar_parms.get("UsePhilipsFloatNotDisplayScaling", None) == 0,
                    ]):
                nifti_img: nib.Nifti1Image = image.load_img(nifti_file)
                nifti_data: np.ndarray = image.get_data(nifti_img)
                RI = json_sidecar_parms["PhilipsRescaleIntercept"]
                RS = json_sidecar_parms["PhilipsRescaleSlope"]
                SS = json_sidecar_parms["PhilipsScaleSlope"]
                new_nifti_data = (nifti_data + (RI / RS)) / SS
                new_nifti = nib.Nifti1Image(dataobj=new_nifti_data, header=nifti_img.header, affine=nifti_img.affine)
                del nifti_img, nifti_data, RI, RS, SS, new_nifti_data
                new_nifti.to_filename(nifti_file)
                del new_nifti
                json_sidecar_parms["UsePhilipsFloatNotDisplayScaling"] = 1

            # Another possibility: Array values were incorrectly converted to Display Values rather than Philips
            # Floating Point
            elif all(["PhilipsRescaleSlope" in json_sidecar_parms,
                      "PhilipsRescaleIntercept" in json_sidecar_parms,
                      "PhilipsScaleSlope" in json_sidecar_parms,
                      "PhilipsRWVSlope" in json_sidecar_parms,
                      json_sidecar_parms.get("UsePhilipsFloatNotDisplayScaling", None) == 1,
                      ]):
                nifti_img: nib.Nifti1Image = image.load_img(nifti_file)
                nifti_data: np.ndarray = image.get_data(nifti_img)
                RS, SS = json_sidecar_parms["PhilipsRescaleSlope"], json_sidecar_parms["PhilipsScaleSlope"]
                new_nifti_data = nifti_data / (RS * SS)
                new_nifti = nib.Nifti1Image(dataobj=new_nifti_data, header=nifti_img.header, affine=nifti_img.affine)
                del nifti_img, nifti_data, RS, SS, new_nifti_data
                new_nifti.to_filename(nifti_file)
                del new_nifti
        elif manufac is None:
            print(f"Encountered error in update_json_sidecar for file {json_file}\n"
                  f"The mandatory Manufacturer key was not present in the json")
            return False, None
        else:
            pass

        with open(json_file, 'w') as w:
            json.dump(json_sidecar_parms, w, indent=3)

    except FileNotFoundError:
        return False, None

    return True, json_sidecar_parms


def asldcm2bids_onedir(dcm_dir: Union[Path, str], config: dict, legacy_mode: bool = False):
    """
    The main function for most ASL-processing steps centered around processing a single dicom directory immediately
    preceding the dicom files.
    :param dcm_dir: the absolute path to the dicom directory
    :param config: the configuration file containing the essential parameters for importing
    :param legacy_mode: whether to adjust file naming to the old ExploreASL format or not.
    :returns: 3 items: whether the operation was a success; the most recent step performed, and the import summary (or
    None if the import for this directory failed)
    """
    # Get the subject, visit, and scan for that directory
    # True, subject, visit, run_dst_name, scan_dst_name
    if isinstance(dcm_dir, str):
        dcm_dir = Path(dcm_dir)
    successful_run, subject, visit, run, scan = get_structure_components(dcm_dir=dcm_dir, config=config)
    if not successful_run:
        print(f"FAILURE ENCOUNTERED AT GETTING THE DIRECTORY STRUCTURE COMPONENTS STEP")
        return False, f"SUBJECT: {subject}; " \
                      f"SCAN: {scan}; " \
                      f"ERROR: Failed at getting directory structure components", None

    # Get the manufacturer tag
    # successful_run, manufacturer = get_manufacturer(dcm_dir=dcm_dir)
    # if not successful_run:
    #     print(f"FAILURE ENCOUNTERED AT GETTING THE MANUFACTURER STEP")
    #     return False, f"SUBJECT: {subject}; " \
    #                   f"SCAN: {scan}; " \
    #                   f"ERROR: Failed retrieving the manufacturer value", None

    # Retrieve the additional DICOM parameters and include the Philips rescale slope indicator
    successful_run, additional_dcm_parameters = get_additional_dicom_parms(dcm_dir=dcm_dir)
    if not successful_run:
        print(f"FAILURE ENCOUNTERED AT THE ADDITIONAL DCM PARMS RETRIEVAL STEP")
        return False, f"SUBJECT: {subject}; " \
                      f"SCAN: {scan}; " \
                      f"ERROR: Failed at retrieving additional dicom parameters", None

    # Generate the directories for dumping dcm2niix output
    successful_run, temp_dst_dir = get_tempdst_dirname(raw_dir=Path(config["RawDir"]),
                                                       subject=subject,
                                                       visit=visit,
                                                       run=run,
                                                       scan=scan,
                                                       legacy_mode=legacy_mode)
    if not successful_run:
        print(f"FAILURE ENCOUNTERED AT THE TEMP FOLDER GENERATION STEP")
        return False, f"SUBJECT: {subject}; " \
                      f"SCAN: {scan}; " \
                      f"ERROR: Failed at Temp folder generation", None

    # Run the main program
    successful_run = run_dcm2niix(temp_dir=temp_dst_dir,
                                  dcm_dir=dcm_dir,
                                  subject=subject,
                                  run=run,
                                  visit=visit,
                                  scan=scan)
    if not successful_run:
        print(f"FAILURE ENCOUNTERED AT THE DCM2NIIX STEP")
        return False, f"SUBJECT: {subject}; " \
                      f"SCAN: {scan}; " \
                      f"ERROR: Failed at DCM2NIIX conversion", None

    # Clean the niftis in the TEMP directory
    (successful_run,
     nifti_parmameters,
     nifti_filepath,
     json_filepath) = clean_niftis_in_temp(temp_dir=temp_dst_dir,
                                           add_parms=additional_dcm_parameters,
                                           subject=subject,
                                           run=run,
                                           visit=visit,
                                           scan=scan,
                                           legacy_mode=legacy_mode)
    if not successful_run:
        print(f"FAILURE ENCOUNTERED AT CLEANING THE NIFTIs IN THE TEMP FOLDER")
        return False, f"SUBJECT: {subject}; " \
                      f"SCAN: {scan}; " \
                      f"ERROR: Failed at post-conversion nifti cleanup", None

    # Minor preprocessing step; overwrite the number of slices based on the nifti shape
    if additional_dcm_parameters["NumberOfSlices"] is None:
        additional_dcm_parameters["NumberOfSlices"] = nifti_parmameters["nz"]

    # For ASL4D and M0 scans, the JSON sidecar from dcm2niix must include additional parameters
    successful_run, json_parameters = update_json_sidecar(json_file=json_filepath,
                                                          nifti_file=nifti_filepath,
                                                          dcm_parms=additional_dcm_parameters)

    if not successful_run:
        print(f"FAILURE ENCOUNTERED AT CORRECTION THE JSON SIDECAR STEP")
        return False, f"SUBJECT: {subject}; " \
                      f"SCAN: {scan}; " \
                      f"ERROR: Failed at updating json sidecar with additional parms", None

    # If everything was a success, prepare the dict for creating the import summary
    import_summary = {}
    import_summary.update(additional_dcm_parameters)
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
    Given a list of individual summaries of each subject/visit/scan, this function will bring all those givens
    together into a single dataframe for easy viewing
    :param import_summaries: a list of dicts, with each dict being the parameters of that subject-visit-scan
    :param config: the import configuration file generated by the GUI to help locate the analysis directory
    """
    analysis_dir = Path(config["RawDir"]).parent / "analysis"
    try:
        df = pd.concat([pd.Series(import_summary) for import_summary in import_summaries], axis=1, sort=True).T
    except ValueError as concat_error:
        print(concat_error)
        return

    df["dt"] = df["RepetitionTime"]
    appropriate_ordering = ['subject', 'visit', 'run', 'scan', 'dx', 'dy', 'dz', 'dt', 'nx', 'ny', 'nz', 'nt',
                            "RepetitionTime", "EchoTime", "NumberOfAverages", "RescaleSlope", "RescaleIntercept",
                            "MRScaleSlope", "AcquisitionTime",
                            "AcquisitionMatrix", "TotalReadoutTime", "EffectiveEchoSpacing"]
    df = df.reindex(columns=appropriate_ordering)
    df = df.sort_values(by=["scan", "subject", "visit", "run"]).reset_index(drop=True)
    print(df)
    now_str = datetime.now().strftime("%a-%b-%d-%Y %H-%M-%S")
    try:
        df.to_csv(analysis_dir / f"import_summary_{now_str}.tsv", sep='\t', index=False, na_rep='n/a')
    except PermissionError:
        df.to_csv(analysis_dir / f"import_summary_backup_{now_str}.tsv", sep='\t', index=False, na_rep='n/a')


def bids_m0_followup(analysis_dir: Path):
    """
    In a BIDS import, this function will run through the imported dataset and adjust any BIDS-standard fields that
    should be present in the m0scan.json sidecar, such as "IntendedFor"
    :param analysis_dir: the absolute path to the analysis directory
    """
    m0_jsons = peekable(analysis_dir.rglob("*_m0scan.json"))
    if not m0_jsons:
        print("bids_m0_followup could not find any _m0scan.json files")
        return
    for m0_json in m0_jsons:
        asl_json = m0_json.with_name("_asl.json")
        asl_nifti = m0_json.with_name("_asl.nii")

        # Ensure that the asl json sidecar and nifti images actually exist adjacent to the m0scan.json
        if asl_json.exists() and asl_nifti.exists():
            # BIDS standard: the "IntendedFor" filepath must be relative to the subject (exclusive)
            # and contain forward slashes
            truncated_asl_nifti = str(asl_nifti).replace(str(analysis_dir), "").replace("\\", "/")
            by_parts = truncated_asl_nifti.split(sep='/')
            truncated_asl_nifti = "/".join(by_parts[2:])

            with open(m0_json) as m0_json_reader:
                m0_parms = json.load(m0_json_reader)
            m0_parms["IntendedFor"] = truncated_asl_nifti
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
