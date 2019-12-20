import modules.globalenv as globalenv
import os
import sys
import json
from fnmatch import fnmatch
import shutil
import pandas as pd
import logging

"""
Changelog:
  QCDC 1.2.0 20190619 [JK]
    - in descriptor json changed entry sub-path to sub_path for Matlab variable name compatibility
"""

logger = logging.getLogger('root')

def show_args_details(args):
    logger.info("""Run QC Data collector using:
    Descriptor: {}
    Workingpath: {}""".format(args.descriptor, args.workingpath))


def check_args(args):
    if not os.path.isabs(args.descriptor):
        logger.error("Error descriptor {} is not an absolute path!".format(args.descriptor))
        sys.exit()

    if not os.path.exists(args.descriptor):
        logger.error("Error descriptor {} file does not exist!".format(args.descriptor))
        sys.exit()

    if not args.descriptor.endswith('.json'):
        logger.error("Error {} file is not a JSON file!".format(args.descriptor))
        sys.exit()

    if not os.path.isabs(args.workingpath):
        logger.error("Error workingpath {} is not an absolute path!".format(args.workingpath))
        sys.exit()

def read_csv(filepath, header=None, sep=None):
    df = pd.read_csv(filepath, header=header, sep=sep, index_col=False, engine='python')
    return df

def read_json(filepath):
    json_config_file = open(filepath)
    return json.loads(json_config_file.read())


# generate_abs_filepath generate an absolute filepath
# if the generated path is a pattern for multiple files only a valid abs_filepath will be returned
def generate_abs_filepath_and_match(workingpath, subpath, filename):
    if not workingpath.endswith("/"):
        workingpath = workingpath + '/'

    if not subpath == '' and not subpath.endswith("/"):
        subpath = subpath + '/'

    abs_file_path_match = workingpath + subpath + filename
    if os.path.exists(abs_file_path_match):
        return abs_file_path_match
    else:
        abs_file_path = get_file_list_matching_pattern(workingpath + subpath, filename)
        return abs_file_path


def get_file_list_matching_pattern(root_path, pattern):
    for path, subdirs, files in os.walk(root_path):
        for name in files:
            if fnmatch(name, pattern):
                return os.path.join(path, name)

    logger.warning("WARNING: util.get_file_list_matching_pattern {} not matching any files".format(root_path+pattern))
    return False

def create_pack_environment():
    if not os.path.exists(globalenv.pack_output):
        os.mkdir(globalenv.pack_output)
    else:
        shutil.rmtree(globalenv.pack_output)
        os.mkdir(globalenv.pack_output)


def write_dict_to_json(abs_path_file, dict):
    dict_json = json.dumps(dict, sort_keys=True, indent=4)
    f = open(abs_path_file, "w")
    f.write(dict_json)
    f.close()



def check_dicom_wrapper(dicom_wrapper):
    if 'sub_path' not in dicom_wrapper:
        logger.error("Error: Mandatory field 'sub_path' does not exist in dicom_wrapper: \n"+json.dumps(dicom_wrapper,sort_keys=True, indent=4))
        return False

    if 'filename' not in dicom_wrapper:
        logger.error("Error: Mandatory field 'filename' does not exist in dicom_wrapper: \n"+json.dumps(dicom_wrapper,sort_keys=True, indent=4))
        return False

    if 'filepath_wadqc_placeholder' not in dicom_wrapper:
        logger.error("Error: Mandatory field 'filepath_wadqc_placeholder' does not exist in dicom_wrapper: \n"+json.dumps(dicom_wrapper,sort_keys=True, indent=4))
        return False

    return True


def copy_dicom_wrapper(dicom_wrapper):
    if not check_dicom_wrapper(dicom_wrapper):
        return False

    abs_filepath = generate_abs_filepath_and_match(globalenv.args.workingpath, dicom_wrapper['sub_path'], dicom_wrapper['filename'])

    if not abs_filepath:
        return False

    copied_dicom = globalenv.pack_output + abs_filepath.split('/')[-1]
    shutil.copyfile(abs_filepath, copied_dicom)

    return copied_dicom


def get_nifti_file(workingpath, nii_filename, qc_item):
    abs_filepath = generate_abs_filepath_and_match(workingpath, qc_item['sub_path'], nii_filename)

    if not abs_filepath and nii_filename.endswith('.nii'):
        nii_filename = nii_filename + '.gz'

    abs_filepath = generate_abs_filepath_and_match(workingpath, qc_item['sub_path'], nii_filename)

    if not abs_filepath and nii_filename.endswith('.gz'):
        nii_filename = os.path.splitext(nii_filename)[0]

    abs_filepath = generate_abs_filepath_and_match(workingpath, qc_item['sub_path'], nii_filename)

    if not abs_filepath:
        logger.error("Error in qc_item {}: No nifti file for has been found: {}\n".format(json.dumps(qc_item, sort_keys=True,
                                                                                                             indent=4)))
        return False

    return abs_filepath