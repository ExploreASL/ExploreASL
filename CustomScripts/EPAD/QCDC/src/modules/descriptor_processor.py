import modules.globalenv as globalenv
import modules.fs_working_env_management as fswem
import os
import json
from shutil import copyfile
import nibabel as nib
import random, string
import logging

logger = logging.getLogger('root')

"""
Changelog:
  QCDC 1.1.0 20190115 [JK]
    - added completeness of report as a result itself. Items is_complete and missing_items
  QCDC 1.2.0 20190619 [JK]
    - in descriptor json changed entry sub-path to sub_path for Matlab variable name compatibility
    - bug fix in descriptor_processor.py: crash if process_qc_item_file did to match a file
    - max length of string import in WAD-QC is 100. Limit missing_items_string length to 90
  QCDC 1.2.1 20190621 [JK]
    - modified send bash script: DICOM UID is no longer changed
    - temporary fix for Azure cloud setup: disabled DICOM send
    - python 2 compatibility fix

"""


# Max length of string import in WAD-QC is 100. Limit missing_items_string length to 90, if shortened add '...'
wad_qc_max_string_length = 90


def load_descriptor(filepath):
    descriptor = fswem.read_json(filepath)

    if not descriptor["dicom_meta"]["dicom_wrapper"]['filepath_wadqc_placeholder'].endswith('/'):
        descriptor["dicom_meta"]["dicom_wrapper"]['filepath_wadqc_placeholder'] + '/'

    return descriptor


##################################################
# Checks functions ------------------------------
##################################################


def check_csv_item_fields(item):
    if 'column' not in item:
        logger.error("Error: Mandatory field 'column' does not exist in qc_item for type 'csv': \n" + json.dumps(item, sort_keys=True, indent=4))
        return False
    if 'row' not in item:
        logger.error("Error: Mandatory field 'row' does not exist in qc_item for type 'csv': \n" + json.dumps(item, sort_keys=True, indent=4))
        return False
    return True



def check_items_fields(item):
    if 'type' not in item:
        logger.error("Error: Mandatory field 'type' does not exist in qc_item: \n"+json.dumps(item,sort_keys=True, indent=4))
        return False

    if item['type'] not in globalenv.available_types:
        logger.error("Error: Type field in qc_item is not valid: \n"+json.dumps(item,sort_keys=True, indent=4))
        logger.error("Type field in qc_item please insert a valid type: "+json.dumps(globalenv.available_types))
        return False

    if 'filename' not in item:
        logger.error("Error: Mandatory field 'filename' does not exist in qc_item: \n" + json.dumps(item, sort_keys=True, indent=4))
        return False

    if 'sub_path' not in item:
        logger.error("Error: Mandatory field 'sub_path' does not exist in qc_item: \n" + json.dumps(item, sort_keys=True, indent=4))
        return False

    if item['type'] == 'csv':
        return check_csv_item_fields(item)

    return True


def check_file_qc_item_exists(workingpath, qc_item):
    filename = qc_item['filename']
    subpath = qc_item['sub_path']

    abs_file_path = fswem.generate_abs_filepath_and_match(workingpath, subpath, filename)

    if (not abs_file_path) or (not os.path.exists(abs_file_path)):
        return False

    return True

##################################################
# Process functions ------------------------------
##################################################

def process_qc_item_csv(qc_item_key, qc_item):
    header = None
    if 'header' in qc_item and qc_item['header']:
        header = 0

    sep = None
    if 'separator' in qc_item:
        sep = qc_item['separator']

    df = fswem.read_csv(fswem.generate_abs_filepath_and_match(globalenv.args.workingpath, qc_item['sub_path'], qc_item['filename']), header=header, sep=sep)

    try:
        val = df.iloc[qc_item['row'], qc_item['column']]
    except IndexError:
        logger.error("Error: in qc_item {} index out of bound for row {} or column {} : {}".format(qc_item_key, qc_item['row'], qc_item['column'], json.dumps(qc_item,sort_keys=True, indent=4)))
        return False

    del df
    return str(val)


def process_qc_item_json(qc_item_key, qc_item):
    json_abs_filepath = fswem.generate_abs_filepath_and_match(globalenv.args.workingpath, qc_item['sub_path'], qc_item['filename'])

    if 'child' not in qc_item:
        logger.error("Error in qc_item {}: Mandatory field 'child' for type json in qc_item: {}\n".format(qc_item_key, json.dumps(qc_item, sort_keys=True,
                                                                                         indent=4)))
        return False

    try:
        json_file = open(json_abs_filepath)
        json_str = json_file.read()
        json_data = json.loads(json_str)
    except json.decoder.JSONDecodeError:
        logger.error("Error in qc_item {}: parsing file {} ".format(qc_item_key, json_file))
        return False

    child_list = qc_item['child'].split('/')

    json_data_tmp = json_data

    for child in child_list:
        if child not in json_data_tmp:
            logger.error("Error in qc_item {}: traversing json {} child not in json {}".format(qc_item_key, child, json.dumps(json_data_tmp,sort_keys=True, indent=4)))
            return False
        json_data_tmp = json_data_tmp[child]

    if type(json_data_tmp) is dict:
        logger.error("Error in qc_item {}: traversing json {} child in json is an object has to be a value {}".format(qc_item_key, child,
                                                                                    json.dumps(json_data_tmp,
                                                                                               sort_keys=True,
                                                                                               indent=4)))
        return False

    return str(json_data_tmp)


def process_qc_item_file(qc_item_key, workingpath, qc_item, dicom_wrapper):
    abs_filepath = fswem.generate_abs_filepath_and_match(workingpath, qc_item['sub_path'], qc_item['filename'])

    if not abs_filepath:
        return False

    copyfile(abs_filepath,globalenv.pack_output+abs_filepath.split('/')[-1])

    return dicom_wrapper['filepath_wadqc_placeholder'] + abs_filepath.split('/')[-1]


def process_qc_item_nii_img(qc_item_key, workingpath, qc_item, dicom_wrapper):
    if 'fsl_slicer_options' not in qc_item:
        logger.error("Error in qc_item {}: Mandatory field 'fsl_slicer_option' for type nii.img in qc_item: {}\n".format(qc_item_key, json.dumps(qc_item, sort_keys=True,
                                                                                         indent=4)))
        return False

    nii_filename = qc_item['filename']

    abs_filepath = fswem.get_nifti_file(workingpath, nii_filename, qc_item)

    if not abs_filepath:
        return False

    filename_img = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(16)) + ".png"

    filepath_img = globalenv.pack_output + filename_img

    command = "slicer {} {} {}".format(abs_filepath,qc_item['fsl_slicer_options'], filepath_img)
    os.system(command)

    copyfile(filepath_img,globalenv.pack_output+abs_filepath.split('/')[-1])

    return dicom_wrapper['filepath_wadqc_placeholder'] + filename_img


def process_qc_item_nii_hdr(qc_item_key, workingpath, qc_item, dicom_wrapper):
    if 'field_name' not in qc_item:
        logger.error("Error in qc_item {}: Mandatory field 'field_name' for type nii.hdr in qc_item: {}\n".format(qc_item_key, json.dumps(qc_item, sort_keys=True,
                                                                                         indent=4)))
        return False

    nii_filename = qc_item['filename']

    abs_filepath = fswem.get_nifti_file(workingpath, nii_filename, qc_item)

    if not abs_filepath:
        return False

    img = nib.load(abs_filepath)
    return img.header[qc_item['field_name']].item(0)


def process_descriptor_items(workingpath, qc_items, dicom_wrapper):
    result_json_array = []
    flag_error = False
    missing_items_string = ""

    for qc_item_key in qc_items:
        qc_item = qc_items[qc_item_key]
        if not check_items_fields(qc_item):
            flag_error = True
            missing_items_string = missing_items_string + " " + qc_item_key
            continue

        # check abs_filepath is existing
        if not check_file_qc_item_exists(workingpath, qc_item):
            abs_filepath = fswem.generate_abs_filepath_and_match(workingpath, qc_item['sub_path'], qc_item['filename'])
            logger.error("Error in qc_item {}: file not exist: {}".format(qc_item_key,abs_filepath))
            flag_error = True
            missing_items_string = missing_items_string + " " + qc_item_key
            continue

        if 'category' not in qc_item:
            logger.error("Error in qc_item {}: Mandatory field 'category' not in qc_item: {}\n".format(qc_item_key,
                                                                                                       json.dumps(
                                                                                                           qc_item,
                                                                                                           sort_keys=True,                                                                                                           indent=4)))
            missing_items_string = missing_items_string + " " + qc_item_key
            continue

        result_item = {}
        result_item['name'] = qc_item_key
        result_item['category'] = qc_item['category']

        if qc_item['type'] == 'csv':
            val = process_qc_item_csv(qc_item_key, qc_item)
            if not val:
                flag_error = True
                missing_items_string = missing_items_string + " " + qc_item_key
                continue
            result_item['val'] = val
        elif qc_item['type'] == 'json':
            val = process_qc_item_json(qc_item_key, qc_item)
            if not val:
                flag_error = True
                missing_items_string = missing_items_string + " " + qc_item_key
                continue
            result_item['val'] = val
        elif qc_item['type'] == 'file':
            val = process_qc_item_file(qc_item_key, workingpath, qc_item, dicom_wrapper)
            if not val:
                flag_error = True
                missing_items_string = missing_items_string + " " + qc_item_key
                continue
            result_item['val'] = val
        elif qc_item['type'] == 'nii.img':
            val = process_qc_item_nii_img(qc_item_key, workingpath, qc_item, dicom_wrapper)
            if not val:
                flag_error = True
                missing_items_string = missing_items_string + " " + qc_item_key
                continue
            result_item['val'] = val
        elif qc_item['type'] == 'nii.hdr':
            val = process_qc_item_nii_hdr(qc_item_key, workingpath, qc_item, dicom_wrapper)
            if not val:
                flag_error = True
                missing_items_string = missing_items_string + " " + qc_item_key
                continue
            result_item['val'] = val

        result_json_array.append(result_item)

    # completeness of results is a result by itself
    result_json_array = add_result(result_json_array, "is_complete", "float", float(not flag_error))

    if flag_error:
        logger.warning("WARNING! Not all the qc_item in descriptor have been detected")

        # include list of missing items as result
        # limit string length or it will produce error upon import in WAD-QC (v2.0.2)
        if len(missing_items_string) <= wad_qc_max_string_length:
            result_json_array = add_result(result_json_array, "missing_items", "string", missing_items_string)
        else:
            result_json_array = add_result(result_json_array, "missing_items", "string", missing_items_string[:wad_qc_max_string_length]+'...')
    else:
        logger.info("Success, all qc_items in descriptor file have been detected")
    


    return result_json_array


def add_descriptor_log(result_json_array, dicom_wrapper):
    result_item = {}
    result_item['name'] = 'qc_data_collector.log'
    result_item['category'] = 'object'
    result_item['val'] = dicom_wrapper['filepath_wadqc_placeholder'] + 'qc_data_collector.log'

    result_json_array.append(result_item)
    return result_json_array
    
def add_result(result_json_array, name, category, value):
    result_item = {}
    result_item['name'] = name
    result_item['category'] = category
    result_item['val'] = value

    result_json_array.append(result_item)
    return result_json_array
