#! /usr/bin/env python
#system requirements:
# pandas: pip3 install pandas
# nibabel: pip3 install nibabel
# fsl scripts slicer

import argparse
import modules.globalenv as globalenv
import modules.fs_working_env_management as fswem
import modules.descriptor_processor as dp
import json
import sys
import modules.wadqc as wadqc
import modules.log as log
import os

"""
Changelog:
  QCDC 1.1.0 20190115 [JK]
    - added completeness of report as a result itself. Items is_complete and missing_items
  QCDC 1.2.0 20190619 [JK]
    - in descriptor json changed entry sub-path to sub_path and ae-title to ae_title for Matlab variable name compatibility
    - bug fix in descriptor_processor.py: crash if process_qc_item_file did to match a file
  QCDC 1.2.1 20190621 [JK]
    - modified send bash script: DICOM UID is no longer changed
    - temporary fix for Azure cloud setup: disabled DICOM send
    - python 2 compatibility fix

"""


parser = argparse.ArgumentParser()
parser.add_argument("descriptor", help="descriptor JSON absolute file path")
parser.add_argument("workingpath", help="where the QC Data Collector needs to be executed")

args = parser.parse_args()

globalenv.args = args

# Fix to allow parallel processing, with multiple unique qcdc_output folders
FileName = os.path.basename(args.descriptor) # get filename
FileName = os.path.splitext(FileName)[0] # remove extension

globalenv.pack_output = args.workingpath + 'qcdc_output_' + FileName + '/' # add FileName

fswem.create_pack_environment()

logger = log.setup_custom_logger('root', globalenv.pack_output)

fswem.show_args_details(args)
fswem.check_args(args)

descriptor = dp.load_descriptor(args.descriptor)

copied_dicom = fswem.copy_dicom_wrapper(descriptor['dicom_meta']['dicom_wrapper'])

if not copied_dicom:
    logger.critical("Critical Error! Wrapper dicom file is missing or not working with this dicom_wrapper: \n" + json.dumps(
        descriptor["dicom_meta"]["dicom_wrapper"], sort_keys=True, indent=4))
    sys.exit(1)

results_json = dp.process_descriptor_items(args.workingpath, descriptor['qc_items'], descriptor["dicom_meta"]["dicom_wrapper"])

results_json = dp.add_descriptor_log(results_json, descriptor["dicom_meta"]["dicom_wrapper"])

fswem.write_dict_to_json(globalenv.pack_output+'results.json', results_json)

# zip all the file inside the workingpath except the *.dcm
logger.info(wadqc.create_zip_package(globalenv.pack_output))

# send the dicom wrapper to WADQC
# JK v1.2.1 disable DICOM send
#logger.info(wadqc.send(descriptor['wad_qc_server'], copied_dicom))
logger.info(wadqc.create_dcm_only_but_do_not_send(descriptor['wad_qc_server'], copied_dicom))


#if is_sent:
#    logger.info("SUCCESS, patient {} sent to WAD QC".format(copied_dicom))
#else:
#    logger.error("ERROR, patient {} has not been sent to WAD QC".format(copied_dicom))
