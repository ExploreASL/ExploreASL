import os
import logging

"""
Changelog:
  QCDC 1.2.0 20190619 [JK]
    - in descriptor json changed entry sub-path to sub_path and ae-title to ae_title for Matlab variable name compatibility

"""


logger = logging.getLogger('root')

def send(wad_qc_server, dicom):
    os.chdir(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../bash'))
    command = "./sendwadqc.sh {} {} {} {}".format(dicom, wad_qc_server['ip_address'],  wad_qc_server['port'],  wad_qc_server['ae_title'])
    result = os.popen(command).read()
    return result

def create_dcm_only_but_do_not_send(wad_qc_server, dicom):
    os.chdir(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../bash'))
    command = "./create_dcm_only_wadqc.sh {} {} {} {}".format(dicom, wad_qc_server['ip_address'],  wad_qc_server['port'],  wad_qc_server['ae_title'])
    result = os.popen(command).read()
    return result

def create_zip_package(pack_output):
    os.chdir(pack_output)
    command = "zip -r -j results.zip * -x 004-0001.dcm "
    result = os.system(command)
    return result
