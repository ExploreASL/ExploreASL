#!/bin/bash
#system requirements:
# uuencode / uudecode available
#sudo apt-get install sharutils
#
# Changelog:
# QCDC 1.2.1 20190621 [JK]
#  - modified send bash script: DICOM UID is no longer changed
#  - temporary fix for Azure cloud setup: disabled DICOM send (in file create_dcm_only_wadqc.sh)
#


if [[ "$1" != "" ]]; then
    dcmfile="$1"
else
    dcmfile=.
fi

if [[ "$2" != "" ]]; then
    wad_qc_server_ip="$2"
else
    wad_qc_server_ip=.
fi

if [[ "$3" != "" ]]; then
    wad_qc_server_port="$3"
else
    wad_qc_server_port=.
fi

if [[ "$4" != "" ]]; then
    wad_qc_server_aetitle="$4"
else
    wad_qc_server_aetitle=.
fi

if [[ "$wad_qc_server_aetitle" == "." ]]; then
    echo "sendwadqc parameter is missing"
    exit 0
fi

zipfile="${dcmfile%/*}/results.zip"
uuefile="${dcmfile%/*}/results.uue"

echo $zipfile
echo $uuefile

# uuencode the zip file for transport in DICOM file
uuencode -m "$zipfile" "results.zip" > "$uuefile"

# DICOM field needs even length; patch with newline if odd length
size=$(stat -c '%s' "$uuefile")
if (($size % 2 == 1))
then
    #echo "patch odd-length UUE file"
    echo "" >> "$uuefile"
fi

# pack in DICOM field
#dcmodify -i "5001,0010=PYWAD_RESULTS" -if "5001,1001="$uuefile "$dcmfile"
# v1.2.1 do not generate new instance UID, remove dcmodify option -gin
#dcmodify -gin -if "TextValue="$uuefile "$dcmfile"
dcmodify -if "TextValue="$uuefile "$dcmfile"

echo $dcmfile

# v1.2.1 temp fix do not send DICOM file to WAD-QC
#storescu -aec "$wad_qc_server_aetitle" "$wad_qc_server_ip" "$wad_qc_server_port" "$dcmfile"

