------------------------------------------------------------------------------------------------------------------------------------
Software Name: Quality Check Data Collector (QCDC)
Author: Gaspare Scherma
(C) G Scherma, S Ingala, AM Wink, HJ Mutsaerts, J Kuijer
License: This software is distributed under GPL license
------------------------------------------------------------------------------------------------------------------------------------

This python module gathers all the information from single subject analysis file csv, json, nifti (header fields, images) and file,
wrapping all elements together into a dicom collector and it submits all this information to the WAD-QC application.

Run:
    python qc_data_collector.py <descriptor_path> <patient_folder>

The input arguments for the scripts are:
-   <descriptor_path>: each json descriptor file is going to describe a single type of analysis, the list of elements you want gather and the WADQC server information.
-   <patient_folder>: or working path is patient folder when you want to run the QCDC

A single run of qc_data_collector will process a single patient analysis

------------------------------------------------------------------------------------------------------------------------------------
Descriptor
------------------------------------------------------------------------------------------------------------------------------------

The descriptor is the main resource needed to generate a patient quality check and it can be located anywhere on the file system.

NB: all the filename fields inside each object have to be intended as a matching expressio. They can be populated with an exact filename indication, with a matching expression case
(for example WMH_LST*.csv will match any file with WMH_LST prefix and .csv extension). In the case multiple files are selected with this match, a random one will be selected.

This is divided in 3 sections
    - dicom_meta: needed to select the dicom file used for wrapping all the information and push it to WAD-QC
      sub-path: under the <patient_folder> a subpath where to find the file to match
      filename: file matching expression
      filepath_wadqc_placeholder:

      example:
      {
        "dicom_meta":{
            "dicom_wrapper":{
                "sub-path": "",
                "filename":"*.dcm",
                "filepath_wadqc_placeholder":"path_placeholder/"
            }
      }

    - qc_item: Is the object containing of all items to gather all around the patient files:
               the type of items can be: csv,json,file, nii.hdr, nii.img. all fields described are required.

               csv:
               "Determinant_Orientation_T1w":{
                  "type": "csv",                                        type of object
                  "sub-path": "",                                       sub-path where to use the matching expression
                  "filename": "CheckOrientation_RigidRegT1.csv",        filename matching expression
                  "column": 0,                                          column number where to get the value
                  "row": 0,                                             row number where to get the value
                  "header": true,                                       set to true is the csv file contains the header
                  "category": "float"                                   category needed for result.json in WADQC
                },

               json:
               "LR_flip_ASL(L)":{
                  "type": "json",                                       type of object
                  "sub-path": "dartel/",                                sub-path where to use the matching expression
                  "filename": "QC_*.json",                              filename matching expression
                  "child":"ASL/LR_flip_ASL",                            path of object tree where to find the value
                  "category": "string"                                  category needed for result.json in WADQC
                },

               file:
               "file test":{
                  "type": "file",                                       type of object
                  "sub-path": "dartel/Tissue_Volume",                   sub-path where to use the matching expression
                  "filename": "Tissue_*.csv",                           filename matching expression
                  "category": "object"                                  category needed for result.json in WADQC
               },

               nii.hdr
               "qform_code gmslice":{
                  "type": "nii.hdr",                                    type of object
                  "sub-path": "dartel/",                                sub-path where to use the matching expression
                  "filename": "GMSlice*.nii.gz",                        filename matching expression
                  "category": "string",                                 category needed for result.json in WADQC
                  "field_name": "qform_code"                            nifti header filed name to select
               },

               nii.img
               "gmslice img":{
                  "type": "nii.img",                                    type of object
                  "sub-path": "dartel/",                                sub-path where to use the matching expression
                  "filename": "GMSlice*.nii.gz",                        filename matching expression
                  "category": "string",                                 category needed for result.json in WADQC
                  "fsl_slicers_options": "-A 1200"                      fsl slicer command options
               }


    -wad_qc_server: WADQC server configuration
      "wad_qc_server":{
        "ip_address": "localhost",
        "port":"11112",
        "ae-title": "dummy"
      }



Completed descriptor example:


{
  "dicom_meta":{
    "dicom_wrapper":{
        "sub-path": "",
        "filename":"*.dcm",
        "filepath_wadqc_placeholder":"path_placeholder/"
    }
  },
  "qc_items": {
    "Determinant_Orientation_T1w":{
      "type": "csv",
      "sub-path": "",
      "filename": "CheckOrientation_RigidRegT1.csv",
      "column": 0,
      "row": 0,
      "header": true,
      "category": "float"
    },
    "Total_Lesion_Volume_WMH(L)":{
      "type": "csv",
      "sub-path": "dartel/Tissue_Volume",
      "filename": "WMH_LST*.csv",
      "column": 0,
      "row": 0,
      "header": true,
      "category": "string"
    },
    "LR_flip_ASL(L)":{
      "type": "json",
      "sub-path": "dartel/",
      "filename": "QC_*.json",
      "child":"ASL/LR_flip_ASL",
      "category": "string"
    },
    "VoxelSize X":{
      "type": "json",
      "sub-path": "dartel/",
      "filename": "QC_*.json",
      "child":"ASL/VoxelSize/Z",
      "category": "string"
    },
    "file test":{
      "type": "file",
      "sub-path": "dartel/Tissue_Volume",
      "filename": "Tissue_*.csv",
      "category": "object"
    },
    "nii test_header":{
      "type": "nii.hdr",
      "sub-path": "dartel/",
      "filename": "GMSlice*.nii.gz",
      "category": "string",
      "field_name": "qform_code"
    }
  },
  "wad_qc_server":{
    "ip_address": "localhost",
    "port":"11112",
    "ae-title": "dummy"
  }
}