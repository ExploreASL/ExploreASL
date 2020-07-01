

## List of ExploreASL Functions

----
### xASL_adm_CatchNumbersFromString

#### Function
```matlab
function [OutputNumber] = xASL_adm_CatchNumbersFromString(InputString)
```

#### Description
...

----
### xASL_adm_CheckFileCount

#### Function
```matlab
function [result, files] = xASL_adm_CheckFileCount(path, expr, mincount, failifmissing)
```

#### Description
Checks the given PATH for files corresponding to the SPM_SELECT regular expression EXPR. Returns if the number of files is equal to or higher than MINCOUNT. If FAILIFMISSING is true and not enough files, then throw and error. If everything goes ok and second output argument is specified then return also the list of files.

-----
### xASL_adm_CheckPermissions

#### Function
```matlab
function [FilesList, FilesExeList, FoldersList] = xASL_adm_CheckPermissions(InputPath, FilesExecutable)
```

#### Description
This function does a recursive search through the root folder & makes a list of the attributes of all files and folders.
It tries to reset the attributes to what we desire, which is by default:

* 664 for files (meaning only reading & writing for users & group, & read-only for others) 
* 775 for folders (meaning reading, writing & opening for current user & current group, & for others only reading & opening) 

For executable files we also want 775. Note that the permission to 'execute a folder' means opening them

----
### xASL_adm_CheckSPM

#### Function
```matlab
function [spm_path, spm_version] = xASL_adm_CheckSPM(modality, proposed_spm_path, check_mode)
```

#### Description
Checks if the spm function exists and if the reported version matches our development version (SPM8 or SPM12). If the spm toolbox is not available yet, it will try the PROPOSED_SPM_PATH (if specified) or the user selected directory and add it to PATH. The function will fail if SPM cannot be found or if detecting an unsupported version.

----
### xASL_adm_CleanUpBeforeRerun

#### Function
```matlab
function xASL_adm_CleanUpBeforeRerun(AnalysisDir, iModule, bRemoveWMH, bAllSubjects, SubjectID, SessionID)
```

#### Description
This function (partly) reverts previous ExploreASL runs, deleting derivatives, while keeping raw data intact. if bAllSubjects==true, then all subjects and all module derivatives will be removed. This function performs the following steps:

1. If a Population folder doesn't exist yet but dartel does, rename it
2. Remove whole-study data files in AnalysisDir if bAllSubjects
3. Remove lock files/folders for reprocessing
4. Restore backupped _ORI (original) files
5. Delete native space CAT12 temporary folders (always, independent of iModule)
6. Remove native space files for iModule
7. Remove standard space files for iModule
8. Remove population module files
9. Remove or clean up stored x-struct & QC file -> THIS HAS NO SESSION SUPPORT YET

----
### xASL_adm_CompareDataSets

#### Function
```matlab
function [RMS] = xASL_adm_CompareDataSets(RefAnalysisRoot,SourceAnalysisRoot,x,type,mutexState)
```

#### Description
Compares two ExploreASL datasets for reproducibility.

----
### xASL_adm_CompareLists.m

#### Function
```matlab
function [NewList] = xASL_adm_CompareLists(list1, list2)
```

#### Description
Compare 2 single dimension lists.

----
### xASL_adm_ConvertDate2Nr.m

#### Function
```matlab
function [Nr DayInYear] = xASL_adm_ConvertDate2Nr(TempDate)
```

#### Description
Converts date to number input mmdd -> output mm (with days in fractions/floating point). Inverse from ConvertNrDate.

----
### xASL_adm_ConvertNr2Time.m

#### Function
```matlab
function Time = xASL_adm_ConvertNr2Time(Nr)
```

#### Description
Converts number to time input hh (with minutes in fractions/floating point) -> output hhmm. Inverse from xASL_adm_ConvertTime2Nr.

----
### xASL_adm_ConvertSubjSess2Subj_Sess.m

#### Function
```matlab
function [iSubj iSess] = xASL_adm_ConvertSubjSess2Subj_Sess(nSessions, iSubjSess)
```

#### Description
Converts combined SubjectSession index to subject & session indices. Useful for data lists in ExploreASL.

----
### xASL_adm_ConvertTime2Nr.m

#### Function
```matlab
function Nr = xASL_adm_ConvertTime2Nr(Time)
```

#### Description
Converts time to number input hhmm -> output hh (with minutes in fractions/floating point). Inverse from xASL_adm_ConvertNr2Time.

----
### xASL_adm_CopyMoveFileList.m

#### Function
```matlab
function [List] = xASL_adm_CopyMoveFileList(OriDir, DstDir, StrRegExp, bMove, bDir, bRecursive, bOverwrite, bVerbose)
```

#### Description
Moves a file to a file, a file to a directory, or a directory to a directory. It keeps the initial extensions, no unzipping or zipping after the move. But it makes sure that only one of .nii and .nii.gz exists in the destination directory. Useful to split a large database.

----
### xASL_adm_CorrectName.m

#### Function
```matlab
function strOut = xASL_adm_CorrectName(strIn, bOption, strExclude)
```

#### Description
Finds and replaces all non-word characters either by empty space or by an underscore. Optionally leaves in few selected special characters. Note that if '\_' is excluded from replacement, but option 2 is on, then underscores are replaced anyway.

----
### xASL_adm_CreateCSVfile.m

#### Function
```matlab
function xASL_adm_CreateCSVfile(CSVfilename,CSVdata)
```

#### Description
Creates a CSV file that can be opened with excel from your data.

----
### xASL_adm_CreateFileReport.m

#### Function
```matlab
function x = xASL_adm_CreateFileReport(x, bHasFLAIR, bHasMoCo, bHasM0, bHasLongitudinal)
```

#### Description
Prints a summary of created files or the individual modules (i.e. Structural, Longiutudinal & ASL modules).
Provides a quick check to see what has been skipped, an whether all files are present.

This script iterates across Native space 1) subject and 2) session files, resampled 3) subject and 4) session files, 5) Lock files and 6) QC Figure files.

For all we perform a A) count of the files present, summarized in fileReportSummary.csv, and we B) list the missing files in
"Missing\*.csv" files

**PM:** simplify/optimize this code, to make filename variable changing, search within subject-directories, etc. Combine the parts searching for missing & summarizing count

----
### xASL_adm_DefineASLResolution.m

#### Function
```matlab
function x = xASL_adm_DefineASLResolution(x)
```

#### Description
...

----
### xASL_adm_DeleteFilePair.m

#### Function
```matlab
function filepaths = xASL_adm_DeleteFilePair(path, varargin)
```

#### Description
Delete the file given in PATH, and also deletes files with the same name, but with extension given in EXT1, and potentially also EXT2, EXT3... 

----
### xASL_adm_Dicom2Parms.m

#### Function
```matlab
function [parms, pathDcmDictOut] = xASL_adm_Dicom2Parms(imPar, inp, parmsfile, dcmExtFilter, bUseDCMTK, pathDcmDictIn)
```

#### Description
The function goes through the INP files, reads the DICOM or PAR/REC files and parses their headers.
It extracts the DICOM parameters important for ASL, makes sure they are in the correct format, if missing then replaces with default value, it also checks if the parameters are consistent across DICOM files for a single sequence.

----
### xASL_adm_FindByRegExp.m

#### Function
```matlab
function [tree, optionalTokens] = xASL_adm_FindByRegExp(root, dirSpecs, varargin)
```

#### Description
Recursively find files in the root directory according to the dirSpecs.

----
### xASL_adm_FindStrIndex.m

#### Function
```matlab
function INDEX = xASL_adm_FindStrIndex(ARRAY, STRING)
```

#### Description
Similar to find, but then for a cell array filled with strings. Only takes 4 dimensions.

----
### xASL_adm_GetFsList.m

#### Function
```matlab
function  RES = xASL_adm_GetFsList(strDirectory, strRegEx, bGetDirNames, bExcludeHidden, bIgnoreCase, nRequired)
```

#### Description
List files or directories from a given path. And optionally uses regular expressions to filter the result with options to exclude hidden files, ignore case, and set a minimal requirement on the number of results. Sorts the results at the end.

----
### xASL_adm_GetNumFromStr.m

#### Function
```matlab
function num = xASL_adm_GetNumFromStr(str)
```

#### Description
Obtains single number from string. CAVE there should only be one number!

----
### xASL_adm_GetPhilipsScaling.m

#### Function
```matlab
function scaleFactor = xASL_adm_GetPhilipsScaling(parms,header)
```

#### Description
This script provides the correct scaling factors for a NIfTI file. It checks the header of the NIfTI that normally has the same scaling as RescaleSlope in DICOM, it checks if dcm2nii (by the info in JSON) has already converted the scale slopes to floating point. And if not, the derive the correct scaling factor to be applied.

----
### xASL_adm_GetUserName.m

#### Function
```matlab
function UserName = xASL_adm_GetUserName()
```

#### Description
...

----
### xASL_adm_Hex2Num.m

#### Function
```matlab
function outNum = xASL_adm_Hex2Num(inStr, type, endian)
```

#### Description
Takes a hexadecimal string and converts it to number. Works also when the string contains escape characters, and for single-floats and for a little and big endian. If containing 8 and less characters than treat as float, if more than as double.

----
### xASL_adm_LesionResliceList.m

#### Function
```matlab
function [INname, OUTname] = xASL_wrp_LesionResliceList(x,bLesion_T1,bLesion_FLAIR,bROI_T1,bROI_FLAIR)
```

#### Description
Creates list of structural image paths to reslice.

----
### xASL_adm_Load4DMemMapping.m

#### Function
```matlab
function LoadFile = xASL_adm_Load4DMemMapping(x, WhichModality)
```

#### Description
Part of ExploreASL analysis module. Loads data & maps it to memory mapping file on disc, if not done before.

----
### xASL_adm_LoadParms.m

#### Function
```matlab
function [Parms, x, Oldx] = xASL_adm_LoadParms(ParmsPath, x, bVerbose)
```

#### Description
This function loads the internal memory x struct, any legacy \*\_parms.mat sidecar, any \*.json BIDS sidecar, to use scan-specific parameters for image processing/quantification. Also, per BIDS inheritance, any x.S.SetsID parameters (from participants.tsv) are loaded as well. This function performs the following steps:

1. Load .mat parameter file
2. Load JSON file
3. Deal with warnings
4. Find fields with scan-specific data in x.S.Sets, and use this if possible (per BIDS inheritance)
5. Sync Parms.\* with x.(Q.)\* (overwrite x/x.Q)
6. Fix M0 parameter if not set

----
### xASL_adm_LoadX.m

#### Function
```matlab
function [x, IsLoaded] = xASL_adm_LoadX(x, Path_xASL, bOverwrite)
```

#### Description
This function loads x.Output & x.Output_im struct fields from the x.mat on the hard drive & adds them to the current x struct located in memory. If it didnt exist in the x.mat, it will set IsLoaded to false, which can be catched externally & a warning issued if managed so in the calling function. If it didnt exist in the memory x struct, or bOverwrite was requested, the contents of x.mat will be loaded to the memory x struct.

----
### xASL_adm_OrderFields.m

#### Function
```matlab
function outStruct = xASL_adm_OrderFields(inStruct,orderStruct)
```

#### Description
Order fields in the structure inStruct to match orderStruct, unmatching fields in inStruct are copied as they are at the end, unmatching fields in orderStruct are ignored. This is just a cosmetic change and no values are edited.

----
### xASL_adm_OtherListSPM.m

#### Function
```matlab
function [OtherListSPM, OtherListOut] = xASL_adm_OtherListSPM(OtherList, bList4D)
```

#### Description
Takes care of the others list for registration functions. 

**bPadComma1:** is to add the ,1 to the end of the pathstring, which SPM uses to assign the first image of a 4D image array (OPTIONAL, DEFAULT = true)

**bList4D:** boolean, true for listing multiple 4D volumes separately in the list (OPTIONAL, DEFAULT=true).

----
### xASL_adm_Par2Parms.m

#### Function
```matlab
function parms = xASL_adm_Par2Parms(pathPar, pathParms, bRecreate)
```

#### Description
Opens the Philips type PAR file. Reads the relevant DICOM headers and saves them to .MAT file. Only recreates an existing file if bRecreate option is set to TRUE.

----
### xASL_adm_ParReadHeader.m

#### Function
```matlab
function info = xASL_adm_ParReadHeader(filename)
```

#### Description
Function for reading the header of a Philips Par / Rec  MR V4.\* file.

----
### xASL_adm_Remove_1_SPM.m

#### Function
```matlab
function [OtherList] = xASL_adm_Remove_1_SPM(OtherList)
```

#### Description
Remove, 1 at end of OtherLists, if exists. These are appended in CoregInit, OldNormalizeWrapper etc, since this should allow 4rd dim (e.g. as in ASL4D).

----
### xASL_adm_ReplaceSymbols.m

#### Function
```matlab
function [x] = xASL_adm_ResetVisualizationSlices(x)
```

#### Description
Removes any predefined slices that should be visualized, allowing to show the default slices. Comes in handy when different pipeline visualization parts are repeated.


----
### xASL_adm_ResetVisualizationSlices.m

#### Function
```matlab
function [x] = xASL_adm_ResetVisualizationSlices(x)
```

#### Description
Removes any predefined slices that should be visualized, allowing to show the default slices. Comes in handy when different pipeline visualization parts are repeated.

----
### xASL_adm_SaveJSON.m

#### Function
```matlab
function xASL_adm_SaveJSON(data, jsonFileName)
```

#### Description
Saves the values in the structure 'data' to a file in JSON forma

----
### xASL_adm_uiGetInput.m

#### Function
```matlab
function [Parms] = xASL_adm_uiGetInput(Parms)
```

#### Description
Checks whether input fields are present, or requests them.

----
### xASL_adm_UnzipOrCopy.m

#### Function
```matlab
function unpackedFiles = xASL_adm_UnzipOrCopy(srcDir, wildCard, destDir, bOverwrite)
```

#### Description
Unpacks (or copy if unpacked) one or more files matching the regular expression.

----
### xASL_adm_Voxel2RealWorldCoordinates.m

#### Function
```matlab
function [X Y Z] = xASL_adm_Voxel2RealWorldCoordinates(X,Y,Z,VoxelSize)
```

#### Description
Converts MNI coordinates from voxel coordinates/indices.
Assumes X Y Z = LR LeftRight AP AnteriorPosterior IS InferiorSuperior.
VoxelSize should be [1 3]-sized input.

----
### xASL_adm_ZipFileList.m

#### Function
```matlab
function filepaths = xASL_adm_ZipFileList(strDirectory, strRegExp, bRecurse, bUseGzip, nRequired, bDelete)
```

#### Description
Zip the files that match regular expression STRREGEXP in the given directory STRDIRECTORY.
Zips recursively if specified in BRECURSE. Zips all files unless the number is specified by NREQUIRED, if the number is not met, then does not zip anything and throws an error.

----
### xASL_bids_Add2ParticipantsTSV.m

#### Function
```matlab
function xASL_bids_Add2ParticipantsTSV(DataIn, DataName, x, bOverwrite)
```

#### Description
This function adds metadata/statistical variables to the participants.tsv in the root/analysis folder, by the following steps. 
This function will iterate over Data provided at DataIn and fill them in the participants.tsv, overwriting if allowed.
Empty data is filled in as 'n/a', and the first column "participants_id" is sorted for participants.

----
### xASL_bids_Dicom2JSON.m

#### Function
```matlab
function [parms, pathDcmDictOut] = xASL_bids_Dicom2JSON(imPar, inp, PathJSON, dcmExtFilter, bUseDCMTK, pathDcmDictIn)
```

#### Description
The function goes through the INP files, reads the DICOM or PAR/REC files and parses their headers.
It extracts the DICOM parameters important for ASL, makes sure they are in the correct format, if missing then replaces with default value, it also checks if the parameters are consistent across DICOM files for a single sequence.

----
### xASL_bids_InsertJSONFields.m

#### Function
```matlab
function [ChildJSON] = xASL_bids_InsertJSONFields(ParentJSON, ChildJSON, Fields2Skip)
```

#### Description
This function takes all parameters from the "parent" JSON & moves them into the "child" JSON.
In case of co-existence of a field with different values, then the value in the child JSON will prevail, per BIDS inheritance.

----
### xASL_bids_parms2BIDS.m

#### Function
```matlab
function outParms = xASL_bids_parms2BIDS(inXasl, inBids, bOutBids, bPriorityBids)
```

#### Description
This functions takes two parameter structures and merges them. At the same time, renames all fields according to the output type (note that only some fields have two standardised names different between the two formats.

In case of duplicities, takes the field value from the preferred format. 
Also takes into account that the units in BIDS are s, but in xASL ms.
This function performs the following steps:

1. Define field names that need to be convert/renamed/merged
2. Convert XASL fields to the output format (BIDS or XASL)
3. Convert BIDS fields to the output format (BIDS or XASL)
4. Merge the BIDS and XASL fields, convert field values

----
### xASL_bids_PARREC2JSON.m

#### Function
```matlab
function parms = xASL_bids_PARREC2JSON(pathPar, PathJSON)
```

#### Description
Opens the Philips type PAR file. Reads the relevant DICOM header fields and adds them to the .json sidecar file.

----
### xASL_fsl_RunFSL.m

#### Function
```matlab
function [x, Result1] = xASL_fsl_RunFSL(FSLCommand, x, OutputZipping, NicenessValue, bVerbose)
```

#### Description
This function runs an FSL command from ExploreASL:

1. Checking the FSL dir
2. Manage CUDA/CPU parallelization (currently disabled, WIP)
3. Setting up FSL environment
4. Running the command

Supports .nii & .nii.gz, Linux, MacOS & Windows (WSL).

----
### xASL_fsl_SetFSLdir.m

#### Function
```matlab
function [FSLdir, x, RootWSLdir] = xASL_fsl_SetFSLdir(x, bAutomaticallyDetectFSL)
```

#### Description
This function finds the FSLdir & puts it out, also in x.FSLdir to allow repeating this function without having to repeat searching.
If the FSLdir & RootFSLDir are already defined in x.FSLdir & x.RootFSLDir, this function is skipped.
Supports Linux, MacOS & Windows (WSL), & several different default installation folders for different Linux distributions.

----
### xASL_fsl_TopUp.m

#### Function
```matlab
function [bSuccess] = xASL_fsl_TopUp(InDir, ScanType, x, OutputPath)
```

#### Description
This function runs FSL TopUp. It assumes that there are 2 TopUp images, i.e. 1 blip up & 1 blip down.

0. Admin: manage ScanType, NIfTI paths, create TopUp parameter file for image to apply TopUp to & for the TopUp NIfTIs, delete files from previous run, define the image with the same acquisition parameters as TopUp (does the image we apply TopUp to, have the Blip up or down?)
1. Register images to image that we apply TopUp to (registration between blip up/down images is performed by TopUp)
2. Run TopUp estimate (i.e. estimate the geometric distortion field from B0 NIfTI & parameters file), this takes quite long. Also has a x.Quality=0 option that is very fast but inaccurate, to try out this pipeline part. Before TopUp, NaNs (e.g. from resampling) are removed from the images TopUp is run with default settings
3. Apply TopUp

----
### xASL_import_json.m

#### Function
```matlab
function [x] = xASL_import_json(DataParFile)
```

#### Description
This function reads in a DATA_PAR file and creates the x structure. The name of the DATA_PAR file is given as a string or character array. The output is the x structure.

If the DATA_PAR file is the dataset_description.json file of the BIDS standard, the x structure is created according to BIDS.

----
### xASL_imwrite.m

#### Function
```matlab
function [ImOut] = xASL_imwrite(ImIn, PathOut, ColorMap, bRescale)
```

#### Description
This functions takes an input image matrix, interpolates it to HD resolution (1920x1080) for visibility, and saves the image as jpg. This function avoids the graphic interface of Matlab, for running from CLI. 
**Careful:** this function overwrites any existing PathOut.

----
### xASL_im_AddIM2QC.m

#### Function
```matlab
function [x] = xASL_im_AddIM2QC(x,parms)
```

#### Description
Checks which images already are loaded, and  adds new image.

----
### xASL_im_BilateralFilter.m

#### Function
```matlab
function [ovol] = xASL_im_BilateralFilter(volIM, mask, VoxelSize, x)
```

#### Description
This function runs a spatial lowpass temporally highpass filter, and removes outliers within this signal, and adapts the time-series accordingly

----
### xASL_im_CenterOfMass.m

#### Function
```matlab
function xASL_im_CenterOfMass(PathNIfTI, OtherList, AllowedDistance)
```

#### Description
This function estimates the center of mass of the image matrix, and if this is too far off the current orientation matrix center, the center will be reset.
This fixes any incorrect orientation outputted by the scanner.
The realignment is only applied when any of the X/Y/Z dimensions have a higher offset than AllowedDistance

----
### xASL_im_CleanupWMHnoise.m

#### Function
```matlab
function xASL_im_CleanupWMHnoise(InputPath, OutputPath, MinLesionVolume, pThresh)
```

#### Description
Threshold white matter lesions, acknowledging the fact that they may be confluent with subresolution connection through a dilation. This part is executed conservatively, as FLAIR hyperintensities inside the GM can be erroneously segmented as WMH, and should not be lesion-filled (otherwise these cannot be fixed later in the Structural module).

Note that LST lesion filling expects a probability map, doesnt work nicely with binary mask.

----
### xASL_im_ClipExtremes.m

#### Function
```matlab
function [NewIM] = xASL_im_ClipExtremes(InputIm, ThreshHigh, ThreshLow, bVerbose)
```

#### Description
Clips image to given percentile. The percentile is found using non-zeros sorted intensities, so both isfinite & non-zeros.

----
### xASL_im_colorbar.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_Column2IM.m

#### Function
```matlab
function [ImageOut] = xASL_im_Column2IM(ColumnIn, BrainMask)
```

#### Description
This function "decompresses" an image matrix (or multiple matrices) from a single-dimensional column, by reconstructing the image matrix from the voxel positions within the BrainMask. 

**NB:** Important to use the same BrainMask as used for converting the image matrix to the column!

**See also:** xASL_im_IM2Column.m

The mask mostly used for xASL_im_IM2Column is x.WBmask, which completely engulfes pGM, pWM & pCSF.

----
### xASL_im_CompareNIfTIResolutionXYZ.m

#### Function
```matlab
function [IsEqualResolution] = xASL_im_CompareNIfTIResolutionXYZ(PathNIfTI1, PathNIfTI2)
```

#### Description
This function checks whether the X, Y and Z resolution of a NIfTI with any number of dimensions is equal. It rounds for 2 floating points, for both NIfTIs, to ensure that the same precision is compared.

----
### xASL_im_ComputeDice.m

#### Function
```matlab
function DiceCoeff = xASL_im_ComputeDice(imA, imB)
```

#### Description
Calculate Dice coefficient of image overlap.

----
### xASL_im_CreateASLDeformationField.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_CreatePseudoCBF.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_CreateSliceGradient.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_CreateVisualFig.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_CreateVisualLongReg.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_CropParmsAcquire.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_CropParmsApply.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_DecomposeAffineTransformation.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_DetermineFlip.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_DilateErodeFull.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_DilateErodeSeparable.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_DilateErodeSphere.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_dilateROI.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_DummyOrientationNIfTI.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_EstimateResolution.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_Flip.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_FlipOrientation.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_FlipOrientation2.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_IM2Column.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_joinColormap.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_JointHist.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_Lesion2CAT.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_Lesion2Mask.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_M0ErodeSmoothExtrapolate.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_MaskNegativeVascularSignal.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_MaskPeakVascularSignal.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_Modulation.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_NormalizeLabelingTerritories.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_OverlapT1_ASL.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_PCA.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_PreSmooth.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_ProcessM0Conventional.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_ProjectLabelsOverData.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_PVCbspline.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_PVCkernel.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_ResampleLinearFair.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_RescaleLinear.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_RestoreOrientation.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_rotate.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_SkullStrip.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_Smooth3D.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_TileImages.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_TransformData2View.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_Upsample.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_VisualizeROIs.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_VisualQC_TopUp.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_im_ZeroEdges.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_init_ConvertM2JSON.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_init_DefaultEffectiveResolution.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_init_DefineStudyData.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_init_FileSystem.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_init_InitializeMutex.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_init_LabelColors.mat

#### Function
```matlab
...
```

#### Description
...

----
### xASL_init_LoadMetadata.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_init_LongitudinalRegistration.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_init_VisualizationSettings.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_io_CreateNifti.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_io_dcm2nii.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_io_DcmtkRead.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_io_MakeNifti4DICOM.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_io_PairwiseSubtraction.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_io_ReadTheDicom.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_io_SplitASL_M0.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_Iteration.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_num2str.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_AsymmetryIndex.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_CAT12_IQR.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_CollectParameters.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_CollectQC_ASL.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_CollectQC_func.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_CollectQC_Structural.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_CollectSoftwareVersions.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_CompareTemplate.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_ComputeFoVCoverage.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_ComputeNiftiOrientation.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_CreatePDF.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_FA_Outliers.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_ObtainQCCategoriesFromJPG.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_PCPStructural.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_PrintOrientation.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_TanimotoCoeff.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_temporalSNR.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_WADQCDC.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_qc_WADQC_GenerateDescriptor.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_quant_AgeSex2Hct.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_quant_FEAST.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_quant_GetControlLabelOrder.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_quant_Hct2BloodT1.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_quant_M0.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_quant_SinglePLD.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_spm_affine.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_spm_BiasfieldCorrection.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_spm_coreg.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_spm_deface.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_spm_deformations.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_AtlasForStats.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_ComputeDifferCoV.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_ComputeMean.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_ComputeSpatialCoV.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_EqualVariancesTest.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_fcdf.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_GetROIstatistics.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_MadNan.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_MeanSSIM.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_MultipleLinReg.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_PrintStats.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_PSNR.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_QuantileNan.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_RobustMean.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_ShapiroWilk.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_StdNan.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_SumNan.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_tcdf.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_ticdf.m

#### Function
```matlab
...
```

#### Description
...

----
### xASL_stat_ttest.m

#### Function
```matlab
function [H,P,CI,stats] = xASL_stat_ttest(X,M,alpha,tail,dim)
```

#### Description
Performs a t-test that the distribution of the input data X has a mean different from 0 (or from a given mean M, or that the distributions X and Y have different means). A normal distribution of the data with an unknown variance is assumed.

----
### xASL_stat_ttest2.m

#### Function
```matlab
function [H,P,CI,stats] = xASL_stat_ttest2(X,Y,alpha,tail,vartype,dim)
```

#### Description
Performs a unpaired t-test that the distribution of the input data X has a mean different from that of Y.  A normal distribution of the data with an unknown variance is assumed.

----
### xASL_stat_VarNan.m

#### Function
```matlab
function y = xASL_stat_VarNan(x,w,dim)
```

#### Description
Calculates variance of values in X while ignoring NaNs.

----
### xASL_str2num.m

#### Function
```matlab
function [DataOut] = xASL_str2num(DataIn)
```

#### Description
str2num wrapper, replacing 'n/a' with NaN (BIDS convention) and converting only strings to numbers. Also allows inputting cells.


