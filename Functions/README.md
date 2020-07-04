

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
### xASL_vis_Imwrite.m

#### Function
```matlab
function [ImOut] = xASL_vis_Imwrite(ImIn, PathOut, ColorMap, bRescale)
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
function axesBar = xASL_io_colorbar(ColorMap,TickLabels)
```

#### Description
Display a color bar, but replace it with a colorbar with given 'ticks' and colormap.

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
function xASL_im_CreateASLDeformationField(x, bOverwrite, EstimatedResolution, PathLowResNIfTI)
```

#### Description
This function smooths a transformation flow field to a lower resolution. Usually, we use a high resolution anatomical image (e.g. 3D T1w) to obtain the flowfields from native space to standard space, and apply these to the lower resolution ASL images. Because of the resolution differences, the flowfields need to be downsampled/smoothed, to avoid deformation effects that are crispier than the functional image that is investigated. This function performs the following steps:

1. Obtain resolutions
2. Fill NaNs at edges y_T1.nii flowfield to prevent interpolation artifact
3. Smooth flowfield
4. Fill NaNs at edges y_ASL.nii

Note that if the resolution of ASL is not significantly (i.e. >0.5 mm in any dimension) lower than T1w, the y_T1.nii is copied to y_ASL.nii

----
### xASL_im_CreatePseudoCBF.m

#### Function
```matlab
function xASL_im_CreatePseudoCBF(x, spatialCoV)
```

#### Description
This function creates a pseudo-CBF image from mean CBF template, arterial transit time (ATT) bias field & vascular artifacts, weighted through spatial CoV.
The first part of this code puts templates in the native space and creates a pseudoCBF image from a downsampled pGM & pWM tissue (PseudoTissue). The latter is used for registration but also as reference for the template registration, to speed this up.
The second part of this code computes a pseudoCBF image based on the pseudoTissue & the CBF templates of CBF, ATT biasfield and vascular peaks, based on spatial CoV.
This submodule performs the following steps:

1. Create the pseudoTissue CBF reference image, if it doesnt exist already
2. Create the native space copies of ASL templates, if they dont exist already
3. Spatial CoV input argument check
4. Load native space copies of templates
5. Create pseudoTissue from segmentations, mix this with the mean CBF template depending on spatial CoV
6. Create pseudoCBF reference image used for CBF-based registration
7. Scale mean_PWI_Clipped source image to the same range as PseudoCBF

----
### xASL_im_CreateSliceGradient.m

#### Function
```matlab
function xASL_im_CreateSliceGradient(x)
```

#### Description

When a 2D readout is used with ASL, post-label delay and hence T1 decay will be dependent on slice timing.
Therefore, quantification part needs slice reference to quantify per slice and correct for effective post-label delay differences

This function uses exact same ASL matrix changes that occurred due to registration to MNI, motion correction and registration to T1.

1.    Create slice gradient in same space as input file
2.    Reslice slice gradient to MNI (using existing ASL matrix changes from e.g. registration to MNI, motion correction, registration to GM)
3.    Creating average slice gradient

----
### xASL_im_CreateVisualFig.m

#### Function
```matlab
function [ImOut, FileName] = xASL_im_CreateVisualFig(x, ImIn, DirOut, IntScale, NamePrefix, ColorMap, bClip, MaskIn, bWhite, MaxWindow, bTransparancy)
```

#### Description
This function creates a visualization Figure by merging flexibly rearranging NIfTI slices, input matrix or path, managing colormaps for different merged image layers. Current use is for visual QC figures and overview in papers.
Function is structured as:

1. Admin, deal with input arguments
2. Process image layers separately
    * xASL_vis_TransformData2View: Reshapes image data into visualization figure
    * xASL_im_ClipExtremes: Clips image to given percentile also we scale for peak intensity, we make sure that there is no visible clipping/distortion
    * Convert to colors, using any input colormaps
3. Combine image layers, using input argument IntScale
4. Print figure

This function assumes that the first image is a grascale background image (e.g. for transparancy reasons), if there are multiple images.

----
### xASL_vis_CropParmsAcquire.m

#### Function
```matlab
function [xmin xmax ymin ymax] = xASL_vis_CropParmsAcquire(temp_image)
```

#### Description
Goes from outside to inside to acquire crop settings. Works with grayscale images (2 dimensions per slice).image position information (2D matrix) should be first 2 dimensions. Could include colordimension later on.

----
### xASL_vis_CropParmsApply.m

#### Function
```matlab
function ImageOut = xASL_vis_CropParmsApply(ImageIn,CropParameters,Xmax,Ymin,Ymax)
```

#### Description
This function crops 2D image matrices.

----
### xASL_im_DecomposeAffineTransformation.m

#### Function
```matlab
function [M, P] = xASL_im_DecomposeAffineTransformation(Mtransformation)
```

#### Description
This function splits a transformation matrix into individual components, which can be useful to guide the SPM reslicing.
The components are the same as in spm_(i)matrix.m, except for the shearing: these are included in the rotations, and the 90 degree rotations, these are separated.

Reason for the separation of the 90 degree rotations, is that these indicate if orientations (transversal, coronal & sagittal) have been switched in the NIfTI.

This can be useful to correct for any erroneous 90degree rotations from registration, or to put two images in the same orientation order or voxelsize without applying their subtle realignment (e.g. for manipulating registration references).

----
### xASL_im_DetermineFlip.m

#### Function
```matlab
function [QCstruct] = xASL_im_DetermineFlip(x,iS,PathOrientationResults,QCstruct)
```

#### Description
Check determinants, should be the same before & after registration, otherwise a left-right flip is applied. This is not visible, but detrimental for image analysis/stats.

----
### xASL_im_DilateErodeFull.m

#### Function
```matlab
function new_mask = xASL_im_DilateErodeFull(mask,type,kernel)
```

#### Description
Runs dilation or erosion on a binary mask in full three dimensions.
It uses its own dilate_erode function and crops the image so that it contains only the mask

----
### xASL_im_DilateErodeSeparable.m

#### Function
```matlab
function new_mask = xASL_im_DilateErodeSeparable(mask, type, kernel_x, kernel_y, kernel_z)
```

#### Description
Runs dilation or erosion on a binary mask separably in three dimensions.
It uses its own dilate_erode function and crops the image so that it contains only the mask

----
### xASL_im_DilateErodeSphere.m

#### Function
```matlab
function el = xASL_im_DilateErodeSphere(R)
```

#### Description
3D structuring element (binary) sphere.

----
### xASL_im_dilateROI.m

#### Function
```matlab
function xASL_im_dilateROI(PathIn, PathTemp)
```

#### Description
...

----
### xASL_im_DummyOrientationNIfTI.m

#### Function
```matlab
function xASL_im_DummyOrientationNIfTI(PathSrc, PathRef, PathDummyOut, bApplyRotationMinor, bApplyRotation90, bApplyZoom, bApplyTranslation)

```

#### Description
This function creates a dummy image as reference for xASL_spm_reslice, allowing to only apply specific parts of the transformation between the two images. E.g. only the rotation, or only the zooming.
This can be useful to correct for any erroneous rotations from registration, or to put two images in the same space without applying their realignment. This function performs the following steps:

1. Load orientations & calculate transformation
2. Calculate the desired transformation
3. Calculate new orientation matrix
4. Calculate the new image size
5. Save the dummy NIfTI

----
### xASL_im_EstimateResolution.m

#### Function
```matlab
function [ resFWHM, resSigma,resErr,imSmo,imMask] = xASL_im_EstimateResolution(imCBF,imGM,imWM,imMaskOrig,PSFtype,maxIter)
```

#### Description
NB: everything in this code is in voxels, not in mm.

----
### xASL_im_Flip.m

#### Function
```matlab
function [MatrixOut] = xASL_im_Flip(MatrixIn, varargin)
```

#### Description
Backwards compatibility for flipping left-right in standard space (NB: this can be different than in native space!).

----
### xASL_im_IM2Column.m

#### Function
```matlab
function [ColumnOut] = xASL_im_IM2Column(ImageIn, BrainMask, ApplyShiftDim)
```

#### Description
QC Converts an image matrix to a single-dimensional column to save memory space & computation time.

----
### xASL_vis_joinColormap.m

#### Function
```matlab
function joinedColormap = xASL_io_joinColormap(blackMiddle,colormap1,colormap2)
```

#### Description
Take two colormaps and put them together, creating a colormap with 256 values input colormaps 1 and 2 are stacked as colormap1 first, then colormap2.
They can be of anylength [N1,3] and [N2,3], joinedColormap has size [256,3].

----
### xASL_im_JointHist.m

#### Function
```matlab
function imHist = xASL_im_JointHist(imA,imB,imMask,minA,maxA,minB,maxB,nBins)
```

#### Description
It calculates a joint histogram of two images of any dimensions over a mask of the same size. The boundaries and number of bins can either be given or min and max values are used. Values outside of the bins are counted to the first/last bin.

----
### xASL_im_Lesion2CAT.m

#### Function
```matlab
function LesionPathOut = xASL_im_Lesion2CAT(PathIn)
```

#### Description
For all lesion masks in the anatomical directory, load them, merge them and save them for the CAT segmentation.

----
### xASL_im_Lesion2Mask.m

#### Function
```matlab
function LesionIM = xASL_im_Lesion2Mask(LesionPath, T1path, pGMpath, pWMpath, x)
```

#### Description
For a standard space lesion mask (or map), this stores the lesion mask, and in additional its perimask (15 mm) and contralateral mask, as 2nd and 3rd volumes.
It plots the masks on a T1 image, and masks the new masks with the subjects' brainmask (pGM+pWM).

----
### xASL_im_M0ErodeSmoothExtrapolate.m

#### Function
```matlab
function [ImOut] = xASL_im_M0ErodeSmoothExtrapolate(ImIn, x)
```

#### Description
This function erodes, smooths & extrapolates M0 in standard space.
It assumes that the M0 image is in standard space & that the GM & WM probability maps are aligned. Here, we mask the M0, to remove high CSF signal and low extracranial signal, enabling us to smooth the image without inserting wrong signal. See also the ExploreASL manuscript for a more extensive explanation. This function performs the following steps:

* Mask 1) Load segmentations, create structural mask
* Mask 2) Create intensity-based mask to remove extracranial signal
* Mask 3) Erode the combined masks
* Mask 4) Determine any odd borders
* 5) Smoothing
* 6) Extrapolating only
* 7) Scale back to the GM M0
* 8) Print visual QC figure

A visual QC figure is created, showing the M0 image processing steps for a single transversal slice (slice 53 in 1.5 mm MNI standard space).

OutputFile = fullfile(x.D.M0regASLdir, ['M0_im_proc_' x.P.SubjectID '.jpg']);
The original M0 image (a) is masked with a (pGM+pWM)>50% mask (b) eroded with a two-voxel sphere to limit the influence of the ventricular and extracranial signal (c) and thresholded to exclude significantly high (i.e. median + 3 \* mean absolute deviation (MAD)) border region values (d) This masked M0 image is smoothed with a 16x16x16 mm full- width-half-maximum Gaussian filter (Mutsaerts et al., 2018) (e) after which the signal at the border is smoothly extrapolated until the full image is filled (f).
             

----
### xASL_im_MaskNegativeVascularSignal.m

#### Function
```matlab
function [NegativeMask, TreatedCBF] = xASL_im_MaskNegativeVascularSignal(x, IsSpace)
```

#### Description
This function segments clusters with significant negative ASL signal. This can be tricky as there is also the negative tail of Gaussian noise from the ASL subtraction. The image feature we use here, is that negative vascular signal will be a relatively large region with significant median negative value, whereas noise will be regions with relatively small negative signal.
Negative signal from wrong background suppression timing (e.g. in the first slice with 2D EPI) can be masked out with this as well.
The procedure works as follows:

1. Obtain mask of negative voxels within pGM>0.5 mask
2. Obtain distribution of subzero clusters
3. Define the negative threshold
4. Create mask by thresholding whole image

Note that the definition of the threshold is obtained within the GM only, but that this threshold is applied to the full image.

----
### xASL_im_MaskPeakVascularSignal.m

#### Function
```matlab
function [MaskIM, TreatedCBF] = xASL_im_MaskPeakVascularSignal(PathCBF, Path_M0, bClip, ClipThresholdValue, CompressionRate)
```

#### Description
This function searches for an acceptable high threshold as definition of high intra-vascular ASL signal.
It also allows to compress the values here (when
bClip==true). Compression retains some variability, but limits their outlying influence on statistics. The procedure works as follows:

1. Segment GM on ASL image at 80th percentile of CBF image distribution
2. For PWI & CBF images, select voxels higher than median + ClipThresholdValue MAD Vascular artifacts will have a high intensity in both images, whereas errors by division by M0 will only have a high intensity on the M0 image, and high values due to a biasfield will only be noticeable on the PWI image
3. Combine the two created masks
4. Obtain average clipping value from selected voxels from the combined masks
5. Apply compression if requested. If not, output image will have NaNs for intra-vascular voxels.

Note that the definition of the threshold is obtained within the GM only, but that this threshold is applied to the full image to also remove extracranial extreme values

----
### xASL_im_Modulation.m

#### Function
```matlab
function xASL_im_Modulation(x)
```

#### Description
Combines the transformations to create Jacobians, & multiplies the standard space segmentations with these to create volumetric images for volumetric analyses.

----
### xASL_im_NormalizeLabelingTerritories.m

#### Function
```matlab
function image_out = xASL_im_NormalizeLabelingTerritories(imageIN, GMmask, x)
```

#### Description
Normalizes per perfusion territory mask should be GM mask.

----
### xASL_vis_OverlapT1_ASL.m

#### Function
```matlab
function xASL_vis_OverlapT1_ASL(x, ASL)
```

#### Description
Part of ExploreASL. Shows spatial agreement ASL and probability maps.

----
### xASL_im_PCA.m

#### Function
```matlab
function [pc, score, eigenvalues, tsquare, loadings, Xmean] = xASL_im_PCA(dataIn)
```

#### Description
...

----
### xASL_im_PreSmooth.m

#### Function
```matlab
function pathOut = xASL_im_PreSmooth(pathRef, pathSrc, pathSmo, resRef, resSrc, srcAffinePath, bInvAffine)
```

#### Description
It assumes that the FWHM is equal to voxel size, unless the real resolution is given. Then takes into account the voxel sizes and orientation difference between the volumes, but performs the smoothing according to the given real resolution (it is possible to supply the resolution for just one image) - this should be helpful primarily when the. It creates a Gaussian covariance matrix for the reference, transforms it according to the difference between the two images a produces the Gaussian covariance matrix of the pre-smoothing in the source space. Then it performs the smoothing.
The following steps are performed:

1. Obtain the voxel size
2. Skip this function if reference resolution is equal to, or lower than source resolution
3. Deal with affine transformation
4. Obtain the transformation matrix from the Reference to the Source space
5. Apply the smoothing filter on the source image(s)
6. Save the smoothed image

----
### xASL_im_ProcessM0Conventional.m

#### Function
```matlab
function [Corr_M0] = xASL_im_ProcessM0Conventional(ImIn, x)
```

#### Description
This function uses the conventional M0 masking, and only a little smoothing, following what Philips uses for its 3D GRASE. Advantages of the newer M0 processing in ExploreASL are the lack of use of M0 threshold-based masking, the removal of high CSF values and higher SNR for ASL division.

----
### xASL_im_ProjectLabelsOverData.m

#### Function
```matlab
function OutputIM = xASL_im_ProjectLabelsOverData(DataIM, LabelIM, x, ScaleFactorData, ScaleFactorLabel)
```

#### Description
This script projects labels over an image, but works only in 2D. Make sure to make a 2D image from a 3D or 4D image using xASL_vis_TransformData2View.m can be used in combination with xASL_vis_Imwrite.m

----
### xASL_im_PVCbspline.m

#### Function
```matlab
function [imPVEC,imCBFrec,imResidual,FWHM] = xASL_im_PVCbspline(imCBF, imPV, bsplineNum)
```

#### Description
PVEc correction of ASL data using prior GM-,WM-partial volume maps.
Follows the principles of the PVEc algorithm by I. Asllani (MRM, 2008).
The PV-corrected CBF_GM and CBF_WM maps are approximated using an uniformly sampled cubic B-splines.
Free for research use without guarantee. 

----
### xASL_im_PVCkernel.m

#### Function
```matlab
function [imPVEC,imCBFrec,imResidual] = xASL_im_PVCkernel(imCBF, imPV,kernel,mode)
```

#### Description
PVEc correction of ASL data using prior GM-,WM-partial volume maps.
Follows the principles of the PVEc algorithm by I. Asllani (MRM, 2008).
Free for research use without guarantee. If used in a study or publication. Please, acknowledge the author.

----
### xASL_im_ResampleLinearFair.m

#### Function
```matlab
function [output_res]=xASL_im_ResampleLinearFair(im_input,newsize)
```

#### Description
Downsample (or upsample, works similarly) old_res image to low_res image, trilinear.

**NB:** new_res should fit exactly integer fold in old_res.

**NB:** All dimensions of new_res should have equal size.


----
### xASL_im_RestoreOrientation.m

#### Function
```matlab
function xASL_im_RestoreOrientation(PathNIfTI)
```

#### Description
This function reverts the NIfTI header orientation matrix to the original orientation from the scanner/dcm2nii conversion.

----
### xASL_im_rotate.m

#### Function
```matlab
function rotated = xASL_im_rotate(im, angle)
```

#### Description
Simple rotation of the first two dimension of a ND image by 0, 90, 180, 270 degrees.

----
### xASL_im_SkullStrip.m

#### Function
```matlab
function xASL_im_SkullStrip(InPath, PathMNIMask, x, OutPath)
```

#### Description
Creates skull-stripped T1w image based on MNI -> native space registration from segmentation.

----
### xASL_im_Smooth3D.m

#### Function
```matlab
function [imSmo,imGaussX,imGaussY,imGaussZ] = xASL_im_Smooth3D(sigma, imIn, PSFtype)
```

#### Description
...

----
### xASL_vis_TileImages.m

#### Function
```matlab
function [ImOut] = xASL_vis_TileImages(ImIn, nColumns)
```

#### Description
Merges selected slices (3D) into one single 2D picture.

Plots all slices in one figure with specified rows and columns, aiming for a square tile.

**PM:** can be extended to multiple slices.

----
### xASL_vis_TransformData2View.m

#### Function
```matlab
function FigureOut = xASL_vis_TransformData2View(ImagesIn, x)
```

#### Description
This function changes the dimensionality and reshapes the input images in such a way that they are nicely tiled in a mosaic for visualization purposes.
Reshaping a series of images with this function can be useful for visualization of SPM/voxel-based analyses.

----
### xASL_im_Upsample.m

#### Function
```matlab
function xASL_im_Upsample(PathOrig, PathDest, NewVoxelSize, LeaveEmpty, PaddingDim, Kernel)
```

#### Description
Upsamples an ASL image, without changing the orientation matrix, which can be used e.g. for PVEc in higher resolution but same space.

----
### xASL_vis_VisualizeROIs.m

#### Function
```matlab
function xASL_vis_VisualizeROIs(x, ROI_T1_list, ROI_FLAIR_list)
```

#### Description
Creates for each subject  a JPEG image containing the original T1w, WMH_SEGM and T1w after lesion-filling.

----
### xASL_vis_VisualQC_TopUp.m

#### Function
```matlab
function [MeanAI_PreTopUp_Perc, MeanAI_PostTopUp_Perc] = xASL_vis_VisualQC_TopUp(PathPopB0, PathPopUnwarped, x, iSubject, CheckDir)
```

#### Description
This function creates a Figure that showes the effect of TopUp with 6 images with axial slices: the NormPE, RevPE and their difference image in colorscale, and this before (upper row) & after (lower row) TopUp.

----
### xASL_im_ZeroEdges.m

#### Function
```matlab
function [IM] = xASL_im_ZeroEdges(IM, EdgeThicknessPerc)
```

#### Description
Resampling can sometimes give some strange errors near image edges. These should be NaNs, but sometimes can be zeros or ones, or even weird numbers. For resampling, NaNs should be set to 0 (this is done in another function) as they can influence the resampling (depending on the transformation matrix). To be sure that the edges are nicely fixed, this function sets a border at the image matrix edges to zero.

----
### xASL_init_ConvertM2JSON.m

#### Function
```matlab
function [PathJSON] = xASL_init_ConvertM2JSON(PathM, bOverwrite)
```

#### Description
This function converts and replaces the legacy data parameter m-format by a JSON file. A DataPar.m was the settings/parameter file, specific to a dataset to be processed by ExploreASL, now replaced to JSON by BIDS. Note that the deployed/compiled version of ExploreASL requires the JSON file, this function should not be compiled along. This function performs the following steps:

1. Run the m-file to load parameters
2. Escape characters that are illegal in JSON
3. Write the JSON

----
### xASL_init_DefaultEffectiveResolution.m

#### Function
```matlab
function [EffectiveResolution] = xASL_init_DefaultEffectiveResolution(PathASL, x)
```

#### Description
This ExploreASL module provides an educated guess on  the effective spatial resolution of ASL. This may depend on the combination of acquisition PSF, reconstruction filter, head motion. Note that the derived/processed images may have a lower effective resolution because of smoothing effects from interpolation. Note that this remains an educated-guess, the actual FWHM may still differ, especially for 3D GRASE sequences, where e.g. the choice of number of segments can affect the smoothness.
This function conducts the following steps:

1. Educated-guess FWHM
2. Attempt accounting for in-plane interpolation in reconstruction
3. Calculate and report effective spatial resolution

----
### xASL_init_DefineStudyData.m

#### Function
```matlab
function [x] = xASL_init_DefineStudyData(x)
```

#### Description
This initialization wrapper initializes the parameters for this pipeline run, i.e. subjects, sessions (runs), timepoints (visits), exclusions, sites, cohorts etc.

Note that ASL sessions are defined here as what BIDS calls "runs".

The "longitudinal_Registration functions here manage different TimePoints, which is what BIDS calls "visits".
With different structural scans, from the same participant. This is managed by subject name suffixes \_1 \_2 \_n, and can be used for comparing visits in the population module, or running SPM's longitudinal within-subject registration

Parallelization is allowed here by calling ExploreASL different times, where it divides the subjects/images for processing across the nWorkers, using iWorker as the reference for the part that the current ExploreASL call will process. This requires having a Matlab license that can be started multiple times on a server, or alternatively running the ExploreASL compilation, and doesn't require the Matlab parallel toolbox.

This function exists from the following parts:

1. Manage included & excluded subjects
2. Create dummy defaults (exclusion list, ASL sessions)
3. Create list of total baseline & follow-up subjects, before exclusions
4. Create TimePoint data-lists
5. Manage exclusions
6. Add sessions as statistical variable, if they exist
7. Parallelization: If running parallel, select cases for this worker
8. Add Longitudinal TimePoints (different T1 volumes) as covariate/set, after excluding
9. Load & add statistical variables
10. Make x.S.SetsOptions horizontal if vertical by transposing
11. Create single site dummy, if there were no sites specified
12. Check for time points sets
13. Create list of baseline & follow-up subjects (i.e. after exclusion)
14. Check what excluded from which TimePoints

----
### xASL_init_FileSystem.m

#### Function
```matlab
function [x] = xASL_init_FileSystem(x)
```

#### Description
This function initializes the file system used throughout ExploreASL, for processing a single dataset/scan.

It is repeated for each scan, and runs the following parts:

1. Create folders
2. Subject/session definitions
3. Add prefixes & suffixes
4. Add Subject-specific prefixes
5. Add sidecars

----
### xASL_init_InitializeMutex.m

#### Function
```matlab
function [x] = xASL_init_InitializeMutex(x, ModuleName)
```

#### Description
This function initializes the mutex/lock system of ExploreASL for a module. Mutex (for mutual exclusion) is a synchronization mechanism for enforcing limits of access to data (here a module for a single scan) to allow parallelization. It also allows stopping and continuing of ExploreASL. This function runs the following steps:

1. Lock folder management
2. Initialize mutex object

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
function [x] = xASL_init_LoadMetadata(x)
```

#### Description
This function loads all metadata used in the study, either statistical covariates (age, MMSE) or groups to compare between (site, sequence, cohort), or parameters to be used in quantification/image processing

These parameters should be provided in .mat files in the root analysis folder. Each .mat file should contain a single type of metadata, and the filename should equal the variable name. Metadata content should be a cell array with subjects as first column and metadata as last column. Sessions (runs) can be included as second column.

Metadata can be in any string or numerical format.

participants.tsv is now added per BIDS. It looks for metadata in participants.tsv first, before going through the mat-files.

This function iterates through the following steps for each variable:

1. Admin (what nOptions do we call ordinal, convert subject numeric to string, remove white spaces from data)
2. Get unique list of data options & check for missing data
3. Deal with data format (correct NaNs, deal with numeric vs strings)
4. Distinguish continous data (e.g. age) or ordinal data (groups to compare, e.g. cohort)
5. Check if data is complete for all subjects
6. Include complete data in x.S.SETS

----
### xASL_init_LongitudinalRegistration.m

#### Function
```matlab
function [SubjectNlist, TimePoint, IsSubject, SubjectID_FirstVolume] = xASL_init_LongitudinalRegistration(x)
```

#### Description
This function initializes the longitudinal registration for ExploreASL, which implements the SPM longitudinal registration.

This function recognizes and defines visits (i.e. multiple scans per subject) in the initialization of ExploreASL, based on the suffixes \_1 \_2 \_n in the subject names (identified as foldernames).

Specifically, this function is called in the registration modules LongReg and DARTEL, the first carrying out within-subject registration and the second between-subject registration, based on the first time point only.
For the first function, we specify here a list of visits/timepoints that should be registered longitudinally, for the second function we specify a list of first visits only, as the between-subject registration in ExploreASL is based on the first scan (as opposed to the average subject's scan).
This function runs the following steps:

1. Get TimePoint-list (list of visits)
2. Find subject IDs

----
### xASL_init_VisualizationSettings.m

#### Function
```matlab
function [x] = xASL_init_VisualizationSettings(x)
```

#### Description
This function defines several visualization settings are used throughout ExploreASL's pipeline and tools, assuming a [121 145 121] matrix with 1.5 mm isotropic resolution in MNI space

INCLUDING:

* Slices*    - defines which transversal slices to show by default
* TransCrop  - defines default cropping settings for all orientations
* ColorMaps  - Generate several colormaps
* VoxelSize  - default voxel-size in standard space/MNI as used by ExploreASL: 1.5 mm
* skull      - brainmask for skullstripping

----
### xASL_io_CreateNifti.m

#### Function
```matlab
function xASL_io_CreateNifti(pathNewNifti, imNew, resMat, nBits, bGZip)
```

#### Description
This function creates a new NIfTI file, using the SPM "nifti" functionality, with the parameters specified as input arguments. This function performs the following steps:

1. Initialize NIfTI
2. Choose datatype (bit resolution)
3. Create scale slopes
4. Create orientation matrix
5. Write the new NIfTI, image matrix & scale slopes
6. Zip and deal with zipping (.nii vs. .nii.gz)

----
### xASL_io_dcm2nii.m

#### Function
```matlab
function [niifiles, ScanNameOut, usedinput, msg] = xASL_io_dcm2nii(inpath, destdir, series_name, varargin)
```

#### Description
Convert DICOM NIfTI/BIDS format using the dcm2nii command line utility.

----
### xASL_io_DcmtkRead.m

#### Function
```matlab
function header = xASL_io_DcmtkRead(filepath, bPixel)
```

#### Description
Calls the MEX function that uses DCMTK library to read the DICOM header.
To change which parameters are read and their names - the MEX file needs to be edited. This function also corrects formating of certain parameters.

----
### xASL_io_MakeNifti4DICOM.m

#### Function
```matlab
function xASL_io_MakeNifti4DICOM(PathIn, x, DataType, OrientationPath, ResamplePath)
```

#### Description
This function converts a NIfTI file to one that is ready to convert to DICOM for PACS visualization purposes:

For scaling/visualization:

1. Remove peak signal
2. Remove valley signal
3. Remove NaNs
4. Rescale to 12 bit integers
5. Save NIfTI. We also zip the NIfTI as this NIfTI won't be opened by ExploreASL
6. Manage scale slope/datatype
7. Apply original orientation
8. Zip NIfTI

----
### xASL_io_PairwiseSubtraction.m

#### Function
```matlab
function xASL_io_PairwiseSubtraction(InputFile,outputPath,do_mask,switch_sign)
```

#### Description
Subtracts controls from labels and takes mean.
Creates new perfusion-weighted delta_M file, prefaced with 's'.
Converts into single precision floating point values (32 bit), removes scale slope.
Only runs if ['s' input_file_ASL] doesn't exist.
Remember to consider motion correction/ SPM realign (isotropically).
Alternative to this function is robust fit (Camille Maumet).

----
### xASL_io_ReadTheDicom.m

#### Function
```matlab
function [Info] = xASL_io_ReadTheDicom(bUseDCMTK, DicomPath)
```

#### Description
This function tries to read a DICOM and throws a warning if it fails to.

----
### xASL_io_SplitASL_M0.m

#### Function
```matlab
function xASL_io_SplitASL_M0(InPath,iM0)
```

#### Description
This function splits ASL4D & M0 if they were in the same sequence.
If dcm2niiX has already splitted the ASL4D NIfTI, this is reconstructed first.
If no M0 exists, or only ASL splitting is desired, leave iM0 empty ([])

Vendor product sequence examples:

* GE 3D spiral sometimes puts the M0 at the last volume of the series -> iM0 = [2];
* Philips 3D GRASE puts the M0 as control-label volume pair -> iM0 = [1 2];
* Siemens 3D GRASE puts the M0 as the first volume -> iM0 = 1;

----
### xASL_Iteration.m

#### Function
```matlab
function [bAborted, xOut] = xASL_Iteration(x, moduleName, dryRun, stopAfterErrors, SelectedSubjects)
```

#### Description
Parses the settings and runs the DatabaseLoop sub-function.

----
### xASL_num2str.m

#### Function
```matlab
function [DataOut] = xASL_num2str(DataIn, f)
```

#### Description
When the provided data is numeric, this function will convert the number to string/characters, rewriting NaN into 'n/a' (BIDS convention) but otherwise preserving the Matlab builtin functionality, also for the second argument "f".
If non-numeric data is provided, it is bypassed (avoiding any issues "num2str" will have with non-numeric data).
See builtin num2str for more details.

----
### xASL_qc_AsymmetryIndex.m

#### Function
```matlab
function [AI_perc] = xASL_qc_AsymmetryIndex(ImageIn)
```

#### Description
Extract voxel-wise asymmetry index for QC purposes.

----
### xASL_qc_CAT12_IQR.m

#### Function
```matlab
function [QA_Output] = xASL_qc_CAT12_IQR(InputImage, InputC1, InputC2, InputC3, bFLAIR)
```

#### Description
Prepare and run CAT12s QC parameters (also for other images).

----
### xASL_qc_CollectParameters.m

#### Function
```matlab
function x = xASL_qc_CollectParameters(x, iSubject, ScanType)
```

#### Description
This function collects QC parameters for a module.

----
### xASL_qc_CollectQC_ASL.m

#### Function
```matlab
function [x] = xASL_qc_CollectQC_ASL(x, iSubject)
```

#### Description
This functions collects QC parameters for the ASL module.
These are stored in x.Output.ASL:

* ID - SubjectName
* ASL_LR_flip_YesNo - Checks whether any image processing changed the left-right orientation by checking whether the determinant differs between nii.mat & nii.mat0 SPM realign (too much motion is suspicious)
* MotionMean_mm    - mean motion
* MotionExcl_Perc  - percentage of excluded outliers
* MotionMax_mm     - max motion
* MotionSD_mm      - SD motion

ASL quantification (strange average CBF, or strange GM-WM contrast). ASL acquisition parameters (should be fairly consistent over subjects/scans):

* TE - echo time
* TR - repetition time
* RescaleSlope - Philips
* Scaleslope - Philips
* Matrix X Y Z - matrix size
* Matrix Z - number of slices
* VoxelSize X Y - in plane resolution
* VoxelSize Z - slice thickness
* RigidBody2Anat_mm - Net Displacement Vector (RMS) from ASL to T1w image (mm) from registration

----
### xASL_qc_CollectQC_func.m

#### Function
```matlab
function [x] = xASL_qc_CollectQC_func(x, iSubject)
```

#### Description
This functions collects QC parameters for the func module. These are stored in x.Output.func:

* ID - SubjectName
* func_LR_flip_YesNo - Checks whether any image processing changed the left-right orientation by checking whether the determinant differs between nii.mat & nii.mat0 SPM realign (too much motion is suspicious)
* MotionMean_mm    - mean motion
* MotionExcl_Perc  - percentage of excluded outliers
* MotionMax_mm     - max motion
* MotionSD_mm      - SD motion

Func quantification (strange average CBF, or strange GM-WM contrast):

* CBF_GM_Median_mL100gmin - median GM CBF
* CBF_WM_Median_mL100gmin - median WM CBF
* SpatialCoV_GM_Perc      - GM spatial CoV
* SpatialCoV_WM_Perc      - WM spatial CoV
* CBF_GM_WM_Ratio         - GM-WM CBF ratio

Func acquisition parameters (should be fairly consistent over subjects/scans):

* TE - echo time
* TR - repetition time
* RescaleSlope - Philips
* Scaleslope - Philips
* Matrix X Y Z - matrix size
* Matrix Z - number of slices
* VoxelSize X Y - in plane resolution
* VoxelSize Z - slice thickness
* RigidBody2Anat_mm - Net Displacement Vector (RMS) from func to T1w image (mm) from registration

----
### xASL_qc_CollectQC_Structural.m

#### Function
```matlab
function [x] = xASL_qc_CollectQC_Structural(x, iSubject)
```

#### Description
This functions collects QC parameters for the structural module. These are stored in x.Output.Structural:

* ID - SubjectName
* T1w_LR_flip_YesNo - Checks whether any image processing changed the left-right orientation by checking whether the determinant differs between nii.mat & nii.mat0

LST output:

* WMH_vol_mL        - WMH volume
* WMH_n             - WMH number of lesions

CAT12 output:

* T1w_IQR_Perc      - CAT12 IQR quality parameter for T1w

Volumetric: GM_vol_mL, WM_vol_mL, CSF_vol_mL, ICV_vol_mL, GM_ICV_Ratio

----
### xASL_qc_CollectSoftwareVersions.m

#### Function
```matlab
function [x] = xASL_qc_CollectSoftwareVersions(x)
```

#### Description
This functions collects software versions for matlab, SPM, CAT, LST & ExploreASL. These are stored in x.Output.Software.

----
### xASL_qc_CompareTemplate.m

#### Function
```matlab
function [QC] = xASL_qc_CompareTemplate(x, ScanTypePrefix, iSubjectSession)
```

#### Description
This function computes several advanced template-based QC parameters:

* RMSE_Perc        - Root Mean Square Error between image and template (%)
* nRMSE_Perc       - Same but then normalized
* AI_Perc          - Asymmetry Index between image and template (%)
* Mean_SSIM_Perc   - mean structural similarity index -> xASL_stat_MeanSSIM.m
* PeakSNR_Ratio    - peak signal-to-noise ratio -> xASL_stat_PSNR.m

----
### xASL_qc_ComputeFoVCoverage.m

#### Function
```matlab
function [CoveragePerc] = xASL_qc_ComputeFoVCoverage(InputPath, x)
```

#### Description
This function computes the intersection/overlap between brainmask on field-of-view (FoV) of low resolution image (native space) & the same brainmask with expanded FoV. It uses the pGM+pWM+pCSF as brainmask.
This assumes that the structural reference image has full brain coverage, and was properly segmented into GM, WM and CSF. Also, we assume that the InputPath contains a single 3D volume.

----
### xASL_qc_ComputeNiftiOrientation.m

#### Function
```matlab
function [Struct] = xASL_qc_ComputeNiftiOrientation(x, PathNIfTI, Struct)
```

#### Description
...

----
### xASL_qc_CreatePDF.m

#### Function
```matlab
function xASL_qc_CreatePDF(x, DoSubject)
```

#### Description
This function iterates over all values in x.Output and all images in x.Output_im, and prints them in a PDF file.
x.Output & x.Output_im should contain the QC/result output of all ExploreASL pipeline steps.

Further code explanation:

Below, using the Matlab & SPM Figure tools we create an image, which is then printed to a PDF file.

* fg = the main Figure handle
* ax = "axes" handles, these are objects containing either 1) text or 2) images, with fg as "parent" (1) & (2) images have ax as "parent".

Positions are calculated in such a way that 4 categories can be printed, which will be the first 4 fields found in x.Output then allowing 8 single slice images, and 15 text lines (name & value columns).

----
### xASL_qc_FA_Outliers.m

#### Function
```matlab
function [FA_Outliers_mL] = xASL_qc_FA_Outliers(InputFA)
```

#### Description
Extract the number of FA outliers, i.e. values of FA above 1 or below 0, from a FA image.

----
### xASL_qc_ObtainQCCategoriesFromJPG.m

#### Function
```matlab
function xASL_qc_ObtainQCCategoriesFromJPG(x)
```

#### Description
This function obtains QC categories as covariant/set, based on the JPGs in //Population/ASLCheck. These are initially sorted by spatial CoV, and should be visually checked & put in the correct folder.

----
### xASL_qc_PCPStructural.m

#### Function
```matlab
function [anatQA] = xASL_qc_PCPStructural(PathT1, Pathc1T1, Pathc2T1, x, PopPathT1)
```

#### Description
This function computes several anatomical QC parameters as proposed in SPM Univariate Plus:

* WM_ref_vol_mL    - volume of the WM reference region (mL)
* WMref_vol_Perc   - same but as percentage of total WM volume
* SNR_GM           - GM signal-to-Noise Ratio (SNR), ie the mean intensity within GM divided by SD of WM reference region. Higher = better.
* CNR_GM_WM        - GM-WM Contrast-to-Noise Ratio (CNR), i.e. the mean of GM - mean of WM divided by the SD of the WM reference region. Higher = better.
* FBER_WMref_Ratio - Foreground to Background Energy Ratio (FBER), i.e. the variance of voxels within the brain (in pGM+pWM mask) divided by the variance of voxels in the WM reference region. Higher = better.
* EFC_bits         - Shannon Entropy Focus Criterion (EFC), i.e. the entropy of voxel intensities proportional to the maximum possibly entropy for a similarly sized image. Indicates ghosting and head motion-induced blurring. Lower = better.
* Mean_AI_Perc     - mean relative voxel-wise absolute Asymmetry Index (AI) within the brain (pGM+pWM mask) (%)
* SD               - same but SD (%)

----
### xASL_qc_PrintOrientation.m

#### Function
```matlab
function xASL_qc_PrintOrientation(DIR, reg_exp_INPUT,OUTPUT_DIR,Name)
```

#### Description
Check orientation of niftis, useful to detect accidental left-right flips (all other flips will be visible). Translations, rotations or shears are not to be worried about, only negative zooms. This can be detected by negative determinants.
So orientation parameters and determinants should be similar across all scans from single scanner/coil, and registration should not give negative determinant.

----
### xASL_qc_TanimotoCoeff.m

#### Function
```matlab
function TC = xASL_qc_TanimotoCoeff(Image1, Image2, imMask, type, bClip, bSmooth)
```

#### Description
Compares images Image1 and Image2 within the mask imMask. TYPE specifies the input data type.

----
### xASL_qc_temporalSNR.m

#### Function
```matlab
function tSNR = xASL_qc_temporalSNR(pathIm4D,pathImTissueProb)
```

#### Description
This function computes several temporal SNR QC parameters as proposed in SPM Univariate Plus:

* tSNR.tSNR_GM_Ratio      : mean GM signal / std GM over time 
* tSNR.tSNR.tSNR_WM_Ratio : mean WM signal / std WM over time
* tSNR.tSNR.tSNR_CSF_Ratio: mean CSF signal / std CSF over time
* tSNR.tSNR_WMref_Ratio   : mean signal/std over time in eroded deep WM
* tSNR.tSNR_GMWM_Ratio    : mean (GM+WM) signal / sqrt(std(GM+WM)^2+std(WMref)^2)
* tSNR.tSNR_GMWM_WMref_Ratio: mean (GM+WM) signal / std WMref over time
* tSNR.tSNR_Physio2Thermal_Ratio: sqrt((tSNR(GM+WM)/tSNR_GMWM_WMref_Ratio))^2-1)
* tSNR.tSNR_Slope_Corr:

Differences to the SPM U+ suggestion: 
* Eroded WM is used for estimating background noise
* Brainmask is determined in the same way as the structural anatQC,
* CSF is determined from the pGM&pWM maps;

----
### xASL_qc_WADQCDC.m

#### Function
```matlab
function xASL_qc_WADQCDC(x, iSubject, ScanType)
```

#### Description
This QC function runs WAD-QC specific Python script to zip QC information & incorporate this into a DICOM field for analysis on the WAD-QC server, by the following:

1. Define QCDC script: this is the Python script written by Gaspare, edited by Joost Kuijer & copied to the EPAD CustomScripts folder of ExploreASL
2. Python installation location is checked, with several known locations, for several servers. If we cannot find it, the QCDC is not ran
3. Previous QCDC results are cleaned. QCDC stores all its results in a separate folder (Something like 2 layers up from the current folder, here referred to as QCDCDir = [x.D.ROOT 'qcdc_output']) from these result files, only the filled DICOM file is interesting, all the rest are copies of the QC results that we embedded into the DICOM
4. Run QCDC (if Python installation detected). The following files need to be set as executable:
    * ('QCDC', 'src', 'qc_data_collector.py')
    * ('QCDC', 'src', 'bash', 'create_dcm_only_wadqc.sh')
    * ('QCDC', 'src', 'bash', 'sendwadqc.sh')
5. Clean up new QCDC results (as above) & move the filled. DICOM to ['qcdc_' DummyFile] within the current ScanType folder.
6. Sending the DICOM to the WAD-QC server using storescu

----
### xASL_qc_WADQC_GenerateDescriptor.m

#### Function
```matlab
function xASL_qc_WADQC_GenerateDescriptor(x, iSubject, ScanTypeIs)
```

#### Description
This QC function generates a JSON descriptor for Gaspare' QCDC script, by the following steps:

1. include information about where to find the dummy DICOM (i.e. placeholder DICOM)
2. For ExploreASL' QC fields (as passed through in x.Output), here we note all these QC fields for each ScanType, as the x.Output should have been collected equally in the QC file 'QC_collection_SubjectName.json' by function xASL_qc_CollectParameters
3. Subfunction xASL_qc_WADQC_images - Includes visual standard space QC images, by searching them on prescribed paths within the Population folder (where currently all derivatives reside)
4. Insert the PDF report; this PDF report is subject-specific, not scan-specific. For completeness it is added to each QCDC descriptor
5. Add WAD-QC server details (i.e. IP address etc)
6. Save the Descriptor JSON file.

----
### xASL_quant_AgeSex2Hct.m

#### Function
```matlab
function Hct = xASL_quant_AgeSex2Hct(age, gender)
```

#### Description
Function to return estimated Hct, based on age and gender. Enter age in years, and for gender: 0=female, 1=male pass NaN is either not known.


----
### xASL_quant_FEAST.m

#### Function
```matlab
function xASL_quant_FEAST(x)
```

#### Description
Computation FEAST-based transit times (uses images that were not vasculary treated) if there is a crushed & non-crushed scan, then transit times will be computed by division of these scans, provided sessions are exactly named as defined below.

----
### xASL_quant_GetControlLabelOrder.m

#### Function
```matlab
function [ControlIm, LabelIm, OrderContLabl] = xASL_quant_GetControlLabelOrder(ASLTimeSeries)
```

#### Description
This function automatically checks (and corrects if required) the control and label order of ASL timeseries based on the larger signal in control volumes.
It supposes that data is acquired in pairs.

----
### xASL_quant_Hct2BloodT1.m

#### Function
```matlab
function BloodT1 = xASL_quant_Hct2BloodT1(Hematocrit, Y, B0, bVerbose)
```

#### Description
This function converts hematocrit to blood T1, according to calculations defined by Patrick Hales. With courtesy and thanks!
Note that we assume a venous O2 saturation of 68% (Yv=0.68)

This function performs the following steps:

1. Check fraction vs percentage hematocrit & Y, should be between 0 and 1
2. Specify defaults (Hb, Fe)
3. Perform calculation
4. Convert s to ms
5. Print what we did

----
### xASL_quant_M0.m

#### Function
```matlab
function [M0IM] = xASL_quant_M0(M0IM, x)
```

#### Description
This function quantifies the M0, except for the difference in voxel size between the M0 and ASL source data (which is scaled in xASL_wrp_ProcessM0.m).

----
### xASL_quant_SinglePLD.m

#### Function
```matlab
function [ScaleImage, CBF] = xASL_quant_SinglePLD(PWI, M0_im, SliceGradient, x)
```

#### Description
This script performs a multi-step quantification, by initializing a ScaleImage that travels through this script & gets changed by the following quantification factors:

1. PLD scalefactor (gradient if 2D multi-slice) (if x.ApplyQuantification(3))
2. Label decay scale factor for single (blood T1) - or dual-compartment (blood+tissue T1) model, CASL or PASL. Single-compartment model: Alsop MRM 2014.
 Dual-compartment model: Wang MRM 2002: Gevers JMRI 2012 (if x.ApplyQuantification(3))
3. Scaling to physiological units \[ml/gr/ms =>ml/100gr/min =>(60,000 ms=>min)(1 gr=>100gr)\] (if x.ApplyQuantification(3))
4. Vendor-specific scalefactor (if x.ApplyQuantification(1) -> future move to dcm2niiX stage)
5. Divide PWI/M0 (if x.ApplyQuantification(5))
6. Print parameters used

Note that the output always goes to the CBF image (in the future this could go to different stages, e.g. dcm2niiX or PWI stage)

----
### xASL_spm_affine.m

#### Function
```matlab
function xASL_spm_affine(srcPath, refPath, fwhmSrc, fwhmRef, otherList, bDCT)
```

#### Description
This SPM wrapper runs SPM's old normalize-estimate function, which calculates the affine transformation (i.e. linear + zooming and shearing) that is required to align the source image with the reference image. By default it does not estimate the low-degree Discrete Cosine Transform (DCT) to have a simple affine transformation  but this can be enabled in this wrapper. Also note that this affine transformation uses a correlation cost function, hence it requires the source and reference images to have similar contrasts and resolution - or provide the resolution estimates so the smoothing can be done.
This function does not change the orientation header by default, but it allows to change those of others through the otherList. If bDCT is used or no otherList given, it stores its estimated affine transformation as orientation difference matrix in a file with the same path but \_sn.mat extension. For the provided smoothing FWHM, note that smoothnesses combine with Pythagoras' rule (i.e. square summing).

----
### xASL_spm_BiasfieldCorrection.m

#### Function
```matlab
function xASL_spm_BiasfieldCorrection(PathIn, SPMdir, Quality, PathMask, PathOut)
```

#### Description
This function is a wrapper around the SPM "old segment" function, for biasfield removal. It is tested for M0 and mean control images. It conducts the following steps:

1. Create implicit mask
2. Define SPM 'old segmentation' settings
3. Run SPM 'old segmentation'
4. Delete temporary files
5. Rename temporary SPM file into output file

----
### xASL_spm_coreg.m

#### Function
```matlab
function xASL_spm_coreg(refPath, srcPath, OtherList, x, sep, FastReg)
```

#### Description
This SPM wrapper runs SPMs coregister-estimate function, which calculates the 6 parameter rigid-body transformation (a.k.a. linear) that is required to align the source image with the reference image. This 6 parameter transformation (i.e. 3 XYZ translations and 3 rotations) is applied to the orientation header of the source NIfTI image, and also to the images provided in OtherList (optional).

Note that this SPM registration function uses a normalized mutual information (NMI) by default, enabling registration between two images with different contrast. Note that this algorithm will use the first volume of a multi-volume NIfTI

----
### xASL_spm_deface.m

#### Function
```matlab
function xASL_spm_deface(PathIn, bReplace)
```

#### Description
This function removes the face from an anatomical NIfTI image, e.g. T1w or FLAIR, for disidentification/privacy purposes. When this script is run after the ExploreASL structural module, it does a pretty good job even for 2D images.
However, note that this can always fail, strip part of the brain, or change the output of pipelines. So best not to compare results from defaced and non-defaced images.
Also, note that defacing makes it difficult to ensure that the FLAIR and T1w are from the same subject.

----
### xASL_spm_deformations.m

#### Function
```matlab
function xASL_spm_deformations(x, PathIn, PathOut, Interpolation, InverseSpace, AffineTrans, DeformationPath, VoxelSize)
```

#### Description
This ExploreASL wrapper manages the SPM deformation tool. It takes multiple (ExploreASL pipeline) transformations and combines/concatenates them into a single transformation prior to applying it to the input images. 
This allows to apply multiple transformations with a single interpolation, avoiding propagation of undesired interpolation effects. Mainly used to get native space images into standard space, or vice versa.
Best to combine as many files as possible within this function, since the deformation calculation (which is the most computation intensive part) needs to be performed once for multi-file resampling.

----
### xASL_stat_AtlasForStats.m

#### Function
```matlab
function [x] = xASL_stat_AtlasForStats(x)
```

#### Description
This function loads atlases, checks them, and their ROI names, for later use as ROI definition in xASL_stat_GetROIstatistics.
Note that the atlases should be integer values, or different 4rd dimensions (i.e. multiple images), that are mutually exclusive. This function takes the following steps:

1. Load atlas ROI names. There should be a TSV sidecar to the atlas NIfTI file, as explained above.
2. deal with memory mapping
3. Resample atlas 50 1.5 mm^3 MNI
4. Converted atlas with integers to 4D binary image
5. Convert/compress masks into Columns
6. Print atlas overview image

----
### xASL_stat_ComputeDifferCoV.m

#### Function
```matlab
function diffCoV = xASL_stat_ComputeDifferCoV(imCBF,imMask,bPVC,imGM,imWM,b3D)
```

#### Description
It calculates the spatial DiffCoV value on finite part of imCBF. Optionally a mask IMMASK is provide, and PVC is done for bPVC==2 using imGM and imWM masks and constructing pseudoCoV from pseudoCBF image. For bPVC\~=2, imGM and imWM are ignored. It is calculated in 2D or assuming also 3D edges based on B3D. Calculate derivate spatial CoV, by summing up differences in CBF between neighbors. The derivative uses Sobels filter.

----
### xASL_stat_ComputeMean.m

#### Function
```matlab
function [CBF_GM, CBF_WM] = xASL_stat_ComputeMean(imCBF, imMask, nMinSize, bPVC, imGM, imWM)
```

#### Description
It behaves in a similar way as VAR.

----
### xASL_stat_ComputeSpatialCoV.m

#### Function
```matlab
function sCov = xASL_stat_ComputeSpatialCoV(imCBF,imMask,nMinSize,bPVC,imGM,imWM)
```

#### Description
It calculates the spatial CoV value on finite part of imCBF. Optionally a mask IMMASK is provide, ROIs of size < NMINSIZE are ignored, and PVC is done for bPVC==2 using imGM and imWM masks and constructing pseudoCoV from pseudoCBF image. For bPVC\~=2, imGM and imWM are ignored.

----
### xASL_stat_EqualVariancesTest.m

#### Function
```matlab
function [resTest, P] = xASL_stat_EqualVariancesTest(X, alpha, type)
```

#### Description
Brown-Forsythe or Levene's test for equality of variances. The response variable is transformed (yij = abs(xij - median(xj)) for Brown-Forsythe and yij = abs(xij - mean(xj)) for Levene's test). And then runs a one-way ANOVA F-test to check if the variances are equal.

----
### xASL_stat_fcdf.m

#### Function
```matlab
function p = xASL_stat_fcdf(F,M,N)
```

#### Description
Calculates the cumulative distribution function of the F-distribution for degrees of freedom M, N at value F.

----
### xASL_stat_GetROIstatistics.m

#### Function
```matlab
function [x] = xASL_stat_GetROIstatistics(x)
```

#### Description
This function computes mean and spatial CoV for each ROI, in a \[1.5 1.5 1.5\] mm MNI space, with several ASL-specific adaptions:

1. Skip ROI masks that are smaller than 1 mL as this would be too noisy for ASL (skipped when x.S.IsASL==false)
2. Expand each ROI mask such that it has sufficient WM content (skipped when IsASL==false)
3. Create for each ROI mask a left, right and bilateral copy
4. Iterate over all subjects:
    * Load partial volume maps
    * Correct for WMH SEGM -> IS THIS STILL REQUIRED???
    * Load data
    * Show ROIs projected on ASL image
    * Actual data computations
        * Partial Volume Correction (PVC) options:
        * PVC==0 -> perform masking only, no regression
        * PVC==1 -> single compartment regression, for pGM
        * PVC==2 -> dual compartment regression for pGM & pWM (i.e. normal PVC)
        * Here we mask out susceptibility artifacts (including outside FoV) for all ASL computations, and also mask out vascular artifacts for computing the mean.
    * While other artifacts/FoV can be masked out on population level (i.e. >0.95 subjects need to have a valid signal in a certain voxel), vascular artifacts differ too much in their location between subjects, so we mask this on subject-level.

Note that the words "mask" and "ROI" are used interchangeably throughout this function, where they can have a different or the same meaning

**PM:** WE COULD CHANGE THIS, INTO MASK BEING USED TO EXCLUDE VOXELS AND ROI FOR INCLUDING VOXELS

----
### xASL_stat_MadNan.m

#### Function
```matlab
function y = xASL_stat_MadNan(x,flag,dim)
```

#### Description
Calculates a Median/Mean Absolute deviation, but ignoring NaNs in the calculation.

xASL_stat_MadNan(X) or xASL_stat_MadNan(X,0) computes xASL_stat_MeanNan(ABS(X-xASL_stat_MeanNan(X))

xASL_stat_MadNan(X,1) computes xASL_stat_MedianNan(ABS(X-xASL_st_MedianNan(X)).

----
### xASL_stat_MeanSSIM.m

#### Function
```matlab
function mssim=xASL_stat_MeanSSIM(imRef,imSrc,dynRange)
```

#### Description
Calculates the similarity index according to Want et al. 

----
### xASL_stat_MultipleLinReg.m

#### Function
```matlab
function [b,CI,pval,stats] = xASL_stat_MultipleLinReg(X,Y,bIntercept)
```

#### Description
Performs a multiple linear regression Y=b \* X+a and provides the intercept and regression coefficients beta including their significance and confidence intervals. It calculates additionally the goodness of the fit.

----
### xASL_stat_PrintStats.m

#### Function
```matlab
function [x] = xASL_stat_PrintStats(x)
```

#### Description
This function prints an overview of statistics from data that were acquired per ROI, in a TSV file. It starts by printing covariates (called "Sets"). Rows will be subjects/sessions, columns will be the sets and

ROI-statistics:

1. First remove previous TSV-file, if already existed printing to a TSV file can be tricky if it is opened by Excel. Make sure to close previous versions first, otherwise this part will crash.
2. Print overview of sets to TSV as explained above. Uses subfunction xASL\_stat\_CreateLegend to put legends. Aim is to create a single TSV file that has a proper overview of the data, & is self-explanatory to those reading/using it.
3. Define number of ASL sessions, force to 1 in case of TT or volume metrics
4. Print the overview

----
### xASL_stat_PSNR.m

#### Function
```matlab
function PSNR=xASL_stat_PSNR(imRef,imSrc)
```

#### Description
Calculates the PSNR, needs two input arguments - 3D images of the same size. 
Uses 95% percentile instead of MAX.

----
### xASL_stat_QuantileNan.m

#### Function
```matlab
function y = xASL_stat_QuantileNan(x,quant,dim)
```

#### Description
Calculates a quantile, but ignoring NaNs in the calculation.

----
### xASL_stat_RobustMean.m

#### Function
```matlab
function [NotOutliers, iOutliers] = xASL_stat_RobustMean(IM, ParameterFunction)
```

#### Description
This function detects outlier images, that can be used to create a robust average, e.g. for template or biasfield creation. This is based either on the sum-of-squares with the mean image (SoS), or on the average relative asymmetry index (AI). Images that are median+/-3 mad off are defined as outliers. MAD = median/mean absolute difference.

----
### xASL_stat_ShapiroWilk.m

#### Function
```matlab
function [H, P, W] = xASL_stat_ShapiroWilk(x, alpha)
```

#### Description
Performs the statistical test of normality - null hypothesis is that the sample is from normal distribution with unspecified mean and variance. Based on the sample kurtosis it performs either Shapiro-Wilk (for platykurtic) or Shapiro-Francia test (for leptokurtic).

----
### xASL_stat_StdNan.m

#### Function
```matlab
function y = xASL_stat_StdNan(varargin)
```

#### Description
It behaves in a similar way as VAR - it directly passes all arguments to xASL_stat_VarNan.

----
### xASL_stat_SumNan.m

#### Function
```matlab
function y = xASL_stat_SumNan(x,dim)
```

#### Description
It uses the function SUM, but it sets all the NaNs to zero before calling it.

----
### xASL_stat_tcdf.m

#### Function
```matlab
function F = xASL_stat_tcdf(T,nu)
```

#### Description
Calculates the cumulative distribution function of the Student's t-distribution for degrees of freedom nu at value T.

----
### xASL_stat_ticdf.m

#### Function
```matlab
function T = xASL_stat_ticdf(P,nu)
```

#### Description
Calculates the inverse of cumulative distribution function of the Student's t-distribution for degrees of freedom nu at value P.

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


