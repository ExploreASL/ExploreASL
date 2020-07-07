----
### xASL\_Iteration.m

#### Format

```matlab
[bAborted, xOut] = xASL_Iteration(x, moduleName[, dryRun, stopAfterErrors])
```

#### Description
Parses the settings and runs the DatabaseLoop sub-function.

##Administration

----
### xASL\_adm\_CatchNumbersFromString.m

#### Format

```matlab
[OutputNumber] = xASL_adm_CatchNumbersFromString(InputString)
```

#### Description
...


----
### xASL\_adm\_CheckFileCount.m

#### Format

```matlab
[result, files] = xASL_adm_CheckFileCount(path, expr[, mincount, failifmissing])
[result]        = xASL_adm_CheckFileCount(...)
```

#### Description
Checks the given **PATH** for files corresponding to the **SPM\_SELECT** regular expression **EXPR**.
Returns if the number of files is equal to or higher than **MINCOUNT**. If **FAILIFMISSING** is true
and not enough files, then throw and error. If everything goes ok and second output argument is specified
then return also the list of files.

----
### xASL\_adm\_CheckPermissions.m

#### Format

```matlab
[FilesList, FilesExeList, FoldersList] = xASL_adm_CheckPermissions(InputPath[, FilesExecutable])
```

#### Description
This function does a recursive search through the root
folder & makes a list of the attributes of all files and folders.
It tries to reset the attributes to what we desire, which is by default:

- 664 for files (meaning only reading & writing for users & group, & read-only for others)
- 775 for folders (meaning reading, writing & opening for current user & current group, & for others only reading & opening)

For executable files we also want 775.
Note that the permission to 'execute a folder' means opening them.

- DataOK checks data permissions.
- ExeOK checks executable permissions.
- DataOK also includes executable permissions for folders.
- This runs recursively (but currently skips the contents of the root-folder) .


----
### xASL\_adm\_CheckSPM.m

#### Format

```matlab
[spm_path, spm_version] = xASL_adm_CheckSPM([modality, proposed_spm_path, check_mode])
[spm_path]              = xASL_adm_CheckSPM(...)
xASL_adm_CheckSPM(...)
```

#### Description
Checks if the spm function exists and if the reported version matches our development
version (**SPM8** or **SPM12**). If the spm toolbox is not available yet, it will try the
**PROPOSED\_SPM\_PATH** (if specified) or the user selected directory and add it to **PATH**.
The function will fail if **SPM** cannot be found or if detecting an unsupported version.

----
### xASL\_adm\_CleanUpBeforeRerun.m

#### Format

```matlab
xASL_adm_CleanupBeforeCompleteRerun(AnalysisDir, iModule, bRemoveWMH, bAllSubjects, SubjectID)
```

#### Description
This function (partly) reverts previous ExploreASL runs,
deleting derivatives, while keeping raw data intact.
if bAllSubjects==true, then all subjects and all module
derivatives will be removed. This function performs the
following steps:

1. If a Population folder doesn't exist yet but dartel does, rename it
2. Remove whole-study data files in AnalysisDir if bAllSubjects
3. Remove lock files/folders for reprocessing
4. Restore backupped \_ORI (original) files
5. Delete native space CAT12 temporary folders (always, independent of iModule)
6. Remove native space files for iModule
7. Remove standard space files for iModule
8. Remove population module files
9. Remove or clean up stored x-struct & QC file -> THIS HAS NO SESSION SUPPORT YET

NB: still need to add xASL\_module\_func & xASL\_module\_dwi for EPAD


----
### xASL\_adm\_CompareDataSets.m

#### Format

```matlab
[RMS] = xASL_adm_CompareDataSets(RefAnalysisRoot,SourceAnalysisRoot,x,type,mutexState)
```

#### Description
Compare data sets is used to ...

- type 0: Only save
- type 1: Save and evaluate
- type 2: Only evaluate



----
### xASL\_adm\_CompareLists.m

#### Format

```matlab
[NewList] = xASL_adm_CompareLists(list1, list2)
```

#### Description
...



----
### xASL\_adm\_ConvertDate2Nr.m

#### Format

```matlab
[Nr DayInYear] = xASL_adm_ConvertDate2Nr(TempDate)
```

#### Description
Converts date to number input mmdd -> output mm (with days in fractions/floating point).
Inverse from ConvertNrDate.



----
### xASL\_adm\_ConvertNr2Time.m

#### Format

```matlab
Time = xASL_adm_ConvertNr2Time(Nr)
```

#### Description
Converts number to time input hh (with minutes in fractions/floating point) -> output hhmm.
Inverse from xASL\_adm\_ConvertTime2Nr.



----
### xASL\_adm\_ConvertSubjSess2Subj\_Sess.m

#### Format

```matlab
[iSubj iSess] = xASL_adm_ConvertSubjSess2Subj_Sess(nSessions, iSubjSess)
```

#### Description
Converts combined SubjectSession index to subject & session
indices. Useful for data lists in ExploreASL.



----
### xASL\_adm\_ConvertTime2Nr.m

#### Format

```matlab
Nr = xASL_adm_ConvertTime2Nr(Time)
```

#### Description
Converts time to number input hhmm -> output hh (with
minutes in fractions/floating point).
Inverse from xASL\_adm\_ConvertNr2Time.



----
### xASL\_adm\_CopyMoveFileList.m

#### Format

```matlab
[List] = xASL_adm_CopyMoveFileList(OriDir, DstDir, StrRegExp, bMove[, bDir, bRecursive, bOverwrite, bVerbose])
```

#### Description
Moves a file to a file, a file to a directory, or a directory to a directory.
It keeps the initial extensions, no unzipping or zipping after the move.
But it makes sure that only one of **.nii** and **.nii.gz** exists in the destination directory.
Useful to split a large database.


----
### xASL\_adm\_CorrectName.m

#### Format

```matlab
strOut = xASL_adm_CorrectName(strIn[, bOption, strExclude])
```

#### Description
Finds and replaces all non-word characters either by empty space or by an underscore.
Optionally leaves in few selected special characters. Note that if '\_' is excluded from
replacement, but option 2 is on, then underscores are replaced anyway.


----
### xASL\_adm\_CreateCSVfile.m

#### Format

```matlab
xASL_adm_CreateCSVfile(CSVfilename,CSVdata)
```

#### Description
Creates a CSV file that can be opened with excel from
your data.



----
### xASL\_adm\_CreateFileReport.m

#### Format

```matlab
x = xASL_adm_CreateFileReport(x, bHasFLAIR, bHasMoCo, bHasM0, bHasLongitudinal)
```

#### Description
Prints a summary of created files or the individual modules
(i.e. Structural, Longiutudinal & ASL modules). Provides a quick check to
see what has been skipped, an whether all files are present.

This script iterates across:
Native space 1) subject and 2) session files,
Resampled 3) subject and 4) session files,
5) Lock files and 6) QC Figure files.

For all we perform a:

- A) Count of the files present, summarized in FileReportSummary.csv
- B) List of the missing files in "Missing\*.csv" files

PM: Simplify/optimize this code, to make filename variable changing,
search within subject-directories, etc. Combine the parts searching for
missing & summarizing count.



----
### xASL\_adm\_DefineASLResolution.m

#### Format

```matlab
x = xASL_adm_DefineASLResolution(x)
```

#### Description
...



----
### xASL\_adm\_DefineASLSequence.m

#### Format

```matlab
[x] = xASL_adm_DefineASLSequence(x)
```

#### Description
This ExploreASL function tries to check what ASL sequence is
being processed, if this was not already defined in x.Sequence.
It does so by checking known combinations of readout dimensionality
(x.readout\_dim) and vendor, knowing the product sequences of the vendors.


----
### xASL\_adm\_DeleteFilePair.m

#### Format

```matlab
filepaths = xASL_adm_DeleteFilePair(path, ext1[, ext2 [, ext3 ...]])
xASL_adm_DeleteFilePair(path, ext1[, ext2 [, ext3 ...]])
```

#### Description
Delete the file given in PATH, and also deletes files with the same name, but with extension
given in EXT1, and potentially also EXT2, EXT3...

----
### xASL\_adm\_Dicom2Parms.m

#### Format

```matlab
[parms pathDcmDictOut] = xASL_adm_Dicom2Parms(inp[, parmsfile, dcmExtFilter, bUseDCMTK, pathDcmDictIn])
```

#### Description
The function goes through the **INP** files, reads the **DICOM** or **PAR/REC** files and parses their headers.
It extracts the **DICOM** parameters important for ASL, makes sure they are in the correct format, if missing then
replaces with default value, it also checks if the parameters are consistent across **DICOM** files for a single sequence.


----
### xASL\_adm\_FindByRegExp.m

#### Format

```matlab
xasl_adm_FindByRegExp(root, dirSpecs[, varargin])
```

#### Description
Recursively find files in the root directory according to the dirSpecs.



----
### xASL\_adm\_FindStrIndex.m

#### Format

```matlab
INDEX = xASL_adm_FindStrIndex(ARRAY, STRING)
```

#### Description
Similar to find, but then for a cell array filled with strings.
Only takes 4 dimensions.



----
### xASL\_adm\_GetFsList.m

#### Format

```matlab
RES = xASL_adm_GetFsList([strDirectory, strRegEx, bGetDirNames, bExcludeHidden, bIgnoreCase, nRequired])
```

#### Description
List files or directories from a given path. And optionally uses regular expressions to filter the result
with options to exclude hidden files, ignore case, and set a minimal requirement on the number of results.
Sorts the results at the end.

----
### xASL\_adm\_GetNumFromStr.m

#### Format

```matlab
num = xASL_adm_GetNumFromStr(str)
```

#### Description
Obtains single number from string.
**CAVE** there should only be one number!



----
### xASL\_adm\_GetPhilipsScaling.m

#### Format

```matlab
scaleFactor = xASL_adm_GetPhilipsScaling(pathParmsMat,pathNifti)
```

#### Description
This script provides the correct scaling factors for a NIfTI file. It checks the header of the NIfTI
that normally has the same scaling as RescaleSlope in DICOM, it checks if dcm2nii (by the info in JSON)
has already converted the scale slopes to floating point. And if not, the derive the correct
scaling factor to be applied.


----
### xASL\_adm\_GetUserName.m

#### Format

```matlab
UserName = xASL_adm_GetUserName()
```

#### Description
...



----
### xASL\_adm\_Hex2Num.m

#### Format

```matlab
outNum = xASL_adm_hex2num(inStr)
```

#### Description
Takes a hexadecimal string and converts it to number. Works
also when the string contains escape characters, and for single-floats and
for a little and big endian. If containing 8 and less
characters than treat as float, if more than as double.


----
### xASL\_adm\_LesionResliceList.m

#### Format

```matlab
[INname, OUTname] = xASL_wrp_LesionResliceList(x,bLesion_T1,bLesion_FLAIR,bROI_T1,bROI_FLAIR)
```

#### Description
Creates list of structural image paths to reslice.



----
### xASL\_adm\_Load4DMemMapping.m

#### Format

```matlab
LoadFile = xASL_adm_Load4DMemMapping(x, WhichModality)
```

#### Description
Part of ExploreASL analysis module.
Loads data & maps it to memory mapping file on disc, if not done before.



----
### xASL\_adm\_LoadParms.m

#### Format

```matlab
[Parms, x, Oldx] = xASL_adm_LoadParms(ParmsPath[, x, bVerbose])
```

#### Description
This function loads the internal memory x struct, any
legacy \*\_parms.mat sidecar, any \*.json BIDS sidecar, to use scan-specific
parameters for image processing/quantification. Also, per BIDS
inheritance, any x.S.SetsID parameters (from participants.tsv) are loaded
as well. This function performs the following steps:

1. Load .mat parameter file
2. Load JSON file
3. Deal with warnings
4. Find fields with scan-specific data in x.S.Sets, and use this if possible (per BIDS inheritance)
5. Sync Parms.\* with x.(Q.)\* (overwrite x/x.Q)
6. Fix M0 parameter if not set



----
### xASL\_adm\_LoadX.m

#### Format

```matlab
[x[, IsLoaded]] = xASL_adm_LoadX(x, Path_xASL[, bOverwrite])
```

#### Description
This function loads x.Output & x.Output\_im struct fields
from the x.mat on the hard drive & adds them to the current x struct
located in memory. If it didnt exist in the x.mat, it will
set IsLoaded to false, which can be catched externally & a warning issued if managed so
in the calling function. If it didnt exist in the memory x
struct, or bOverwrite was requested, the contents of x.mat
will be loaded to the memory x struct


----
### xASL\_adm\_OrderFields.m

#### Format

```matlab
outStruct = xASL_adm_OrderFields(inStruct,orderStruct)
```

#### Description
Order fields in the structure **inStruct** to match
{{orderStruct}}, unmatching fields in inStruct are copied as
they are at the end, unmatching fields in **orderStruct** are
ignored. This is just a cosmetic change and no values are
edited.



----
### xASL\_adm\_OtherListSPM.m

#### Format

```matlab
[OtherListSPM, OtherListOut] = xASL_adm_OtherListSPM(OtherList, bList4D)
```

#### Description
bPadComma1 is to add the ,1 to the end of the pathstring, which SPM uses
to assign the first image of a 4D image array (OPTIONAL, DEFAULT = true)
bList4D: boolean, true for listing multiple 4D volumes separately in the
list (OPTIONAL, DEFAULT=true).



----
### xASL\_adm\_Par2Parms.m

#### Format

```matlab
parms = xASL_adm_Par2Parms(pathPar, pathParms[, bRecreate])
```

#### Description
Opens the Philips type PAR file. Reads the relevant DICOM headers and saves them to .MAT file.
Only recreates an existing file if bRecreate option is set to TRUE.


----
### xASL\_adm\_ParReadHeader.m

#### Format

```matlab
info =xASL_adm_ParReadHeader(filename)
```

#### Description
Function for reading the header of a Philips Par /
Rec  MR V4.\* file.




----
### xASL\_adm\_Remove\_1\_SPM.m

#### Format

```matlab
[OtherList] = xASL_adm_Remove_1_SPM(OtherList)
```

#### Description
Remove ,1 at end of OtherLists, if exists.
These are appended in CoregInit, OldNormalizeWrapper etc,
since this should allow 4rd dim (e.g. as in ASL4D).



----
### xASL\_adm\_ReplaceSymbols.m

#### Format

```matlab
strOut = xASL_adm_ReplaceSymbols(strIn, symbolTable[, bracketLeft, bracketRight])
```

#### Description
It takes the STRIN on input, then looks for symbols between BRACKETLEFT and BRACKETRIGHT and replaces these symbols in
in the string by the values provided in the SYMBOLTABLE as SYMBOLTABLE.SYMBOL, SYMBOLTABLE.D.SYMBOL, or SYMBOLTABLE.P.SYMBOL


----
### xASL\_adm\_ResetVisualizationSlices.m

#### Format

```matlab
[x] = xASL_adm_ResetVisualizationSlices(x)
```

#### Description
Removes any predefined slices that should be visualized,
allowing to show the default slices. Comes in handy when different
pipeline visualization parts are repeated.



----
### xASL\_adm\_SaveJSON.m

#### Format

```matlab
xASL_adm_SaveJSON(data, jsonFileName)
```

#### Description
Saves the values in the structure 'data' to a file in JSON format.



----
### xASL\_adm\_UnzipOrCopy.m

#### Format

```matlab
unpackedFiles = xASL_adm_UnzipOrCopy(srcDir, wildCard, destDir [, bOverwrite])
```

#### Description
This is a simple wrapper function to (g)unzip one or more files to the specified destination
directory. Existing files or directories will not be overwritten, unless forced with bOverwrite.
A regular file-copy will be used if the source files don't have gz or zip filename extensions.


----
### xASL\_adm\_Voxel2RealWorldCoordinates.m

#### Format

```matlab
[X Y Z] = xASL_adm_Voxel2RealWorldCoordinates(X,Y,Z,VoxelSize)
```

#### Description
Converts MNI coordinates from voxel coordinates/indices.
Assumes X Y Z = LR LeftRight AP AnteriorPosterior IS InferiorSuperior.
VoxelSize should be [1 3]-sized input.


----
### xASL\_adm\_ZipFileList.m

#### Format

```matlab
filepaths = xASL_adm_ZipFileList(strDirectory, strRegExp[, bRecurse, bUseGzip, nRequired])
xASL_adm_ZipFileList(strDirectory, strRegExp[, bRecurse, bUseGzip, nRequired])
```

#### Description
Zip the files that match regular expression STRREGEXP in the given directory STRDIRECTORY.
Zips recursively if specified in BRECURSE. Zips all files unless the number is specified
by NREQUIRED, if the number is not met, then does not zip anything and throws an error.

----
### xASL\_adm\_uiGetInput.m

#### Format

```matlab
[Parms] = xASL_adm_uiGetInput(Parms)
```

#### Description
Checks whether input fields are present, or requests them.



##BIDS

----
### xASL\_bids\_Add2ParticipantsTSV.m

#### Format

```matlab
xASL_bids_Add2ParticipantsTSV(DataIn, DataName, x, bOverwrite)
```

#### Description
This function adds metadata/statistical variables to the
participants.tsv in the root/analysis folder, by the following steps.
This function will iterate over Data provided at DataIn and fill them
in the participants.tsv, overwriting if allowed.
Empty data is filled in as 'n/a', and the first column "participants\_id"
is sorted for participants.

This function runs the following steps:

1. Admin - Validate that there are not too many columns
2. Admin - Detect nSubjectsSessions
3. Admin - Load pre-existing participants.tsv or create one
4. Admin - Get column number of data
5. Add data to CellArray
6. Sort rows on subjects
7. Fill empty cells
8. Write data to participants.tsv


----
### xASL\_bids\_Dicom2JSON.m

#### Format

```matlab
[parms pathDcmDictOut] = xASL_bids_Dicom2Parms(inp[, parmsfile, dcmExtFilter, bUseDCMTK, pathDcmDictIn])
```

#### Description
The function goes through the INP files, reads the DICOM or PAR/REC files and parses their headers.
It extracts the DICOM parameters important for ASL, makes sure they are in the correct format, if missing then
replaces with default value, it also checks if the parameters are consistent across DICOM files for a single sequence.




----
### xASL\_bids\_InsertJSONFields.m

#### Format

```matlab
[ChildJSON] = xASL_bids_InsertJSONFields(ParentJSON, ChildJSON[, Fields2Skip])
```

#### Description
This function takes all parameters from the "parent" JSON & moves them into the "child" JSON.
In case of co-existence of a field with different values,
then the value in the child JSON will prevail, per BIDS inheritance.

This function runs the following steps:

1. Load JSON or parms.mat (legacy), if the inputs were paths
2. Insert the fields
3. Save a new JSON file (if ChildJSON was a path)



----
### xASL\_bids\_PARREC2JSON.m

#### Format

```matlab
parms = xASL_adm_Par2Parms(pathPar, PathJSON)
```

#### Description
Opens the Philips type PAR file. Reads the relevant DICOM header fields and adds them to the .json sidecar file.


----
### xASL\_bids\_parms2BIDS.m

#### Format

```matlab
outBids = xASL_bids_parms2BIDS(inParms[, inBids, bOutBids, priorityBids])
```

#### Description
This functions takes two parameter structures and merges them. At the same time, renames all fields
according to the output type (note that only some fields have two standardised names different between the two formats.
In case of duplicities, takes the field value from the preferred format.
Also takes into account that the units in BIDS are s, but in xASL ms.
This function performs the following steps:

1. Define field names that need to be convert/renamed/merged
2. Convert XASL fields to the output format (BIDS or XASL)
3. Convert BIDS fields to the output format (BIDS or XASL)
4. Merge the BIDS and XASL fields, convert field values


##FSL

----
### xASL\_fsl\_RunFSL.m

#### Format

```matlab
[x] = xASL_adm_RunFSL(FSLCommand, x[, OutputZipping])
```

#### Description
This function runs an FSL command from ExploreASL:

1. Checking the FSL dir
2. Manage CUDA/CPU parallelization (currently disabled, WIP)
3. Setting up FSL environment
4. Running the command

Supports .nii & .nii.gz, Linux, MacOS & Windows (WSL)


----
### xASL\_fsl\_SetFSLdir.m

#### Format

```matlab
[FSLdir[, x, RootWSLdir]] = xASL_adm_SetFSLdir(x, bUseLatestVersion)
```

#### Description
This function finds the FSLdir & puts it out, also in
x.FSLdir to allow repeating this function without having to repeat
searching.
If the FSLdir & RootFSLDir are already defined in x.FSLdir & x.RootFSLDir, this function
is skipped.
Supports Linux, MacOS & Windows (WSL), & several different
default installation folders for different Linux
distributions


----
### xASL\_fsl\_TopUp.m

#### Format

```matlab
xASL_fsl_TopUp(InDir[, ScanType], x)
```

#### Description
This function runs FSL TopUp. It assumes that there are 2
TopUp images, i.e. 1 blip up & 1 blip down.

0. Admin: manage ScanType, NIfTI paths, create TopUp
parameter file for image to apply TopUp to & for the TopUp NIfTIs,
delete files from previous run, define the image with the
same acquisition parameters as TopUp (does the image
we apply TopUp to, have the Blip up or down?)
1. Register images to image that we apply TopUp to
(registration between blip up/down images is performed by
TopUp)
2. Run TopUp estimate (i.e. estimate the geometric distortion field from B0 NIfTI &
parameters file), this takes quite long. Also has a x.Quality=0 option that is very fast
but inaccurate, to try out this pipeline part. Before
TopUp, NaNs (e.g. from resampling) are removed from the images
TopUp is run with default settings
3. Apply TopUp


##Imaging

----
### xASL\_im\_BilateralFilter.m

#### Format

```matlab
[ovol] = xASL_im_BilateralFilter(volIM, mask, VoxelSize, x)
```

#### Description
This function runs a spatial lowpass temporally
highpass filter, and removes outliers within this signal, and adapts the
time-series accordingly.



----
### xASL\_im\_CenterOfMass.m

#### Format

```matlab
xASL_im_CenterOfMass(PathNIfTI, OtherList, AllowedDistance)
```

#### Description
This function estimates the center of mass of the image
matrix, and if this is too far off the current orientation
matrix center, the center will be reset.
This fixes any incorrect orientation outputted by the
scanner.
The realignment is only applied when any of the X/Y/Z
dimensions have a higher offset than AllowedDistance.




----
### xASL\_im\_CleanupWMHnoise.m

#### Format

```matlab
xASL_im_CleanupWMHnoise(InputPath, OutputPath, MinLesionVolume, pThresh)
```

#### Description
Threshold white matter lesions,
acknowledging the fact that they may be confluent with subresolution connection
through a dilation. This part is executed conservatively, as FLAIR hyperintensities
inside the GM can be erroneously segmented as WMH, and should not be lesion-filled
(otherwise these cannot be fixed later in the Structural module).

Note that LST lesion filling expects a probability map, doesnt work nicely with binary mask



----
### xASL\_im\_ClipExtremes.m

#### Format

```matlab
[NewIM] = xASL_im_ClipExtremes(InputIm, ThreshHigh, ThreshLow, bVerbose)
```

#### Description
Clips image to given percentile. The percentile is found
using non-zeros sorted intensities, so both isfinite & non-zeros.



----
### xASL\_im\_Column2IM.m

#### Format

```matlab
[ImageOut] = xASL_im_Column2IM(ColumnIn, BrainMask)
```

#### Description
This function "decompresses" an image matrix (or multiple matrices)
from a single-dimensional column, by reconstructing the image matrix
from the voxel positions within the BrainMask.
NB: Important to use the same BrainMask as used for converting the
image matrix to the column!
See also: xASL\_im\_IM2Column.m

The mask mostly used for xASL\_im\_IM2Column is x.WBmask, which completely
engulfes pGM, pWM & pCSF.


----
### xASL\_im\_CompareNIfTIResolutionXYZ.m

#### Format

```matlab
[IsEqualResolution] = xASL_im_CompareNIfTIResolutionXYZ(PathNIfTI1, PathNIfTI2)
```

#### Description
This function checks whether the X, Y and Z resolution of a
NIfTI with any number of dimensions is equal. It rounds for 2 floating
points, for both NIfTIs, to ensure that the same precision is compared.


----
### xASL\_im\_ComputeDice.m

#### Format

```matlab
DiceCoeff = xASL_im_ComputeDice(imA, imB)
```

#### Description
Calculate Dice coefficient of image overlap.



----
### xASL\_im\_CreateASLDeformationField.m

#### Format

```matlab
xASL_im_CreateASLDeformationField(x, bOverwrite, EstimatedResolution)
```

#### Description
This function smooths a transformation flow field to a lower
resolution. Usually, we use a high resolution anatomical
image (e.g. **3D T1w**) to obtain the flowfields from native
space to standard space, and apply these to the lower
resolution ASL images. Because of the resolution
differences, the flowfields need to be downsampled/smoothed,
to avoid deformation effects that are crispier than the
functional image that is investigated. This function
performs the following steps:

1. Obtain resolutions
2. Fill NaNs at edges y\_T1.nii flowfield to prevent interpolation artifact
3. Smooth flowfield
4. Fill NaNs at edges y\_ASL.nii

Note that if the resolution of ASL is not significantly (i.e. >0.5 mm in
any dimension) lower than T1w, the y\_T1.nii is copied to y\_ASL.nii
--------------------------------------------------------------------------------------------------------------------

----
### xASL\_im\_CreatePseudoCBF.m

#### Format

```matlab
xASL_im_CreatePseudoCBF(x, spatialCoV)
```

#### Description
This function creates a pseudo-CBF image from mean CBF template,
arterial transit time (ATT) bias field & vascular artifacts, weighted through spatial CoV
The first part of this code puts templates in the native space and
creates a pseudoCBF image from a downsampled pGM & pWM tissue (PseudoTissue). The latter
is used for registration but also as reference for the template
registration, to speed this up.
The second part of this code computes a pseudoCBF image based on the
pseudoTissue & the CBF templates of CBF, ATT biasfield and vascular peaks, based on spatial CoV.

This submodule performs the following steps:

1. Create the pseudoTissue CBF reference image, if it doesnt exist already
2. Create the native space copies of ASL templates, if they dont exist already
3. Spatial CoV input argument check
4. Load native space copies of templates
5. Create pseudoTissue from segmentations, mix this with the mean CBF template depending on spatial CoV
6. Create pseudoCBF reference image used for CBF-based registration
7. Scale mean\_PWI\_Clipped source image to the same range as PseudoCBF


----
### xASL\_im\_CreateSliceGradient.m

#### Format

```matlab
xASL_im_CreateSliceGradient(x)
```

#### Description
1. Create slice gradient in same space as input file
2. Reslice slice gradient to MNI (using existing ASL matrix changes from e.g. registration to MNI, motion correction, registration to GM)
3. Creating average slice gradient



----
### xASL\_im\_DecomposeAffineTransformation.m

#### Format

```matlab
[M, P] = xASL_im_DecomposeAffineTransformation(Mtransformation)
```

#### Description
This function splits a transformation matrix into individual
components, which can be useful to guide the SPM reslicing.
The components are the same as in spm\_(i)matrix.m, except for
the shearing: these are included in the rotations, and
the 90 degree rotations, these are separated.

Reason for the separation of the 90 degree rotations, is
that these indicate if orientations (transversal, coronal &
sagittal) have been switched in the NIfTI.

This can be useful to correct for any erroneous 90degree rotations from registration,
or to put two images in the same orientation order or voxelsize without applying
their subtle realignment (e.g. for manipulating registration references)

THEORY 90 degree rotations:
Any rotation will always swap the other dims (X rotation swaps Y/Z, Y
rotation swaps X/Z etc.) because they are perpendicular (haaks)

Dims X Y Z care for LR, AP and IS translation.
- X-rotation will rotate the transverse slice (LR <-> AP)
swapping Y (coronal) & Z (saggital)
- Y-rotation will rotate the coronal slice (IS <-> LR) slice,
swapping X (transversal) and Z (sagittal)
- Z-rotation will rotate the sagittal slice (AP <-> IS)
swapping X (transversal) and Y (sagittal)

E.g., MPRAGE is acquired in sagittal slices, and ASL/fMRI/BOLD in
transversal slices. This is an Y rotation (you look into the coronal
plane, rotate this, which will swap the sagittal slices into transversal)

This function performs the following steps:

0. Start with an output P and M struct
1. Obtain translations
2. Obtain zoom
3. Obtain 90degree rotations
4. Obtain subtle rotations & shearing
5. Check the rounding errors of the decomposition


----
### xASL\_im\_DetermineFlip.m

#### Format

```matlab
[QCstruct] = xASL_im_DetermineFlip(x,iS,PathOrientationResults,QCstruct)
```

#### Description
Check determinants, should be the same
before & after registration, otherwise a left-right flip is applied
This is not visible, but detrimental for image analysis/stats.



----
### xASL\_im\_DilateErodeFull.m

#### Format

```matlab
new_mask = xASL_im_DilateErodeFull(mask,type,kernel)
```

#### Description
Runs dilation or erosion on a binary mask in full three dimensions
It uses its own dilate\_erode function and crops the image so that it
contains only the mask.

Works only with odd sized kernels.



----
### xASL\_im\_DilateErodeSeparable.m

#### Format

```matlab
new_mask = xASL_im_DilateErodeSeparable(mask,type,kernel_x,kernel_y,kernel_z)
```

#### Description
Runs dilation or erosion on a binary mask separably in three dimensions
It uses its own dilate\_erode function and crops the image so that it
contains only the mask.

Works only with odd sized kernels



----
### xASL\_im\_DilateErodeSphere.m

#### Format

```matlab
el = xASL_im_DilateErodeSphere(R)
```

#### Description
3D structuring element (binary) sphere.



----
### xASL\_im\_DummyOrientationNIfTI.m

#### Format

```matlab
xASL_im_DummyOrientationNIfTI(PathSrc, PathRef, PathDummyOut[, bApplyRotationMinor, bApplyRotation90, bApplyZoom, bApplyTranslation])
```

#### Description
This function creates a dummy image as reference for xASL\_spm\_reslice,
allowing to only apply specific parts of the transformation between the
two images. E.g. only the rotation, or only the zooming.
This can be useful to correct for any erroneous rotations from registration,
or to put two images in the same space without applying their
realignment. This function performs the following steps:

1. Load orientations & calculate transformation
2. Calculate the desired transformation
3. Calculate new orientation matrix
4. Calculate the new image size
5. Save the dummy NIfTI



----
### xASL\_im\_EstimateResolution.m

#### Format

```matlab
[resFWHM, resSigma,resErr,imSmo,imMask] = xASL_im_EstimateResolution(imCBF,imGM,imWM,imMaskOrig,PSFtype,maxIter)
```

#### Description
NB: everything in this code is in voxels, not in mm



----
### xASL\_im\_Flip.m

#### Format

```matlab
[MatrixOut] = xASL_im_Flip(MatrixIn, varargin)
```

#### Description
Backwards compatibility for flipping left-right in standard
space (NB: this can be different than in native space!).



----
### xASL\_im\_IM2Column.m

#### Format

```matlab
[ColumnOut] = xASL_im_IM2Column(ImageIn, BrainMask[, ApplyShiftDim])
```

#### Description
This function "compresses" an image matrix (or multiple matrices)
for optimization of memory and CPU resources. The output column only includes
voxels that lie within the BrainMask. This excludes extracranial
zero-information voxels from computations and memory use.
NB: Important to use the same BrainMask for converting the
column back to an image matrix!
See also: xASL\_im\_Column2IM.m

The mask mostly used for xASL\_im\_IM2Column is x.WBmask, which completely
engulfes pGM, pWM & pCSF


----
### xASL\_im\_JointHist.m

#### Format

```matlab
imHist = xASL_im_JointHist(imA,imB[,imMask,minA,maxA,minB,maxB,nBins])
```

#### Description
It calculates a joint histogram of two images of any dimensions over a mask of the same size.
The boundaries and number of bins can either be given or min and max values are used. Values
outside of the bins are counted to the first/last bin.

----
### xASL\_im\_Lesion2CAT.m

#### Format

```matlab
LesionPathOut = xASL_im_Lesion2CAT(PathIn)
```

#### Description
For all lesion masks in the anatomical directory, load
them, merge them and save them for the CAT segmentation.
If there are no lesions found, the images are untouched.



----
### xASL\_im\_Lesion2Mask.m

#### Format

```matlab
LesionIM = xASL_im_Lesion2Mask(LesionPath, T1path, pGMpath, pWMpath, x)
```

#### Description
For a standard space lesion mask (or map), this stores
the lesion mask, and in additional its perimask (15 mm) and contralateral
mask, as 2nd and 3rd volumes.
It plots the masks on a T1 image, and masks the new masks with the
subjects' brainmask (pGM+pWM).



----
### xASL\_im\_M0ErodeSmoothExtrapolate.m

#### Format

```matlab
[ImOut] = xASL_im_M0ErodeSmoothExtrapolate(ImIn, x)
```

#### Description
This function erodes, smooths & extrapolates M0 in standard space.
It assumes that the M0 image is in standard space & that the GM & WM probability maps
are aligned. Here, we mask the M0, to remove high CSF signal and low extracranial signal,
enabling us to smooth the image without inserting wrong signal. See also
the ExploreASL manuscript for a more extensive explanation. This function
performs the following steps:

1. Mask: Load segmentations, create structural mask
2. Mask: Create intensity-based mask to remove extracranial signal
3. Mask: Erode the combined masks
4. Mask: Determine any odd borders
5. Smoothing
6. Extrapolating only
7. Scale back to the GM M0
8. Print visual QC figure

A visual QC figure is created, showing the M0 image processing steps for a single transversal slice (slice 53 in 1.5 mm MNI standard space)
OutputFile = fullfile(x.D.M0regASLdir,['M0\_im\_proc\_' x.P.SubjectID '.jpg']);
The original M0 image (a) is masked with a (pGM+pWM)>50% mask (b)
eroded with a two-voxel sphere to limit the influence of the ventricular and extracranial signal (c)
and thresholded to exclude significantly high (i.e. median + 3\*mean absolute deviation (MAD)) border region values (d)
This masked M0 image is smoothed with a 16x16x16 mm full- width-half-maximum Gaussian filter (Mutsaerts et al., 2018) (e)
after which the signal at the border is smoothly extrapolated until the full image is filled (f).
Whereas the masking avoids mixing with cerebrospinal fluid or extracranial signal, the extrapolation avoids M0 division artifacts


----
### xASL\_im\_MaskNegativeVascularSignal.m

#### Format

```matlab
[NegativeMask, TreatedCBF] = xASL_quant_DetectNegativeVascularSignal(x)
```

#### Description
This function segments clusters with significant negative
ASL signal. This can be tricky as there is also the negative tail of Gaussian noise
from the ASL subtraction. The image feature we use here, is that negative
vascular signal will be a relatively large region with
significant median negative value, whereas noise will be
regions with relatively small negative signal.
Negative signal from wrong background suppression timing
(e.g. in the first slice with 2D EPI) can be masked out with
this as well.
The procedure works as follows:

1) Obtain mask of negative voxels within pGM>0.5 mask
2) Obtain distribution of subzero clusters
3) Define the negative threshold
4) Create mask by thresholding whole image

Note that the definition of the threshold is obtained within
the GM only, but that this threshold is applied to the full image


----
### xASL\_im\_MaskPeakVascularSignal.m

#### Format

```matlab
[MaskIM, CBF] = xASL_quant_VascularContrast(PathCBF, Path_M0, CompressionRate, ClipThresholdValue, bClip)
```

#### Description
This function searches for an acceptable high
threshold as definition of high intra-vascular ASL signal.
It also allows to compress the values here (when
bClip==true). Compression retains some variability, but limits their outlying influence
on statistics.
The procedure works as follows:

1. Segment **GM** on **ASL** image at 80th percentile of **CBF** image distribution
2. For PWI & CBF images, select voxels higher than median + ClipThresholdValue **MAD**
Vascular artifacts will have a high intensity in both images, whereas errors by division by M0 will only have a high
intensity on the M0 image, and high values due to a biasfield will only be noticeable on the PWI image
3. Combine the two created masks
4. Obtain average clipping value from selected voxels from the combined masks
5. Apply compression if requested. If not, output image will
have NaNs for intra-vascular voxels.

Note that the definition of the threshold is obtained within
the GM only, but that this threshold is applied to the full
image to also remove extracranial extreme values.


----
### xASL\_im\_Modulation.m

#### Format

```matlab
xASL_im_Modulation(x)
```

#### Description
Combines the transformations to create Jacobians, &
multiplies the standard space segmentations with these to create volumetric
images for volumetric analyses.



----
### xASL\_im\_NormalizeLabelingTerritories.m

#### Format

```matlab
image_out = xASL_im_NormalizeLabelingTerritories( imageIN, GMmask, x)
```

#### Description
Normalizes per perfusion territory mask should be GM mask.



----
### xASL\_im\_PCA.m

#### Format

```matlab
[pc, score, eigenvalues, tsquare, loadings, Xmean] = xASL_im_PCA(dataIn)
```

#### Description
...



----
### xASL\_im\_PVCbspline.m

#### Format

```matlab
[imPVEC,imCBFrec,imResidual,FWHM] = xASL_im_PVCbspline(imCBF,imPV,bsplineNum)


```

#### Description
PVEc correction of ASL data using prior GM-,WM-partial volume maps.
Follows the principles of the PVEc algorithm by I. Asllani (MRM, 2008).
The PV-corrected CBF\_GM and CBF\_WM maps are approximated using an
uniformly sampled cubic B-splines.
Free for research use without guarantee.
Created by Jan Petr, j.petr@hzdr.de



----
### xASL\_im\_PVCkernel.m

#### Format

```matlab
[imPVEC,imCBFrec,imResidual] = xASL_im_PVCkernel(imCBF, imPV,kernel,mode)


```

#### Description
PVEc correction of ASL data using prior GM-,WM-partial volume maps.
Follows the principles of the PVEc algorithm by I. Asllani (MRM, 2008).
Free for research use without guarantee. If used in a study or
publication. Please, acknowledge the author.
Created by Jan Petr, j.petr@hzdr.de



----
### xASL\_im\_PreSmooth.m

#### Format

```matlab
pathOut = xASL_im_PreSmooth(pathRef,pathSrc[,pathSmo,resRef,resSrc,srcAffinePath, bInvAffine])
```

#### Description
It assumes that the FWHM is equal to voxel size, unless the real resolution is given.
Then takes into account the voxel sizes and orientation difference between the volumes, but
performs the smoothing according to the given real resolution (it is possible to supply the
resolution for just one image) - this should be helpful primarily when the
It creates a Gaussian covariance matrix for the reference, transforms it according to
the difference between the two images a produces the Gaussian covariance matrix
of the pre-smoothing in the source space. Then it performs the smoothing.

The following steps are performed:
1) Obtain the voxel size
2) Skip this function if reference resolution is equal to, or lower than source resolution
3) Deal with affine transformation
4) Obtain the transformation matrix from the Reference to the Source space
5) Apply the smoothing filter on the source image(s)
6) Save the smoothed image


----
### xASL\_im\_ProcessM0Conventional.m

#### Format

```matlab
[Corr_M0] = xASL_im_ProcessM0Conventional(ImIn, x)
```

#### Description
This function uses the conventional M0 masking,
and only a little smoothing, following what Philips uses for its 3D
{{GRASE}}. Advantages of the newer M0 processing in ExploreASL are the lack
of use of **M0** threshold-based masking, the removal of high CSF values and
higher **SNR** for **ASL** division.



----
### xASL\_im\_ProjectLabelsOverData.m

#### Format

```matlab
OutputIM = xASL_im_ProjectLabelsOverData(DataIM,LabelIM,x,ScaleFactorData,ScaleFactorLabel)
```

#### Description
This script projects labels over an image,
but works only in 2D. Make sure to make a 2D image from a 3D or 4D image
using xASL\_vis\_TransformData2View.m
can be used in combination with xASL\_vis\_Imwrite.m



----
### xASL\_im\_ResampleLinearFair.m

#### Format

```matlab
[output_res]=xASL_im_ResampleLinearFair(im_input,newsize)
```

#### Description
Downsample (or upsample, works similarly) old\_res image to
low\_res image, trilinear.

**NB:** new\_res should fit exactly integer fold in old\_res

**NB:** all dimensions of new\_res should have equal size



----
### xASL\_im\_RestoreOrientation.m

#### Format

```matlab
xASL_im_RestoreOrientation(PathNIfTI)
```

#### Description
This function reverts the NIfTI header orientation matrix
to the original orientation from the scanner/dcm2nii conversion.



----
### xASL\_im\_SkullStrip.m

#### Format

```matlab
xASL_im_SkullStrip(InPath, PathMNIMask, x, OutPath)
```

#### Description
Creates skull-stripped **T1w** image based on **MNI** -> native
space registration from segmentation.



----
### xASL\_im\_Smooth3D.m

#### Format

```matlab
[imSmo,imGaussX,imGaussY,imGaussZ] = xASL_im_Smooth3D(sigma,imIn,PSFtype)
```

#### Description
...



----
### xASL\_im\_Upsample.m

#### Format

```matlab
xASL_im_Upsample(PathOrig, PathDest, NewVoxelSize, LeaveEmpty, PaddingDim, Kernel)
```

#### Description
Upsamples an ASL image, without changing the orientation
matrix, which can be used e.g. for PVEc in higher
resolution but same space.



----
### xASL\_im\_ZeroEdges.m

#### Format

```matlab
[IM] = xASL_im_ZeroEdges(IM[, EdgeThicknessPerc])
```

#### Description
Resampling can sometimes give some strange errors near image edges. These should be NaNs,
but sometimes can be zeros or ones, or even weird numbers. For resampling, NaNs should be set to 0 (this is done
in another function) as they can influence the resampling (depending on the transformation matrix). To be sure
that the edges are nicely fixed, this function sets a border at the image matrix edges to zero.


----
### xASL\_im\_dilateROI.m

#### Format

```matlab
xASL_im_dilateROI(PathIn, PathTemp)
```

#### Description
...



----
### xASL\_im\_rotate.m

#### Format

```matlab
rotated = xASL_im_rotate(im, angle)
```

#### Description
Simple rotation of the first two dimension of a ND image by
0, 90, 180, 270 degrees.



----
### xASL\_import\_json.m

#### Format

```matlab
[x] = xASL_import_json(DataParFile)
```

#### Description
This function reads in a DATA\_PAR file and creates the x
structure. The name of the DATA\_PAR file is given as a string or
character array. The output is the x structure.

If the DATA\_PAR file is the dataset\_description.json file of the BIDS
standard, the x structure is created according to BIDS.




##Initialization

----
### xASL\_init\_ConvertM2JSON.m

#### Format

```matlab
[PathJSON] = xASL_init_ConvertM2JSON(PathM)
```

#### Description
This function converts and replaces the legacy data parameter m-format
by a JSON file. A DataPar.m was the settings/parameter file, specific to a dataset to be
processed by ExploreASL, now replaced to JSON by BIDS.
Note that the deployed/compiled version of ExploreASL
requires the JSON file, this function should not be compiled
along. This function performs the following steps:

1) Run the m-file to load parameters
2) Escape characters that are illegal in JSON
3) Write the JSON


----
### xASL\_init\_DefaultEffectiveResolution.m

#### Format

```matlab
[EffectiveResolution] = xASL_init_DefaultEffectiveResolution(PathASL, x)
```

#### Description
This ExploreASL module provides an educated guess on
the effective spatial resolution of ASL. This may depend on the
combination of acquisition PSF, reconstruction filter, head motion.
Note that the derived/processed images may have a lower effective
resolution because of smoothing effects from interpolation. Note that
this remains an educated-guess, the actual FWHM may still differ,
especially for 3D GRASE sequences, where e.g. the choice of number of
segments can affect the smoothness.

This function conducts the following steps:
1) Educated-guess FWHM
2) Attempt accounting for in-plane interpolation in reconstruction
3) Calculate and report effective spatial resolution


----
### xASL\_init\_DefineStudyData.m

#### Format

```matlab
[x] = xASL_init_DefineStudyData(x)
```

#### Description
This initialization wrapper initializes the parameters for
this pipeline run, i.e. subjects, sessions (runs), timepoints (visits),
exclusions, sites, cohorts etc.

Note that ASL sessions are defined here as what BIDS calls "runs".

The "longitudinal\_Registration functions here manage different
TimePoints, which is what BIDS calls "visits".
With different structural scans, from the same participant. This is
managed by subject name suffixes \_1 \_2 \_n, and can be used for comparing
visits in the population module, or running SPM's longitudinal within-subject
registration.

Parallelization is allowed here by calling ExploreASL different times,
where it divides the subjects/images for processing across the nWorkers,
using iWorker as the reference for the part that the current ExploreASL
call will process. This requires having a Matlab license that can be
started multiple times on a server, or alternatively running the
ExploreASL compilation, and doesn't require the Matlab parallel toolbox.

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
### xASL\_init\_FileSystem.m

#### Format

```matlab
[x] = xASL_init_FileSystem(x)
```

#### Description
This function initializes the file system used throughout ExploreASL, for processing a single dataset/scan.
It is repeated for each scan, and runs the following parts:
1) Create folders
2) Subject/session definitions
3) Add prefixes & suffixes
4) Add Subject-specific prefixes
5) Add sidecars


----
### xASL\_init\_InitializeMutex.m

#### Format

```matlab
[x] = xASL_init_InitializeMutex(x, ModuleName)
```

#### Description
This function initializes the mutex/lock system of
ExploreASL for a module. Mutex (for mutual exclusion) is a
synchronization mechanism for enforcing limits of access to data (here a
module for a single scan) to allow parallelization. It also allows
stopping and continuing of ExploreASL. This function runs the following
steps:
1) Lock folder management
2) Initialize mutex object


----
### xASL\_init\_LoadMetadata.m

#### Format

```matlab
[x] = xASL_init_LoadMetadata(x)
```

#### Description
This function loads all metadata used in the study, either statistical
covariates (age, MMSE) or groups to compare between (site, sequence,
cohort), or parameters to be used in quantification/image processing

These parameters should be provided in .mat files in the root analysis
folder. Each .mat file should contain a single type of metadata, and
the filename should equal the variable name.
Metadata content should be a cell array
with subjects as first column and metadata as last column.
Sessions (runs) can be included as second column.

Metadata can be in any string or numerical format.

participants.tsv is now added per BIDS. It looks for metadata in participants.tsv first,
before going through the mat-files

This function iterates through the following steps for each variable:

1) Admin (what nOptions do we call ordinal, convert subject numeric to
string, remove white spaces from data)
2) Get unique list of data options & check for missing data
3) Deal with data format (correct NaNs, deal with numeric vs strings)
4) Distinguish continous data (e.g. age) or ordinal data (groups to compare, e.g. cohort)
5) Check if data is complete for all subjects
6) Include complete data in x.S.SETS



----
### xASL\_init\_LongitudinalRegistration.m

#### Format

```matlab
[SubjectNlist, TimePoint, IsSubject, SubjectID_FirstVolume] = xASL_init_LongitudinalRegistration(x)
```

#### Description
This function initializes the longitudinal registration for ExploreASL,
which implements the SPM longitudinal registration.

This function recognizes and defines visits (i.e. multiple scans per
subject) in the initialization of ExploreASL, based on the suffixes \_1 \_2
\_n in the subject names (identified as foldernames).

Specifically, this function is called in the registration modules LongReg and DARTEL,
the first carrying out within-subject registration and
the second between-subject registration, based on the first time point
only.
For the first function, we specify here a list of visits/timepoints that
should be registered longitudinally, for the second function we specify a
list of first visits only, as the between-subject registration in
ExploreASL is based on the first scan (as opposed to the average
subject's scan).
This function runs the following steps:
1) Get TimePoint-list (list of visits)
2) Find subject IDs


----
### xASL\_init\_VisualizationSettings.m

#### Format

```matlab
[x] = xASL_init_VisualizationSettings(x)
```

#### Description
This function defines several visualization settings are
used throughout ExploreASL's pipeline and tools, assuming a [121 145 121]
matrix with 1.5 mm isotropic resolution in MNI space.


##Input and Output

----
### xASL\_io\_CreateNifti.m

#### Format

```matlab
xASL_io_CreateNifti(pathNewNifti, imNew, resMat, nBits, bGZip)
```

#### Description
This function creates a new NIfTI file, using the SPM "nifti" functionality, with the parameters
specified as input arguments. This function performs the
following steps:

1) Initialize NIfTI
2) Choose datatype (bit resolution)
3) Create scale slopes
4) Create orientation matrix
5) Write the new NIfTI, image matrix & scale slopes
6) Zip and deal with zipping (.nii vs. .nii.gz)


----
### xASL\_io\_DcmtkRead.m

#### Format

```matlab
header = xASL_io_DcmtkRead(filepath, bPixel)
```

#### Description
SHORT  Reads DICOM headers using DCMTK

FORMAT: header = xASL\_io\_DcmtkRead(filepath, bPixel)

INPUT:
filepath (string) - full path to the DICOM file
bPixel (bool) - read pixel data, default 0
OUTPUT:
header (structure) - structure containing parsed DICOM header


Calls the MEX function that uses DCMTK library to read the DICOM header.
To change which parameters are read and their names - the MEX file needs to be edited.
This function also corrects formating of certain parameters.



----
### xASL\_io\_MakeNifti4DICOM.m

#### Format

```matlab
xASL_io_MakeNifti4DICOM(PathIn, x)
```

#### Description
This function converts a NIfTI file to one that is ready to convert to DICOM for
PACS visualization purposes:
For scaling/visualization:
1) Remove peak signal
2) Remove valley signal
3) Remove NaNs
4) Rescale to 12 bit integers
5) Save NIfTI. We also zip the NIfTI as this NIfTI won't be opened by ExploreASL
6) Manage scale slope/datatype
7) Apply original orientation
8) Zip NIfTI


----
### xASL\_io\_PairwiseSubtraction.m

#### Format

```matlab
xASL_io_PairwiseSubtraction(InputFile,outputPath,do_mask,switch_sign)
```

#### Description
Subtracts controls from labels and takes mean.
Creates new perfusion-weighted delta\_M file, prefaced with 's'.
Converts into single precision floating point values (32 bit), removes scale slope.
Only runs if ['s' input\_file\_ASL] doesn't exist.
Remember to consider motion correction/ SPM realign (isotropically).
Alternative to this function is robust fit (Camille Maumet).



----
### xASL\_io\_ReadTheDicom.m

#### Format

```matlab
[Info] = xASL_io_ReadTheDicom(bUseDCMTK, DicomPath)
```

#### Description
This function tries to read a DICOM and throws a warning if it fails to


----
### xASL\_io\_SplitASL\_M0.m

#### Format

```matlab
xASL_io_SplitASL_M0(InPath,iM0)
```

#### Description
This function splits ASL4D & M0 if they were in the same sequence.
If dcm2niiX has already splitted the ASL4D NIfTI, this is reconstructed first.
If no M0 exists, or only ASL splitting is desired, leave iM0 empty ([])

Vendor product sequence examples:
GE 3D spiral sometimes puts the M0 at the last volume of the series -> iM0 = [2];
Philips 3D GRASE puts the M0 as control-label volume pair -> iM0 = [1 2];
Siemens 3D GRASE puts the M0 as the first volume -> iM0 = 1;


----
### xASL\_io\_dcm2nii.m

#### Format

```matlab
[niifiles, ScanNameOut, usedinput, msg] = xASL_io_dcm2nii(inpath, destdir, series_name, varargin)
```

#### Description
Convert DICOM NIfTI/BIDS format using the dcm2nii command line utility.



----
### xASL\_num2str.m

#### Format

```matlab
[DataOut] = xASL_num2str(DataIn[, f])
```

#### Description
when the provided data is numeric, this function will convert the number to string/characters,
rewriting NaN into 'n/a' (BIDS convention) but otherwise preserving the Matlab builtin functionality, also for the second argument "f".
If non-numeric data is provided, it is bypassed (avoiding any issues "num2str" will have with non-numeric data).
See builtin num2str for more details


##QC

----
### xASL\_qc\_AsymmetryIndex.m

#### Format

```matlab
[AI_perc] = xASL_qc_AsymmetryIndex(ImageIn)
```

#### Description
Extract voxel-wise asymmetry index for QC purposes.



----
### xASL\_qc\_CAT12\_IQR.m

#### Format

```matlab
[QA_Output] = xASL_qc_CAT12_IQR(InputImage, InputC1, InputC2, InputC3, bFLAIR)
```

#### Description
...



----
### xASL\_qc\_CollectParameters.m

#### Format

```matlab
x = xASL_qc_CollectParameters(x, iSubject, ScanType, CollectQCFunction)
```

#### Description
This function collects QC parameters for a module


----
### xASL\_qc\_CollectQC\_ASL.m

#### Format

```matlab
[x] = xASL_qc_CollectQC_ASL(x, iSubject)
```

#### Description
This functions collects QC parameters for the ASL module
These are stored in x.Output.ASL:
ID - SubjectName
ASL\_LR\_flip\_YesNo - Checks whether any image processing changed the left-right orientation
by checking whether the determinant differs between nii.mat & nii.mat0
SPM realign (too much motion is suspicious)
MotionMean\_mm    - mean motion
MotionExcl\_Perc  - percentage of excluded outliers
MotionMax\_mm     - max motion
MotionSD\_mm      - SD motion
ASL quantification (strange average CBF, or strange GM-WM contrast)
ASL acquisition parameters (should be fairly consistent over subjects/scans):
TE - echo time
TR - repetition time
RescaleSlope - Philips
Scaleslope - Philips
Matrix X Y Z - matrix size
Matrix Z - number of slices
VoxelSize X Y - in plane resolution
VoxelSize Z - slice thickness
RigidBody2Anat\_mm - Net Displacement Vector (RMS) from ASL to T1w image (mm) from registration


----
### xASL\_qc\_CollectQC\_Structural.m

#### Format

```matlab
[x] = xASL_qc_CollectQC_Structural(x, iSubject)
```

#### Description
This functions collects QC parameters for the structural module
These are stored in x.Output.Structural:
ID - SubjectName
T1w\_LR\_flip\_YesNo - Checks whether any image processing changed the left-right orientation
by checking whether the determinant differs between nii.mat & nii.mat0
LST output:
WMH\_vol\_mL        - WMH volume
WMH\_n             - WMH number of lesions
CAT12 output:
T1w\_IQR\_Perc      - CAT12 IQR quality parameter for T1w
volumetric: GM\_vol\_mL, WM\_vol\_mL, CSF\_vol\_mL, ICV\_vol\_mL, GM\_ICV\_Ratio


----
### xASL\_qc\_CollectQC\_func.m

#### Format

```matlab
[x] = xASL_qc_CollectQC_func(x, iSubject)
```

#### Description
This functions collects QC parameters for the func module
These are stored in x.Output.func:
ID - SubjectName
func\_LR\_flip\_YesNo - Checks whether any image processing changed the left-right orientation
by checking whether the determinant differs between nii.mat & nii.mat0
SPM realign (too much motion is suspicious)
MotionMean\_mm    - mean motion
MotionExcl\_Perc  - percentage of excluded outliers
MotionMax\_mm     - max motion
MotionSD\_mm      - SD motion
func quantification (strange average CBF, or strange GM-WM contrast)
CBF\_GM\_Median\_mL100gmin - median GM CBF
CBF\_WM\_Median\_mL100gmin - median WM CBF
SpatialCoV\_GM\_Perc      - GM spatial CoV
SpatialCoV\_WM\_Perc      - WM spatial CoV
CBF\_GM\_WM\_Ratio         - GM-WM CBF ratio
func acquisition parameters (should be fairly consistent over subjects/scans):
TE - echo time
TR - repetition time
RescaleSlope - Philips
Scaleslope - Philips
Matrix X Y Z - matrix size
Matrix Z - number of slices
VoxelSize X Y - in plane resolution
VoxelSize Z - slice thickness
RigidBody2Anat\_mm - Net Displacement Vector (RMS) from func to T1w image (mm) from registration


----
### xASL\_qc\_CollectSoftwareVersions.m

#### Format

```matlab
[x] = xASL_qc_CollectSoftwareVersions(x)
```

#### Description
This functions collects software versions for matlab, SPM, CAT, LST & ExploreASL
These are stored in x.Output.Software


----
### xASL\_qc\_CompareTemplate.m

#### Format

```matlab
[QC] = xASL_qc_CompareTemplate(x, ModPrefix, iSubjectSession)
```

#### Description
This function computes several advanced template-based QC parameters:
RMSE\_Perc        - Root Mean Square Error between image and template (%)
nRMSE\_Perc       - Same but then normalized
AI\_Perc          - Asymmetry Index between image and template (%)
Mean\_SSIM\_Perc   - mean structural similarity index -> xASL\_stat\_MeanSSIM.m
PeakSNR\_Ratio    - peak signal-to-noise ratio -> xASL\_stat\_PSNR.m


----
### xASL\_qc\_ComputeFoVCoverage.m

#### Format

```matlab
[CoveragePerc] = xASL_qc_ComputeFoVCoverage(InputPath, x)
```

#### Description
This function computes the intersection/overlap between
brainmask on field-of-view (FoV) of low resolution image
(native space) & the same brainmask with expanded FoV.
It uses the pGM+pWM+pCSF as brainmask
This assumes that the structural reference image has full brain coverage,
and was properly segmented into GM, WM and CSF
Also, we assume that the InputPath contains a single 3D volume


----
### xASL\_qc\_ComputeNiftiOrientation.m

#### Format

```matlab
[Struct] = xASL_qc_ComputeNiftiOrientation(x, PathNIfTI, Struct)
```

#### Description
...



----
### xASL\_qc\_CreatePDF.m

#### Format

```matlab
xASL_qc_CreatePDF(x[, DoSubject])
```

#### Description
This function iterates over all values in x.Output and all
images in x.Output\_im, and prints them in a PDF file.
x.Output & x.Output\_im should contain the QC/result output
of all ExploreASL pipeline steps.

Further code explanation:
Below, using the Matlab & SPM Figure tools we create an image, which is
then printed to a PDF file
fg = the main Figure handle
ax = "axes" handles, these are objects containing either 1) text or 2)
images, with fg as "parent" (1) & (2) images have ax as "parent"
Positions are calculated in such a way that 4 categories can be printed,
which will be the first 4 fields found in x.Output
then allowing 8 single slice images, and 15 text lines (name & value
columns)


----
### xASL\_qc\_FA\_Outliers.m

#### Format

```matlab
[FA_Outliers_mL] = xASL_qc_FA_Outliers(InputFA)
```

#### Description
Extract the number of FA outliers, i.e. values of FA
above 1 or below 0, from a FA image.




----
### xASL\_qc\_ObtainQCCategoriesFromJPG.m

#### Format

```matlab
xASL_qc_ObtainQCCategoriesFromJPG(x)
```

#### Description
This function obtains QC categories as covariant/set,
based on the JPGs in //Population/ASLCheck. These are initially sorted by
spatial CoV, and should be visually checked & put in the correct folder.
-------------------------------------------------------------------------------------------------------------------------

----
### xASL\_qc\_PCPStructural.m

#### Format

```matlab
[anatQA] = xASL_qc_PCPStructural(PathT1, Pathc1T1, Pathc2T1, x, PopPathT1)
```

#### Description
This function computes several anatomical QC parameters as proposed in SPM Univariate Plus:
WM\_ref\_vol\_mL    - volume of the WM reference region (mL)
WMref\_vol\_Perc   - same but as percentage of total WM volume
SNR\_GM           - GM signal-to-Noise Ratio (SNR), ie the mean intensity within GM divided
by SD of WM reference region. Higher = better.
CNR\_GM\_WM        - GM-WM Contrast-to-Noise Ratio (CNR), i.e. the mean of GM - mean of WM
divided by the SD of the WM reference region. Higher = better.
FBER\_WMref\_Ratio - Foreground to Background Energy Ratio (FBER), i.e. the variance of voxels within the brain (in pGM+pWM mask)
divided by the variance of voxels in the WM reference region. Higher = better.
EFC\_bits         - Shannon Entropy Focus Criterion (EFC), i.e. the entropy of voxel intensities proportional to the maximum
possibly entropy for a similarly sized image. Indicates ghosting and head motion-induced blurring. Lower = better.
Mean\_AI\_Perc     - mean relative voxel-wise absolute Asymmetry Index (AI) within the brain (pGM+pWM mask) (%)
SD               - same but SD (%)

REFERENCES:
Preprocessed Connectome Project Quality Assurance Protocol (QAP):
http://preprocessed-connectomes-project.org/quality-assessment-protocol/
http://ieeexplore.ieee.org/document/650886/


----
### xASL\_qc\_PrintOrientation.m

#### Format

```matlab
xASL_qc_PrintOrientation(DIR, reg_exp_INPUT,OUTPUT_DIR,Name);
```

#### Description
Check orientation of niftis, useful to detect
accidental left-right flips (all other flips will be visible).
translations, rotations or shears are not to be worried about,
only negative zooms. This can be detected by negative determinants.
So orientation parameters and determinants should be similar across
all scans from single scanner/coil, and registration should not
give negative determinant.



----
### xASL\_qc\_TanimotoCoeff.m

#### Format

```matlab
TC = xASL_qc_TanimotoCoeff(Image1, Image2[, imMask, type])
```

#### Description
Compares images Image1 and Image2 within the mask imMask. TYPE specifies the input data type.

RATIONALE:   Note that the Tanimoto Coefficient is a measure of image
overlap/intersection, similar to the Dice coefficient. With
the option type 3, this is a fuzzy coefficient, which
doesn't require to convert the two images to a binary mask.
The TC can be interpreted as a stringent Kappa, ranging from
0 (completely dissimilar) to 100% (identical images).
Assuming that perfect registration should not lead to
identical images but still retain physiological differences,
TC>70% can be regarded as excellent image agreement. The TC
will be overestimated when smoothing, but this may lead to more stable artifact detection.


----
### xASL\_qc\_WADQCDC.m

#### Format

```matlab
xASL_qc_WADQCDC(x, iSubject[, ScanType])
```

#### Description
This QC function runs WAD-QC specific Python script to zip QC information &
incorporate this into a DICOM field for analysis on the
WAD-QC server, by the following:

1. Define QCDC script: this is the Python script written by
Gaspare, edited by Joost Kuijer & copied to the EPAD
CustomScripts folder of ExploreASL
2. Python installation location is checked, with several known
locations, for several servers. If we cannot find it,
the QCDC is not ran
3. Previous QCDC results are cleaned. QCDC stores all its
results in a separate folder (Something like 2 layers up from the current
folder, here referred to as QCDCDir = [x.D.ROOT 'qcdc\_output'])
from these result files, only the filled DICOM file is
interesting, all the rest are copies of the QC results
that we embedded into the DICOM
4. Run QCDC (if Python installation detected)
The following files need to be set as executable:
('QCDC', 'src', 'qc\_data\_collector.py')
('QCDC', 'src', 'bash', 'create\_dcm\_only\_wadqc.sh')
('QCDC', 'src', 'bash', 'sendwadqc.sh')
5. Clean up new QCDC results (as above) & move the filled
DICOM to ['qcdc\_' DummyFile] within the current ScanType
folder
6. Sending the DICOM to the WAD-QC server using storescu


----
### xASL\_qc\_WADQC\_GenerateDescriptor.m

#### Format

```matlab
xASL_qc_WADQC_GenerateDescriptor(x, iSubject)
```

#### Description
This QC function generates a JSON descriptor for Gaspare'
QCDC script, by the following steps:

a) include information about where to find the dummy DICOM (i.e. placeholder DICOM)
b) For ExploreASL' QC fields (as passed through in
x.Output), here we note all these QC fields for each
ScanType, as the x.Output should have been collected
equally in the QC file 'QC\_collection\_SubjectName.json'
by function xASL\_qc\_CollectParameters
c) Subfunction xASL\_qc\_WADQC\_images - Includes visual standard space QC
images, by searching them on prescribed paths within the
Population folder (where currently all derivatives reside)
d) Insert the PDF report; this PDF report is
subject-specific, not scan-specific. For completeness it
is added to each QCDC descriptor
e) Add WAD-QC server details (i.e. IP address etc)
f) Save the Descriptor JSON file.


----
### xASL\_qc\_temporalSNR.m

#### Format

```matlab
tSNR = xASL_qc_temporalSNR(pathIm4D,pathImTissueProb)
```

#### Description
This function computes several temporal SNR QC parameters as proposed in SPM Univariate Plus:
tSNR.tSNR\_GM\_Ratio      : mean GM signal / std GM over time
tSNR.tSNR.tSNR\_WM\_Ratio : mean WM signal / std WM over time
tSNR.tSNR.tSNR\_CSF\_Ratio: mean CSF signal / std CSF over time
tSNR.tSNR\_WMref\_Ratio   : mean signal/std over time in eroded deep WM
tSNR.tSNR\_GMWM\_Ratio    : mean (GM+WM) signal / sqrt(std(GM+WM)^2+std(WMref)^2)
tSNR.tSNR\_GMWM\_WMref\_Ratio: mean (GM+WM) signal / std WMref over time
tSNR.tSNR\_Physio2Thermal\_Ratio: sqrt((tSNR(GM+WM)/tSNR\_GMWM\_WMref\_Ratio))^2-1)
tSNR.tSNR\_Slope\_Corr:
Differences to the SPM U+ suggestion:
- eroded WM is used for estimating background noise
- Brainmask is determined in the same way as the structural anatQC,
- CSF is determined from the pGM&pWM maps;

REFERENCES:
1) Thomas Liu (2016). Noise contributions to the fMRI signal: An overview NeuroImage, 343, 141-151
http://dx.doi.org/10.1016/j.neuroimage.2016.09.008
2) Cesar Caballero-Gaudes and Richard C. Reynolds (2016). Methods For Cleaning The BOLD fMRI Signal. NeuroImage, 154,128-149
3) Lawrence Wald and Jonathan R Polimeni (2016). Impacting the effect of fMRI noise through
hardware and acquisition choices ??? Implications for controlling false positive rates. NeuroImage, 154,15-22
4) SPM Utility + toolbox. Cyril Pernet. https://osf.io/wn3h8/


##Quantization

----
### xASL\_quant\_AgeSex2Hct.m

#### Format

```matlab
[Hematocrit] = xASL_quant_AgeSex2Hct([age, sex])
```

#### Description
This function estimates a participants hematocrit, based on
literature-based values for age and sex. It performs the following steps:

1. Warning unknown age/sex
2. Imputing unknown age/sex
3. Define hematocrit per age for unknown sex
4. Define Hematocrit per age for males
5. Define Hematocrit per age for females



----
### xASL\_quant\_FEAST.m

#### Format

```matlab
xASL_quant_FEAST(x)
```

#### Description
This function quantifies ATT using the FEAST equations,
using crushed and non-crushed sessions, of which the ratio is
proportional to ATT.
Note that the order of sessions should be 1) crushed 2) non-crushed

This function runs the following steps:
1. Skip this function if no FEAST data available
2. Admin
3. Load data & correct for timing differences (PLD etc)
4. Smooth and clip CBF maps & FEAST ratio
5. Compute TT maps


----
### xASL\_quant\_GetControlLabelOrder.m

#### Format

```matlab
[ControlIM, LabelIM, OrderContLabl] = xASL_quant_GetControlLabelOrder(FramesIn, x)
```

#### Description
This function automatically checks (and corrects if required)
the control and label order of ASL timeseries
based on the larger signal in control volumes.
It supposes that data is acquired in pairs.


----
### xASL\_quant\_Hct2BloodT1.m

#### Format

```matlab
BloodT1 = xASL_quant_Hct2BloodT1(Hematocrit, Y, B0, bVerbose)
```

#### Description
This function converts hematocrit to blood T1, according to
calculations defined by Patrick Hales. With courtesy and thanks!
Note that we assume a venous O2 saturation of 68% (Yv=0.68)

This function performs the following steps:
1) Check fraction vs percentage hematocrit & Y, should be between 0 and 1
2) Specify defaults (Hb, Fe)
3) Perform calculation
4) Convert s to ms
5) Print what we did
--------------------------------------------------------------------------------------------------------------

----
### xASL\_quant\_M0.m

#### Format

```matlab
[M0IM] = xASL_quant_M0(M0IM, x)
```

#### Description
This function quantifies the M0, except for the difference in voxel size
between the M0 and ASL source data (which is scaled in
xASL\_wrp\_ProcessM0.m). This function runs the following steps:

1. Correct scale slopes, if Philips
2. Skip M0 quantification if ~x.ApplyQuantification(4)
3. Set TR specifically for GE
4. Check for correct TR values
5. Quantify the M0, either for single 3D volume or slice-wise
6. Apply custom scalefactor if requested (x.M0\_GMScaleFactor)



----
### xASL\_quant\_SinglePLD.m

#### Format

```matlab
[ScaleImage[, CBF]] = xASL_quant_SinglePLD(PWI, M0_im, SliceGradient, x)
```

#### Description
This script performs a multi-step quantification, by
initializing a ScaleImage that travels through this script & gets changed by the following quantification
factors:
1)    PLD scalefactor (gradient if 2D multi-slice) (if x.ApplyQuantification(3))
2)    Label decay scale factor for single (blood T1) - or dual-compartment (blood+tissue T1) model, CASL or PASL
Single-compartment model: Alsop MRM 2014
Dual-compartment model: Wang MRM 2002: Gevers JMRI 2012 (if x.ApplyQuantification(3))
3)    Scaling to physiological units [ml/gr/ms =>ml/100gr/min =>(60,000 ms=>min)(1 gr=>100gr)]
(if x.ApplyQuantification(3))
4)    Vendor-specific scalefactor (if x.ApplyQuantification(1) -> future move to dcm2niiX stage)
Finally, we:
5)    Divide PWI/M0 (if x.ApplyQuantification(5))
6)    Print parameters used
Note that the output always goes to the CBF image (in the
future this could go to different stages, e.g. dcm2niiX or
PWI stage)


##SPM

----
### xASL\_spm\_BiasfieldCorrection.m

#### Format

```matlab
xASL_spm_BiasfieldCorrection(PathIn, SPMdir, Quality, MaskName, PathOut)
```

#### Description
This function is a wrapper around the SPM "old segment"
function, for biasfield removal. It is tested for M0 and mean control
images. It conducts the following steps:

1) Create implicit mask
2) Define SPM 'old segmentation' settings
3) Run SPM 'old segmentation'
4) Delete temporary files
5) Rename temporary SPM file into output file


----
### xASL\_spm\_affine.m

#### Format

```matlab
xASL_spm_affine(srcPath, refPath, fwhmSrc, fwhmRef[,otherList, bDCT])
```

#### Description
This SPM wrapper runs SPM's old normalize-estimate function, which calculates the affine transformation (i.e. linear + zooming and shearing) that is required to
align the source image with the reference image. By default it does not estimate the low-degree Discrete Cosine Transform (DCT) to have a simple affine transformation
but this can be enabled in this wrapper. Also note that this affine transformation uses a correlation cost function, hence it requires the source and reference images
to have similar contrasts and resolution - or provide the resolution estimates so the smoothing can be done.
This function does not change the orientation header by default, but it allows to change those of others through the otherList. If bDCT is used or no otherList given,
it stores its estimated affine transformation as orientation difference matrix in a file with the same path but \_sn.mat extension.
For the provided smoothing FWHM, note that smoothnesses combine with Pythagoras' rule (i.e. square summing).



----
### xASL\_spm\_coreg.m

#### Format

```matlab
xASL_spm_coreg(refPath, srcPath[, OtherList, x, sep, FastReg])
```

#### Description
This SPM wrapper runs SPMs coregister-estimate function, which calculates the 6 parameter rigid-body transformation (a.k.a. linear) that is required to
align the source image with the reference image. This 6 parameter transformation (i.e. 3 XYZ translations and 3 rotations) is applied to
the orientation header of the source NIfTI image, and also to the images provided in OtherList (optional).
Note that this SPM registration function uses a normalized mutual information (NMI) by default, enabling registration between two images
with different contrast.
Note that this algorithm will use the first volume of a multi-volume NIfTI



----
### xASL\_spm\_deface.m

#### Format

```matlab
xASL_spm_deface(PathIn, bReplace)
```

#### Description
This function removes the face from an anatomical NIfTI
image, e.g. T1w or FLAIR, for disidentification/privacy purposes.
When this script is run after the ExploreASL structural
module, it does a pretty good job even for 2D images.
However, note that this can always fail, strip part of
the brain, or change the output of pipelines. So best not
to compare results from defaced and non-defaced images.
Also, note that defacing makes it difficult to ensure that
the FLAIR and T1w are from the same subject.




----
### xASL\_spm\_deformations.m

#### Format

```matlab
xASL_spm_deformations([x,] PathIn, PathOut[, Interpolation, InverseSpace, AffineTrans, DeformationPath])
```

#### Description
This ExploreASL wrapper manages the SPM deformation tool.
It takes multiple (ExploreASL pipeline) transformations and combines/concatenates them
into a single transformation prior to applying it to the input images.
This allows to apply multiple transformations with a single interpolation, avoiding
propagation of undesired interpolation effects. Mainly used to get native
space images into standard space, or vice versa.
Best to combine as many files as possible within this function, since the
deformation calculation (which is the most computation intensive part) needs to be performed once for multi-file resampling



##Statistics

----
### xASL\_stat\_AtlasForStats.m

#### Format

```matlab
[x] = xASL_stat_AtlasForStats(x)
```

#### Description
This function loads atlases, checks them, and
their ROI names, for later use as ROI definition in xASL\_stat\_GetROIstatistics
Note that the atlases should be integer values, or different 4rd
dimensions (i.e. multiple images), that are mutually
exclusive. This function takes the following steps:
1) Load atlas ROI names
There should be a TSV sidecar to the atlas NIfTI file, as
explained above.
2) deal with memory mapping
3) Resample atlas 50 1.5 mm^3 MNI
4) Converted atlas with integers to 4D binary image
5) Convert/compress masks into Columns
6) Print atlas overview image


----
### xASL\_stat\_ComputeDifferCoV.m

#### Format

```matlab
diffCoV = xASL_stat_ComputeDifferCoV(imCBF)
diffCoV = xASL_stat_ComputeDifferCoV(imCBF,imMask)
diffCoV = xASL_stat_ComputeDifferCoV(imCBF,imMask,bPVC,imGM,imWM)
diffCoV = xASL_stat_ComputeDifferCoV(imCBF,imMask,bPVC,imGM,imWM,b3D)
```

#### Description
It calculates the spatial DiffCoV value on finite part of imCBF. Optionally a mask IMMASK is provide,
and PVC is done for bPVC==2 using imGM and imWM masks and constructing
pseudoCoV from pseudoCBF image. For bPVC~=2, imGM and imWM are ignored. It is calculated in 2D or assuming also 3D edges based on B3D.
Calculate derivate spatial CoV, by summing up differences in CBF between neighbors.
The derivative uses Sobels filter.


----
### xASL\_stat\_ComputeMean.m

#### Format

```matlab
[CBF_GM CBF_WM] = xASL_stat_ComputeMean(imCBF [,imMask,nMinSize,bPVC,imGM,imWM])
```

#### Description
It behaves in a similar way as VAR.


----
### xASL\_stat\_ComputeSpatialCoV.m

#### Format

```matlab
sCov = xASL_stat_ComputeSpatialCoV(imCBF)
sCov = xASL_stat_ComputeSpatialCoV(imCBF,imMask)
sCov = xASL_stat_ComputeSpatialCoV(imCBF,imMask,nMinSize)
sCov = xASL_stat_ComputeSpatialCoV(imCBF,imMask,nMinSize,bPVC,imGM,imWM)
```

#### Description
It calculates the spatial CoV value on finite part of imCBF. Optionally a mask IMMASK is provide,
ROIs of size < NMINSIZE are ignored, and PVC is done for bPVC==2 using imGM and imWM masks and constructing
pseudoCoV from pseudoCBF image. For bPVC~=2, imGM and imWM are ignored


----
### xASL\_stat\_EqualVariancesTest.m

#### Format

```matlab
[resTest, P] = xASL_stat_EqualVariancesTest(X[, alpha, type])
```

#### Description
Brown-Forsythe or Levene's test for equality of variances. The response variable is
transformed (yij = abs(xij - median(xj)) for Brown-Forsythe and yij = abs(xij - mean(xj))
for Levene's test). And then runs a one-way ANOVA F-test to check if the variances are equal.


----
### xASL\_stat\_GetROIstatistics.m

#### Format

```matlab
[x] = xASL_stat_GetROIstatistics(x)
```

#### Description
This function computes mean and spatial CoV for each ROI,
in a [1.5 1.5 1.5] mm MNI space,
with several ASL-specific adaptions:

1. Skip ROI masks that are smaller than 1 mL
as this would be too noisy for ASL (skipped when x.S.IsASL==false)
2. Expand each ROI mask such that it has sufficient WM
content (skipped when IsASL==false)
3. Create for each ROI mask a left, right and bilateral copy
4. Iterate over all subjects:
- a) Load partial volume maps
- b) Correct for WMH SEGM -> IS THIS STILL REQUIRED???
- c) Load data
- d) Show ROIs projected on ASL image
- e) Actual data computations
Partial Volume Correction (PVC) options:
PVC==0 -> perform masking only, no regression
PVC==1 -> single compartment regression, for pGM
PVC==2 -> dual compartment regression for pGM & pWM (i.e. normal
PVC)
Here we mask out susceptibility artifacts (including
outside FoV) for all ASL computations, and also mask
out vascular artifacts for computing the mean.

While other artifacts/FoV can be masked out on population
level (i.e. >0.95 subjects need to have a valid signal in a
certain voxel), vascular artifacts differ too much in their
location between subjects, so we mask this on subject-level.

Note that the words "mask" and "ROI" are used
interchangeably throughout this function, where they can
have a different or the same meaning
PM: WE COULD CHANGE THIS, INTO MASK BEING USED TO EXCLUDE
VOXELS AND ROI FOR INCLUDING VOXELS



----
### xASL\_stat\_MadNan.m

#### Format

```matlab
y = xASL_stat_MadNan(x[,flag, dim])
```

#### Description
Calculates a Median/Mean Absolute deviation, but ignoring NaNs in the calculation.
xASL\_stat\_MadNan(X) or xASL\_stat\_MadNan(X,0) computes xASL\_stat\_MeanNan(ABS(X-xASL\_stat\_MeanNan(X))
xASL\_stat\_MadNan(X,1) computes xASL\_stat\_MedianNan(ABS(X-xASL\_st\_MedianNan(X)).



----
### xASL\_stat\_MeanSSIM.m

#### Format

```matlab
mssim=xASL_stat_MeanSSIM(imRef,imSrc[,dynRange])
```

#### Description
Calculates the similarity index according to Want et al.


----
### xASL\_stat\_MultipleLinReg.m

#### Format

```matlab
[b,CI,pval,stats] = xASL_stat_MultipleLinReg(X,Y[,bIntercept])
```

#### Description

Performs a multiple linear regression Y=b\*X+a and provides the intercept and regression coefficients beta
including their significance and confidence intervals. It calculates additionally the goodness of the fit.


----
### xASL\_stat\_PSNR.m

#### Format

```matlab
PSNR=xASL_stat_PSNR(imRef,imSrc)
```

#### Description
Calculates the PSNR, needs two input arguments - 3D images of the same size.
Uses 95% percentile instead of MAX.




----
### xASL\_stat\_PrintStats.m

#### Format

```matlab
[x] = xASL_stat_PrintStats(x)
```

#### Description
This function prints an overview of statistics from
data that were acquired per ROI, in a TSV file. It starts by
printing covariates (called "Sets"). Rows will be
subjects/sessions, columns will be the sets and
ROI-statistics:
1) First remove previous TSV-file, if already existed
printing to a TSV file can be tricky if it is opened by
Excel. Make sure to close previous versions first,
otherwise this part will crash.
2) Print overview of sets to TSV
as explained above. Uses subfunction
xASL\_stat\_CreateLegend to put legends. Aim is to create a
single TSV file that has a proper overview of the data,
& is self-explanatory to those reading/using it.
3) Define number of ASL sessions, force to 1 in case of TT or volume metrics
4) Print the overview


----
### xASL\_stat\_QuantileNan.m

#### Format

```matlab
y = xASL_stat_QuantileNan(x[,quant, dim])
```

#### Description
Calculates a quantile, but ignoring NaNs in the calculation


----
### xASL\_stat\_RobustMean.m

#### Format

```matlab
[NoOutliers, iOutliers, ThresholdDeviation] = xASL_stat_RobustMean(IM, ParameterFunction)
```

#### Description
This function detects outlier images, that can be used to create
a robust average, e.g. for template or biasfield creation. This is based either on the sum-of-squares
with the mean image (SoS), or on the average relative asymmetry index (AI). Images that are
median+/-3 mad off are defined as outliers. MAD = median/mean absolute difference


----
### xASL\_stat\_ShapiroWilk.m

#### Format

```matlab
[H, P, W] = xASL_stat_ShapiroWilk(x[, alpha])
```

#### Description
Performs the statistical test of normality - null hypothesis is that the sample is from normal
distribution with unspecified mean and variance. Based on the sample kurtosis it performs either
Shapiro-Wilk (for platykurtic) or Shapiro-Francia test (for leptokurtic).

----
### xASL\_stat\_StdNan.m

#### Format

```matlab
y = xASL_stat_StdNan(x[,w,dim])
```

#### Description
It behaves in a similar way as VAR - it directly passes all arguments to xASL\_stat\_VarNan.


----
### xASL\_stat\_SumNan.m

#### Format

```matlab
y = xASL_stat_SumNan(x[,dim])
```

#### Description
It uses the function SUM, but it sets all the NaNs to zero before calling it.


----
### xASL\_stat\_VarNan.m

#### Format

```matlab
y = xASL_stat_VarNan(x[,w,dim])
```

#### Description
It behaves in a similar way as VAR.


----
### xASL\_stat\_fcdf.m

#### Format

```matlab
F = xASL_stat_fcdf(F,M,N)
```

#### Description
Calculates the cumulative distribution function of the F-distribution for degrees of freedom M,N at value F.


----
### xASL\_stat\_tcdf.m

#### Format

```matlab
F = xASL_stat_tcdf(T,nu)
```

#### Description
Calculates the cumulative distribution function of the Student's t-distribution for degrees of freedom nu at value T.


----
### xASL\_stat\_ticdf.m

#### Format

```matlab
T = xASL_stat_ticdf(P,nu)
```

#### Description
Calculates the inverse of cumulative distribution function of the Student's t-distribution for degrees of freedom nu at value P.


----
### xASL\_stat\_ttest.m

#### Format

```matlab
[H,P,CI,stats] = xASL_stat_ttest(X[,M,alpha,tail,dim])
```

#### Description
Performs a t-test that the distribution of the input data X has a mean different from 0 (or from a
given mean M, or that the distributions X and Y have different means). A normal distribution of the data
with an unknown variance is assumed.


----
### xASL\_stat\_ttest2.m

#### Format

```matlab
[H,P,CI,stats] = xASL_stat_ttest2(X,Y[,alpha,tail,vartype,dim])
```

#### Description
Performs a unpaired t-test that the distribution of the input data X has a mean different from that of Y.
A normal distribution of the data with an unknown variance is assumed.


----
### xASL\_str2num.m

#### Format

```matlab
[DataOut] = xASL_str2num(DataIn)
```

#### Description
str2num wrapper, replacing 'n/a' with NaN (BIDS convention)
and converting only strings to numbers. Also allows inputting cells.



##Visualization

----
### xASL\_vis\_AddIM2QC.m

#### Format

```matlab
[x] = xASL_vis_AddIM2QC(x,parms);
```

#### Description
Checks which images already are loaded, and  adds new image.



----
### xASL\_vis\_CreateVisualFig.m

#### Format

```matlab
[ImOut, FileName] = xASL_vis_CreateVisualFig(x, ImIn, DirOut, IntScale, NamePrefix, ColorMap, bClip)
```

#### Description
This function creates a visualization Figure by merging flexibly rearranging NIfTI slices, input matrix or
path, managing colormaps for different merged image layers. Current use is for visual QC figures and overview in papers.
Function is structured as:

1. Admin, deal with input arguments
2. Process image layers separately
\* xASL\_im\_TransformData2View: Reshapes image data into visualization figure
\* xASL\_im\_ClipExtremes: Clips image to given percentile
also we scale for peak intensity, we make sure that there is no
visible clipping/distortion
\* Convert to colors, using any input colormaps
3. combine image layers, using input argument IntScale
4. print figure

This function assumes that the first image is a grayscale background
image (e.g. for transparancy reasons), if there are multiple
images


----
### xASL\_vis\_CropParmsAcquire.m

#### Format

```matlab
[xmin xmax ymin ymax] = xASL_vis_CropParmsAcquire(temp_image)
```

#### Description
Goes from outside to inside to acquire crop settings.
Works with grayscale images (2 dimensions per slice).
Image position information (2D matrix) should be first
2 dimensions. Could include colordimension later on.



----
### xASL\_vis\_CropParmsApply.m

#### Format

```matlab
ImageOut = xASL_vis_CropParmsApply(ImageIn,CropParameters)
```

#### Description
This function crops 2D image matrices.




----
### xASL\_vis\_Imwrite.m

#### Format

```matlab
[ImOut] = xASL_vis_Imwrite(ImIn, PathOut[, ColorMap, bRescale])
```

#### Description
This functions takes an input image matrix, interpolates it
to HD resolution (1920x1080) for visibility, and saves the image as jpg.
This function avoids the graphic interface of Matlab, for running from CLI
Careful: this function overwrites any existing PathOut.




----
### xASL\_vis\_OverlapT1\_ASL.m

#### Format

```matlab
xASL_vis_OverlapT1_ASL( x, ASL)
```

#### Description
Part of ExploreASL.
Shows spatial agreement ASL and probability maps.



----
### xASL\_vis\_TileImages.m

#### Format

```matlab
...
```

#### Description
Merges selected slices (3D) into one single 2D picture.
Plots all slices in one figure with specified rows and
columns, aiming for a square tile.

PM: can be extended to multiple slices



----
### xASL\_vis\_TransformData2View.m

#### Format

```matlab
FigureOut = xASL_vis_TransformData2View(ImagesIn, x)
```

#### Description
This function changes the dimensionality and reshapes the input images
in such a way that they are nicely tiled in a mosaic for visualization purposes.
Reshaping a series of images with this function can be useful for
visualization of SPM/voxel-based analyses.


----
### xASL\_vis\_VisualQC\_TopUp.m

#### Format

```matlab
[MeanAI_PreTopUp_Perc, MeanAI_PostTopUp_Perc] = xASL_vis_VisualQC_TopUp(PathPopB0, PathPopUnwarped, x, iSubject, CheckDir)
```

#### Description
This function creates a Figure that showes the effect of TopUp
with 6 images with axial slices: the NormPE, RevPE and
their difference image in colorscale, and this before (upper
row) & after (lower row) TopUp.



----
### xASL\_vis\_VisualizeROIs.m

#### Format

```matlab
xASL_vis_VisualizeROIs(x, ROI_T1_list, ROI_FLAIR_list)
```

#### Description
Creates for each subject a JPEG image containing
the original T1w, WMH\_SEGM and T1w after lesion-filling.



