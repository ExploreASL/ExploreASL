# SPM xASL Functions

----
### xASL\_Copy.m

#### Format

```matlab
xASL_Copy(SrxASL_SysCopyath, DstPath[, bOverwrite, bVerbose)])
```

#### Description
Copies a file to a file or a directory to a directory. For a file, it zips it in the end if the destination path contains nii.gz.
It also makes sure that only one of .nii and .nii.gz exists
in the destination directory. It is faster than the default
Matlab function, and runs on multiple OSes.

NB: This function calls xASL\_SysMove for the actual copying.
Run xASL\_SysMove instead of xASL\_Move if you don't want the
.nii|.nii.gz management/checking


----
### xASL\_Move.m

#### Format

```matlab
xASL_Move(SrcPath, DstPath[, bOverwrite, bVerbose])
```

#### Description
Moves a file to a file, a file to a directory, or a directory to a directory. It keeps the initial extensions, no unzipping or zipping
after the move. But it makes sure that only one of .nii and .nii.gz exists in the destination directory.
Bypass inefficient matlab stuff on linux and windows, but
can only move on the same file system.

NB: This function calls xASL\_SysMove for the actual moving.
Run xASL\_SysMove instead of xASL\_Move if you don't want the
.nii|.nii.gz management/checking


----
### xASL\_SysCopy.m

#### Format

```matlab
xASL_SysCopy(SrcPath, DstPath, bOverwrite, bVerbose)
```

#### Description
Copies a file to a file or a directory to a directory. Bypass inefficient matlab stuff on linux and windows,
but can only move on the same file system.

----
### xASL\_SysMove.m

#### Format

```matlab
xASL_SysMove(SrcPath, DstPath[, bForce, bSourceCheck])
```

#### Description
Moves a file to a file, a file to a directory, or a directory to a directory. SBypass inefficient matlab stuff on linux and windows, but can only move on same file system!

----
### xASL\_TrackProgress.m

#### Format

```matlab
xASL_TrackProgress(iCurrent[, iMax])
```

#### Description
Counts the percentage of the work done and display on the screen. Either iCurrent of iMax are given. Or only the percentages are given.

----
### xASL\_adm\_ConvertSeconds2TimeString.m

#### Format

```matlab
TimeString = xASL_adm_ConvertSeconds2TimeString(Seconds)
```

#### Description
Converts number to time
input hh (with minutes in fractions/floating point) -> output hhmm
Inverse from xASL\_adm\_ConvertTime2Nr.


----
### xASL\_adm\_ConvertSlash.m

#### Format

```matlab
[newString] = xASL_adm_ConvertSlash( StringOriginal,ForceUnix)
```

#### Description
Converts Windows forward slashes to backward slashes
Prevents confusion file separation & regular expression forward slashes
in Windows.

----
### xASL\_adm\_CreateDir.m

#### Format

```matlab
status = xASL_adm_CreateDir(strPath)                  create all missing subdirs
status = xASL_adm_CreateDir(strPath, strBranch)       create strBranch (if missing) under existing strPath
status = xASL_adm_CreateDir(strPath, nMaxNewDirs)     impose limit on the number of new directories
```

#### Description
Recursively creates missing directories at the given path or for given subdirectories, with an option
to limit the number of newly created directories.

----
### xASL\_adm\_DeleteFileList.m

#### Format

```matlab
filepaths = xASL_adm_DeleteFileList(strDirectory, strRegEx[, bRecurse, nRequired])
xASL_adm_DeleteFileList(strDirectory, strRegEx[, bRecurse, nRequired])
```

#### Description
Delete the files that match regular expression STRREGEXP in the given directory STRDIRECTORY.
Deletes recursively if specified in BRECURSE. Deletes all files unless the number is specified
by NREQUIRED, if the number is not met, then does not delete anything and throws an error.

----
### xASL\_adm\_GetFileList.m

#### Format

```matlab
filepaths = xASL_adm_GetFileList(strDirectory[, strRegEx, mode, nRequired, bGetDirNames])
```

#### Description
List files or directories from a given path. And optionally uses regular expressions to filter the result
with option to set a minimal requirement on the number of results.

----
### xASL\_adm\_GzipNifti.m

#### Format

```matlab
pathOut = xASL_adm_GzipNifti(pathIn [,bOverwrite])
```

#### Description
Take the input file, zips it, overwriting any existing zipped file and return the path of the zipped file.


----
### xASL\_adm\_ManageMoCoMat.m

#### Format

```matlab
xASL_adm_ManageMoCoMat(PathIn)
```

#### Description
This function manages the orientation matrices that SPM puts
in an external .mat sidecar file if there are more than 1
volumes. The first volume should be equal to the orientation
header of the NIfTI, if not, we assume that the NIfTI
header is correct.
This function performs several checks & corrects if
necessary, combined with throwing a warning:
A) the nVolumes in .mat & .nii image should be equal, if not, delete sidecar
B) .mat should have more than one volume, if not delete sidecar
C) If there are illegal numbers in the diagonal of the .mat
orientation matrices (here only checked for zeros or non
finite values) then the .mat is removed
D) If this is true for the first volume only, the .mat is
retained but the first volume orientation is overwritten
with a zero matrix


----
### xASL\_adm\_UnixPath.m

#### Format

```matlab
[PathIs] = xASL_adm_UnixPath(PathIs)
```

#### Description
This function performs the following steps to convert a path to a path that is compatible with the Unix-filesystem
as used in e.g. Linux/MacOS/Windows Subsystem for Linux (WSL).

1. Skip this function without Unix-filesystem
2. Trim whitespace
3. Selectively convert forward to backward slashes (ignore already escaped whitespace)
4. Escape characters and residual whitespaces (ignore already escaped whitespaces)
5. If WSL: add mounting prefix


----
### xASL\_adm\_UnzipNifti.m

#### Format

```matlab
pathOut = xASL_adm_UnzipNifti(pathIn[, bOverwrite])
```

#### Description
Takes the input file, unzips if needed, delete the zipped file and return the path to the unzipped file.
If the input is already unzipped, then does nothing, but returns the original filename - so it
can be run just to be sure a file is unzipped without much overhead.
Returns error if more than one file is in the archive, if the filename does not exist, is a directory etc.
If there's a NII and NII.GZ already existing, then return error, or just overwrite in case overwrite is set to 1


----
### xASL\_adm\_ZipFileNameHandling.m

#### Format

```matlab
[srcOut, dstOut] = xASL_adm_ZipFileNameHandling(srcIn, dstIn)
```

#### Description
Adjusts the source and destination filenames of a nifti file to reflect if NII or NII.GZ exist on the input.
If either .nii or .nii.gz is corrupt, it automatically deletes the corrupt one and keeps the healthy one,
while reporting a warning. This happens when you restart the pipeline after it crashed, if it crashed while unzipping.


----
### xASL\_bids\_csv2tsvReadWrite.m

#### Format

```matlab
[PathTSV, CellContents] = xASL_bids_csv2tsvReadWrite(PathIn[, bDeleteCSV, bWriteTSV])
```

#### Description
This function PathIn and loads it, also trying CSV or TSV
extensions if these exist. It outputs the contents to a cell array. If a
CSV file exists but not a TSV file, it converts and replaces the CSV to
TSV file, per BIDS. This function has the following parts:

1) Read the CSV or TSV file
2) Write the TSV file (if requested)
3) Delete the CSV file (if requested)


----
### xASL\_csvRead.m

#### Format

```matlab
[CellContents] = xASL_csvRead(PathCSV)
```

#### Description
This function loads a comma-separated value (csv) file - which
is the format that BIDS prefers - and outputs it to a cell array.


----
### xASL\_csvWrite.m

#### Format

```matlab
xASL_csvWrite(InputCell, PathCSV, bOverwrite)
```

#### Description
Rudimentary function, please use xASL\_tsvWrite instead.
For usage, type help xASL\_tsvWrite.
This function will still work though.


----
### xASL\_delete.m

#### Format

```matlab
xASL_delete(InputPath)
```

#### Description
Delete the file in the given path. If a NIFTI file with extension '.nii' or '.nii.gz' is given,
Then delete both the .nii and .nii.gz files.

----
### xASL\_exist.m

#### Format

```matlab
xASL_exist(PathIn[,Type])
```

#### Description
Check if the given path exists, wrapper around the Matlab
exist function, to allow checking for either .nii or .nii.gz
Otherwise, exist is used normally.


----
### xASL\_fileparts.m

#### Format

```matlab
[Fpath, Ffile, Fext] = xASL_fileparts(InputPath)
```

#### Description
Returns the path, file name, and file extension for InputPath using the fileparts.m function.
If a file ending at nii.gz is given, then the whole nii.gz is returned as the extension.
Does not verify the existence of the file, or existence of .nii or .nii.gz

----
### xASL\_im\_ConvertMap2Mask.m

#### Format

```matlab
[IMout] = xASL_im_ConvertMap2Mask(IMin)
```

#### Description
Provides a robust way of conversion of
a continuous map to a binary mask, which can be used for lesions, ROIs,
or tissue probability maps. Based on the assumption that a map should
be thresholded at 50% to form a map, which is often the case for
automatic segmentations.


----
### xASL\_im\_DistanceTransform.m

#### Format

```matlab
[dist, x, y, z] = xASL_im_DistanceTransform(im)
```

#### Description
Calculates the distance transform in a binary image
Uses Borgefors Chamfers computation of Euclidean distance in 3D using a
5x5x5 window.


----
### xASL\_im\_FillNaNs.m

#### Format

```matlab
xASL_im_FillNaNs(InputPath[, UseMethod, bQuality])
```

#### Description
This function fills any NaNs in an image. In SPM, any voxels
outside the boundary box/field of view are filled by NaNs
when resampling. These NaNs can confuse some algorithms,
hence it doesn't hurt replacing them in some cases (e.g. for
flowfields). Also, smoothing restricted in a mask is done in
ExploreASL with the function xASL\_im\_ndnanfilter, after
first setting all voxels outside the mask to NaN. In this
case, this functon can be useful to extrapolate the smoothed
image to avoid any division artifact near brain edges (e.g.
for reducing the M0 image to a smooth biasfield).
This function performs the following 3 steps:
1) Load image
2) Replace NaNs
3) Save image
--------------------------------------------------------------------------------------------------------------------

----
### xASL\_im\_LesionRemoval4CAT.m

#### Format

```matlab
[Ycls, LesionImOut] = xASL_im_LesionRemoval4CAT(Ycls, PathIn)
```

#### Description
For all lesion masks in the anatomical directory, remove
them from the current segmentations.


----
### xASL\_im\_ResampleIM.m

#### Format

```matlab
[imOut] = xASL_im_ResampleIM(imIn,matIn,matOut,dimOut,intMethod)
```

#### Description
Resamples an image using Matlab interp3 function.


----
### xASL\_im\_SaveOriginal4CAT.m

#### Format

```matlab
xASL_im_SaveOriginal4CAT(Ycls, PathIn)
```

#### Description
Save the segmentation before lesion masking.


----
### xASL\_im\_conv3Dsep.m

#### Format

```matlab
[imConv] = xASL_mex_conv3Dsep(im,kX,[kY,kZ])
```

#### Description
3D separable convolution with a supplied kernel
It converts the results to double
Returned is the convoluted image
The wrapper makes sure that kX are Nx1 format, removes nan, and removes
excessive zeros at the ends.


----
### xASL\_im\_ndnanfilter.m

#### Format

```matlab
[Y,fil] = xASL_im_ndnanfilter(X,filterType,F,WNAN)
```

#### Description

This function applies a 3-dimensional convolution of X with given kernel.
NaNs elements are taken into account (ignored).

By default, edges are not padded and one-sided filter is used at the image edges.

Notes:
\* Accepts empty value for any input. When X is empty, the program can
be used as a N-dimensional window generator.
\* NaNs elements surrounded by no-NaNs elements (which will depend on
window width) are the ones that will be interpolated. The others
are leaved untouched.
\* When WNAN=2, the programs acts like an NAN-interpolat/GAP-filling,
leaving untouched the no-NaNs elements but the filtering is
perfomed anyway. I recommend the default behaviour (WNAN=0) in order
to keep the filtered data in the workspace, and then use the code
at the end of this function to get/remove the interpolated NaNs
\* To achieve similar results as ndnanfilter previously, use same F
as with the 'rect' filter.
\* Note that the FWHM of Gaussian is given in VOXELS, not in mm
\* For the Gaussian filter, use (previous N, new FWHM)
N= 1 ~ FWHM 0.94
N= 2 ~ FWHM 1.885
N= 4 ~ FWHM 3.76
N= 6 ~ FWHM 5.652
N= 8 ~ FWHM 7.536
N=10 ~ FWHM 9.42
N=12 ~ FWHM 11.3
N=16 ~ FWHM 15.07
N=20 ~ FWHM 18.84
N=10/2.355 ~ FWHM 4
Basically divide by 1.06


----
### xASL\_io\_Nifti2Im.m

#### Format

```matlab
imOut = xASL_io_Nifti2Im(niftiIn [, ImageSize])
```

#### Description
This function loads a NIfTI image matrix with flexible input
(as explained under INPUT: niftiIn). It does the following.
1) Try to load a NIfTI
2) If NIfTI successfully loaded, try to load the NIfTI image
3) If the above didnt work, try to create a dummy image
4) Convert to single precision data format
5) Also able to load NIfTI as .nii.mat format



----
### xASL\_io\_ReadNifti.m

#### Format

```matlab
[NiftiObject, pathIn] = xASL_io_ReadNifti(pathIn)
```

#### Description
Read Nifti file given by the path. Return the NII object. And also return the actual path to the loaded
Nifti if by any reason the name changed during the function runtime (e.g. unzipping).


----
### xASL\_io\_SaveNifti.m

#### Format

```matlab
xASL_io_SaveNifti(pathOrigNifti, pathNewNifti, imNew[, nBits, bGZip, newMat])
```

#### Description
It loads the pathOrigNifti, takes all the parameters from it, and creates a new Nifti file with
these parameters, but new image matrix from imNew. It saves the result in pathNewNifti.


----
### xASL\_round.m

#### Format

```matlab
[OutputN] = xASL_round(InputN[, PrecisionN])
```

#### Description
Recent Matlab versions support a second input that specifies that number of decimals to round at,
but earlier Matlab versions do not support this. For backward compatibility, use this wrapper instead of round.


----
### xASL\_spm\_admin.m

#### Format

```matlab
[InPath] = xASL_spm_admin(InPath, bPadComma1)
```

#### Description
n/a


----
### xASL\_spm\_reslice.m

#### Format

```matlab
xASL_spm_reslice(refPath, srcPath[, srcAffinePath, bInvAffine, bQuality, NewName, InterpolationPar])
```

#### Description
This wrapper runs SPM's reslice function (a.k.a. coregister: reslice) which resamples a source image into the space of a reference image,
taking into account any orientation differences between the two images that are defined in the orientation matrix in the NIfTI header.
When the source image contains multiple volumes, they are all resampled.
The source image will get the same orientation matrix as the reference image,
as it is now in the same space as the reference image. This can be useful when two images are to be compared voxel-by-voxel, e.g. when
overlaying a CBF image over a structural image, or when wanting to use a mask from a different space. When after running this function, the
reference and source images are not in alignment, they need to be registered first (i.e. xASL\_spm\_register).
Resampling/reslicing needs an interpolation method to know from which voxels of the source image, a voxel in the new image will be computed.
A simplistic explanation is that this determines the number of surrounding neighborhood voxels it uses to determine the new voxel value.
The example syntax below would reslice/resample the CBF image to the T1w space (assuming the Affine was done to register CBF to T1w)
It also works with the external .mat file of the source file that has the same name as the source file. It also can optionally take a \_sn.mat containing the
affine transformation information.


----
### xASL\_spm\_smooth.m

#### Format

```matlab
xASL_spm_smooth(pathIn, fwhmSmooth[, pathNew])
```

#### Description
This SPM wrapper runs SPM's smooth function, which spatially smooths the input image with a Gaussian kernel.
In the case of multiple volumes (i.e. a 4D NIfTI), each 3D volume is spatially smoothed separately.
Note that smoothnesses combine with Pythagoras' rule (i.e. sum quadratically)



----
### xASL\_stat\_MeanNan.m

#### Format

```matlab
y = xASL_stat_MeanNan(x[,dim])
```

#### Description
It calculates the sum using the SUM functions and divides by the number of values but ignoring NaNs.


----
### xASL\_stat\_MedianNan.m

#### Format

```matlab
y = xASL_stat_MedianNan(x[,dim])
```

#### Description
It calculates the MEDIAN along the given dimension, but it sets all the NaNs to zero before calling it.


----
### xASL\_tsvRead.m

#### Format

```matlab
[CellContents] = xASL_tsvRead(PathTSV[, bStruct])
```

#### Description
This function loads a tab-separated value (TSV) file - which
is the format that BIDS prefers - and outputs it to a cell array.


----
### xASL\_tsvWrite.m

#### Format

```matlab
xASL_tsvWrite(InputCell, PathTSV[, bOverwrite, bCSV])
```

#### Description
This function loads a cell array and prints it to a
tab-separated value (TSV) file, which is the format that BIDS prefers.


----
### xASL\_wrp\_DARTELSaveIntermedTrans.m

#### Format

```matlab
xASL_wrp_DARTELSaveIntermedTrans(Yy,u,odim,rdim,idim,Mar,mat,M0,M1,fname,it)
```

#### Description
n/a


----
### xASL\_wrp\_GSSaveIntermedTrans.m

#### Format

```matlab
xASL_wrp_GSSaveIntermedTrans(y,idim,odim,rdim,M0,M1,R,M1t,M1r,fname,it)
```

#### Description
n/a


