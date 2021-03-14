# ExploreASL v1.5.0

----
## Versions included software
Versions included & used third-party tools (see /External/README_SPM.txt):

SPM12 7219 
CAT12 r1615 
LST 2.0.15 

----
## Feature improvements (still backward compatible)

* #39: Create PV-corrected GM & WM CBF maps in native space
* #56,#410: Add Mindboggle atlas to ExploreASL and restructure general atlas access in population module
* #283: `xASL_stat_GetROIstatistics` provide more feedback on missing images
* #299: Move `CustomScripts` with study-specific scripts to a separate repository
* #302: Remove server calls in CAT12 functions
* #313: Move GUI to a separate repository
* #351: T2 and T1c files are now also aligned to the T1w and outputted to standard space
* #354: Added an option `x.DummyScanPositionInASL4D` that removes marked dummy scans when splitting ASL to ASL+M0+dummy
* #356,#396,#397: Internally restructure SliceTime allowing ExploreASL now to work with multi-band 2D EPI as well or any other SliceTime order

----
## Work in progress

### ASL-BIDS
* #163,#189,#357,#373: Conversion from DICOM to BIDS
* #334,#382: Import of PAR-REC to BIDS
* #343: Add separate M0 option to mTrial import

### Compilation/stand-alone version
* #335: All input arguments can be passed in the deployed mode
* #380: Enable advanced input parsing for xASL compiled

----

## Bug Fixes

* #228: Fix CAT12 warnings with non-existent field cm
* #272: Fix errors in JSON import of ASL sidecars
* #273,#285,#291,#363: Minor fixes in input parameter administration
* #276,#280,#288,#329,#400: Fix error in reading TSV files with unclear number of columns
* #282: Population module is run serially in otherwise parallel mode
* #292: `xASL_qc_SortBySpatialCoV` now use all subjects without skipping
* #305: `xASL_adm_UnixPath`: bug with Windows+WSL
* #306,#378: Fix JSON reading with 'i' interpreted as a complex number
* #309: Fix non-linear registration of T1w to standard space with too high T1w values
* #312: `xASL_stat_GetROIstatistics` fix skipping of actual ROI extraction
* #325: `xASL_adm_CleanUpBeforeRerun` delete files correctly
* #339: Fix JSON reading of special characters
* #399: Fix special characters in Windows filenames
* #405: Fix range-check error in Background Suppression timing calculation
* #406: Fix `xASL_stat_MedianNan` for all-NaN input
* #408,#409: Skip missing fields in CAT during reports in compiled ExploreASL

----

## Documentation
* #7: Create README files in subfolders and added to interactive documentation 
* #279, #345: Move documentation to a separate repository
* #300,#318: Improve ExploreASL tutorial
* #355: Documentation improvements regarding input parameters

## Testing
* #193,#350,#398: Testing DICOM to BIDS conversion against a reference
* #294: Implement initial unit testing framework
* #326: Parse warnings/errors from all log in all subdirectories
* #369: Unit testing of `xASL_test_getLogContent`
* #371,#376,#404: Testing script for the DICOM->BIDS->Legacy conversion and processing

# ExploreASL v1.4.0

----
## Versions included software
Versions included & used third-party tools (see /External/README_SPM.txt):

SPM12 7219 
CAT12 r1615 
LST 2.0.15 

----
## Feature improvements (still backward compatible)
* #20 Implement BASIL -> changed order quantification/masking, and always save resampled PWI4D for quantification, facilitating BASIL
* #35 Calculation background suppression efficiency for pseudo-M0 -> in case of missing separate M0 images but still using background suppressed mean control as pseudo-M0, a single correction value (3D) or slice-wise correction value (2D) are applied to the pseudo-M0 image/slices
* #131 ExploreASL_GUI beta-testing enhancements set 1 -> Aesthetic improvements in certain modules, fixed incorrect removal of Philips-related json-sidecar fields in DCM2BIDS / Import module. Correction of ASL image flickering bug and ability for the user to subset without having to reload the data. Added ability to clarify data type of variables without the need to reload data. Added auto-select / auto-complete functionality in the ParmsMaker module as soon as the user indicates an analysis directory (i.e auto-completion of SliceReadoutTime for pCASL)
* #166 Update ADNI import in ExploreASL_Import
* #179 Finish creation average maps for CICERO -> improvements xASL_wrp_CreatePopulationTemplates
* #190 Hammers atlas option CAT12 -> restored the original Hammers atlas option for a colleague (note it's license though!)
* #225 Create DataPar option for running SPM12 longitudinal registration
* #241 Add warning when loading data without x output structure 
* #243 Shorten SPM initialization time -> removed configuration loading of unused toolboxes (check if the SPM DICOM import module is still needed)

----
## Work in progress


### ASL-BIDS
* #193 Comparing BIDS folders (for testing purpose)
* #226 Add new DICOM tags to DCMTK import

----

## Bug Fixes
* #191 WMH warning when no FLAIR analyzed -> this warning is now removed if no FLAIR was present in the data
* #248 Temporary fix native space processing -> in Population module
* #252 Population modules analysis masks - minor errors
* #267 Error in reading JSONs from EPAD

----

## Documentation
* #196 All contributors -> all contributors are now automatically added to the main README.md
* #217 Documentation improvements
* #219 Add user to documentation

----


# ExploreASL v1.3.0

----
## Versions included software
Versions included & used third-party tools (see /External/README_SPM.txt):

SPM12 7219 
CAT12 r1615 
LST 2.0.15 

----
## Feature improvements (still backward compatible)

* #59: Assign "weights" to status files, allowing the external Python ExploreASL GUI to provide a better estimate of the progress
* #71 Remove custom lesions `(Lesion\_(FLAIR|T1)\_\\d\\.nii)` from `WMH_SEGM.nii`
* #76: Improve functionality of ExploreASL through Command Line Interface (CLI), i.e. without using the Graphical User Interface (GUI) of Matlab
* #123 Create status files also for skipped processing parts: this is mainly the case for running the structural module without a FLAIR scan. Having all status files helps third-party tools such as the Python ExploreASL GUI to know that processing has succesfully completed (duplicate issues #137 and #129)
* The same was done for the ASL realignment status file, in case realignment is skipped for a 3D scan
* #145: Improve .nii(.gz) management in xASL_spm_deformations: allow either .nii or .nii.gz as input, treat them equally, and when .nii.gz is provided as output path, zip the resulting deformed image

----
## Work in progress
* #32 Docker integration
* #55 ExploreASL GUI, written in Python
*  #106 `xASL_im_SplitImageLabels`: Allow splitting labels, and warping them to standard space.
This is part of a continuous development on creating average flow territory templates and figures.
* #162 Remove `bNativeSpaceProcessing` from TestDataSet for now, return this when `bNativeSpaceProcess` is made more modular

### ASL-BIDS
* #82 Avoid 4D files with nT==1, which is not allowed in the BIDS validator
### Compilation/stand-alone version
* #88 `xASL_SysMove` error in Windows when a path includes whitespaces ' '

----

## Bug Fixes

* #85 Improvement ApplyQuantification
* #99 Improve loading of metadata (`xASL_str2num` & `xASL_init_LoadMetaData`)
* #105 In case of missing data, fill `x.S.SUBJECTID` and `S.DAT` data for the last subject/session
* #119 Create status file for last subject, if it's processing is skipped
* #120 `xASL_im_CreateAnalysisMask` in native space mode when in parallel execution
* #138 Fix structure TestDataJSON & its JSON files
* #139 Ensure that x fields are not case sensitive, by e.g. replacing `strcmp` by `strcmpi`
* #141 Solve conflicts between develop and master, these were minor edits that weren't implemented in develop yet
* #143 ensure that VBA masks are also created for a 3D spiral sequence (this was not created yet as the susceptibility masks were missing for this sequence)
* #148 Syncing the ROI/lesion processing of T1w & FLAIR
* #151 Minor bugfixes for TopUp
* #177 Fix `xASL_adm_UnzipNifti` & `xASL_io_SaveNifti` when path is incomplete
* * Skip warnings for small populations (in the start of `xASL_module_Population`)
* Fix regular expression in `xASL_init_LoadMetaData`
* `iRow` counting fix in `xASL_bids_Add2ParticipantsTSV`
----

## Documentation
#7 Documentation/revamp `xASL_im_ClipExtremes`
#159 Ensure that all sequence notations use underscores instead of whitespace, e.g. `3D_spiral` instead of `3D spiral`


## Testing
* #86 `xASL_qc_TestExploreASL`: improve Table creation
* #112 Save Tanimoto Coefficient (i.e. a fuzzy overlap/Dice score) of the final ASL-T1w registration
* #128 Improved one internal test dataset
* #130 `xASL_qc_TestExploreASL`: Complete functionality by allowing Windows parallelization & testing the Windows ExploreASL compilation. Also added unit testing framework in the same issue.


# ExploreASL v1.2.2

----

## Bug Fixes
* #119 xASL_wrp_LinearReg_T1w2MNI.m: ROI .nii files correctly aligned with T1

----

# ExploreASL v1.2.1

----

## Bug Fixes
* xASL_qc_TestExploreASL: Remove locked folders if rerun
* #90 xASL_io_Nifti2Im: manage the detection of odd scaling
* #93 xASL_SysMove: diz illegal symbols Windows
* #102 xASL_wrp_CreatePopulationTemplates: minor bugfix
* #104 Fix creation ResultsTable for TestCases (also #86)
* #114 Acquiring Matlab version doesn't crash anymore in deployed mode
* #115 Fixing NaNs problem in M0 mask computation
* #118 xASL_adm_GzipAllFiles doesn't crash anymore in Windows

* #116 ExploreASL testing fixes:
** Edit header of saved .tsv-file
** Cosmetic changes\nUnzip before SPM
** Clear variables before loading .mat
** Replace spaces in headers with underscores (which happened automatically apparently upon saving)
** Remove SPM cellstring 

---
## ASL-BIDS-related bug fixes
* #96 Fix order of magnitude in JSON sidecars
* #109 No warning if SliceReadoutTime not provided in DataPar
* #110 xASL_bids_parms2BIDS.m now deals correctly with vectors and SliceReadoutTime = 'shortestTR': created function for this: xASL_quant_SliceReadoutTime_Shortest_TR

----

# ExploreASL v1.2.0

----
## Versions included software
Versions included & used third-party tools (see /External/README_SPM.txt):

SPM12 7219 
CAT12 r1615 
LST 2.0.15

----

## Major feature improvements (still backward compatible)

* Add user flexibility for creating average maps, allowing flipping
* Provide a lesion or ROI mask, to be used not only for cost function masking but also as standard space ROI for ROI-analysis
  This use is now easier, by simply adding Lesion_FLAIR_1.nii or Lesion_T1_2.nii etc, and visualization improved. These masks are now automatically created   (where lesion can be any other ROI):
>> 1. Intralesional
>> 2. Perilesional, pGM+pWM
>> 3. Hemisphere (ipsilateral to lesion)
>> 4. Contralateral (i)
>> 5. Contralateral (ii)
>> 6. Contralateral (iii)
* Option x.S.bMasking added, allowing specifying masking separately for:
>> 1. bSusceptibilityMask
>> 2. bVascularMask
>> 3. subject-wise bGMMask (e.g. the pGM>0.7)
>> 4. brainmasking when loading for lower memory usage
* Affine registration improved & Discrete Cosine Transform (DCT) non-linear registration option added, including an option with partial volume correction built-in for improved DCT-based registration

----

## Bug Fixes

* Allow zipping in Unix-based OS without JavaVirtualMachine 
* Quantification issue with Philips scale slopes
* DCTMK fix, import parameters
* Use _xASL_adm_UnixPath()_ for Unix system calls, for correct path usage (e.g. for spaces that need escaping)
* Double escaping of backslashes in converting .m to .json for DataPar file - subject-regexp
* Compilation path error
* Correctly concatenate numbers when _xASL_num2str_
* Improvements _xASL_adm_LoadParms_ for converting parameters ASL flavors to BIDS/ExploreASL internally

----

## Features

* New startup option for starting ExploreASL, loading data, without processing data
* Shortcut _ExploreASL_ for _ExploreASL_Master_
* Now we have _xASL_csvWrite_, _xASL_csvRead_, _xASL_tsvRead_, _xASL_tsvWrite_
* _xASL_io_Nifti2Im_ now detects erroneously extreme high scaling (potential import issue with Philips RescaleSlope)
  issues a warning and/or tries to fix automatically for FLAIR/T1w images

---
## Work in progress
* Docu Crawler for automatic documentation creation 
* ExploreASL_Import for different ASL flavors
* BIDS implementation import

----

## Documentation
 
* create new prefix for visualization functions (_xASL_vis\_\*_ instead of _xASL_im\_\*_)
* revamp quantification functions for better headers, comments etc

----

# ExploreASL v1.1.3

----

## Bug Fixes

* hotfix minor bug in running the import using DCMTK without the Matlab Image Processing Toolbox #30

# ExploreASL v1.1.2

----

## Bug Fixes

* hotfix minor bug in loading NIfTIs containing lesion masks in CAT12 #28

----

# ExploreASL v1.1.1

----

## Bug Fixes

* hotfix minor bug in creating participants.tsv #23

----


# ExploreASL v1.1.0


----
## Versions included software
Versions included & used third-party tools (see /External/README_SPM.txt):

SPM12 7219 
CAT12 r1615 
LST 2.0.15

----

## Bug Fixes

* Bug fixes and overall code improvements related to the BIDS import workflow (#11)
* Registration with poor CBF contrast will not issue an error anymore but correctly switch to control-T1w registration only (#17)
* Unexisting x.Sequence field fixed, now an appropriate warning is issued and this field is defined automatically by ```xASL_adm_DefineASLSequence.m``` (#16)
----

## Features
* Quantification can now be fully disabled by: ```x.ApplyQuantification = [0 0 0 0 0];``` (#14)
* Insert option to disable M0-ASL registration (#13)
* Update of CAT (Computational Anatomy Toolbox) from version 12.5 to 12.7 (#2)

---
## Work in progress
* Minor improvements of custom scripts for BBB-ASL and BIDS (#8)
* Minor improvements regarding unit testing of ExploreASL (#10)
* Additional warnings for ExploreASL users (#12)
----

## Documentation
 
* Recent changes include the improvement of the documentation within the ExploreASL structure using markdown files and the introduction of a new documentation repository (#7)
* Some function headers were added for increased understandability (#19). These can be viewed in Matlab by: ```help ExploreASL_Master``` where you can replace ExploreASL_Master by the actual function name


----

# ExploreASL v1.0.0

  ----
This is the first release version.

----
## Versions included software
Versions included & used third-party tools (see /External/README_SPM.txt):

SPM12 7219 
CAT12 r1363
LST 2.0.15
