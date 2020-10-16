# ExploreASL v1.3.0

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
