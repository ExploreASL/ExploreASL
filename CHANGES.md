
# ExploreASL Change Log

---
## ExploreASL v1.7.0

### Versions included software

Versions included & used third-party tools (see `/External/README_SPM.txt`):

* SPM12 7219 
* CAT12 r1615 
* LST 2.0.15 

### Feature improvements

* Issue #455: Automatically compare results of TestDataSets with a saved reference results
* Issues #480, #623, #649, #661: Restructure `x`: `x.opts` (input arguments and their derivatives), `x.dir` (directories), `x.settings` (mostly booleans for pipeline settings), `x.dataset` (dataset related fields), `x.external`, ...
* Issue #572: Restructure JSON handling during NiFTI to BIDS import
* Issue #580: Add parsing of Gold Standard Phantoms **ASL-DRO**
* Issues #588, #612: `ExploreASL` reads folders and automatically searches for **sourceStructure**, **studyPar** and **dataPar** JSON files
* Issue #600: Put participants.tsv to the derivatives folder during import to legacy ExploreASL format
* Issue #602: Remove option for cloning the NIfTI output after import as BIDS directory is a read-only archive
* Issue #603: Give ExploreASL version in JSON files after BIDS to Legacy conversion
* Issue #631: Remove repeated warnings
* Issue #632: Add comparison script for untouched NIfTI comparison
* Issue #643: `bids.layout`: avoid printing the same warning repetitively in case multiple scans in a data set have the same issue
* Issue #656: Improve warnings (data loading)


### Bug fixes

* Issue #583: Proper testing of flavors using ExploreASL_Master
* Issue #584: Print the subject name depending on the existence of its definition in `x.SUBJECT` to avoid crashes for error reporting in the population module
* Issue #586: Avoid crashing `xASL_adm_GetPopulationSessions` if no sessions are found
* Issue #591: **MultiTE** import puts TE before PLD in the time series and corrects the JSON output
* Issue #618: Add session name to all M0Check and ASLCheck QC files in the Population folder
* Issue #620: `xASL_adm_GzipAllFiles`: Allow spaces in an input path for macOS/Linux
* Issue #625: Fix bug related to session format
* Issue #627: Remove a BIDS fiels and BIDS2Legacy should crash and show you why it crashed
* Issue #628: Fix parsing sessions and runs for converting rawdata to derivatives
* Issue #630: Move creation of population folder
* Issue #646: Improve BIDS warnings
* Issue #652: `xASL_vis_CreateVisualFig`: allow empty overlays
* Issue #659: `xASL_stat_PrintStats`: Fix visits bug (legacy format)
* Issue #655: `xASL_adm_GetPopulationSessions` gave incorrect warnings
* Issue #666: Warning when multiple `dataPar*.json` or `studyPar*.json` or `sourcestructure*.json` are present
* Issue #670: Fix warnings and behavior of `ExploreASL_Initialize`

### Other improvements

* Issue #465: Add projects to acknowledgments
* Issue #615: Add change log to documentation
* Issue #637: Restyle ExploreASL change log

----
## ExploreASL v1.6.2


### Bug Fixes

* Issue #589 Fix scaling issues in JSONs in TestDataSet/derivatives

---
## ExploreASL v1.6.1

### Bug Fixes

* Issue #578 Fix incorrect path searching by providing error if no .json-file is inputted

---
## ExploreASL v1.6.0

### Versions included software

Versions included & used third-party tools (see `/External/README_SPM.txt`):

* SPM12 7219 
* CAT12 r1615 
* LST 2.0.15 

### Feature improvements (still backward compatible)

* Issues #62, #266, #449: JSON i/o unify and use `spm_jsonread` and `spm_jsonwrite` for all operations
* Issue #349: Improve screenprint of the current subjects/sessions/modules by xASL_Iteration
* Issue #384: Add method to export image matrices as structured points in VTK format
* Issue #523: Update ExploreASL_ImportConfig to work with JSONs and add alert for users
* Issue #538: Function to user-replace label values in atlas NIfTI, written for the CICERO study
* Issue #559: Add option `x.MakeNIfTI4DICOM` to create CBF optimized for DICOM creation/PACS export

### Work in progress

#### BASIL
* Issue #20: Data pre-processing prepared for BASIL
* Issue #391: Add single-PLD model for BASIL

#### ASL-BIDS

* Issues #290,#483,#484: Initial version of the ASL-BIDS import workflow
* Issue #353: Correct conversion BIDS->Legacy for M0 with reversed PE direction
* Issues #394,#514,#545: Improve modularity of the ASL-BIDS import module
* Issue #411: Delete temporary folders in DICOM->BIDS conversion
* Issue #421: Use `ImageType` DICOM field to detect scan order in GE in DICOM->BIDS import
* Issue #426: Reading PLD from DICOM for GE in import to BIDS
* Issue #479: ASL-BIDS import for Hadamard encoded FME sequences

#### DRO and QASPER
* Issues #361,#443: Import and import test of DRO
* Issues #467: Improve script to generate ASL-BIDS version of ASL DRO v2.2.0

#### Compilation/stand-alone version

### Bug Fixes

* Issue #184: Skip PVC in Population statistics, when this does not make sense for a given ROI
* Issue #262: Improve GZIP on windows
* Issue #341: Reduce extensive usage of ExploreASL CLI progress tracker
* Issue #368: Fix problem with writing trailing zeros in real numbers in spm_jsonwrite
* Issue #387: Remove graphical waitbar interface from xASL_im_ResampleLinearFair
* Issue #415: xASL_Copy: Allow recursive copying of directories
* Issue #418: Verifies the SliceTiming parameter if timing difference is consistent between slices
* Issue #424,#454: Remove `string` and `contains` functions to ensure Matlab compatibility for R2016
* Issue #430: Fix `xASL_adm_ReplaceSymbols` crash when trying to replace sub-strings in PhoenixProtocol field
* Issues #433,#474,#542: Splitting of ASL and M0 - fix on rerun, split metadata, backup aslcontext.tsv
* Issue #451: Clean unused code and cleaned the development directory
* Issue #466: Warning if the equal sign is used in JSON files instead of colon
* Issue #475: Fix error with studies that have a special character in their name
* Issue #477: Change vascular territory atlases file extensions to .nii.gz
* Issue #496: Fixed session numbers in population module
* Issue #502: Fix smoothing of 4D NIfTIs
* Issue #505: Allow Token Aliases in import to be row vectors in import
* Issue #518: ExploreASL logo and verbose overview shown only once in a single pipeline run
* Issue #520: Stop import workflow from creating Population folder in root directory
* Issue #543: Minor fix of xASL_num2str behavior
* Issue #563: Minor design fix for ExploreASL dataset initialization

### Documentation
* Issues #403,#423,#457: Improved inline comments and headers
* Issue #452: Provide descriptions of available atlas options
* Issue #463: Remove remaining markdown file to a separate Documentation repository
* Issues #486,#489: Introduce document templates for GitHub QC Workflow
* Issue #499: Make heading in documentation work in dark mode as well
* Issue #515: Added tutorials to documentation
* Issue #536: Create a first version of the QC walkthrough document in markdown

### Testing
* Issue #156: Make internal error messages more specific by providing subject information
* Issue #352: Improve parsing of errors and warnings from log files
* Issue #395,#416: Improve testing of BIDS import
* Issue #517: Improve unit testing scripts
* Issue #529: TestDataSet is in BIDS derivatives format
* Issue #570: Release testing and minor documentation improvements


----
## ExploreASL v1.5.1

### Bug Fixes

* Issue #439 Fix population module error by correctly renaming `MNI_Structural.*` files

----
## ExploreASL v1.5.0

### Versions included software

Versions included & used third-party tools (see `/External/README_SPM.txt`):

* SPM12 7219 
* CAT12 r1615 
* LST 2.0.15 

### Feature improvements (still backward compatible)

* Issue #39: Create PV-corrected GM & WM CBF maps in native space
* Issue #56,#410: Add Mindboggle atlas to ExploreASL and restructure general atlas access in population module
* Issue #283: `xASL_stat_GetROIstatistics` provide more feedback on missing images
* Issue #299: Move `CustomScripts` with study-specific scripts to a separate repository
* Issue #302: Remove server calls in CAT12 functions
* Issue #313: Move GUI to a separate repository
* Issue #351: T2 and T1c files are now also aligned to the T1w and outputted to standard space
* Issue #354: Added an option `x.DummyScanPositionInASL4D` that removes marked dummy scans when splitting ASL to ASL+M0+dummy
* Issues #356,#396,#397: Internally restructure SliceTime allowing ExploreASL now to work with multi-band 2D EPI as well or any other SliceTime order

### Work in progress

#### ASL-BIDS

* Issues #163,#189,#357,#373: Conversion from DICOM to BIDS
* Issues #334,#382: Import of PAR-REC to BIDS
* Issue #343: Add separate M0 option to mTrial import

#### Compilation/stand-alone version
* Issue #335: All input arguments can be passed in the deployed mode
* Issue #380: Enable advanced input parsing for xASL compiled

### Bug Fixes

* Issue #228: Fix CAT12 warnings with non-existent field cm
* Issue #272: Fix errors in JSON import of ASL sidecars
* Issues #273,#285,#291,#363: Minor fixes in input parameter administration
* Issues #276,#280,#288,#329,#400: Fix error in reading TSV files with unclear number of columns
* Issue #282: Population module is run serially in otherwise parallel mode
* Issue #292: `xASL_qc_SortBySpatialCoV` now use all subjects without skipping
* Issue #305: `xASL_adm_UnixPath`: bug with Windows+WSL
* Issues #306,#378: Fix JSON reading with 'i' interpreted as a complex number
* Issue #309: Fix non-linear registration of T1w to standard space with too high T1w values
* Issue #312: `xASL_stat_GetROIstatistics` fix skipping of actual ROI extraction
* Issue #325: `xASL_adm_CleanUpBeforeRerun` delete files correctly
* Issue #339: Fix JSON reading of special characters
* Issue #399: Fix special characters in Windows filenames
* Issue #405: Fix range-check error in Background Suppression timing calculation
* Issue #406: Fix `xASL_stat_MedianNan` for all-NaN input
* Issues #408,#409: Skip missing fields in CAT during reports in compiled ExploreASL

### Documentation
* Issue #7: Create README files in subfolders and added to interactive documentation 
* Issues #279, #345: Move documentation to a separate repository
* Issues #300,#318: Improve ExploreASL tutorial
* Issue #355: Documentation improvements regarding input parameters

### Testing
* Issues #193,#350,#398: Testing DICOM to BIDS conversion against a reference
* Issue #294: Implement initial unit testing framework
* Issue #326: Parse warnings/errors from all log in all subdirectories
* Issue #369: Unit testing of `xASL_test_getLogContent`
* Issues #371,#376,#404: Testing script for the DICOM->BIDS->Legacy conversion and processing

----
## ExploreASL v1.4.0

### Versions included software
Versions included & used third-party tools (see /External/README_SPM.txt):

* SPM12 7219 
* CAT12 r1615 
* LST 2.0.15 

### Feature improvements (still backward compatible)

* Issue #20 Implement BASIL -> changed order quantification/masking, and always save resampled PWI4D for quantification, facilitating BASIL
* Issue #35 Calculation background suppression efficiency for pseudo-M0 -> in case of missing separate M0 images but still using background suppressed mean control as pseudo-M0, a single correction value (3D) or slice-wise correction value (2D) are applied to the pseudo-M0 image/slices
* Issue #131 ExploreASL_GUI beta-testing enhancements set 1 -> Aesthetic improvements in certain modules, fixed incorrect removal of Philips-related json-sidecar fields in DCM2BIDS / Import module. Correction of ASL image flickering bug and ability for the user to subset without having to reload the data. Added ability to clarify data type of variables without the need to reload data. Added auto-select / auto-complete functionality in the ParmsMaker module as soon as the user indicates an analysis directory (i.e auto-completion of `SliceReadoutTime` for pCASL)
* Issue #166 Update ADNI import in ExploreASL_Import
* Issue #179 Finish creation average maps for CICERO -> improvements xASL_wrp_CreatePopulationTemplates
* Issue #190 Hammers atlas option CAT12 -> restored the original Hammers atlas option for a colleague (note it's license though!)
* Issue #225 Create DataPar option for running SPM12 longitudinal registration
* Issue #241 Add warning when loading data without x output structure 
* Issue #243 Shorten SPM initialization time -> removed configuration loading of unused toolboxes (check if the SPM DICOM import module is still needed)

### Work in progress

#### ASL-BIDS

* Issue #193 Comparing BIDS folders (for testing purpose)
* Issue #226 Add new DICOM tags to DCMTK import

### Bug Fixes

* Issue #191 WMH warning when no FLAIR analyzed -> this warning is now removed if no FLAIR was present in the data
* Issue #248 Temporary fix native space processing -> in Population module
* Issue #252 Population modules analysis masks - minor errors
* Issue #267 Error in reading JSONs from EPAD

### Documentation
* Issue #196 All contributors -> all contributors are now automatically added to the main README.md
* Issue #217 Documentation improvements
* Issue #219 Add user to documentation

----
## ExploreASL v1.3.0

### Versions included software
Versions included & used third-party tools (see /External/README_SPM.txt):

* SPM12 7219 
* CAT12 r1615 
* LST 2.0.15 

### Feature improvements (still backward compatible)

* Issue #59: Assign "weights" to status files, allowing the external Python ExploreASL GUI to provide a better estimate of the progress
* Issue #71 Remove custom lesions `(Lesion\_(FLAIR|T1)\_\\d\\.nii)` from `WMH_SEGM.nii`
* Issue #76: Improve functionality of ExploreASL through Command Line Interface (CLI), i.e. without using the Graphical User Interface (GUI) of Matlab
* Issue #123 Create status files also for skipped processing parts: this is mainly the case for running the structural module without a FLAIR scan. Having all status files helps third-party tools such as the Python ExploreASL GUI to know that processing has succesfully completed (duplicate issues #137 and #129)
* The same was done for the ASL realignment status file, in case realignment is skipped for a 3D scan
* Issue #145: Improve `.nii(.gz)` management in `xASL_spm_deformations`: allow either .nii or .nii.gz as input, treat them equally, and when `.nii.gz` is provided as output path, zip the resulting deformed image


### Work in progress

* Issue #32 Docker integration
* Issue #55 ExploreASL GUI, written in Python
* Issue #106 `xASL_im_SplitImageLabels`: Allow splitting labels, and warping them to standard space.
This is part of a continuous development on creating average flow territory templates and figures.
* Issue #162 Remove `bNativeSpaceProcessing` from TestDataSet for now, return this when `bNativeSpaceProcess` is made more modular

#### ASL-BIDS

* Issue #82 Avoid 4D files with nT==1, which is not allowed in the BIDS validator

#### Compilation/stand-alone version

* Issue #88 `xASL_SysMove` error in Windows when a path includes whitespaces ' '

### Bug Fixes

* Issue #85 Improvement ApplyQuantification
* Issue #99 Improve loading of metadata (`xASL_str2num` & `xASL_init_LoadMetaData`)
* Issue #105 In case of missing data, fill `x.S.SUBJECTID` and `S.DAT` data for the last subject/session
* Issue #119 Create status file for last subject, if it's processing is skipped
* Issue #120 `xASL_im_CreateAnalysisMask` in native space mode when in parallel execution
* Issue #138 Fix structure TestDataJSON & its JSON files
* Issue #139 Ensure that x fields are not case sensitive, by e.g. replacing `strcmp` by `strcmpi`
* Issue #141 Solve conflicts between develop and master, these were minor edits that weren't implemented in develop yet
* Issue #143 ensure that VBA masks are also created for a 3D spiral sequence (this was not created yet as the susceptibility masks were missing for this sequence)
* Issue #148 Syncing the ROI/lesion processing of T1w & FLAIR
* Issue #151 Minor bugfixes for TopUp
* Issue #177 Fix `xASL_adm_UnzipNifti` & `xASL_io_SaveNifti` when path is incomplete
    * Skip warnings for small populations (in the start of `xASL_module_Population`)
    * Fix regular expression in `xASL_init_LoadMetaData`
    * `iRow` counting fix in `xASL_bids_Add2ParticipantsTSV`


### Documentation

* Issue #7 Documentation/revamp `xASL_im_ClipExtremes`
* Issue #159 Ensure that all sequence notations use underscores instead of whitespace, e.g. `3D_spiral` instead of `3D spiral`


### Testing

* Issue #86 `xASL_qc_TestExploreASL`: improve Table creation
* Issue #112 Save Tanimoto Coefficient (i.e. a fuzzy overlap/Dice score) of the final ASL-T1w registration
* Issue #128 Improved one internal test dataset
* Issue #130 `xASL_qc_TestExploreASL`: Complete functionality by allowing Windows parallelization & testing the Windows ExploreASL compilation. Also added unit testing framework in the same issue.


----
## ExploreASL v1.2.2


### Bug Fixes

* Issue #119 xASL_wrp_LinearReg_T1w2MNI.m: ROI .nii files correctly aligned with T1

----
## ExploreASL v1.2.1

### Bug Fixes

* xASL_qc_TestExploreASL: Remove locked folders if rerun
* Issue #90 `xASL_io_Nifti2Im`: manage the detection of odd scaling
* Issue #93 `xASL_SysMove`: diz illegal symbols Windows
* Issue #102 `xASL_wrp_CreatePopulationTemplates`: minor bugfix
* Issue #104 Fix creation ResultsTable for TestCases (also #86)
* Issue #114 Acquiring Matlab version doesn't crash anymore in deployed mode
* Issue #115 Fixing NaNs problem in M0 mask computation
* Issue #118 `xASL_adm_GzipAllFiles` doesn't crash anymore in Windows

* Issue #116 ExploreASL testing fixes:
    * Edit header of saved .tsv-file
    * Cosmetic changes\nUnzip before SPM
    * Clear variables before loading .mat
    * Replace spaces in headers with underscores (which happened automatically apparently upon saving)
    * Remove SPM cellstring 

### ASL-BIDS-related bug fixes

* Issue #96 Fix order of magnitude in JSON sidecars
* Issue #109 No warning if SliceReadoutTime not provided in DataPar
* Issue #110 `xASL_bids_parms2BIDS.m` now deals correctly with vectors and `SliceReadoutTime = 'shortestTR'`: created function for this: `xASL_quant_SliceReadoutTime_Shortest_TR`

----
## ExploreASL v1.2.0

### Versions included software

Versions included & used third-party tools (see /External/README_SPM.txt):

* SPM12 7219 
* CAT12 r1615 
* LST 2.0.15

### Major feature improvements (still backward compatible)

* Add user flexibility for creating average maps, allowing flipping
* Provide a lesion or ROI mask, to be used not only for cost function masking but also as standard space ROI for ROI-analysis
  This use is now easier, by simply adding `Lesion_FLAIR_1.nii` or `Lesion_T1_2.nii` etc, and visualization improved. These masks are now automatically created   (where lesion can be any other ROI):

1. Intralesional
2. Perilesional, pGM+pWM
3. Hemisphere (ipsilateral to lesion)
4. Contralateral (i)
5. Contralateral (ii)
6. Contralateral (iii)

* Option x.S.bMasking added, allowing specifying masking separately for:

1. bSusceptibilityMask
2. bVascularMask
3. subject-wise bGMMask (e.g. the `pGM>0.7`)
4. brainmasking when loading for lower memory usage

* Affine registration improved & Discrete Cosine Transform (DCT) non-linear registration option added, including an option with partial volume correction built-in for improved DCT-based registration


### Bug Fixes

* Allow zipping in Unix-based OS without JavaVirtualMachine 
* Quantification issue with Philips scale slopes
* DCTMK fix, import parameters
* Use `xASL_adm_UnixPath()` for Unix system calls, for correct path usage (e.g. for spaces that need escaping)
* Double escaping of backslashes in converting `.m` to `.json` for DataPar file - subject-regexp
* Compilation path error
* Correctly concatenate numbers when _xASL_num2str_
* Improvements `xASL_adm_LoadParms` for converting parameters ASL flavors to BIDS/ExploreASL internally

### Features

* New startup option for starting ExploreASL, loading data, without processing data
* Shortcut `ExploreASL` for `ExploreASL_Master`
* Now we have `xASL_csvWrite, xASL_csvRead, xASL_tsvRead, xASL_tsvWrite`
* `xASL_io_Nifti2Im` now detects erroneously extreme high scaling (potential import issue with Philips RescaleSlope)
  issues a warning and/or tries to fix automatically for FLAIR/T1w images


### Work in progress
* Docu Crawler for automatic documentation creation 
* `ExploreASL_Import` for different ASL flavors
* BIDS implementation import


### Documentation
 
* create new prefix for visualization functions (`xASL_vis\_\*` instead of `xASL_im\_\*`)
* revamp quantification functions for better headers, comments etc

----
## ExploreASL v1.1.3

### Bug Fixes

* hotfix minor bug in running the import using DCMTK without the Matlab Image Processing Toolbox #30

----
## ExploreASL v1.1.2

### Bug Fixes

* hotfix minor bug in loading NIfTIs containing lesion masks in CAT12 #28

----
## ExploreASL v1.1.1

----

### Bug Fixes

* hotfix minor bug in creating participants.tsv #23

----
## ExploreASL v1.1.0

### Versions included software
Versions included & used third-party tools (see `/External/README_SPM.txt`):

* SPM12 7219 
* CAT12 r1615 
* LST 2.0.15

### Bug Fixes

* Bug fixes and overall code improvements related to the BIDS import workflow (#11)
* Registration with poor CBF contrast will not issue an error anymore but correctly switch to control-T1w registration only (#17)
* Unexisting x.Sequence field fixed, now an appropriate warning is issued and this field is defined automatically by `xASL_adm_DefineASLSequence.m` (#16)


### Features

* Quantification can now be fully disabled by: ```x.ApplyQuantification = [0 0 0 0 0];``` (#14)
* Insert option to disable M0-ASL registration (#13)
* Update of CAT (Computational Anatomy Toolbox) from version 12.5 to 12.7 (#2)


### Work in progress

* Minor improvements of custom scripts for BBB-ASL and BIDS (#8)
* Minor improvements regarding unit testing of ExploreASL (#10)
* Additional warnings for ExploreASL users (#12)

### Documentation
 
* Recent changes include the improvement of the documentation within the ExploreASL structure using markdown files and the introduction of a new documentation repository (#7)
* Some function headers were added for increased understandability (#19). These can be viewed in Matlab by: `help ExploreASL_Master` where you can replace ExploreASL_Master by the actual function name

----
## ExploreASL v1.0.0

This is the first release version.

### Versions included software
Versions included & used third-party tools (see /External/README_SPM.txt):

* SPM12 7219 
* CAT12 r1363
* LST 2.0.15


