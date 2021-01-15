# Modules

----
## 1. Import Module

----
### ExploreASL\_Import

#### Format

```matlab
ExploreASL_Import(imPar[, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bRunDCM2NII, bClone2Source])
```

#### Description

Import batch T1, FLAIR, DWI, fMRI, M0, ASL data from dicom 2 NIfTI.
Uses dcm2niiX for the conversion, and additionally collects important DICOM header data
and puts them in \_parms.mat and .json sidecars to be used with the ExploreASL pipeline.
This function takes any folder input, but the folder input should be
specified in the imPar definition in the ExploreASL\_ImportConfig.m (later
to be converted to e.g. a CSV file). Follow the steps below, for study "MyStudy" located on "//MyDisk":

1. Make sure you have your DICOM data. Export them from XNAT, download them, or whatsoever
Create a root folder with study ID name, and put the DICOMs in any structure in the raw folder within the study ID root folder
Examples:
imPar.StudyID: MyStudy
StudyRoot folder: //MyDisk/MyStudy
Raw folder containing DICOMs: //MyDisk/MyStudy/raw
2. Make sure that your DICOM data has any structure that can be retrieved
from the folder and/or file names. This function doesn't yet read the DICOM headers
For a quick and dirty (but actually slow) function that converts a
DICOM folder/file structure into readable format, first run
ConvertDicomFolderStructure\_CarefulSlow.m. This will read each DICOM
individually, and put it in a folder with the name identical to the
DICOMs SeriesName/ProtocolName.
3. Once you have all DICOMs in folderstructure with identifyable names
inside //MyDisk/MyStudy/raw, set up the folderstructure in
ExploreASL\_ImportConfig.m. This setup uses the SPM form of regular
expressions, which can be daunting at first, but are very flexible.
Easiest is to study other examples, before creating your own.
For this example, let's say we have //MyDisk/MyStudy/raw/ScanType/SubjectName
because we downloaded our data from XNAT, ordered per ScanType first,
and then per subject.

BRIEF EXPLANATION:
Let's suppose we don't have sessions (only a single structural and functional scan per subject)
The names of our scans comes out of XNAT as '3D\_FLAIR\_eyesClosed', 'T1w\_MPRAGE' and 'PCASL\_10\_min'

and the subject names are 'MyStudy001' .. 'MyStudy002' .. etc.

- imPar.folderHierarchy   - contains a a cell array of regular expressions, with each cell specifying a directory layer/level
the parts within brackets () tell the script that this is a token (i.e. subject, session, ScanType)
Examples:
imPar.folderHierarchy = {'^(3D\_FLAIR|T1w|PCASL).*', '^(Sub-\d{3})$'};
here we say that there are two folder layers '', separated by comma ,
where the names between brackets are used to define what is what.
^ means that the foldername has to start with the following, $ means that the previous has to be the end of the foldername
.* means anything, anylength, \d{3} means three digits
- imPar.tokenOrdering     - defines which tokens are captured by the brackets () in imPar.folderHierarchy: position 1==subject, 2==visit, 3==session, 4==ScanType
Examples:
imPar.tokenOrdering = [2 3 0 1]; stating that subject is the 2nd token, visit is the 3rd token, session has no token (i.e. no session) and ScanType is the 1st token
- imPar.tokenVisitAliases - cell array that defines the aliases for the Visits, i.e. it tells the script which scans are which timepoint/visit.
Similar as explained below for ScanAliases.
First column contains the names that are
recognized in raw DICOM folders for visits,
second column how it is named in NIfTI
structure (should be \_1 \_2 \_3 etc).
Examples:
imPar.tokenVisitAliases = {'Screening','\_1'; 'Month\_12','\_2'; 'Month\_24','\_3'; 'Month\_36','\_4'; 'Month\_48','\_5'};
Note that if you specify tokenVisitAliases, the folders will receive
the indices (e.g. \_1 \_2 \_3), or even \_1 only with a single Visit). If you don't specify
them, they will not get this postfix.
- imPar.tokenScanAliases  - cell array that defines the aliases for the ScanTypes, i.e. it tells the script which scans are which ScanType.
First column should contain regular expression corresponding with the matching criteria in imPar.folderHierarchy
whereas the second column contains the
alias. Following valid aliases exist:
'T1' 'FLAIR' 'ASL4D' 'M0' 'ASL4D\_RevPE' 'func' 'func\_NormPE' 'func\_RevPE' 'dwi' 'dwi\_RevPE' 'DSC4D'
Examples:
imPar.tokenScanAliases = {'^3D\_FLAIR$', 'FLAIR'; '^T1w$', 'T1'; '^PCASL$', 'ASL4D'};
- imPar.tokenSessionAliases-same as tokenScanAliases but for sessions
Examples:
imPar.tokenSessionAliases = {}; as we don't have sessions
- imPar.bMatchDirectories - true if the last layer is a folder, false if the last layer is a filename (as e.g. with PAR/REC, enhanced DICOMs)


----
## 2. Structural Module

----
### xASL\_module\_Structural

#### Format

```matlab
[result, x] = xASL_module_Structural(x)
```

#### Description
This first ExploreASL module processes the structural
images, i.e. high-resolution T1w and FLAIR (if present), on an individual (i.e. subject-to-subject) basis.
If a FLAIR is present, this is processed first to obtain a WMH mask to fill the hypointense lesions on the T1w,
before segmenting the T1w. For the T1w segmentation this module uses CAT12
by default but if this fails it falls back to SPM after trying to
optimize the image contrast. This module has the following steps/submodules/wrappers:

- 010\_LinearReg\_T1w2MNI         - Ensure the alignment of subjects' anterior commissure (AC) with the AC in MNI & apply this to all images
- 020\_LinearReg\_FLAIR2T1w       - Align the FLAIR (if present) with T1w
- 030\_FLAIR\_BiasfieldCorrection - Perform a biasfield correction (if not performed  by LST in following steps)
- 040\_LST\_Segment\_FLAIR\_WMH     - Segment WMH lesions on FLAIR (if present)
- 050\_LST\_T1w\_LesionFilling\_WMH - Use WMH segmentation to fill lesions on T1w
- 060\_Segment\_T1w               - Tissue segmentation on T1w
- 070\_CleanUpWMH\_SEGM           - Extra WMH cleanup of some over- and under-segmentation
- 080\_Resample2StandardSpace    - Clone all images to standard space
- 090\_GetVolumetrics            - Obtain whole-brain volumes of GM, WM, CSF, WMH
- 100\_VisualQC                  - Obtain QC parameters & save QC Figures
- 110\_DoWADQCDC                 - QC for WAD-QC DICOM server (OPTIONAL)


----
## 3. ASL Module

----
### xASL\_module\_ASL

#### Format

```matlab
[result, x] = xASL_module_ASL(x)
```

#### Description
This ExploreASL module processes the ASL
images, i.e. ASL4D, M0, etc (if present), on an individual (i.e. session-to-session, where session indicates BIDS "run") basis.
Both 2D and 3D options are automatially chosen, as well as processing of time-series (if available), such as motion correction and outlier
exclusion. This module has the following submodules/wrappers:

- 010\_TopUpASL          - FSL TopUp geometric distortion correction (if M0 images with reversed phase-encoding are present)
- 020\_RealignASL        - If time-series are present, motion correction and outlier exclusion (ENABLE)
- 030\_RegisterASL       - Registration of ASL to T1w anatomical images (if lacking, to MNI images)
- 040\_ResliceASL        - Resample ASL images to standard space
- 050\_PreparePV         - Create partial volume images in ASL space with ASL resolution
- 060\_ProcessM0         - M0 image processing
- 070\_Quantification    - CBF quantification
- 080\_CreateAnalysisMask- Create mask using FoV, vascular outliers & susceptibility atlas
- 090\_VisualQC\_ASL      - Generate QC parameters & images
- 100\_WADQC             - QC for WAD-QC DICOM server (OPTIONAL)


----
## 4. Population Module

----
### xASL\_module\_Population

#### Format

```matlab
[result, x] = xASL_module_Population(x)
```

#### Description
This ExploreASL module processes all available images on the
group level. It assumes that all images were adequately processed in the
previous modules. It will perform the following group-wise processing and
checks:

- 010\_CreatePopulationTemplates - Create population average images, to compare scanners, cohorts etc without physiological variance
- 020\_CreateAnalysisMask        - Generate a group-level mask by combining individuals masks, for ROI-based analysis & VBA
- 030\_CreateBiasfield           - When there are multiple scanners, create scanner-specific biasfields (uses Site.mat for this)
- 040\_GetDICOMStatistics        - Create TSV file with overview of DICOM parameters
- 050\_GetVolumeStatistics       - Create TSV file with overview of volumetric parameters
- 060\_GetMotionStatistics       - Create TSV file with overview of motion parameters
- 065\_GetRegistrationStatistics - Create TSV file with overview of the registration statistics
- 070\_GetROIstatistics          - Create TSV file with overview of regional values (e.g. qCBF, mean control, pGM etc)
- 080\_SortBySpatialCoV          - Sort ASL\_Check QC images by their spatial CoV in quality bins
- 090\_DeleteAndZip              - Delete temporary files and gzip all NIfTIs


