function [x] = xASL_module_Import(studyPath, imParPath, studyParPath, bRunSubmodules, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source, x)
%xASL_module_Import Imports the DICOM or PAR/REC source data to NIFTIs in ASL-BIDS format
%
% FORMAT: xASL_module_Import(studyPath, imParPath, studyParPath, bRunSubmodules, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source, x)
%
% INPUT:
%   studyPath           - path to the study directory containing the 'sourcedata' directory with the DICOM files (REQUIRED)
%   imParPath           - path to the JSON file with structure with import parameters, output of ExploreASL_ImportConfig.m (originally)
%                         All other input parameters are configured within this function. (OPTIONAL)
%                         The path is optional, but the file has to be there. Either provided as a full path or a filename in the path,
%                         or default names (case-insensitive) sourceStructure.json, ImagePar.json are seeked
%   studyParPath        - path to the JSON file with the BIDS parameters relevant for the whole study. These parameters are used
%                         if they cannot be extracted from the DICOMs automatically. (OPTIONAL)
%                         Looking automatically for file studyPar.json
%   bRunSubmodules      - Specify which of the parts should be run (OPTIONAL, DEFAULT [1 1 0])
%                         [1 0 0] - Run the DICOM to NIFTI conversion
%                         [0 1 0] - Run the NIFTI transformation to the proper ASL-BIDS
%                         [0 0 1] - Run the defacing and full anonymization
%   bCopySingleDicoms   - if true, copies a single DICOM with each NIfTI
%                         dataset/ScanType, that can be used to retrieve missing parameters from
%                         the DICOM header, or as dummy DICOM to dump embed data into (e.g. WAD-QC) (DEFAULT=false)
%   bUseDCMTK           - if true, then use DCMTK, otherwise use DICOMINFO from Matlab (DEFAULT=false)
%   bCheckPermissions   - if true, check whether data permissions are set correctly, before trying to read/copy the files (DEFAULT=false)
%   Clone2Source        - if true, then makes a copy of everything it converted to NIfTI.
%                         Can be useful to have a separate source BIDS structure to store
%                         all source NIfTIs, and to keep the derivatives in the
%                         analysisfolder (OPTIONAL, DEFAULT=false)
%   x                   - if x is provided, initialization of ExploreASL is skipped
%
%
% OUTPUT: n/a
%
% OUTPUT FILES:
%   //AnalysisDir/dcm2niiCatchedErrors.(mat|json) - overview of catched dcm2nii errors, or other errors in this function
%   //AnalysisDir/import_log_StudyID_yyyymmdd_hhmmss.txt - diary log of this function
%   //AnalysisDir/import_summary.csv - hence the name
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
% Import batch `T1`, `T2`, `FLAIR`, `DWI`, `fMRI`, `M0`, `ASL` data from dicom 2 NIfTI in ASL-BIDS format and structure.
% Uses dcm2niiX for the conversion, and additionally collects important DICOM header data
% and puts them in `.json` sidecars to be used with the ExploreASL pipeline.
% This function takes any folder input, but the folder input should be
% specified in the imPar definition. Follow the steps below, for study `"MyStudy"` located on `"//MyDisk"`:
%
% 1. Make sure you have your DICOM data. Export them from XNAT, download them, or whatsoever
%    Create a root folder with study ID name, and put the DICOMs in any structure in the sourcedata folder within the study ID root folder
%    Examples:
%    imPar.StudyID: MyStudy
%    StudyRoot folder: `//MyDisk/MyStudy`
%    sourcedata folder containing DICOMs: `//MyDisk/MyStudy/sourcedata`
% 2. Make sure that your DICOM data has any structure that can be retrieved
%    from the folder and/or file names. This function doesn't yet read the DICOM headers
%    For a quick and dirty (but actually slow) function that converts a
%    DICOM folder/file structure into readable format, first run
%    ConvertDicomFolderStructure_CarefulSlow.m. This will read each DICOM
%    individually, and put it in a folder with the name identical to the
%    DICOMs SeriesName/ProtocolName.
% 3. Once you have all DICOMs in folderstructure with identifyable names
%    inside `//MyDisk/MyStudy/sourcedata`, set up the folderstructure in
%    ExploreASL_ImportConfig.m. This setup uses the SPM form of regular
%    expressions, which can be daunting at first, but are very flexible.
%    Easiest is to study other examples, before creating your own.
%    For this example, let's say we have `//MyDisk/MyStudy/sourcedata/ScanType/SubjectName`
%    because we downloaded our data from XNAT, ordered per ScanType first,
%    and then per subject.
%
%    BRIEF EXPLANATION:
%    Let's suppose we don't have sessions (only a single structural and functional scan per subject)
%    The names of our scans comes out of XNAT as `'3D_FLAIR_eyesClosed'`, `'T1w_MPRAGE'` and `'PCASL_10_min'`
%    and the subject names are `'MyStudy001'` .. `'MyStudy002'` .. etc.
%
%    - imPar.folderHierarchy     - contains a a cell array of regular expressions, with each cell specifying a directory layer/level
%                                  the parts within brackets `()` tell the script that this is a token (i.e. subject, session, ScanType)
%                                  Examples:
%                                  `imPar.folderHierarchy = {'^(3D_FLAIR|T1w|PCASL).*', '^(Sub-\d{3})$'};`
%                                  here we say that there are two folder layers '', separated by comma ,
%                                  where the names between brackets are used to define what is what.
%                                  `^` means that the foldername has to start with the following, $ means that the previous has to be the end of the foldername
%                                  `.*` means anything, anylength, `\d{3}` means three digits
%    - imPar.tokenOrdering       - defines which tokens are captured by the brackets `()` in imPar.folderHierarchy: position `1==subject`, `2==visit`, `3==session`, `4==ScanType`
%                                  Examples:
%                                  `imPar.tokenOrdering = [2 3 0 1];` stating that subject is the 2nd token, visit is the 3rd token, session has no token (i.e. no session) and ScanType is the 1st token
%    - imPar.tokenVisitAliases   - cell array that defines the aliases for the Visits, i.e. it tells the script which scans are which timepoint/visit.
%                                  Similar as explained below for ScanAliases.
%                                  First column contains the names that are
%                                  recognized in sourcedata DICOM folders for visits,
%                                  second column how it is named in NIfTI
%                                  structure (should be _1 _2 _3 etc).
%                                  Examples:
%                                  `imPar.tokenVisitAliases = {'Screening','_1'; 'Month_12','_2'; 'Month_24','_3'; 'Month_36','_4'; 'Month_48','_5'};`
%                                  Note that if you specify tokenVisitAliases, the folders will receive
%                                  the indices (e.g. `_1 _2 _3`), or even `_1` only with a single Visit). If you don't specify
%                                  them, they will not get this postfix.
%    - imPar.tokenScanAliases    - cell array that defines the aliases for the ScanTypes, i.e. it tells the script which scans are which ScanType.
%                                  First column should contain regular expression corresponding with the matching criteria in imPar.folderHierarchy
%                                  whereas the second column contains the
%                                  alias. Following valid aliases exist:
%                                  `'T1'` `'FLAIR'` `'ASL4D'` `'M0'` `'ASL4D_RevPE'` `'func'` `'func_NormPE'` `'func_RevPE'` `'dwi'` `'dwi_RevPE'` `'DSC4D'`
%                                  Examples:
%                                  `imPar.tokenScanAliases = {'^3D_FLAIR$', 'FLAIR'; '^T1w$', 'T1'; '^PCASL$', 'ASL4D'};`
%    - imPar.tokenSessionAliases - same as tokenScanAliases but for sessions
%                                  Examples:
%                                  `imPar.tokenSessionAliases = {}; % as we don't have sessions`
%    - imPar.bMatchDirectories   - true if the last layer is a folder, false if the last layer is a filename (as e.g. with PAR/REC, enhanced DICOMs)
%
% EXAMPLE: xASL_module_Import('//MyDisk/MyStudy');
%          xASL_module_Import('//MyDisk/MyStudy','sourceStructure.json','studyHiQ.json');
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% 1. Initialize the parameters

    % First do the basic parameter admin and initialize the default values
    if nargin < 1 ||  isempty(studyPath)
        error('The studyPath needs to be defined');
    end

    if strcmp(studyPath(end),'\') || strcmp(studyPath(end),'/')
        studyPath = studyPath(1:end-1); % bugfix
    end

    % Check the imagePar input file
    if nargin < 2 || isempty(imParPath)
        % If the path is empty, then try to find sourceStructure.json or sourcestruct.json
        fListImPar = xASL_adm_GetFileList(studyPath,'^(source|Source)(Structure|Struct|structure|struct).json$', 'List', [], 0);
        if length(fListImPar) < 1
            error('Could not find the sourceStructure.json file');
        end
        imParPath = fullfile(studyPath,fListImPar{1});
    else
        [fpath, ~, ~] = fileparts(imParPath);
        if isempty(fpath)
            imParPath = fullfile(studyPath,imParPath);
        end
    end

    % Find the studyPar input file
    if nargin < 3 || isempty(studyParPath)
        % If the path is empty, then try to find studyPar.json
        fListStudyPar = xASL_adm_GetFileList(studyPath,'^(study|Study)(Par|par).json$', 'List', [], 0);
        if length(fListStudyPar) < 1
            warning('Could not find the StudyPar.json file');
        else
            studyParPath = fullfile(studyPath,fListStudyPar{1});
        end
    else
        [fpath, ~, ~] = fileparts(studyParPath);
        if isempty(fpath)
            studyParPath = fullfile(studyPath,studyParPath);
        end
    end

    if nargin<4 || isempty(bRunSubmodules)
        bRunSubmodules = [1 1 0];
    else
        if length(bRunSubmodules) ~= 3
            error('bRunSubmodules must have length 3');
        end
    end

    % By default don't copy DICOMs for anonymization reasons
    if nargin<5 || isempty(bCopySingleDicoms)
        bCopySingleDicoms = false;
    end

    if nargin<6 || isempty(bUseDCMTK)
        % Default set to using DCM-TK
        bUseDCMTK = true;
    elseif ~bUseDCMTK && isempty(which('dicomdict'))
        error('Dicomdict missing, image processing probably not installed, try DCMTK instead');
    end

    if nargin<7 || isempty(bCheckPermissions)
        if isunix
            bCheckPermissions = true;
        else
            bCheckPermissions = false;
        end
    end

    if nargin<8 || isempty(bClone2Source)
        bClone2Source = false;
    end

    % Only initialize ExploreASL if this wasn't initialized before
    if nargin<9 || isempty(x)
        x = ExploreASL_Initialize;
    end
    
    % Write import settings to modules field of x structure
    x.modules.import.settings.bCopySingleDicoms = bCopySingleDicoms;
    x.modules.import.settings.bUseDCMTK = bUseDCMTK;
    x.modules.import.settings.bCheckPermissions = bCheckPermissions;
    x.modules.import.settings.bClone2Source = bClone2Source;
    

    %% 2. Initialize the setup of the dicom2nii conversion
    imPar = xASL_imp_DCM2NII_Initialize(studyPath, imParPath);


    %% 3. Run the DCM2NIIX
    if bRunSubmodules(1)
        xASL_imp_DCM2NII(imPar, x);
    end


    %% 4. Run the NIIX to ASL-BIDS
    if bRunSubmodules(2)
        % Run NII to BIDS
        xASL_imp_NII2BIDS(imPar, studyPath, studyParPath);
        % Update x.DataParPath
        [x] = xASL_imp_Import_UpdateDataParPath(x, studyPath);
    end


    %% 5. Run defacing
    if bRunSubmodules(3)
        xASL_imp_Anonymize(imPar);
    end


end


% Update x.DataParPath to dataset_description.json after NII2BIDS conversion
function [x] = xASL_imp_Import_UpdateDataParPath(x, studyPath)

    % Search for dataset_description.json within the rawdata subfolder
    foundFiles = xASL_adm_GetFileList(fullfile(studyPath,'rawdata'),'dataset_description.json');
    
    % Check if valid dataset_description.json exists within the rawdata folder
    if isempty(foundFiles)
        warning('No valid dataset_description.json found within the rawdata directory...');
    else
        x.DataParPath = foundFiles{1};
    end

end



