function [result, x] = xASL_module_Import(x)
%xASL_module_Import Imports the DICOM or PAR/REC source data to NIFTIs in ASL-BIDS format
%
% FORMAT: [result, x] = xASL_module_Import(x)
%
% INPUT:
%   x                     - ExploreASL x structure. (STRUCT, REQUIRED)
%   x.dir.DatasetRoot     - Path to the study directory containing the 'sourcedata'
%                           directory with the DICOM files. (REQUIRED)
%   x.dir.sourceStructure - Path to the JSON file with structure with import parameters.
%                           All other input parameters are configured in this function. (OPTIONAL)
%                           The path is optional, but the file has to be there. Either provided as a full 
%                           path or a filename in the path, or default names (case-insensitive) sourceStructure.json 
%                           are seeked and then loaded inside x.modules.import.imPar.
%   x.dir.studyPar        - Path to the JSON file with the BIDS parameters relevant for the whole study.
%                           These parameters are used if they cannot be extracted from the DICOMs automatically.
%                           Looking automatically for file studyPar.json. (OPTIONAL)
%   x.opts.bImport        - Specify which of the parts should be run (OPTIONAL, DEFAULT [0 0 0])
%                           [1 0 0] - Run the DICOM to NIFTI conversion
%                           [0 1 0] - Run the NIFTI transformation to the proper ASL-BIDS
%                           [0 0 0] - Run the Defacing module
%   x.modules.import.settings.bCopySingleDicoms
%                         - If true, copies a single DICOM with each NIfTI dataset/ScanType,
%                           that can be used to retrieve missing parameters from the DICOM header,
%                           or as dummy DICOM to dump embed data into (e.g. WAD-QC). (DEFAULT = false)
%   x.modules.import.settings.bUseDCMTK
%                         - If true, then use DCMTK, otherwise use DICOMINFO from Matlab. (DEFAULT = false)
%   x.modules.import.settings.bCheckPermissions   
%                         - If true, check whether data permissions are set correctly, 
%                           before trying to read/copy the files. (DEFAULT=false)
% 
% 
% OUTPUT:
%   x        - ExploreASL x structure
%   result   - True for successful run of this module, false for insuccessful run
% 
% 
% OUTPUT FILES:
%   //temp/dcm2niiCatchedErrors.(mat|json) - overview of catched dcm2nii errors, or other errors in this function
%   //temp/import_log_StudyID_yyyymmdd_hhmmss.txt - diary log of this function
%   //temp/import_summary.csv - hence the name
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
% Import batch `T1`, `T2`, `FLAIR`, `DWI`, `fMRI`, `M0`, `ASL` data from dicom 2 NIfTI in ASL-BIDS format and structure.
% Uses dcm2niiX for the conversion, and additionally collects important DICOM header data
% and puts them in `.json` sidecars to be used with the ExploreASL pipeline.
% This function takes any folder input, but the folder input should be
% specified in the x.modules.import.imPar definition. Follow the steps below, for study `"MyStudy"` located on `"//MyDisk"`:
%
% 1. Make sure you have your DICOM data. Export them from XNAT, download them, or whatsoever
%    Create a root folder with study ID name, and put the DICOMs in any structure in the sourcedata folder within the study ID root folder
%    Examples:
%    x.modules.import.imPar.StudyID: MyStudy
%    Dataset Root folder: `//MyDisk/MyStudy`
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
%    imPar is now contained inside x.modules.import.imPar
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
%    - imPar.bMatchDirectories   - true if the last layer is a folder, false if the last layer is a filename (e.g., \.PAR, \.nii, \.nii\.gz, or \.dcm for enhanced DICOM)
%
% Any (X|Y|Z) expression is referred to as a "captured group"
% All captured groups defined in tokenOrdering, scanAliases, visitAliases, and sessionAliases are "tokens"
% While it is simplest if all captured groups in folderHierarchy represent tokens (for the user and bugfixing), this is not required 
%
% EXAMPLE: [~, x] = xASL_init_Iteration(x,'xASL_module_Import');
% __________________________________
% Copyright 2015-2022 ExploreASL


    %% Import Module
    
    result = false;
    
    % DCM2NIIX and other tools seem to stop the diary logging automatically, here we extract the current 
    % diary file path to make sure that at the beginning of each module the logging is still enabled.
    x.dir.diaryFile = get(0,'DiaryFile');
    
    % Start Mutex
    x = xASL_init_InitializeMutex(x, 'Import');
    fprintf('\n');
    
    % Define lock states
    StateName{1} = '010_DCM2NII';
    StateName{2} = '020_NII2BIDS';
	StateName{3} = '030_DEFACE';
   
    
    %% 0. Initialization
    
    % First do the basic parameter admin and initialize the default values
    if nargin<1 || isempty(x)
        error('The x struct needs to be defined...');
    end
    
    % Initialize x struct
    x = xASL_init_SubStructs(x);
    
    %% 1. Run the DCM2NIIX
    iState = 1;
    if x.opts.bImport(1) && ~x.mutex.HasState(StateName{1})
        x = xASL_wrp_DCM2NII(x);
        x.mutex.AddState(StateName{iState});
    elseif x.opts.bImport(1) && x.mutex.HasState(StateName{1})
        fprintf('DCM2NIIX was run before...   \n');
    end

    
    %% 2. Run the NIfTI to ASL-BIDS
    iState = 2;
    if x.opts.bImport(2) && ~x.mutex.HasState(StateName{2})
        x = xASL_wrp_NII2BIDS(x);
        x.mutex.AddState(StateName{iState});
    elseif x.opts.bImport(2) && x.mutex.HasState(StateName{2})
        fprintf('NIIX to ASL-BIDS was run before...   \n');
	end
    
	%% 3. Run the DCM2NIIX
    iState = 3;
    if x.opts.bImport(3) && ~x.mutex.HasState(StateName{3})
        xASL_wrp_Deface(x);
        x.mutex.AddState(StateName{iState});
    elseif x.opts.bImport(3) && x.mutex.HasState(StateName{3})
        fprintf('DEFACE was run before...   \n');
    end
    
    %% 4. Clean-up
    x = xASL_adm_CleanUpX(x);
    
    % We need to terminate the module correctly
    x.mutex.AddState('999_ready');
    x.mutex.Unlock();
    x.result = true;
    close all;
    
    % Return the results
    result = x.result;
end

