function xASL_module_Import(studyPath, imParPath, studyParPath, bRunSubmodules, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source, x)
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
% Import batch T1, T2, FLAIR, DWI, fMRI, M0, ASL data from dicom 2 NIfTI in ASL-BIDS format and structure.
% Uses dcm2niiX for the conversion, and additionally collects important DICOM header data
% and puts them in .json sidecars to be used with the ExploreASL pipeline.
% This function takes any folder input, but the folder input should be
% specified in the imPar definition. Follow the steps below, for study "MyStudy" located on "//MyDisk":
%
% 1. Make sure you have your DICOM data. Export them from XNAT, download them, or whatsoever
%    Create a root folder with study ID name, and put the DICOMs in any structure in the sourcedata folder within the study ID root folder
%    Example:
%    imPar.StudyID: MyStudy
%    StudyRoot folder: //MyDisk/MyStudy
%    sourcedata folder containing DICOMs: //MyDisk/MyStudy/sourcedata
% 2. Make sure that your DICOM data has any structure that can be retrieved
%    from the folder and/or file names. This function doesn't yet read the DICOM headers
%    For a quick and dirty (but actually slow) function that converts a
%    DICOM folder/file structure into readable format, first run
%    ConvertDicomFolderStructure_CarefulSlow.m. This will read each DICOM
%    individually, and put it in a folder with the name identical to the
%    DICOMs SeriesName/ProtocolName.
% 3. Once you have all DICOMs in folderstructure with identifyable names
%    inside //MyDisk/MyStudy/sourcedata, set up the folderstructure in
%    ExploreASL_ImportConfig.m. This setup uses the SPM form of regular
%    expressions, which can be daunting at first, but are very flexible.
%    Easiest is to study other examples, before creating your own.
%    For this example, let's say we have //MyDisk/MyStudy/sourcedata/ScanType/SubjectName
%    because we downloaded our data from XNAT, ordered per ScanType first,
%    and then per subject.
%
%    BRIEF EXPLANATION:
%    Let's suppose we don't have sessions (only a single structural and functional scan per subject)
%    The names of our scans comes out of XNAT as '3D_FLAIR_eyesClosed', 'T1w_MPRAGE' and 'PCASL_10_min'
%
%    and the subject names are 'MyStudy001' .. 'MyStudy002' .. etc.
%
%    imPar.folderHierarchy   - contains a a cell array of regular expressions, with each cell specifying a directory layer/level
%                              the parts within brackets () tell the script that this is a token (i.e. subject, session, ScanType)
%                              Example:
%                              imPar.folderHierarchy = {'^(3D_FLAIR|T1w|PCASL).*', '^(Sub-\d{3})$'};
%                              here we say that there are two folder layers '', separated by comma ,
%                              where the names between brackets are used to define what is what.
%                              ^ means that the foldername has to start with the following, $ means that the previous has to be the end of the foldername
%                              .* means anything, anylength, \d{3} means three digits
%    imPar.tokenOrdering     - defines which tokens are captured by the brackets () in imPar.folderHierarchy: position 1==subject, 2==visit, 3==session, 4==ScanType
%                              Example:
%                              imPar.tokenOrdering = [2 3 0 1]; stating that subject is the 2nd token, visit is the 3rd token, session has no token (i.e. no session) and ScanType is the 1st token
%    imPar.tokenVisitAliases - cell array that defines the aliases for the Visits, i.e. it tells the script which scans are which timepoint/visit.
%                              Similar as explained below for ScanAliases.
%                              First column contains the names that are
%                              recognized in sourcedata DICOM folders for visits,
%                              second column how it is named in NIfTI
%                              structure (should be _1 _2 _3 etc).
%                              Example:
%                              imPar.tokenVisitAliases = {'Screening','_1'; 'Month_12','_2'; 'Month_24','_3'; 'Month_36','_4'; 'Month_48','_5'};
%                              Note that if you specify tokenVisitAliases, the folders will receive
%                              the indices (e.g. _1 _2 _3), or even _1 only with a single Visit). If you don't specify
%                              them, they will not get this postfix.
%    imPar.tokenScanAliases  - cell array that defines the aliases for the ScanTypes, i.e. it tells the script which scans are which ScanType.
%                              First column should contain regular expression corresponding with the matching criteria in imPar.folderHierarchy
%                              whereas the second column contains the
%                              alias. Following valid aliases exist:
%                              'T1' 'FLAIR' 'ASL4D' 'M0' 'ASL4D_RevPE' 'func' 'func_NormPE' 'func_RevPE' 'dwi' 'dwi_RevPE' 'DSC4D'
%                              Example:
%                              imPar.tokenScanAliases = {'^3D_FLAIR$', 'FLAIR'; '^T1w$', 'T1'; '^PCASL$', 'ASL4D'};
%    imPar.tokenSessionAliases-same as tokenScanAliases but for sessions
%                              Example:
%                              imPar.tokenSessionAliases = {}; as we don't have sessions
%    imPar.bMatchDirectories - true if the last layer is a folder, false if the last layer is a filename (as e.g. with PAR/REC, enhanced DICOMs)
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

if nargin<5 || isempty(bCopySingleDicoms)
    bCopySingleDicoms = false; % by default don't copy dicoms for anonymization reasons
end

if nargin<6 || isempty(bUseDCMTK)
    bUseDCMTK = true; % default set to using DCM-TK
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

if nargin<9 || isempty(x)
    x = ExploreASL_Initialize; % only initialize ExploreASL if this wasnt initialized before
end

%% 2. Initialize the setup of the dicom2nii conversion

% Initialize the imPar field
[fpath, fname, fext] = fileparts(studyPath);

% Load the imPar from the file
imPar = spm_jsonread(imParPath);

% Initialize
imPar = xASL_bids_DCM2NII_Initialize(imPar, fpath, fname, fext);

%% 3. Run the DCM2NIIX
if bRunSubmodules(1)
	xASL_bids_DCM2NII(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bClone2Source,x);
end

%% 4. Run the NIIX to ASL-BIDS
if bRunSubmodules(2)
    xASL_bids_NII2BIDS(studyParPath, imPar);
end


%% 5. Run defacing
if bRunSubmodules(3)
	listSubjects = xASL_adm_GetFileList(imPar.AnalysisRoot,[],false,[],true);
	for iSubject = 1:length(listSubjects)
		
		subjectLabel = xASL_adm_CorrectName(listSubjects{iSubject},2);
		
		% Check if the anatomical directory exists
		if exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat'),'dir')
			% Process all anatomical files
			fAnat = xASL_adm_GetFileList(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat'),'^.+\.nii',false,[]);
			for iAnat = 1:length(fAnat)
				%Unzip the file for SPM
				pathUnzipped = xASL_adm_UnzipNifti(fullfile(imPar.BidsRoot,['sub-' subjectLabel],'anat',fAnat{iAnat}));
				% Remove the face
				xASL_spm_deface(pathUnzipped,true);
				% Zip again
				gzip(pathUnzipped);
				delete(pathUnzipped);
			end
		end
	end
end % End of defacing

% end of import
end


% -----------------------------------------------------------------------------
% Append Nifti Parameters
% -----------------------------------------------------------------------------
function s = AppendNiftiParameters(nii_files)
% This function outputs s=[FileName voxel size XYZ matrix size XYZ]
s = [];

if ischar(nii_files)
    nii_files = {nii_files};
end

for iNii=1:length(nii_files)
    [~, Ffile, Fext] = fileparts(nii_files{iNii});
    s = sprintf(',"%s"', [Ffile Fext]); % filename

    tempnii = xASL_io_ReadNifti(nii_files{iNii});
    s = [s sprintf(',%g', tempnii.hdr.pixdim(2:5) )]; % voxel size XYZ
    s = [s sprintf(',%g', tempnii.hdr.dim(2:5) )]; % matrix size XYZ
end
end

% -----------------------------------------------------------------------------
% Append Parms Parameters
% -----------------------------------------------------------------------------
function [s, FieldNames] = AppendParmsParameters(parms)
% This function outputs s=fields of _parms.mat
s = [];

FieldNames = {'RepetitionTimePreparation' 'RepetitionTime' 'EchoTime' 'NumberOfAverages' 'RescaleSlope' 'RescaleSlopeOriginal'...
    'MRScaleSlope' 'RescaleIntercept' 'AcquisitionTime' 'AcquisitionMatrix' 'TotalReadoutTime'...
    'EffectiveEchoSpacing'};

if ~isempty(parms)
    for iField=1:length(FieldNames)
        if isfield(parms,FieldNames{iField})
            s = [s ',' xASL_num2str(parms.(FieldNames{iField}))];
        else
            s = [s ',n/a'];
        end
    end
end
    
end

% -----------------------------------------------------------------------------
% Catch Errors
% -----------------------------------------------------------------------------
function [dcm2niiCatchedErrors] = CatchErrors(WarningID, WarningMessage, WarningLine, WarningFileName, WarningPath, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar, StackIn)
% Catch reported warnings/errors, print them if verbose, & add them to a structure of warnings/errors to be stored for later QC

    if imPar.bVerbose % print warning if we want verbose
        warning(WarningMessage);
    end

    % Find index of the warning to store
    if isempty(fields(dcm2niiCatchedErrors))
        IndexN = 1;
    else
        IndexN = length(dcm2niiCatchedErrors)+1;
    end

    % store the warning/error
    dcm2niiCatchedErrors(IndexN).scan_name = scan_name;
    dcm2niiCatchedErrors(IndexN).scanpath = scanpath;
    dcm2niiCatchedErrors(IndexN).destdir = destdir;
    dcm2niiCatchedErrors(IndexN).identifier = WarningID;
    dcm2niiCatchedErrors(IndexN).message = WarningMessage;
    dcm2niiCatchedErrors(IndexN).cause = 'n/a';

    if exist('StackIn', 'var')
        dcm2niiCatchedErrors(IndexN).stack = StackIn;
    else
        dcm2niiCatchedErrors(IndexN).stack.file = fullfile(WarningPath, [WarningFileName '.m']);
        dcm2niiCatchedErrors(IndexN).stack.name = WarningFileName;
        dcm2niiCatchedErrors(IndexN).stack.line = WarningLine(end).line;
    end

end





function imPar = xASL_bids_DCM2NII_Initialize(imPar, fpath, fname, fext)

    % Specify paths
    if ~isfield(imPar,'studyID') || isempty(imPar.studyID)
        imPar.studyID = [fname fext];
    end
    if ~isfield(imPar,'AnalysisRoot') || isempty(imPar.AnalysisRoot) 
        imPar.AnalysisRoot = fpath;
    end
    if ~isfield(imPar,'RawRoot') || isempty(imPar.RawRoot)
        imPar.RawRoot = fpath;
    end
    if ~isfield(imPar,'BidsRoot') || isempty(imPar.BidsRoot)
        imPar.BidsRoot = fpath;
    end

    % Finalize the directories
    imPar.RawRoot = fullfile(imPar.RawRoot,imPar.studyID, 'sourcedata'); % default name
    imPar.AnalysisRoot = fullfile(imPar.AnalysisRoot,imPar.studyID,'analysis');
    imPar.BidsRoot = fullfile(imPar.BidsRoot,imPar.studyID,'rawdata');

    % Specify the tokens
    if ~isfield(imPar,'folderHierarchy')
        imPar.folderHierarchy = {}; % must define this per study; use a cell array of regular expressions. One cell per directory level.
    end
    if ~isfield(imPar,'tokenOrdering')
        imPar.tokenOrdering = []; % must match imPar.folderHierarchy: 1==subject, 2=visit, 3==session, 4==scan (if visit or session are omitted, they will be skipped)
    else
        imPar.tokenOrdering = imPar.tokenOrdering(:)';
    end
    if ~isfield(imPar,'tokenScanAliases')
        imPar.tokenScanAliases = [];
    else
        if (size(imPar.tokenScanAliases,2) > 2) || (size(imPar.tokenScanAliases,2) == 1)
            tokenScanAliasesOld = imPar.tokenScanAliases;
            imPar.tokenScanAliases = tokenScanAliasesOld(1:2:end);
            imPar.tokenScanAliases(:,2) = tokenScanAliasesOld(2:2:end);
        end
    end
    if ~isfield(imPar,'tokenVisitAliases')
        imPar.tokenVisitAliases = [];
    else
        if (size(imPar.tokenVisitAliases,2) > 2) || (size(imPar.tokenVisitAliases,2) == 1)
            tokenVisitAliasesOld = imPar.tokenVisitAliases;
            imPar.tokenVisitAliases = tokenVisitAliasesOld(1:2:end);
            imPar.tokenVisitAliases(:,2) = tokenVisitAliasesOld(2:2:end);
        end
    end
    if ~isfield(imPar,'tokenSessionAliases')
        imPar.tokenSessionAliases = [];
    else
        if (size(imPar.tokenSessionAliases,2) > 2) || (size(imPar.tokenSessionAliases,2) == 1)
            tokenSessionAliasesOld = imPar.tokenSessionAliases;
            imPar.tokenSessionAliases =tokenSessionAliasesOld(1:2:end);
            imPar.tokenSessionAliases(:,2) = tokenSessionAliasesOld(2:2:end);
        end
    end
    if ~isfield(imPar,'bMatchDirectories')
        imPar.bMatchDirectories  = false;
    end

    % Specify the additional details of the conversion
    if ~isfield(imPar,'bVerbose') || isempty(imPar.bVerbose)
        imPar.bVerbose = true;
    end
    if ~isfield(imPar,'bOverwrite') || isempty(imPar.bOverwrite)
        imPar.bOverwrite  = false; % NB, the summary file will be recreated anyway and dicom conversion in temp is always done, even if dest. exists
    end
    if ~isfield(imPar,'visitNames') || isempty(imPar.visitNames)
        imPar.visitNames = {};
    end
    if ~isfield(imPar,'nMaxVisits') || isempty(imPar.nMaxVisits)
        imPar.nMaxVisits = 0;
    end
    if ~isfield(imPar,'sessionNames') || isempty(imPar.sessionNames)
        imPar.sessionNames = {};
    end
    if ~isfield(imPar,'nMaxSessions') || isempty(imPar.nMaxSessions)
        imPar.nMaxSessions = 0;
    end

end


