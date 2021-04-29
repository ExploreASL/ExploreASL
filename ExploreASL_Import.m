function ExploreASL_Import(imPar, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bRunDCM2NII, bClone2Source, x)
%ExploreASL_Import Imports the DICOM or PAR/REC raw data to NIFTIs
%
% FORMAT: ExploreASL_Import(imPar[, bCopySingleDicoms, bUseDCMTK, bCheckPermissions, bRunDCM2NII, bClone2Source])
%
% INPUT:
%   imPar               - structure with import parameters, output of ExploreASL_ImportConfig.m
%                         All other input parameters are configured within this function.
%   bCopySingleDicoms   - if true, copies a single DICOM with each NIfTI
%                         dataset/ScanType, that can be used to retrieve missing parameters from
%                         the DICOM header, or as dummy DICOM to dump embed data into (e.g. WAD-QC) (DEFAULT=false)
%   bUseDCMTK           - if true, then use DCMTK, otherwise use DICOMINFO from Matlab (DEFAULT=false)
%   bCheckPermissions   - if true, check whether data permissions are set correctly, before trying to read/copy the files (DEFAULT=false)
%   bRunDCM2NII         - if true, run dcm2niiX. Setting this to false allows to skip dcm2niiX and only create .mat files (DEFAULT=true)
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
% Import batch T1, FLAIR, DWI, fMRI, M0, ASL data from dicom 2 NIfTI.
% Uses dcm2niiX for the conversion, and additionally collects important DICOM header data
% and puts them in _parms.mat and .json sidecars to be used with the ExploreASL pipeline.
% This function takes any folder input, but the folder input should be
% specified in the imPar definition in the ExploreASL_ImportConfig.m (later
% to be converted to e.g. a CSV file). Follow the steps below, for study "MyStudy" located on "//MyDisk":
%
% 1. Make sure you have your DICOM data. Export them from XNAT, download them, or whatsoever
%    Create a root folder with study ID name, and put the DICOMs in any structure in the raw folder within the study ID root folder
%    Examples:
%    imPar.StudyID: MyStudy
%    StudyRoot folder: //MyDisk/MyStudy
%    Raw folder containing DICOMs: //MyDisk/MyStudy/raw
% 2. Make sure that your DICOM data has any structure that can be retrieved
%    from the folder and/or file names. This function doesn't yet read the DICOM headers
%    For a quick and dirty (but actually slow) function that converts a
%    DICOM folder/file structure into readable format, first run
%    ConvertDicomFolderStructure_CarefulSlow.m. This will read each DICOM
%    individually, and put it in a folder with the name identical to the
%    DICOMs SeriesName/ProtocolName.
% 3. Once you have all DICOMs in folderstructure with identifyable names
%    inside //MyDisk/MyStudy/raw, set up the folderstructure in
%    ExploreASL_ImportConfig.m. This setup uses the SPM form of regular
%    expressions, which can be daunting at first, but are very flexible.
%    Easiest is to study other examples, before creating your own.
%    For this example, let's say we have //MyDisk/MyStudy/raw/ScanType/SubjectName
%    because we downloaded our data from XNAT, ordered per ScanType first,
%    and then per subject.
%
%    BRIEF EXPLANATION:
%    Let's suppose we don't have sessions (only a single structural and functional scan per subject)
%    The names of our scans comes out of XNAT as '3D_FLAIR_eyesClosed', 'T1w_MPRAGE' and 'PCASL_10_min'
%
%    and the subject names are 'MyStudy001' .. 'MyStudy002' .. etc.
%
%    - imPar.folderHierarchy   - contains a a cell array of regular expressions, with each cell specifying a directory layer/level
%                                the parts within brackets () tell the script that this is a token (i.e. subject, session, ScanType)
%                                Examples:
%                                imPar.folderHierarchy = {'^(3D_FLAIR|T1w|PCASL).*', '^(Sub-\d{3})$'};
%                                here we say that there are two folder layers '', separated by comma ,
%                                where the names between brackets are used to define what is what.
%                                ^ means that the foldername has to start with the following, $ means that the previous has to be the end of the foldername
%                                .* means anything, anylength, \d{3} means three digits
%    - imPar.tokenOrdering     - defines which tokens are captured by the brackets () in imPar.folderHierarchy: position 1==subject, 2==visit, 3==session, 4==ScanType
%                                Examples:
%                                imPar.tokenOrdering = [2 3 0 1]; stating that subject is the 2nd token, visit is the 3rd token, session has no token (i.e. no session) and ScanType is the 1st token
%    - imPar.tokenVisitAliases - cell array that defines the aliases for the Visits, i.e. it tells the script which scans are which timepoint/visit.
%                                Similar as explained below for ScanAliases.
%                                First column contains the names that are
%                                recognized in raw DICOM folders for visits,
%                                second column how it is named in NIfTI
%                                structure (should be _1 _2 _3 etc).
%                                Examples:
%                                imPar.tokenVisitAliases = {'Screening','_1'; 'Month_12','_2'; 'Month_24','_3'; 'Month_36','_4'; 'Month_48','_5'};
%                                Note that if you specify tokenVisitAliases, the folders will receive
%                                the indices (e.g. _1 _2 _3), or even _1 only with a single Visit). If you don't specify
%                                them, they will not get this postfix.
%    - imPar.tokenScanAliases  - cell array that defines the aliases for the ScanTypes, i.e. it tells the script which scans are which ScanType.
%                                First column should contain regular expression corresponding with the matching criteria in imPar.folderHierarchy
%                                whereas the second column contains the
%                                alias. Following valid aliases exist:
%                                'T1' 'FLAIR' 'ASL4D' 'M0' 'ASL4D_RevPE' 'func' 'func_NormPE' 'func_RevPE' 'dwi' 'dwi_RevPE' 'DSC4D'
%                                Examples:
%                                imPar.tokenScanAliases = {'^3D_FLAIR$', 'FLAIR'; '^T1w$', 'T1'; '^PCASL$', 'ASL4D'};
%    - imPar.tokenSessionAliases-same as tokenScanAliases but for sessions
%                                Examples:
%                                imPar.tokenSessionAliases = {}; as we don't have sessions
%    - imPar.bMatchDirectories - true if the last layer is a folder, false if the last layer is a filename (as e.g. with PAR/REC, enhanced DICOMs)
%
% EXAMPLE: ExploreASL_Import(ExploreASL_ImportConfig('//MyDisk/MyStudy'));
% __________________________________
% Copyright 2015-2020 ExploreASL



%% Check input parameters
dcm2niiCatchedErrors = struct; % initialization
if nargin<2 || isempty(bCopySingleDicoms)
    bCopySingleDicoms = false; % by default don't copy dicoms for anonymization reasons
end
if nargin<5 || isempty(bRunDCM2NII)
    bRunDCM2NII = true;
end
if nargin<3 || isempty(bUseDCMTK)
    bUseDCMTK = true; % default set to using DCM-TK
elseif ~bUseDCMTK && isempty(which('dicomdict'))
    error('Dicomdict missing, image processing probably not installed, try DCMTK instead');
end
if nargin<4 || isempty(bCheckPermissions)
    if isunix
        bCheckPermissions = true;
    else
        bCheckPermissions = false;
    end
end
if nargin<6 || isempty(bClone2Source)
    bClone2Source = false;
end
if nargin<7 || isempty(x)
    x = ExploreASL_Initialize; % only initialize ExploreASL if this wasnt initialized before
end
if nargin<1 || ~isfield(imPar,'studyID')
    error('ExploreASL_Import: Please provide studyID (name of the study folder) as input parameter');
end

dcm2niiDir = fullfile(x.MyPath, 'External', 'MRIcron');
xASL_adm_CheckPermissions(dcm2niiDir, true); % dcm2nii needs to be executable

%%%%%%%%%%%%%%%%%%%%%%
% Initialize defaults
%%%%%%%%%%%%%%%%%%%%%%
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
if ~isfield(imPar,'dcm2nii_version') || isempty(imPar.dcm2nii_version)
    imPar.dcm2nii_version = '20190902'; % OR for PARREC imPar.dcm2nii_version = '20101105'; THIS IS AUTOMATED BELOW
end
if ~isfield(imPar,'dcmExtFilter') || isempty(imPar.dcmExtFilter)
    imPar.dcmExtFilter = '^(.*\.dcm|.*\.img|.*\.IMA|[^.]+|.*\.\d*)$'; % the last one is because some convertors save files without extension, but there would be a dot/period before a bunch of numbers
end
if isempty(imPar.RawRoot)
    error('imPar.RawRoot was empty');
elseif isempty(imPar.AnalysisRoot)
    error('imPar.AnalysisRoot was empty');
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
if ~isfield(imPar,'SkipSubjectIfExists') || isempty(imPar.SkipSubjectIfExists)
    % allows to skip existing subject folders in the analysis folder, when this is set to true,
    % avoiding partly re-importing/converting dcm2niiX when processing has been partly done
    imPar.SkipSubjectIfExists = false;
else
    warning('Skipping existing subjects in analysis folder');
    fprintf('If you want to overwrite, first remove the full subject folder');
end

%% Create the basic folder structure for raw & derivative data
imPar.RawRoot = fullfile(imPar.RawRoot,imPar.studyID, 'raw'); % default name

if ~exist(imPar.RawRoot, 'dir')
    warning(['Couldnt find ' imPar.RawRoot ', trying to find a different folder instead...']);
    
    % find any folder except for analysis, source, raw, derivatives
    % xASL_adm_GetFileList uses regular expressions, to create a nice list of foldernames,
    % with/without FullPath (FPList), with/without recursive (FPListRec)
    % very powerful once you know how these work
    FolderNames = xASL_adm_GetFileList(fullfile(imPar.RawRoot, imPar.studyID), '^(?!(analysis|derivatives|source|raw)).*$', 'FPList', [0 Inf], true); 
    
    if length(FolderNames)==1
        imPar.RawRoot = FolderNames{1};
        fprintf('%s\n', ['Found ' imPar.RawRoot ' as raw folder']);
    else
        error('Couldnt find a raw folder, please rename one, or move other folders');
    end
end

imPar.AnalysisRoot = fullfile(imPar.AnalysisRoot,imPar.studyID,'analysis');
imPar.SourceRoot = fullfile(imPar.AnalysisRoot,imPar.studyID,'source');

xASL_adm_CreateDir(imPar.AnalysisRoot);

if bCheckPermissions
    xASL_adm_CheckPermissions(imPar.RawRoot, false); % don"t need execution permisions
    xASL_adm_CheckPermissions(imPar.AnalysisRoot, false);  % don"t need execution permisions
end

% here we try to fix backwards compatibility, but this may break
if length(imPar.tokenOrdering)==3 % backwards compatibility Visits
	if size(imPar.tokenOrdering,1) > 1
		% Vertical vector
		imPar.tokenOrdering = imPar.tokenOrdering';
	end
	imPar.tokenOrdering = [imPar.tokenOrdering(1) 0 imPar.tokenOrdering(2:3)]; % insert Visits(none)
elseif length(imPar.tokenOrdering)==2 % backwards compatibility Visits & sessions
    imPar.tokenOrdering = [imPar.tokenOrdering(1) 0 0 imPar.tokenOrdering(2)]; % insert sessions & visits(none)
end

% Path to the dictionary to initialize - we need to keep track if the dictionary has been set, because
% Dicominfo can be used despite bUSEDCMTK==1 when DCMTK fails
pathDcmDict = fullfile(x.MyPath,'External','xASL_DICOMLibrary.txt');
if ~bUseDCMTK
    % -----------------------------------------------------------------------------
    % Initialize dicom dictionary by appending private philips stuff to a temporary copy
    % -----------------------------------------------------------------------------
    dicomdict('set', pathDcmDict);
end

fid_summary = -1; % initialize to be able to catch errors and close if valid

% change dcmnii_version for PARREC if needed
if ~isempty(strfind(char(imPar.folderHierarchy(end)),'PAR'))
    imPar.dcm2nii_version = '20101105';
end

% redirect output to a log file
diary_filepath = fullfile(imPar.AnalysisRoot, ['import_log_' imPar.studyID '_' datestr(now,'yyyymmdd_HHMMSS') '.txt']);
diary(diary_filepath);

if imPar.bMatchDirectories
    strLookFor = 'Directories';
else
    strLookFor = 'Files';
end

%%
% -----------------------------------------------------------------------------
% Start with defining the subjects, visits, sessions (i.e. BIDS runs) and scans (i.e. ScanTypes) by listing or typing
% -----------------------------------------------------------------------------

% Recursively scan the directory tree using regular exspressions at each directory level. Use ()-brackets to extract tokens
% that will be used to identify subjects, sessions and scans. In the loop below it is possible to translate the tokens
% to more convenient strings before using them in destination paths.
[matches, tokens] = xASL_adm_FindByRegExp(imPar.RawRoot, imPar.folderHierarchy, 'StripRoot', true, 'Match', strLookFor,'IgnoreCase',true);
if isempty(matches)
    warning('No matching files, skipping');
    return;
elseif imPar.bVerbose
    fprintf('\nMatching files:\n');
    disp(matches);
    fprintf('#=%g\n',length(matches));
end

% Copy the columns into named vectors. This construction allows for arbitrary directory hierarchies.
% Make sure to select the columns that correspond to the folder ordering defined using the regular expressions above.
% Define Subjects

% SUBJECTS
vSubjectIDs = tokens(:,imPar.tokenOrdering(1)); % cell vector with extracted subject IDs (for all sessions and scans)

% VISITS
if imPar.tokenOrdering(2)==0
    % a zero means: no visits applicable
    bUseVisits = false;
    vVisitIDs = cellfun(@(y) '1', vSubjectIDs, 'UniformOutput', false); % each subject has a single visit
    imPar.tokenVisitAliases = {'^1$', '_1'};
else
    bUseVisits = true;
    vVisitIDs = tokens(:,imPar.tokenOrdering(2)); % cell vector with extracted session IDs (for all subjects and scans)
end

% SESSIONS
if imPar.tokenOrdering(3)==0
    % a zero means: no sessions applicable
    bUseSessions = false;
    vSessionIDs = cellfun(@(y) '1', vSubjectIDs, 'UniformOutput', false); % each subject-visit has a single session
    imPar.tokenSessionAliases = {'^1$', 'ASL_1'};
else
    bUseSessions = true;
    vSessionIDs = tokens(:,imPar.tokenOrdering(3)); % cell vector with extracted session IDs (for all subjects and scans)
end

% SCANTYPES
vScanIDs = tokens(:,imPar.tokenOrdering(4)); % cell vector with extracted scan IDs (for all subjects and sessions)

% Convert the vectors to unique & sort sets by the output aliases
subjectIDs  = sort(unique(vSubjectIDs));
nSubjects = length(subjectIDs);

visitIDs  = unique(vVisitIDs);
% sort by output
if length(visitIDs)>1
    for iV=1:length(visitIDs)
        IDrow(iV) = find(cellfun(@(y) strcmp(y,visitIDs{iV}), imPar.tokenVisitAliases(:,1)));
    end
    visitIDs = visitIDs(IDrow);
end
nVisits = length(visitIDs);

sessionIDs  = sort(unique(vSessionIDs));
nSessions = length(sessionIDs);

scanIDs = sort(unique(lower(vScanIDs)));
nScans = length(scanIDs);

% VISIT NAMES
if isempty(imPar.visitNames)
    if isempty(visitIDs)
        imPar.visitNames = cell(nVisits,1);
        for kk=1:nVisits
            imPar.visitNames{kk}=sprintf('ASL_%g',kk);
        end
    else
        imPar.visitNames = visitIDs;
    end
end

% SESSION NAMES
% optionaly we can have human readble session names; by default they are the same as the original tokens in the path
if isempty(imPar.sessionNames)
    if isempty(sessionIDs)
        imPar.sessionNames = cell(nSessions,1);
        for kk=1:nSessions
            imPar.sessionNames{kk}=sprintf('ASL_%g',kk);
        end
    else
        imPar.sessionNames = sessionIDs;
    end
end

% SCAN NAMES
scanNames = scanIDs;

% sanity check for missing elements
if nSubjects==0
    error('No subjects')
end
if nVisits==0
    error('No visits')
end
if nSessions==0
    error('No sessions')
end
if nScans==0
    error('No scans')
end

% preallocate space for (global) counts
globalCounts.converted_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
globalCounts.skipped_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
globalCounts.missing_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans

% define a cell array for storing info for parameter summary file
xASL_adm_CreateDir(imPar.AnalysisRoot);
summary_lines = cell(nSubjects,nVisits,nSessions,nScans);

%% -----------------------------------------------------------------------------
% import subject by subject, visit by visit, session by session, scan by scan
% -----------------------------------------------------------------------------
fprintf('%s\n', 'Running import (i.e. dcm2niiX)');

% Create ID struct
listsIDs.vSubjectIDs = vSubjectIDs;
listsIDs.vVisitIDs = vVisitIDs;
listsIDs.vSessionIDs = vSessionIDs;
listsIDs.vScanIDs = vScanIDs;
listsIDs.subjectIDs = subjectIDs;
listsIDs.visitIDs = visitIDs;
listsIDs.sessionIDs = sessionIDs;
listsIDs.scanIDs = scanIDs;

% Create number of struct
numOf.nSubjects = nSubjects;
numOf.nVisits = nVisits;
numOf.nSessions = nSessions;
numOf.nScans = nScans;

% Create settings struct
settings.bUseVisits = bUseVisits;
settings.bClone2Source = bClone2Source;
settings.bUseDCMTK = bUseDCMTK;
settings.bCopySingleDicoms = bCopySingleDicoms;
settings.bUseSessions = bUseSessions;

% Iterate over subjects
for iSubject=1:nSubjects
    [imPar, summary_lines, PrintDICOMFields, globalCounts, scanNames, dcm2niiCatchedErrors, pathDcmDict] = ...
        xASL_imp_DCM2NII_Subject(x, imPar, listsIDs, numOf, settings, globalCounts, scanNames, iSubject, summary_lines, matches, dcm2niiCatchedErrors, pathDcmDict);
end


%% -----------------------------------------------------------------------------
% create summary file
% -----------------------------------------------------------------------------
xASL_imp_CreateSummaryFile(imPar, numOf, listsIDs, PrintDICOMFields, globalCounts, scanNames, summary_lines, fid_summary);

% cleanup
if ~bUseDCMTK || isempty(pathDcmDict)
    dicomdict('factory');
end
diary('off');

if ~isempty(fields(dcm2niiCatchedErrors))
    fclose all;
    SavePath = fullfile(imPar.AnalysisRoot, 'dcm2niiCatchedErrors.mat');
    SaveJSON = fullfile(imPar.AnalysisRoot, 'dcm2niiCatchedErrors.json');
    xASL_delete(SavePath);
    xASL_delete(SaveJSON);
    save(SavePath,'dcm2niiCatchedErrors');
    spm_jsonwrite(SaveJSON, dcm2niiCatchedErrors);
end

fprintf('\n');

end



