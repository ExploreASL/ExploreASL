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
% 1) Make sure you have your DICOM data. Export them from XNAT, download them, or whatsoever
%    Create a root folder with study ID name, and put the DICOMs in any structure in the raw folder within the study ID root folder
%    Example:
%    imPar.StudyID: MyStudy
%    StudyRoot folder: //MyDisk/MyStudy
%    Raw folder containing DICOMs: //MyDisk/MyStudy/raw
% 2) Make sure that your DICOM data has any structure that can be retrieved
%    from the folder and/or file names. This function doesn't yet read the DICOM headers
%    For a quick and dirty (but actually slow) function that converts a
%    DICOM folder/file structure into readable format, first run
%    ConvertDicomFolderStructure_CarefulSlow.m. This will read each DICOM
%    individually, and put it in a folder with the name identical to the
%    DICOMs SeriesName/ProtocolName.
% 3) Once you have all DICOMs in folderstructure with identifyable names
%    inside //MyDisk/MyStudy/raw, set up the folderstructure in
%    ExploreASL_ImportConfig.m This setup uses the SPM form of regular
%    expressions, which can be daunting at first, but are very flexible.
%    Easiest is to study other examples, before creating your own.
%    For this example, let's say we have //MyDisk/MyStudy/raw/ScanType/SubjectName
%    because we downloaded our data from XNAT, ordered per ScanType first, and then per subject
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
%                              recognized in raw DICOM folders for visits,
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
%    imPar.RawRootModName    - If you don't define this field at all, the algorithm will search for the default 'raw' folder.
%                              You can select 'adaptive' to search for a random subfolder (assuming that there is only one. 
%                              Or you can define a specific folder.
%
% EXAMPLE: ExploreASL_Import(ExploreASL_ImportConfig('//MyDisk/MyStudy'));
% __________________________________
% Copyright 2015-2019 ExploreASL



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
elseif bUseDCMTK && isempty(which('dicomdict'))
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
    x = ExploreASL_Initialize('',0); % only initialize ExploreASL if this wasnt initialized before
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
if ~isfield(imPar,'SkipSubjectIfExists') || isempty(imPar.SkipSubjectIfExists)
    % allows to skip existing subject folders in the analysis folder, when this is set to true,
    % avoiding partly re-importing/converting dcm2niiX when processing has been partly done
    imPar.SkipSubjectIfExists = false;
else
    warning('Skipping existing subjects in analysis folder');
    fprintf('If you want to overwrite, first remove the full subject folder');
end

%% Create the directory for analysis
if isfield(imPar,'RawRootModName')
    if strcmp(imPar.RawRootModName,'adaptive')
        % Search for raw data in archive -> Modified for datasets without a hardcoded "raw" folder
        curArchive = dir(fullfile(imPar.RawRoot,imPar.studyID));
        curArchive = curArchive([curArchive.isdir]'); % Only check directories
        remFolders = strcmp({curArchive.name}','.') | strcmp({curArchive.name}','..') | strcmp({curArchive.name}','analysis');
        curArchive = curArchive(~remFolders);
        if numel(curArchive)==1
            imPar.RawRoot = fullfile(imPar.RawRoot,imPar.studyID,curArchive.name);
        else
            % Fallback: Use 'raw'
            imPar.RawRoot = fullfile(imPar.RawRoot,imPar.studyID,'raw');
        end
    else
        % Option to define the raw folder name manually
        imPar.RawRoot = fullfile(imPar.RawRoot,imPar.studyID,imPar.RawRootModName);
    end
else
    % Default/Fallback solution
    imPar.RawRoot = fullfile(imPar.RawRoot,imPar.studyID,'raw');
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
converted_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
skipped_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans
missing_scans = zeros(nSubjects,nVisits,nSessions,nScans,'uint8'); % keep a count of all individual scans

% define a cell array for storing info for parameter summary file
xASL_adm_CreateDir(imPar.AnalysisRoot);
summary_lines = cell(nSubjects,nVisits,nSessions,nScans);

%% -----------------------------------------------------------------------------
% import subject by subject, visit by visit, session by session, scan by scan
% -----------------------------------------------------------------------------
fprintf('%s\n', 'Running import (i.e. dcm2niiX)');
separatorline = repmat(char('+'),1,80);
for iSubject=1:nSubjects
    subjectID = subjectIDs{iSubject};

    for iVisit=1:nVisits
        visitID = visitIDs{iVisit};

        % convert visit ID to a suitable name
        if size(imPar.tokenVisitAliases,2)==2
            iAlias = find(~cellfun(@isempty,regexp(visitIDs{iVisit},imPar.tokenVisitAliases(:,1),'once')));
            if ~isempty(iAlias)
                imPar.visitNames{iVisit} = imPar.tokenVisitAliases{iAlias,2};
            end
        end

        if bUseVisits % only pad VisitID _1 _2 _3 etc if there are visits specified
            % Multiple visits is defined by the tokenVisitAliases.
            % If this is non-existing, it is set to 1, and if it does exist,
            % it will put the _1 _2 _3 etc in the folder
            % this fix allows to import a single visit from a range of
            % specified visits
            SubjDir = fullfile(imPar.AnalysisRoot, [subjectID imPar.visitNames{iVisit}]);
        % if strcmp(imPar.visitNames{iVisit},'_1') % only pad the visitID _1 _2 _3 etc if there are multiple visits
        else
            SubjDir = fullfile(imPar.AnalysisRoot, subjectID);
        end

        if imPar.SkipSubjectIfExists && exist(SubjDir, 'dir')
            continue; % we found the subject dir (i.e. SubjectVisit), so we skip it
            % this is ignored when imPar.SkipSubjectIfExists is set to
            % false (default)
        end

        fprintf('%s\nImporting subject=%s:   \n',separatorline,[subjectID imPar.visitNames{iVisit}]); % display subject-visit ID

        % loop through all sessions
        for iSession=1:nSessions
            sessionID = sessionIDs{iSession};

            % convert session ID to a suitable name
            if size(imPar.tokenSessionAliases,2)==2
                iAlias = find(~cellfun(@isempty,regexp(sessionID,imPar.tokenSessionAliases(:,1),'once')));
                if ~isempty(iAlias)
                    imPar.sessionNames{iSession} = imPar.tokenSessionAliases{iAlias,2};
                end
            end

            for iScan=1:nScans
                scanID = scanIDs{iScan};
                summary_line = [];
                first_match = [];
                summary_lines{iSubject,iVisit,iSession,iScan} = 'n/a';

                if ~imPar.bVerbose % if not verbose, track % progress
                    CounterN = (iSession-1)*nScans+iScan;
                    CounterT = nSessions*nScans;
                    xASL_TrackProgress(CounterN, CounterT);
                end

                % convert scan ID to a suitable name and set scan-specific parameters
                if size(imPar.tokenScanAliases,2)==2
                    iAlias = find(~cellfun(@isempty,regexp(lower(scanID),lower(imPar.tokenScanAliases(:,1)),'once')));
                    if ~isempty(iAlias)
                        scanNames{iScan} = imPar.tokenScanAliases{iAlias,2};
                    else
                        % keep the original name
                        WarningMessage = ['ExploreASL_Import: Unknown scan ID ' scanID ' found, don"t know what this is'];
                        dcm2niiCatchedErrors = CatchErrors('isempty(iAlias)', WarningMessage, dbstack, mfilename, pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar);
                    end
                end
                scan_name = scanNames{iScan};

                % minimalistic feedback of where we are
                if imPar.bVerbose; fprintf('>>> Subject=%s, visit=%s, session=%s, scan=%s\n',subjectID, visitID, num2str(iSession), scan_name); end

                bOneScanIsEnough = false; % default
                bPutInSessionFolder = true; % by default put in session folder
                switch scan_name
                    case {'ASL4D', 'M0', 'ASL4D_RevPE', 'func_bold'}
                        bPutInSessionFolder = true;
                    case {'T1', 'WMH_SEGM', 'FLAIR'}
                        bPutInSessionFolder = false;
                end

                if ~isempty(strfind(char(imPar.folderHierarchy(end)),'PAR'))
                    imPar.dcm2nii_version = '20101105';
                end

                % now pick the matching one from the folder list
                iMatch = find(strcmp(vSubjectIDs,subjectID) & strcmp(vVisitIDs, xASL_adm_CorrectName(visitID,2,'_')) & strcmp(vSessionIDs,sessionID) & strcmp(lower(vScanIDs),lower(scanID)) ); % only get the matching session
                if isempty(iMatch)
                    % only report as missing if we need a scan for each session (i.e. ASL)
                    if sum(converted_scans(iSubject,iVisit,:,iScan))==0
                        WarningMessage = ['Missing scan: ' [subjectID imPar.visitNames{iVisit}] ', ' num2str(iSession) ', ' scan_name];
                        if imPar.bVerbose; warning(WarningMessage); end
                        missing_scans(iSubject, iVisit, iSession, iScan) = 1;
                    end

                    summary_lines{iSubject, iVisit, iSession, iScan} = summary_line;
                    continue;
                    warning('Dont forget to comment continue here for debugging');
                end

                % determine input and output paths
                bSkipThisOne = false;
                branch = matches{iMatch};
                scanpath = fullfile(imPar.RawRoot,branch);

                if ~isempty(strfind(scanNames{iScan}, 'ASL4D')) || ~isempty(strfind(scanNames{iScan}, 'M0'))
                    session_name = ['ASL_' num2str(iSession)];
                elseif ~isempty(strfind(scanNames{iScan}, 'DSC4D'))
                     session_name = ['DSC_' num2str(iSession)];
                else
                    session_name = [scanNames{iScan} '_' num2str(iSession)]; % Allow multiple ScanTypes for sessions
                end

                if bPutInSessionFolder
                    destdir = fullfile(SubjDir, session_name);
                else % put in subject folder instead of session folder
                    destdir = SubjDir;
                end

                if bOneScanIsEnough && sum(converted_scans(iSubject,iVisit,:,iScan))~=0
                    % one scan is enough, so skip this one if there was already a scan converted of this type (i.e. T1)
                    if imPar.bVerbose; fprintf('Skipping scan: %s, %s, %s\n',[subjectID imPar.visitNames{iVisit}],session_name,scan_name); end
                    bSkipThisOne = true;
                    destdir = []; % just in case
                end

                % start the conversion if this scan should not be skipped
                if bSkipThisOne
                    summary_line = sprintf(',"skipped",,,,,,,,');
                    skipped_scans(iSubject, iVisit, iSession, iScan) = 1;
                else
                    nii_files = {};
                    xASL_adm_CreateDir(destdir);

                    % check if we have a nii(gz) file, or something that needs to be converted (parrec/dicom)
                    if ~exist(scanpath, 'dir') && ~isempty(regexp(lower(scanpath),'(\.nii|\.nii\.gz)$'))
                        % we found a NIfTI file
                        % check if output exists
                        first_match = fullfile(destdir, [scan_name '.nii']);
                        if imPar.bOverwrite || ~xASL_exist(first_match,'file')
                            [~, fname, fext] = fileparts(scanpath);
                            destfile = fullfile(destdir, [fname fext]); % will be renamed later
                            xASL_Copy(scanpath, destfile, imPar.bOverwrite, imPar.bVerbose);
                            % gunzip if required
                            destfile = xASL_adm_UnzipNifti(destfile);
                            xASL_Move(destfile, first_match, imPar.bOverwrite, imPar.bVerbose);
                        end
                        nii_files{1} = first_match;
                    else % we found dicom files
                        % -----------------------------------------------------------------------------
                        % start the conversion. Note that the dicom filter is only in effect when a directory is specified as input.
                        % -----------------------------------------------------------------------------
                        if bRunDCM2NII
                            try
                                [nii_files, scan_name, first_match, MsgDcm2nii] = xASL_io_dcm2nii(scanpath, destdir, scan_name, 'DicomFilter', imPar.dcmExtFilter, 'Verbose', imPar.bVerbose, 'Overwrite', imPar.bOverwrite, 'Version', imPar.dcm2nii_version, 'x', x);

                                % If dcm2nii produced a warning or error, catch this & store it
                                if ~isempty(MsgDcm2nii) && ~isempty(regexp(lower(MsgDcm2nii),'.*(warning|error).*')) % if it contains a warning/error
                                    dcm2niiCatchedErrors = CatchErrors('xASL_io_dcm2nii', MsgDcm2nii, dbstack, ['dcm2nii_' imPar.dcm2nii_version], pwd, scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar);
                                end

                            catch ME
                                dcm2niiCatchedErrors = CatchErrors(ME.identifier, ME.message, [], [], [], scan_name, scanpath, destdir, dcm2niiCatchedErrors, imPar, ME.stack);

                                if imPar.bVerbose; warning(['dcm2nii ' scanpath ' crashed, skipping']); end
                                if imPar.bVerbose; warning('Check whether the scan is complete'); end
                                first_match = xASL_adm_GetFileList(scanpath, ['.*' imPar.dcmExtFilter],'FPList',[0 Inf]);
                                if  ~isempty(first_match); first_match = first_match{1}; end
                            end
                        else
                            first_match = xASL_adm_GetFileList(scanpath, ['.*' imPar.dcmExtFilter],'FPList',[0 Inf]);
                            if  ~isempty(first_match); first_match = first_match{1}; end
                        end
                    end

                    %% Merge NIfTIs if there are multiple for ASL only
                    % check the number of created nifiti files in case of ASL: control and label should be merged as one 4D
                    if length(nii_files)>1 && ~isempty(strfind(scan_name,'ASL4D'))
                        nii_files = merge_2_ASL_nii_files(nii_files, scan_name);
                    end

                    % Extract relevant parameters from nifti header and append to summary file
                    summary_line = AppendNiftiParameters(nii_files);
                    converted_scans(iSubject,iSession,iScan) = 1;
                end

                % extract relevant parameters from dicom header, if not
                % already exists
                % Find JSONpath that is there already
                SavePathJSON = {};
                SavePathJSON{1} = fullfile(destdir, [scan_name '.json']);
                SavePathJSON{2} = fullfile(destdir, [session_name '.json']);
                for iPath=1:length(nii_files)
                    % now we add the path only if it didnt exist already in this list
					tmpNewPath = [nii_files{iPath}(1:end-4) '.json'];
                    if ~max(cellfun(@(y) strcmp(y, tmpNewPath), SavePathJSON))
                        SavePathJSON{end+1} = tmpNewPath;
                    end
                end

                for iPath=1:length(SavePathJSON)
                    if exist(SavePathJSON{iPath}, 'file') && ~isempty(first_match)
                        [~, ~, fext] = fileparts(first_match);
                        if  strcmpi(fext,'.PAR')
                            parms = xASL_adm_Par2Parms(first_match, SavePathJSON{iPath});
                        elseif strcmpi(fext,'.nii')
                            parms = [];
                        elseif imPar.bMatchDirectories
                            Fpath  = fileparts(first_match);
                            [parms, pathDcmDict] = xASL_bids_Dicom2JSON(imPar, Fpath, SavePathJSON{iPath}, imPar.dcmExtFilter, bUseDCMTK, pathDcmDict);
                            clear Fpath Ffile Fext
                        else
                            [parms, pathDcmDict] = xASL_bids_Dicom2JSON(imPar, first_match, SavePathJSON{iPath}, imPar.dcmExtFilter, bUseDCMTK, pathDcmDict);
                        end
                    end
                end

                % correct nifti rescale slope if parms.RescaleSlopeOriginal =~1
                % but nii.dat.scl_slope==1 (this can happen in case of
                % hidden scale slopes in private Philips header,
                % that is dealt with by xASL_bids_Dicom2JSON but not by
                % dcm2niiX

                if ~isempty(nii_files) && exist('parms','var')
                    [TempLine, PrintDICOMFields] = AppendParmsParameters(parms);
                    summary_line = [summary_line TempLine];
                end

                if bClone2Source % make a copy of analysisdir in sourcedir
                    if ~isempty(nii_files)
                        for iFile=1:length(nii_files)
                            % replace 'analysis' by 'source'
                            [iStart, iEnd] = regexp(nii_files{iFile}, 'analysis');
                            DestPath = [nii_files{iFile}(1:iStart-1) 'source' nii_files{iFile}(iEnd+1:end)];
                            xASL_Copy(nii_files{iFile}, DestPath, true);
                            % do the same for other extensions
                            Extensions = {'.json' '_parms.json' '_parms.mat'};
                            for iExt=1:length(Extensions)
                                [Fpath, Ffile] = xASL_fileparts(nii_files{iFile});
                                CopyPath = fullfile(Fpath, [Ffile Extensions{iExt}]);
                                [Fpath, Ffile] = xASL_fileparts(DestPath);
                                DestPath = fullfile(Fpath, [Ffile Extensions{iExt}]);
                                if xASL_exist(CopyPath)
                                    xASL_Copy(CopyPath, DestPath, true);
                                end
                            end
                        end
                    end
                end

                % Copy single dicom as QC placeholder
                if bCopySingleDicoms && ~isempty(first_match)
                    xASL_Copy(first_match, fullfile(destdir, ['DummyDicom_' scan_name '.dcm']), imPar.bOverwrite, imPar.bVerbose);
                end

                % store the summary info so it can be sorted and printed below
                summary_lines{iSubject, iVisit, iSession, iScan} = summary_line;
            end % scansIDs
        end % sessionIDs
    end % visitIDs
end % subjectIDs




%% -----------------------------------------------------------------------------
% create summary file
% -----------------------------------------------------------------------------
summary_filepath = fullfile(imPar.AnalysisRoot, 'import_summary.csv');
fid_summary = fopen(summary_filepath,'wt');
% Print headers for parameters obtained from NIfTI file
fprintf(fid_summary,'subject,visit,session,scan,filename,dx,dy,dz,dt,nx,ny,nz,nt');
% Print headers for parameters obtained from DICOM file
if exist('PrintDICOMFields','var')
    for iField=1:length(PrintDICOMFields)
        fprintf(fid_summary,[',' PrintDICOMFields{iField}]);
    end
end
fprintf(fid_summary,'\n');

for iScan=1:nScans
    for iSubject=1:nSubjects
        for iVisit=1:nVisits
            for iSession=1:nSessions
                if converted_scans(iSubject, iVisit, iSession, iScan) || skipped_scans(iSubject, iVisit, iSession, iScan) || missing_scans(iSubject, iVisit, iSession, iScan)
                    fprintf(fid_summary,'"%s","%s","%s","%s"%s,\n', subjectIDs{iSubject}, visitIDs{iVisit}, imPar.sessionNames{iSession}, scanNames{iScan}, summary_lines{iSubject, iVisit, iSession, iScan});
                end
            end
        end
    end
end
fprintf(fid_summary,'\n');

nMissing = sum(missing_scans(:));
nSkipped = sum(skipped_scans(:));

% report totals
% header first
fprintf(fid_summary,'\n');
fprintf(fid_summary,'\nSubject,nConverted');
fprintf(fid_summary,[',MissingScans (n=' num2str(nMissing) ')']);
fprintf(fid_summary,[',SkippedScans (n=' num2str(nSkipped) ')\n']);
% then subjects row-by-row
for iSubject=1:nSubjects
    for iVisit=1:nVisits
        fprintf(fid_summary,'"%s"', [subjectIDs{iSubject} visitIDs{iVisit}]);
        fprintf(fid_summary,',%d',sum(converted_scans(iSubject,:,:,:)));

        for iSession=1:nSessions
            fprintf(fid_summary,',"');
            fprintf(fid_summary,'%s ',scanNames{logical(missing_scans(iSubject, iVisit, iSession,:))});
            fprintf(fid_summary,'"');
        end

        for iSession=1:nSessions
            fprintf(fid_summary,',"');
            fprintf(fid_summary,'%s ',scanNames{logical(skipped_scans(iSubject, iVisit, iSession,:))});
            fprintf(fid_summary,'"');
        end

        fprintf(fid_summary,'\n');
    end
end

% and a grand total of missing and skipped
if nMissing>0
    fprintf(2,'Number of missing scans: %d\n',nMissing);
end
if nSkipped>0
    fprintf(2,'Number of skipped scans: %d\n',nSkipped);
end
fclose(fid_summary);

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
    xASL_adm_SaveJSON(dcm2niiCatchedErrors, SaveJSON);
end

fprintf('\n');

end

% -----------------------------------------------------------------------------
%
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
%
% -----------------------------------------------------------------------------
function [s, FieldNames] = AppendParmsParameters(parms)
% This function outputs s=fields of _parms.mat
s = [];

FieldNames = {'RepetitionTime' 'EchoTime' 'NumberOfAverages' 'RescaleSlope' 'RescaleSlopeOriginal'...
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
%
% -----------------------------------------------------------------------------
function nii_files = merge_2_ASL_nii_files(nii_files, basename)
% merge_2_ASL_nii_files
% CAVE: deletes original files

if length(nii_files)>1
    fprintf('Warning EXPLOREASL_IMPORT: concatenating multiple NIfTIs & jsons as output from dcm2niiX\n');
    
    % First rename the NIfTI and JSON files to 4 digit format & sort them
    % this avoids 1 10 2 issues
    for iFile=1:length(nii_files)
        [Fpath, Ffile, Fext] = xASL_fileparts(nii_files{iFile});
        [iStart, iEnd] = regexp(Ffile,'\d*$');
        FfileNew = [Ffile(1:iStart-1) sprintf('%04d', str2num(Ffile(iStart:iEnd)))];
        PathNew = fullfile(Fpath, [FfileNew Fext]);
		xASL_Move(nii_files{iFile}, PathNew);
		nii_files{iFile} = PathNew;
		
		% Rename also the JSON if it exists
		PathOldJSON = fullfile(Fpath,[Ffile '.json']);
		PathNewJSON = fullfile(Fpath, [FfileNew '.json']);
		if exist(PathOldJSON,'file')
			xASL_Move(PathOldJSON,PathNewJSON);
		end
    end
    
    nii_files = sort(nii_files);
    
    % Then create image
    % Here we issue a warning for dimensionality>4
	countJSON = 1;
    for iFile=1:length(nii_files)
        tempIM = xASL_io_Nifti2Im(nii_files{iFile});
        if length(size(tempIM))>4
            error('Dimensionality incorrect for this ASL NIfTI file');
		end
		
		bFileOK = 1;
		if iFile==1
			IM = tempIM;
		else
			% Check the file sizes and merge only if of similar size
			sizeFirst = size(IM);
			sizeFirst = sizeFirst(1:3);
			sizeNew   = size(tempIM);
			sizeNew   = sizeNew(1:3);
			if isequal(sizeNew, sizeFirst)
				IM(:,:,:,end+1:end+size(tempIM,4)) = tempIM;
				% Immediately delete the file
				xASL_delete(nii_files{iFile});
			else
				bFileOK = 0;
			end
		end
		% If file was merged, then proceed to delete
		if bFileOK
			[tempFpath, tempFfile, ~] = xASL_fileparts(nii_files{iFile});
			PathOldJSON = fullfile(tempFpath,[tempFfile '.json']);
			if countJSON == 1
				if exist(PathOldJSON,'file')
					% Save the first JSON as ASL4D.json
					countJSON = countJSON + 1;
					xASL_Move(PathOldJSON, fullfile(Fpath, 'ASL4D.json'),1);
				end
			else
				% Delete all the other JSONs
				if exist(PathOldJSON,'file')
					xASL_delete(PathOldJSON);
				end
			end
		end
		
	end
    
	% Save the file with the same header and remove the last remaining file
    NewNIfTI = fullfile(fileparts(nii_files{1}), 'ASL4D.nii');
    xASL_io_SaveNifti(nii_files{1}, NewNIfTI, IM, [], 0);
	xASL_delete(nii_files{1});
    
    fprintf('Corrected dcm2niiX output for\n');
    fprintf('%s\n', NewNIfTI);
    
    nii_files = {NewNIfTI};
end

end


% -----------------------------------------------------------------------------
%
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
