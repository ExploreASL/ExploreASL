function [x] = xASL_bids_BIDS2Legacy(pathStudy, x, bOverwrite)
%xASL_bids_BIDS2Legacy Convert BIDS rawdata to ExploreASL legacy format
%
% FORMAT: [x] = xASL_bids_BIDS2Legacy(pathStudy, x[, bOverwrite])
% 
% INPUT:
%   pathStudy  - path to the study folder containing the BIDS data in rawdata subfolder (REQUIRED)
%   x          - ExploreASL x structure (REQUIRED, STRUCT)
%   bOverwrite - boolean, true for overwriting files (OPTIONAL, DEFAULT = true)
%   
% OUTPUT: 
%   x          - ExploreASL x structure (STRUCT)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function converts BIDS rawdata (in pathStudy/rawdata/) 
% to xASL legacy derivative format (e.g. pathStudy/derivatives/ExploreASL/)
%
% Can be updated step-by-step when ExploreASL's derivative structure moves to BIDS
% NB: ask how Visits/session layer is defined in bids-matlab (should be
% separate layer within subjects, but now isn't?)
%
% This function performs the following steps:
%
% 1. Parse a folder using bids-matlab
% 2. Define SubjectSession
% 3. Define Session
% 4. Parse modality
% - Parse scantype
% - Compile paths for copying
% - Manage sidecars to copy
% - Copy files
% 5. Parse M0
% 6. Clean up
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [x] = xASL_bids_BIDS2Legacy(pathStudy, x, bOverwrite);
% __________________________________
% Copyright 2015-2021 ExploreASL


%% 0. Admin

% Verify the input parameters
if nargin<1 || isempty(pathStudy)
	error('pathStudy is a required parameter.');
end

if nargin<3 || isempty(bOverwrite)
    bOverwrite = 1;
end

% Verify that the rawdata subfolder exists
if ~exist(fullfile(pathStudy,'rawdata'),'dir')
    warning('Invalid folder selected, not containing rawdata folder');
    return;
end

% Creates the derivatives directory
pathLegacy = fullfile(pathStudy, 'derivatives', 'ExploreASL');
if exist(pathLegacy, 'dir') && bOverwrite
    fprintf('The derivatives directory already exists, overwriting...\n');
elseif exist(pathLegacy, 'dir')
    fprintf('%s\n', [pathLegacy ' exists, merging']);
else
    xASL_adm_CreateDir(pathLegacy);
end

% Loads the configuration for file renaming from the BIDS configuration file
bidsPar = xASL_bids_Config();


%% 1. Parse a folder using bids-matlab
BIDS = bids.layout(fullfile(pathStudy,'rawdata'));
nSubjects = numel(BIDS.subjectName);
nVisits = numel(BIDS.sessionName); % this is called sessions in BIDS
% we use this below to see if the legacy subjectname gets _1 as visit suffix or not
[~, studyName] = fileparts(pathStudy);
fprintf('Converting from BIDS to Legacy: %s   \n', studyName);


%% 2. Define SubjectSession
for iSubjSess=1:numel(BIDS.subjects) 
    % Iterate over BIDS.subjects (indices that include both subjects & sessions)
    % so  1 subject  6 session/visits, will give numel(BIDS.subjects)=6
    % and 6 subjects 1 session/visits, will give numel(BIDS.subjects)=6
    xASL_TrackProgress(iSubjSess, nSubjects);
    
    % Subject ID
    SubjectID = BIDS.subjects(iSubjSess).name;
    
    % Subject-wise processing (because xASL_Iteration runs over subjects)
    if ~isempty(regexpi(SubjectID,x.SUBJECT))
    
        % Session ID (Currently, ExploreASL concatenates subject_visit/timepoint in the same folder layer, so we only use SubjectSession)
        SessionID = BIDS.subjects(iSubjSess).session;


        %% 3. Define Session
        iVisit = find(strcmp(BIDS.sessionName, SessionID));
        % remove iteration for iVisit=1 % iterate visit/session in this "BIDS.subjects" (always 1 session per BIDS.subjects)
        % ExploreASL uses visit as a number (e.g. _1 _2 _3 etc)
        if nVisits==1
            pathLegacy_SubjectVisit = fullfile(pathLegacy, SubjectID);
            VisitString = '';
        else
            if isempty(iVisit)
                fprintf('\nEmpty session number, setting session number to 1...\n');
                iVisit = 1;
            end

            pathLegacy_SubjectVisit = fullfile(pathLegacy, [SubjectID '_' xASL_num2str(iVisit)]);
            VisitString = [' visit ' SessionID];
        end
        SubjectVisit = [SubjectID VisitString];
        xASL_adm_CreateDir(pathLegacy_SubjectVisit);


        %% 4. Parse modality
        % Modalities - the BIDS scantypes
        ModalitiesUnique = unique(bidsPar.BIDS2LegacyFolderConfiguration(2:end, 2));
        nModalities = length(ModalitiesUnique);
        xASL_bids_BIDS2Legacy_ParseModality(BIDS, bidsPar, SubjectVisit, iSubjSess, ModalitiesUnique, nModalities, bOverwrite, pathLegacy_SubjectVisit);

    end

end
% Print new line after track progress bar
fprintf('   \n');


% Get directories of current subject
SubjectDirs = xASL_adm_GetFileList(fullfile(pathLegacy), ['^.+' x.SUBJECT '.+$'], [],[],true);


%% 5. Parse M0 of current subject
ListASL4D = [];
if ~isempty(SubjectDirs)
    for iDir=1:numel(SubjectDirs)
        currentASL = xASL_adm_GetFileList(SubjectDirs{iDir}, '^ASL4D\.nii$', 'FPListRec');
        ListASL4D = vertcat(ListASL4D,currentASL);
    end
end
% Parse M0s
if ~isempty(ListASL4D)
    for iList=1:numel(ListASL4D)
        xASL_bids_parseM0(ListASL4D{iList});
        [~, currentNifti] = fileparts(ListASL4D{iList});
        fprintf('M0 parsed for subject %s image %s ...\n', SubjectID, currentNifti);
    end
else
    warning('No ASL4D file found in %s...', pathLegacy);
end


%% 6. Clean up
xASL_imp_BIDS2Legacy_CleanUp(pathStudy);


end



%% Clean-Up subfunction
function xASL_imp_BIDS2Legacy_CleanUp(pathStudy)

    % Start with empty file list
    filesCleanUp = {};
    
    % Search for dcm2niix summary file
    dcm2niixImportLogFile = xASL_adm_GetFileList(pathStudy,'^import_log.+$');
    
    % Merge log file lists
    if ~isempty(dcm2niixImportLogFile)
        filesCleanUp = vertcat(filesCleanUp,dcm2niixImportLogFile);
    end
    
    % Move files to derivatives
    if ~isempty(filesCleanUp)
        for iFile = 1:size(filesCleanUp,1)
            sourceCleanUp = filesCleanUp{iFile};
            [~, fileCleanUp, extCleanUp] = xASL_fileparts(sourceCleanUp);
            destCleanUp = fullfile(pathStudy, 'derivatives', 'ExploreASL', [fileCleanUp extCleanUp]);
            if xASL_exist(sourceCleanUp,'file')
                xASL_Move(sourceCleanUp,destCleanUp);
            end
        end
    end
    
    % Print new line after track progress bar
    fprintf('   \n');

end





