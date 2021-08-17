function [dataPar] = xASL_bids_BIDS2Legacy(pathStudy, x, bOverwrite, dataPar)
%xASL_bids_BIDS2Legacy Convert BIDS rawdata to ExploreASL legacy format
%
% FORMAT: [dataPar] = xASL_bids_BIDS2Legacy(pathStudy, x[, bOverwrite, dataPar])
% 
% INPUT:
%   pathStudy  - path to the study folder containing the BIDS data in rawdata subfolder (REQUIRED)
%   x          - ExploreASL x structure (REQUIRED, STRUCT)
%   bOverwrite - boolean, true for overwriting files (OPTIONAL, DEFAULT = true)
%   dataPar    - dataPar values to be filled in a basic dataPar created with the conversion to legacy 
%                (OPTIONAL, DEFAULT = basic dataPar settings)
%   
% OUTPUT: dataPar - dataPar struct (STRUCT)
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
% 2. Define Subject
% 3. Define SubjectVisit
% 4. Parse modality
% - Parse scantype
% - Compile paths for copying
% - Manage sidecars to copy
% - Copy files
% 5. Parse M0
% 6. Create DataPar.json
% 7. Copy participants.tsv
% 8. Add dataset_description.json
% 9. Clean up
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [dataPar] = xASL_bids_BIDS2Legacy(pathStudy, x, bOverwrite, dataPar);
% __________________________________
% Copyright 2015-2021 ExploreASL


%% 0. Admin

% Verify the input parameters
if nargin<1 || isempty(pathStudy)
	error('pathStudy is a required parameter.');
end

if nargin<2 || isempty(bOverwrite)
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
    warning([pathLegacy ' already exists, overwriting']);
elseif exist(pathLegacy, 'dir')
    fprintf('%s\n', [pathLegacy ' exists, merging']);
else
    xASL_adm_CreateDir(pathLegacy);
end

% Creates a default dataPar 
if nargin < 3 || isempty(dataPar)
	dataPar = struct();
end

% Fills in important information in the dataPar if missing
if ~isfield(dataPar,'x')
    % Add x field
    dataPar.x = struct;
end
% Dataset fields
if ~isfield(dataPar.x,'dataset')
    dataPar.x.dataset = struct;
end
% Add default subject regular expression
if ~isfield(dataPar.x.dataset,'subjectRegexp')
    dataPar.x.dataset.subjectRegexp = '^sub-.*$';
end

% Check for settings fields
if ~isfield(dataPar.x,'settings')
    dataPar.x.settings = struct;
end

% Check for quality field
if ~isfield(dataPar.x.settings,'Quality')
    dataPar.x.settings.Quality = 1;
end

% Check for DELETETEMP field
if ~isfield(dataPar.x.settings,'DELETETEMP')
	dataPar.x.settings.DELETETEMP = 1;
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


%% 2. Define Subject
for iSubjSess=1:numel(BIDS.subjects) % iterate over BIDS.subjects (indices that include both subjects & sessions)
    % so  1 subject  6 session/visits, will give numel(BIDS.subjects)=6
    % and 6 subjects 1 session/visits, will give numel(BIDS.subjects)=6
    xASL_TrackProgress(iSubjSess, nSubjects);
    SubjectID = BIDS.subjects(iSubjSess).name;
    SessionID = BIDS.subjects(iSubjSess).session;
    % Currently, ExploreASL concatenates subject_visit/timepoint in the
    % same folder layer, so we only use SubjectSession
    
    
    %% 3. Define SubjectVisit
    iVisit = find(strcmp(BIDS.sessionName, SessionID));
    % remove iteration for iVisit=1 % iterate visit/session in this "BIDS.subjects" (always 1 session per BIDS.subjects)
    % ExploreASL uses visit as a number (e.g. _1 _2 _3 etc)
    if nVisits==1
        pathLegacy_SubjectVisit = fullfile(pathLegacy, SubjectID);
        VisitString = '';
    else
        pathLegacy_SubjectVisit = fullfile(pathLegacy, [SubjectID '_' xASL_num2str(iVisit)]);
        VisitString = [' visit ' SessionID];
    end
    SubjectVisit = [SubjectID VisitString];
    xASL_adm_CreateDir(pathLegacy_SubjectVisit);
    
    
    %% 4. Parse modality
    % Modalities - the BIDS domains of scantypes
    ModalitiesUnique = unique(bidsPar.BIDS2LegacyFolderConfiguration(2:end, 2));
    nModalities = length(ModalitiesUnique);
    xASL_bids_BIDS2Legacy_ParseModality(BIDS, bidsPar, SubjectVisit, iSubjSess, ModalitiesUnique, nModalities, bOverwrite, pathLegacy_SubjectVisit);
        
end

fprintf('\n');

%% 5. Parse M0
ListASL4D = xASL_adm_GetFileList(pathLegacy, '^ASL4D\.nii$', 'FPListRec');
if ~isempty(ListASL4D)
    for iList=1:numel(ListASL4D)
        xASL_bids_parseM0(ListASL4D{iList});
        [~, currentNifti] = fileparts(ListASL4D{iList});
        fprintf('M0 parsed for %s %s ...\n', studyName, currentNifti);
    end
else
    warning(['No ASL4D file found in ' pathLegacy]);
end

%% 6. Create dataPar.json

% Write DataParFile if it does not exist already
if ~xASL_exist(fullfile(pathLegacy, 'dataPar.json'),'file')
    fprintf('Creating default dataPar.json since file was not found...\n');
    spm_jsonwrite(fullfile(pathLegacy, 'dataPar.json'), dataPar);
end

% Overwrite dataPar.json in x structure
fprintf('Overwriting x.dir.dataPar...\n');

% Add the path to the dataPar.x struct that we return to the Master script
dataPar.x.dir.dataPar = fullfile(pathLegacy, 'dataPar.json');

%% 7. Copy participants.tsv

% Define file names
pathParticipantsTSV = fullfile(pathStudy, 'participants.tsv');
pathParticipantsTSVxASL = fullfile(pathLegacy, 'participants.tsv');

% Backwards compatibility: rename to lowercase
fileListParticipantsTSVold = xASL_adm_GetFileList(pathStudy,'^Participants.tsv$',false);
if ~isempty(fileListParticipantsTSVold) % We added the _temp copy step so that the code works on case insensitive systems like windows as well. Please don't remove that step for backwards compatibility (at least not until release 2.0.0).
    pathParticipantsTSVold = fullfile(pathStudy, 'Participants.tsv');
    pathParticipantsTSVoldTemp = fullfile(pathStudy, 'Participants_temp.tsv');
    xASL_Move(pathParticipantsTSVold,pathParticipantsTSVoldTemp);
    xASL_Move(pathParticipantsTSVoldTemp,pathParticipantsTSV);
end

% Check if participants.tsv exists & copy it to the derivatives
if xASL_exist(pathParticipantsTSV,'file')
    xASL_Copy(pathParticipantsTSV,pathParticipantsTSVxASL);
end

%% 8. Add "GeneratedBy" fields
try
    % Copy dataset_description JSON file
    xASL_Copy(fullfile(pathStudy, 'rawdata', 'dataset_description.json'), fullfile(pathLegacy, 'dataset_description.json'));

    % Get all JSON files in derivatives and add the "GeneratedBy" field
    jsonFiles = xASL_adm_GetFileList(pathLegacy, '.json$', 'FPListRec');

    % Check file list
    if size(jsonFiles,1)>0
        % Iterate over files
        for iFile = 1:size(jsonFiles,1)
            % Check if file should be excluded first (exclude participants.json)
            [~,fileName] = xASL_fileparts(jsonFiles{iFile});
            if ~strcmpi(fileName,'participants')
                thisPath = jsonFiles{iFile};
                xASL_bids_AddGeneratedByField(x, thisPath);
            end
        end
    end
catch ME
    warning('Adding the GeneratedBy fields failed...');
    fprintf('%s\n', ME.message);
end

%% 9. Clean up
try
    filesCleanUp = xASL_adm_GetFileList(pathStudy,'^import_.+$');
    reportFiles = xASL_adm_GetFileList(pathStudy,'^bidsReport.+$');
    if ~isempty(reportFiles)
        filesCleanUp = vertcat(filesCleanUp,reportFiles);
    end
    if ~isempty(filesCleanUp)
        for iFile = 1:size(filesCleanUp,1)
            sourceCleanUp = filesCleanUp{iFile};
            [~, fileCleanUp, extCleanUp] = xASL_fileparts(sourceCleanUp);
            destCleanUp = fullfile(pathStudy, 'derivatives', 'ExploreASL', [fileCleanUp extCleanUp]);
            xASL_Move(sourceCleanUp,destCleanUp);
        end
    end
catch ME
    warning('Clean up failed...');
    fprintf('%s\n', ME.message);
end

end



