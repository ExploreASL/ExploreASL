function [dataPar] = xASL_bids_BIDS2Legacy(pathStudy, x, bOverwrite, dataPar)
%xASL_bids_BIDS2Legacy Convert BIDS rawdata to ExploreASL legacy format
%
% FORMAT: xASL_bids_BIDS2Legacy(pathStudy[, bOverwrite, dataPar])
% 
% INPUT:
%   pathStudy  - path to the study folder containing the BIDS data in rawdata subfolder (REQUIRED)
%   x          - ExploreASL x structure (REQUIRED, STRUCT)
%   bOverwrite - boolean, true for overwriting files (OPTIONAL, DEFAULT = true)
%   dataPar    - dataPar values to be filled in a basic dataPar created with the conversion to legacy 
%                (OPTIONAL, DEFAULT = basic dataPar settings)
%   
% OUTPUT: n/a
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
% 5. Parse scantype
% 6. Compile paths for copying
% 7. Manage sidecars to copy
% 8. Copy files
% 9. Parse M0
% 10. Create DataPar.json
% 11. Copy participants.tsv
% 12. Add dataset_description.json
% 13. Clean up
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_bids_BIDS2xASL('pathMyStudy')
% __________________________________
% Copyright 2015-2021 ExploreASL


%% ------------------------------------------------------------------------------------
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
% Add default subject regular expression
if ~isfield(dataPar.x,'subject_regexp')
    dataPar.x.subject_regexp = '^sub-.*$';
end

if ~isfield(dataPar.x,'Quality')
	dataPar.x.Quality = 1;
end

if ~isfield(dataPar.x,'DELETETEMP')
	dataPar.x.DELETETEMP = 1;
end

% Loads the configuration for file renaming from the BIDS configuration file
bidsPar = xASL_bids_Config();

%% ------------------------------------------------------------------------------------
%% 1. Parse a folder using bids-matlab
BIDS = bids.layout(fullfile(pathStudy,'rawdata'));
nSubjects = numel(BIDS.subjectName);
nVisits = numel(BIDS.sessionName); % this is called sessions in BIDS
% we use this below to see if the legacy subjectname gets _1 as visit suffix or not

fprintf('%s\n', ['Converting from BIDS to Legacy: ' pathStudy]);

%% -----------------------------------------------
%% 2. Define Subject
for iSubjSess=1:numel(BIDS.subjects) % iterate over BIDS.subjects (indices that include both subjects & sessions)
    % so  1 subject  6 session/visits, will give numel(BIDS.subjects)=6
    % and 6 subjects 1 session/visits, will give numel(BIDS.subjects)=6
    xASL_TrackProgress(iSubjSess, nSubjects);
    SubjectID = BIDS.subjects(iSubjSess).name;
    SessionID = BIDS.subjects(iSubjSess).session;
    % Currently, ExploreASL concatenates subject_visit/timepoint in the
    % same folder layer, so we only use SubjectSession

    %% -----------------------------------------------
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

        %% -----------------------------------------------
        %% 4. Parse modality
		% Modalities - the BIDS domains of scantypes
        ModalitiesUnique = unique(bidsPar.BIDS2LegacyFolderConfiguration(2:end, 2));
        nModalities = length(ModalitiesUnique);
        for iModality=1:nModalities % iterate modalities in this Subject/Visit
            ModalityIs = ModalitiesUnique{iModality};
            if isfield(BIDS.subjects(iSubjSess), ModalityIs) && ~isempty(BIDS.subjects(iSubjSess).(ModalityIs))
                ModalityFields = BIDS.subjects(iSubjSess).(ModalityIs);
                nScans = length(ModalityFields);

                % Parse fields of this modality for combinations ScanType & Run, in a reference table
                Reference = {'index' 'ScanType' 'run'};
                for iScan=1:nScans % iterate NIfTIs in this Subject/Visit/Modality
                    Reference{iScan+1, 1} = iScan; % index
                    Reference{iScan+1, 2} = ModalityFields(iScan).type; % ScanType
                    Reference{iScan+1, 3} = xASL_str2num(ModalityFields(iScan).run); % run
                    if isempty(Reference{iScan+1, 3}) || isnan(Reference{iScan+1, 3})
                        Reference{iScan+1, 3} = 1; % default to 1st run
                    end
                end
                Reference(2:end,:) = sortrows(Reference(2:end,:), [2, 3]); % first sort for ScanType then run

                RunsAre = cellfun(@(y) y, Reference(2:end, 3));
                RunsUnique = unique(RunsAre);
                
                %% -----------------------------------------------
                %% 5. Parse scantype
                modalityIndices = find(strcmp(bidsPar.BIDS2LegacyFolderConfiguration(2:end,2), ModalityIs));
                modalityConfiguration = bidsPar.BIDS2LegacyFolderConfiguration(1+modalityIndices, :);
                
                for iType=1:size(modalityConfiguration, 1) % iterate scantypes in this Subject/Visit/Modality
                    TypeIs = modalityConfiguration{iType,1};
                    TypeIndices = cellfun(@(y) strcmp(y, TypeIs), Reference(2:end, 2)); % this are the indices for this ScanType

                    if ~isempty(TypeIndices)
                        for iRun=1:length(RunsUnique) % iterate runs in this Subject/Visit/Modality
                            RunIs = RunsUnique(iRun);
                            RunIndices = RunsAre==RunsUnique(iRun);
                            TypeRunIndex = find(RunIndices & TypeIndices);
                            if length(TypeRunIndex)>1
                                warning(['Multiple NIfTIs found for ' SubjectVisit '_run-' xASL_num2str(RunIs) '_' TypeIs ', using first only']);
                                TypeRunIndex = TypeRunIndex(1);
                            end

                            %% -----------------------------------------------
                            %% 6. Compile paths for copying                          
                            if length(TypeRunIndex)==1 % if this scantype-run combination exists

                                % ModalityIs = current modality (e.g. 'anat' 'perf')
                                % TypeIs = current scantype, e.g. 'asl' 'm0' 't1w'
                                % RunIs = current run (e.g. 1, 2, 3)
                                % TypeRunIndex = index of current scantype & run inside the above created Reference Table
                                
                                % Define folder & filename
                                ConfigIndex = cellfun(@(y) strcmp(TypeIs, y), bidsPar.BIDS2LegacyFolderConfiguration(2:end, 1)); % find scantype index
                                ConfigIndex2 = cellfun(@(y) strcmp(ModalityIs, y), bidsPar.BIDS2LegacyFolderConfiguration(2:end, 2)); % find modality index
                                ConfigIndex = find(ConfigIndex & ConfigIndex2); % Combine them & convert to index
                                
								% xASL legacy subfolder name (empty means no subfolder)
                                FolderIs = bidsPar.BIDS2LegacyFolderConfiguration{ConfigIndex+1,3};
								
								% xASL legacy filename
                                FileIs = bidsPar.BIDS2LegacyFolderConfiguration{ConfigIndex+1,4};
                                
                                % Manage runs
                                PrintRun = false;
								% if xASL legacy requires to specify the first run (e.g. T1.nii.gz vs ASL_1/ASL4D.nii.gz)
                                if RunIs==1 && bidsPar.BIDS2LegacyFolderConfiguration{ConfigIndex+1, 6}==1
                                    % if we need to print first run
                                    PrintRun = true;
                                elseif RunIs>1
                                    PrintRun = true;
                                end
                                
                                if PrintRun==1
									% xASL legacy location of run specification (e.g. T1_2.nii.gz vs ASL_2/ASL4D.nii.gz for file vs folder location)
                                    if strcmp(bidsPar.BIDS2LegacyFolderConfiguration{ConfigIndex+1,5}, 'file')
                                        FileIs = [FileIs '_' xASL_num2str(RunIs)];
                                    elseif strcmp(bidsPar.BIDS2LegacyFolderConfiguration{ConfigIndex+1,5}, 'folder')
                                        FolderIs = [FolderIs '_' xASL_num2str(RunIs)];
                                    end
                                end
                                
                                pathOrig = '';
                                pathDest = '';
                                
                                pathOrig{1} = fullfile(BIDS.subjects(iSubjSess).path, ModalityIs, ModalityFields(TypeRunIndex).filename);
                                [~, ~, Fext] = xASL_fileparts(ModalityFields(TypeRunIndex).filename);
                                pathDest{1} = fullfile(pathLegacy_SubjectVisit, FolderIs, [FileIs Fext]);
                                
                                %% -----------------------------------------------
                                %% 7. Manage sidecars to copy
								% Sidecars definitions are loaded by xASL_bids_Config at function start

                                % Assuming that each .nii has a .json sidecar, do the same for .json (and for
                                % other sidecars only if they exist per bidsPar.sidecarRequired)
                                
                                iCount = 1;
                                for iCar=1:length(bidsPar.sidecarName)
                                    [Fpath, Ffile] = xASL_fileparts(pathOrig{1});
                                    
                                    if ~bidsPar.sidecarSuffixType(iCar)
                                        Ffile = Ffile(1:end-length(TypeIs)-1);
                                    end
                                    TempSidecar = fullfile(Fpath, [Ffile bidsPar.sidecarName{iCar}]);
                                    
                                    if ~strcmp(bidsPar.sidecarTypeSpecific{iCar}, 'no') && ~strcmp(bidsPar.sidecarTypeSpecific{iCar}, TypeIs)
                                        % skip this sidecar (e.g. some asl-specific sidecars for non-asl NIfTIs)
                                    elseif ~exist(TempSidecar, 'file') && bidsPar.sidecarRequired(iCar)
                                            warning([TempSidecar ' missing']);
                                    elseif exist(TempSidecar, 'file')
                                        pathOrig{iCount+1} = TempSidecar;

                                        [Fpath, Ffile] = xASL_fileparts(pathDest{1});
                                        pathDest{iCount+1} = fullfile(Fpath, [Ffile bidsPar.sidecarName{iCar}]);
                                        iCount = iCount+1;
                                    end
                                end
                                    
                                %% -----------------------------------------------
                                %% 8. Copy files
                                for iFile=1:length(pathOrig)
                                    xASL_bids_BIDS2xASL_CopyFile(pathOrig{iFile}, pathDest{iFile}, bOverwrite);
                                end

                            end
                        end % iterate runs in this Subject/Visit/Modality
                    end
                end
            end
        end
    % end
end

fprintf('\n');

%% 9. Parse M0
ListASL4D = xASL_adm_GetFileList(pathLegacy, '^ASL4D\.nii$', 'FPListRec');
if ~isempty(ListASL4D)
    for iList=1:numel(ListASL4D)
        xASL_bids_parseM0(ListASL4D{iList});
    end
    fprintf('%s\n', ['M0 parsed for ' ListASL4D{iList}]);
else
    warning(['No ASL4D file found in ' pathLegacy]);
end

%% 10. Create DataPar.json

% Write DataParFile if it does not exist already
if ~(xASL_exist(fullfile(pathLegacy, 'dataPar.json'),'file')==2)
    fprintf('Creating default dataPar.json since file was not found...\n');
    spm_jsonwrite(fullfile(pathLegacy, 'dataPar.json'), dataPar);
end

% Overwrite dataPar.json in x structure
fprintf('Overwriting x.dir.dataPar...\n');

% Add the path to the dataPar.x struct that we return to the Master script
dataPar.x.dir.dataPar = fullfile(pathLegacy, 'dataPar.json');

%% 11. Copy participants.tsv

% Check if participants.tsv exists
if xASL_exist(fullfile(pathStudy, 'participants.tsv'),'file')==2
    xASL_Copy(fullfile(pathStudy, 'participants.tsv'),fullfile(pathLegacy, 'Participants.tsv'));
elseif xASL_exist(fullfile(pathStudy, 'Participants.tsv'),'file')==2
    xASL_Copy(fullfile(pathStudy, 'Participants.tsv'),fullfile(pathLegacy, 'Participants.tsv'));
end

%% 12. Add "GeneratedBy" fields
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

%% 13. Clean up
try
    filesCleanUp = xASL_adm_GetFileList(pathStudy,'^import_.+$');
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


%% ===========================================================================
function xASL_bids_BIDS2xASL_CopyFile(pathOrig, pathDest, bOverwrite)
%xASL_bids_BIDS2xASL_CopyFile

    % Create folder(s) if didnt exist
    xASL_adm_CreateDir(fileparts(pathDest));

    % check for existance & overwriting
    if ~xASL_exist(pathOrig)
        warning(['Couldnt find ' pathOrig, ' skipping']);
    elseif xASL_exist(pathDest) && ~bOverwrite
        warning([pathDest ' already existed, skipping']);
    else % if didnt already exist, or bOverwrite is true
        xASL_Copy(pathOrig, pathDest, 1);
    end
    
end


%% ===========================================================================
function xASL_bids_AddGeneratedByField(x, pathJSONin, pathJSONout)

    % If only one path is provided, then overwrite original file
    if nargin<3
        pathJSONout = pathJSONin;
    end

    % Check if input JSON exists
    if xASL_exist(pathJSONin,'file')==2
        % Load JSON file
        thisJSON = spm_jsonread(pathJSONin);
        % Add "GeneratedBy" field
        thisJSON.GeneratedBy = struct;
        thisJSON.GeneratedBy.Name = 'ExploreASL';
        thisJSON.GeneratedBy.Version = x.Version;
        spm_jsonwrite(pathJSONout, thisJSON);
    else
        warning('Adding the "GeneratedBy" JSON field failed...');
    end

end


