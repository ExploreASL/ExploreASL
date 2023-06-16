function [x] = xASL_wrp_BIDS2Legacy(x, bOverwrite)
%xASL_wrp_BIDS2Legacy Convert BIDS rawdata to ExploreASL legacy format
%
% FORMAT: [x] = xASL_wrp_BIDS2Legacy(x[, bOverwrite])
% 
% INPUT:
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
% EXAMPLE: [x] = xASL_wrp_BIDS2Legacy(x, bOverwrite);
% __________________________________
% Copyright (c) 2015-2023 ExploreASL


%% 0. Admin

% Verify the input parameters

if nargin<2 || isempty(bOverwrite)
    bOverwrite = 1;
end

% Creates the derivatives directory
if exist(x.dir.xASLDerivatives, 'dir') && bOverwrite
    fprintf('%s\n', [x.dir.xASLDerivatives ' already exists, overwriting any pre-existing derivatives...']);
elseif exist(x.dir.xASLDerivatives, 'dir')
    fprintf('%s\n', [x.dir.xASLDerivatives ' exists, merging with new derivatives']);
else
    xASL_adm_CreateDir(x.dir.xASLDerivatives);
end

% Loads the configuration for file renaming from the BIDS configuration file
bidsPar = xASL_bids_Config();

%% 1. Parse a folder using the output of bids-matlab (was run before this point)

nSubjects = numel(x.modules.bids2legacy.BIDS.subjectName);
nVisits = numel(x.modules.bids2legacy.BIDS.sessionName); % this is called sessions in BIDS
% we use this below to see if the legacy subjectname gets _1 as visit suffix or not
[~, studyName] = fileparts(x.dir.DatasetRoot);
fprintf('Converting from BIDS to Legacy: %s   \n', studyName);

%% 2. Define SubjectSession
for iSubjSess=1:numel(x.modules.bids2legacy.BIDS.subjects) 
    % Iterate over x.modules.bids2legacy.BIDS.subjects (indices that include both subjects & sessions)
    % so  1 subject  6 session/visits, will give numel(x.modules.bids2legacy.BIDS.subjects)=6
    % and 6 subjects 1 session/visits, will give numel(x.modules.bids2legacy.BIDS.subjects)=6
    xASL_TrackProgress(iSubjSess, nSubjects);
    
    % Subject ID
    SubjectID = x.modules.bids2legacy.BIDS.subjects(iSubjSess).name;
    
    % Subject-wise processing (because xASL_Iteration runs over subjects)
    if ~isempty(regexpi(SubjectID, x.SUBJECT, 'once'))
    
        % Session ID (Currently, ExploreASL concatenates subject_visit/timepoint in the same folder layer, so we only use SubjectSession)
        SessionID = x.modules.bids2legacy.BIDS.subjects(iSubjSess).session;

        %% 3. Define Session
        iVisit = find(strcmp(x.modules.bids2legacy.BIDS.sessionName, SessionID));
        
        if isempty(SessionID)
            SessionID = 'unknown'; % default to have at least some output
        end
        
        % Remove iteration for iVisit = 1 -> iterate visit/session in this "x.modules.bids2legacy.BIDS.subjects" (always 1 session per x.modules.bids2legacy.BIDS.subjects)
        % ExploreASL uses visit as a number (e.g. _1 _2 _3 etc)
        if nVisits==1
            pathLegacy_SubjectVisit = fullfile(x.dir.xASLDerivatives, SubjectID);
            VisitString = '';
        else
			if isempty(iVisit)
                fprintf('\nEmpty session number, setting session number to 1...\n');
                pathLegacy_SubjectVisit = fullfile(x.dir.xASLDerivatives, [SubjectID '_' xASL_num2str(iVisit)]);
			else
				% If Visit name is of a form ses-number then use this number otherwise the ID
				if regexpi(SessionID, 'ses-\d+')
					pathLegacy_SubjectVisit = fullfile(x.dir.xASLDerivatives, [SubjectID '_' SessionID(5:end)]);
				else
					error(['Not yet supporting session names as strings, but only as numbers for session ', SessionID]);
				end
			end
            VisitString = [' visit ' SessionID];
        end
        
        SubjectVisit = [SubjectID VisitString];
        
        % Create subject/session directory (to enable reruns for pre-imported or crashed datasets, we need a subject level here/above!)
        xASL_adm_CreateDir(pathLegacy_SubjectVisit);

        %% 4. Parse modality
        % Modalities - the BIDS scantypes
        ModalitiesUnique = unique(bidsPar.BIDS2LegacyFolderConfiguration(2:end, 2));
        nModalities = length(ModalitiesUnique);
        xASL_bids_BIDS2Legacy_ParseModality(x.modules.bids2legacy.BIDS, bidsPar, SubjectVisit, iSubjSess, ModalitiesUnique, nModalities, bOverwrite, pathLegacy_SubjectVisit);

    end

end
% Print new line after track progress bar
fprintf('   \n');


% Get directories of current subject. The BIDS subject name can be prefixed
% with sub- and suffixed with _1 _2 _3 etc for visits, in the legacy format
SubjectDirs = xASL_adm_GetFileList(fullfile(x.dir.xASLDerivatives), ['^(|sub-)' x.SUBJECT '(|_\d*)$'], [],[],true);


%% 5. Parse M0 of current subject
ListASL4D = [];
if ~isempty(SubjectDirs)
    for iDir=1:numel(SubjectDirs)
        currentASL = xASL_adm_GetFileList(SubjectDirs{iDir}, '^ASL4D\.nii$', 'FPListRec');
        ListASL4D = vertcat(ListASL4D, currentASL);
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
    warning(['When parsing M0: no ASL4D.nii found for ' SubjectID 'in ' x.dir.xASLDerivatives '...']);
end


%% 6. Clean up
xASL_bids_BIDS2Legacy_CleanUp(x.dir.DatasetRoot);


end



%% Clean-Up subfunction
function xASL_bids_BIDS2Legacy_CleanUp(pathStudy)

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





