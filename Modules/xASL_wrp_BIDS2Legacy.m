function [x] = xASL_wrp_BIDS2Legacy(x, bOverwrite, bVerbose)
%xASL_wrp_BIDS2Legacy Convert BIDS rawdata to ExploreASL legacy format
%
% FORMAT: [x] = xASL_wrp_BIDS2Legacy(x[, bOverwrite, bVerbose])
% 
% INPUT:
%   x          - ExploreASL x structure (REQUIRED, STRUCT)
%   bOverwrite - boolean, true for overwriting files (OPTIONAL, DEFAULT = true)
%   bVerbose   - boolean, true for verbose output (OPTIONAL, DEFAULT = true)
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
% 1. Parse a folder using bids-matlab output
% 2. Define SubjectSession
% 3. Define Session
% 4. Parse modality
%    - Parse scantype
%    - Compile paths for copying
%    - Manage sidecars to copy
%    - Copy files
% 5. Parse M0
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [x] = xASL_wrp_BIDS2Legacy(x);
% __________________________________
% Copyright (c) 2015-2023 ExploreASL


%% 0. Admin

% Verify the input parameters

if nargin<3 || isempty(bVerbose)
    bVerbose = 1; % default to verbosity
end
if nargin<2 || isempty(bOverwrite)
    bOverwrite = 1; % default to overwrite
end


%% 1. Parse a folder using the output of bids-matlab (was run before this point)
nSubjects = numel(x.modules.bids2legacy.BIDS.subjectName);
nVisits = numel(x.modules.bids2legacy.BIDS.sessionName); % this is called sessions in BIDS
% we use this below to see if the legacy subjectname gets _1 as visit suffix or not
[~, studyName] = fileparts(x.dir.DatasetRoot);

if bVerbose
    fprintf('Converting from BIDS to Legacy: %s   \n', studyName);
end

%% 2. Define Subject
subjectBIDS = ['sub-' x.SUBJECT];

indicesCurrentSubject = find(cellfun( @(y) strcmp(y, subjectBIDS), {x.modules.bids2legacy.BIDS.subjects.name}));


for iSubjSess=indicesCurrentSubject
    % Iterate over x.modules.bids2legacy.BIDS.subjects for the current subject
    % (indices that include both subjects & sessions)
    % so  1 subject 6 session/visits, will give numel(x.modules.bids2legacy.BIDS.subjects)=6
    
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
                if bVerbose
                    fprintf('\nEmpty session number, setting session number to 1...\n');
                    pathLegacy_SubjectVisit = fullfile(x.dir.xASLDerivatives, [SubjectID '_1']);
                end
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
        ModalitiesUnique = unique(x.modules.bids2legacy.bidsPar.BIDS2LegacyFolderConfiguration(2:end, 2));
        nModalities = length(ModalitiesUnique);
        xASL_bids_BIDS2Legacy_ParseModality(x.modules.bids2legacy.BIDS, x.modules.bids2legacy.bidsPar, SubjectVisit, iSubjSess, ModalitiesUnique, nModalities, bOverwrite, pathLegacy_SubjectVisit);

    end
end


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
        [currentFolder, currentNifti] = xASL_fileparts(ListASL4D{iList});
        [~, currentFolder] = fileparts(fileparts(currentFolder));

        % Get session
        sessionIndex = regexp(currentFolder, '_\d+');
        if isempty(sessionIndex)
            SessionID = '1';
        else
            SessionID = num2str(currentFolder(sessionIndex+1:end));
        end

        if bVerbose
            fprintf('%s\n', ['M0 parsed for subject ' SubjectID ' session ' SessionID ': image ' currentNifti]);
        end
    end
else
    if bVerbose
        warning(['When parsing M0: no ASL4D.nii found for ' SubjectID 'in ' x.dir.xASLDerivatives '...']);
    end
end



end