function [result, x] = xASL_module_BIDS2Legacy(x, bOverwrite, bVerbose)
%xASL_module_BIDS2Legacy Convert BIDS rawdata to ExploreASL legacy format
%
% FORMAT: [x] = xASL_module_BIDS2Legacy(x [,bOverwrite, bVerbose]);
%
% INPUT:
%   x             - ExploreASL x structure, containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   bOverwrite    - Boolean indicating whether to overwrite existing files (OPTIONAL, DEFAULT = true)
%   bVerbose   - boolean, true for verbose output (OPTIONAL, DEFAULT = false)
%
% OUTPUT:
%   x             - ExploreASL x structure, containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%  result         - Boolean indicating whether module succeeded
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
% 1. Get subjectID & sessionID from x.modules.bids2legacy.BIDS
% 2. Parse modality
% 3. Parse M0s
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        x = xASL_module_BIDS2Legacy(x);
%
% __________________________________
% Copyright (c) 2015-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    %% 0. Initialization
    if nargin<3 || isempty(bVerbose)
        bVerbose = 0; % default to not have verbose output for each subject (the mutex set up will already do this)
    end
    if nargin<2 || isempty(bOverwrite)
        bOverwrite = 1; % Default is to overwrite existing files
    end
    % Make sure that logging is still active
    if isfield(x.dir,'diaryFile')
        diary(x.dir.diaryFile); % Start diary
    end
    
    [x] = xASL_init_InitializeMutex(x, 'BIDS2Legacy'); % Start Mutex
    result = true; % Default for result
    if x.mutex.bAnyModuleLocked
        % If any module is locked for this subject, we skip this module for
        % this subject
        return;
    end
    if x.mutex.HasState('999_ready')
        bO = false; % no Output, as everything has been done already
    else
        bO = true; % yes Output, not completely done, so we want to know what has or has not been done
    end    
    
    if ~x.mutex.HasState('010_BIDS2LEGACY')
        %% 1. Get subjectID & sessionID from x.modules.bids2legacy.BIDS
        %   !!!!! NOTE THAT SESSION (BIDS) IS A VISIT (LEGACY) HERE, NOT A RUN (BIDS) !!!!!
        iSubjSess = find(strcmp(x.SUBJECTS, x.SUBJECT));
        
        % each index in x.modules.bids2legacy.BIDS.subjects has its unique subject-session combination
        % which is identical to the order in ExploreASL legacy x.SUBJECTS
        % which is defined in xASL_init_BIDS2Legacy, e.g.,
        % x.modules.bids2legacy.BIDS.subjects(1) -> name = sub-10015124 session = ses-1, legacy sub-10015124_1
        % x.modules.bids2legacy.BIDS.subjects(2) -> name = sub-10015124 session = ses-2, legacy sub-10015124_2
        % x.modules.bids2legacy.BIDS.subjects(3) -> name = sub-10015125 session = ses-1, legacy sub-10015125_1
        %
        % Note that a subject-session combination can have multiple runs, e.g.,
        % x.modules.bids2legacy.BIDS.subjects(1).perf(1) -> filename: 'sub-10015124_ses-1_run-1_asl.nii.gz'
        % x.modules.bids2legacy.BIDS.subjects(1).perf(2) -> filename: 'sub-10015124_ses-1_run-2_asl.nii.gz'
        % x.modules.bids2legacy.BIDS.subjects(1).perf(3) -> filename: 'sub-10015124_ses-1_run-1_m0scan.nii.gz'
        % x.modules.bids2legacy.BIDS.subjects(1).perf(4) -> filename: 'sub-10015124_ses-1_run-2_m0scan.nii.gz'  
		%
		% Also note that sessions can now be any name, number, data, or a string and are always treated as a string
        % Subject ID
        SubjectID = x.modules.bids2legacy.BIDS.subjects(iSubjSess).name;
        SessionID = x.modules.bids2legacy.BIDS.subjects(iSubjSess).session;
        if isempty(SessionID)
            SessionID = '1';
        else
            SessionID = SessionID(5:end);
        end
        SubjectSession = [SubjectID '_' SessionID]; % this is the legacy name
        if ~strcmp(SubjectSession, x.SUBJECT)
            error(['xASL_init_BIDS2Legacy ' x.SUBJECT ' should match with xASL_module_BIDS2Legacy ' SubjectSession]);
        else
            %% 2. Parse modality
            
            pathLegacy_SubjectSession = fullfile(x.dir.xASLDerivatives, SubjectSession);
            bFolderExisted = exist(pathLegacy_SubjectSession, 'dir')==7;
            % Backwards compatibility: rename sub-*** to sub-***_1
            if strcmp(SessionID, '1')
                pathLegacy_SubjectSession_OLD = fullfile(x.dir.xASLDerivatives, SubjectID); % it should be [SubjectID '_' SessionID] now
                if ~bFolderExisted && exist(pathLegacy_SubjectSession_OLD, 'dir')==7
                    warning(['Renaming ' pathLegacy_SubjectSession_OLD '->' SubjectSession]);
                    
                    % Rename all files
                    fList = xASL_adm_GetFileList(x.dir.xASLDerivatives, SubjectID, 'FPListRec');
                    for iList=1:length(fList)
                        newPath = replace(fList{iList}, SubjectID, SubjectSession);
                        if isempty(regexpi(newPath, [SubjectSession '_1'])) && ~xASL_exist(newPath, 'file') % avoid pre-existing files (e.g., current BIDS2Legacy lock & log files)
                            xASL_Move(fList{iList}, newPath);
                        end
                    end
                    
                    % Rename all folders
                    fList = xASL_adm_GetFileList(x.dir.xASLDerivatives, [SubjectID '$'], 'FPListRec', [], 1);
                    for iList=1:length(fList)
                        newPath = replace(fList{iList}, SubjectID, SubjectSession);
                        if isempty(regexpi(newPath, [SubjectSession '_1'])) % avoid pre-existing files (e.g., current BIDS2Legacy lock & log files)
                            if xASL_exist(newPath, 'dir')
                                addpath(fullfile(x.opts.MyPath, 'WorkInProgress', 'General'));
                                xASL_MoveDirMerge(fList{iList}, newPath);
                                rmdir(fList{iList});
                            else
                                xASL_Move(fList{iList}, newPath);
                            end
                        end
                    end
                    bFolderExisted = true;
                end
            end
            % Provide warning that we are going to overwrite the folder
            if bFolderExisted
                warning(['Any rawdata copies will be overwritten in the /derivatives folder: ' pathLegacy_SubjectSession]);
                fprintf('%s\n', ['Note that /rawdata itself remains untouched. And we only overwrite the copies of /rawdata in /derivatives, ' ...
                    'any processing results in /derivatives will not be overwritten,' ...
                    'but we merge previously processed results with fresh copies of original files from /rawdata']);
                fprintf('\n%s\n\n', 'Nevertheless, it is always cleanest to remove all derivatives before rerunning ExploreASL');
            end
            % Create subject/session directory (to enable reruns for pre-imported or crashed datasets, we need a subject level here/above!)
            xASL_adm_CreateDir(pathLegacy_SubjectSession);
            % Modalities - the BIDS scantypes
            ModalitiesUnique = unique(x.modules.bids2legacy.bidsPar.BIDS2LegacyFolderConfiguration(2:end, 2));
            nModalities = length(ModalitiesUnique);
            xASL_bids_BIDS2Legacy_ParseModality(x.modules.bids2legacy.BIDS, x.modules.bids2legacy.bidsPar, SubjectSession, iSubjSess, ModalitiesUnique, nModalities, bOverwrite, pathLegacy_SubjectSession);
        end
        %% 3. Parse M0s
        ListASL4D = xASL_adm_GetFileList(pathLegacy_SubjectSession, '^ASL4D\.nii$', 'FPListRec');
        if bVerbose && isempty(ListASL4D)
            warning(['When parsing M0: no ASL4D.nii found for ' SubjectID 'in ' x.dir.xASLDerivatives '...']);
        end
        
        % Parse M0s
        for iList=1:numel(ListASL4D)
            xASL_bids_parseM0(ListASL4D{iList});
            if bVerbose
                fprintf('%s\n', ['M0 parsed for subject ' SubjectID ' session ' SessionID]);
            end
        end
        x.mutex.AddState('010_BIDS2LEGACY'); % Add mutex state
    elseif x.mutex.HasState('010_BIDS2LEGACY')
        if bO; fprintf('%s\n', ['BIDS2Legacy already done, skipping ' x.SUBJECT '.    ']); end
    end
    
    %% Finalize and unlock mutex for this module
    x.mutex.AddState('999_ready'); % Add ready state
    x.mutex.Unlock(); % Unlock mutex
        
end