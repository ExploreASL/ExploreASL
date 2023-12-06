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
% 0. Initialization
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
% 6. Finalize and unlock mutex for this module
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        x = xASL_module_BIDS2Legacy(x);
%
% __________________________________
% Copyright (c) 2015-2023 ExploreASL


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

    
    xASL_adm_BreakString('BIDS to ExploreASL LEGACY CONVERSION'); % Print feedback
	

    if ~x.mutex.HasState('010_BIDS2LEGACY')





        %% 1. Parse a folder using the output of bids-matlab (was run before this point)


        %% 2. Define Subject
        iSubjSess = find(strcmp(x.SUBJECTS, x.SUBJECT));
        
        % each index in x.modules.bids2legacy.BIDS.subjects has its unique subject-visit combination
        % which is identical to the order in ExploreASL legacy x.SUBJECTS
        % which is defined in xASL_init_BIDS2Legacy, e.g.,
        % x.modules.bids2legacy.BIDS.subjects(1) -> name = sub-10015124 session = ses-1, legacy sub-10015124_1
        % x.modules.bids2legacy.BIDS.subjects(2) -> name = sub-10015124 session = ses-2, legacy sub-10015124_2
        % x.modules.bids2legacy.BIDS.subjects(3) -> name = sub-10015125 session = ses-1, legacy sub-10015125_1
        %
        % Note that a subject-visit combination can have multiple runs, e.g.,
        % x.modules.bids2legacy.BIDS.subjects(1).perf(1) -> filename: 'sub-10015124_ses-1_run-1_asl.nii.gz'
        % x.modules.bids2legacy.BIDS.subjects(1).perf(2) -> filename: 'sub-10015124_ses-1_run-2_asl.nii.gz'
        % x.modules.bids2legacy.BIDS.subjects(1).perf(3) -> filename: 'sub-10015124_ses-1_run-1_m0scan.nii.gz'
        % x.modules.bids2legacy.BIDS.subjects(1).perf(4) -> filename: 'sub-10015124_ses-1_run-2_m0scan.nii.gz'            

        % Subject ID
        SubjectID = x.modules.bids2legacy.BIDS.subjects(iSubjSess).name;
        SessionID = x.modules.bids2legacy.BIDS.subjects(iSubjSess).session;
        if isempty(SessionID)
            iVisit = 1; % we default to a single visit
        else
            iVisit = str2num(SessionID(5:end));
        end

        SubjectVisit = [SubjectID '_' xASL_num2str(iVisit)]; % this is the legacy name

        if ~strcmp(SubjectVisit, x.SUBJECT)
            error(['xASL_init_BIDS2Legacy ' x.SUBJECT ' should match with xASL_module_BIDS2Legacy ' SubjectVisit]);
        else

            %% 3. Define Session
            
            pathLegacy_SubjectVisit = fullfile(x.dir.xASLDerivatives, SubjectVisit);
            
            % Create subject/session directory (to enable reruns for pre-imported or crashed datasets, we need a subject level here/above!)
            xASL_adm_CreateDir(pathLegacy_SubjectVisit);

            %% 4. Parse modality
            % Modalities - the BIDS scantypes
            ModalitiesUnique = unique(x.modules.bids2legacy.bidsPar.BIDS2LegacyFolderConfiguration(2:end, 2));
            nModalities = length(ModalitiesUnique);
            xASL_bids_BIDS2Legacy_ParseModality(x.modules.bids2legacy.BIDS, x.modules.bids2legacy.bidsPar, SubjectVisit, iSubjSess, ModalitiesUnique, nModalities, bOverwrite, pathLegacy_SubjectVisit);

        end



        %% 5. Parse M0s
        ListASL4D = xASL_adm_GetFileList(pathLegacy_SubjectVisit, '^ASL4D\.nii$', 'FPListRec');
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
        fprintf('%s\n', ['BIDS2Legacy already done, skipping ' x.SUBJECT '.    ']);
    end

    
    %% 6. Finalize and unlock mutex for this module
    x.mutex.AddState('999_ready'); % Add ready state
    x.mutex.Unlock(); % Unlock mutex
        
end