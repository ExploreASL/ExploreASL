function x = xASL_imp_NII2BIDS_Session(x, imPar, bidsPar, studyPar, listSessions, nameSubjectSession, bidsLabel, iSession)
%xASL_imp_NII2BIDS_Session NII2BIDS conversion for a single session.
%
% FORMAT: x = xASL_imp_NII2BIDS_Session(x, imPar, bidsPar, studyPar, listSessions, nameSubjectSession, bidsLabel, iSession)
% 
% INPUT:
%   x                     - ExploreASL x structure (REQUIRED, STRUCT)
%   imPar                 - JSON file with structure with import parameter (STRUCT, REQUIRED)
%   bidsPar               - Output of xASL_imp_Config (STRUCT, REQUIRED)
%   studyPar              - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
%   listSessions          - list of sessions (CELL ARRAY, REQUIRED, e.g.: {'ASL_1'})
%   nameSubjectSession    - name of the subject (CELL ARRAY, REQUIRED)
%   bidsLabel             - BIDS label (CHAR ARRAY, REQUIRED)
%   iSession              - Session number (INTEGER, REQUIRED)
%
% OUTPUT:
%   x                     - ExploreASL x structure (STRUCT)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: NII2BIDS conversion for a single sessions.
% 
% 1. Make a subject directory
% 2. Iterate over runs
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     x = xASL_imp_NII2BIDS_Session(x, imPar, bidsPar, studyPar, listSessions, nameSubjectSession, bidsLabel, iSession);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Print the session that is being converted
    fprintf('\n====================================== CONVERT SESSION =======================================\n');
    fprintf('Converting subject %s, session %s, ', bidsLabel.subject, bidsLabel.visit);

    %% 1. Make a subject directory
    if length(listSessions)>1 || ~isempty(bidsLabel.visit)
        if ~isempty(bidsLabel.visit)
            sessionLabel = ['ses-' bidsLabel.visit];
        else
            sessionLabel = ['ses-' listSessions{iSession}(5:end)];
        end
        
        xASL_adm_CreateDir(fullfile(imPar.BidsRoot,['sub-' bidsLabel.subject], sessionLabel,'perf'));
        
        inSessionPath = fullfile(imPar.TempRoot, nameSubjectSession, listSessions{iSession});
        outSessionPath = fullfile(imPar.BidsRoot, ['sub-' bidsLabel.subject], sessionLabel);

        % Need to add the underscore so that it doesn't need to be added automatically and can be skipped for empty session
        sessionLabel = ['_' sessionLabel];
    else
        % Session label is skipped
        sessionLabel = '';

        % Only one session - no session labeling
        if ~exist(fullfile(imPar.BidsRoot,['sub-' bidsLabel.subject]),'dir')
            xASL_adm_CreateDir(fullfile(imPar.BidsRoot, ['sub-' bidsLabel.subject], bidsPar.strPerfusion));
        end
        inSessionPath = fullfile(imPar.TempRoot, nameSubjectSession, listSessions{iSession});
        outSessionPath = fullfile(imPar.BidsRoot, ['sub-' bidsLabel.subject]);
    end

    % Check if there are multiple runs per session
    listRuns = xASL_adm_GetFileList(inSessionPath, '^ASL4D_\d.nii+$', false, [], false);

    %% 2. Iterate over runs
    for iRun = 1:(max(length(listRuns),1))
        x = xASL_imp_NII2BIDS_Run(x, imPar, bidsPar, studyPar, [bidsLabel.subject sessionLabel], inSessionPath, outSessionPath, listRuns, iRun, nameSubjectSession);
    end

end

