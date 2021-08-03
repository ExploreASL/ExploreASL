function xASL_imp_NII2BIDS_Session(imPar, bidsPar, studyPar, listSessions, nameSubjectSession, bidsLabel, iSession)
%xASL_imp_NII2BIDS_Session NII2BIDS conversion for a single session.
%
% FORMAT: xASL_imp_NII2BIDS_Session(imPar, bidsPar, studyPar, listSessions, nameSubjectSession, subjectLabel, iSession)
% 
% INPUT:
%   imPar                 - JSON file with structure with import parameter (STRUCT, REQUIRED)
%   bidsPar               - Output of xASL_imp_Config (STRUCT, REQUIRED)
%   studyPar              - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
%   listSessions          - list of sessions (CELL ARRAY, REQUIRED, e.g.: {'ASL_1'})
%   nameSubjectSession    - name of the subject (CELL ARRAY, REQUIRED)
%   subjectLabel          - subject label (CHAR ARRAY, REQUIRED)
%   iSession              - Session number (INTEGER, REQUIRED)
%
% OUTPUT:
%  n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: NII2BIDS conversion for a single sessions.
% 
% 1. Make a subject directory
% 2. Iterate over runs
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Unpack bidsLabel struct
    subjectLabel = bidsLabel.subject;
    visitLabel = bidsLabel.visit;

    %% 1. Make a subject directory
    if length(listSessions)>1 || ~isempty(visitLabel)
        if ~isempty(visitLabel)
            sessionLabel = ['ses-' visitLabel];
        else
            sessionLabel = ['ses-' listSessions{iSession}(5:end)];
        end
        
        xASL_adm_CreateDir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel));
        xASL_adm_CreateDir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel,'perf'));
        
        inSessionPath = fullfile(imPar.TempRoot,nameSubjectSession,listSessions{iSession});
        outSessionPath = fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel);

        % Need to add the underscore so that it doesn't need to be added automatically and can be skipped for empty session
        sessionLabel = ['_' sessionLabel];
    else
        % Session label is skipped
        sessionLabel = '';

        % Only one session - no session labeling
        if ~exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel]),'dir')
            xASL_adm_CreateDir(fullfile(imPar.BidsRoot,['sub-' subjectLabel]));
            xASL_adm_CreateDir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],bidsPar.strPerfusion));
        end
        inSessionPath = fullfile(imPar.TempRoot,nameSubjectSession,listSessions{iSession});
        outSessionPath = fullfile(imPar.BidsRoot,['sub-' subjectLabel]);
    end

    % Check if there are multiple runs per session
    listRuns = xASL_adm_GetFileList(inSessionPath,'^ASL4D_\d.nii+$',false,[],false);

    %% 2. Iterate over runs
    for iRun = 1:(max(length(listRuns),1))
        xASL_imp_NII2BIDS_Run(imPar, bidsPar, studyPar, [subjectLabel sessionLabel], inSessionPath, outSessionPath, listRuns, iRun, nameSubjectSession);
    end

end
