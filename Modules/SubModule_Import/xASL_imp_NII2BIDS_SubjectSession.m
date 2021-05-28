function xASL_imp_NII2BIDS_SubjectSession(imPar, bidsPar, studyPar, listSessions, nameSubject, subjectLabel, iSession)
%xASL_imp_NII2BIDS_SubjectSession NII2BIDS conversion for a single sessions.
%
% FORMAT: xASL_imp_NII2BIDS_SubjectSession(imPar, bidsPar, studyPar, listSessions, nameSubject, subjectLabel, iSession)
% 
% INPUT:
%   imPar          - JSON file with structure with import parameter (STRUCT, REQUIRED)
%   bidsPar        - Output of xASL_imp_Config (STRUCT, REQUIRED)
%   studyPar       - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
%   listSessions   - list of sessions (CELL ARRAY, REQUIRED, e.g.: {'ASL_1'})
%   nameSubject    - name of the subject (CELL ARRAY, REQUIRED)
%   subjectLabel   - subject label (CHAR ARRAY, REQUIRED)
%   iSession       - Session number (INTEGER, REQUIRED)
%
% OUTPUT:
%  n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: NII2BIDS conversion for a single sessions.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3.1 Make a subject directory
    if length(listSessions)>1
        sessionLabel = ['ses-' listSessions{iSession}(5:end)];

        if ~exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel),'dir')
            mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel));
            mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel,'asl'));
        end
        inSessionPath = fullfile(imPar.TempRoot,nameSubject,listSessions{iSession});
        outSessionPath = fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel);

        % Need to add the underscore so that it doesn't need to be added automatically and can be skipped for empty session
        sessionLabel = ['_' sessionLabel];
    else
        % Session label is skipped
        sessionLabel = '';

        % Only one session - no session labeling
        if ~exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel]),'dir')
            mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel]));
            mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],bidsPar.strPerfusion));
        end
        inSessionPath = fullfile(imPar.TempRoot,nameSubject,listSessions{iSession});
        outSessionPath = fullfile(imPar.BidsRoot,['sub-' subjectLabel]);
    end

    % Check if there are multiple runs per session
    listRuns = xASL_adm_GetFileList(inSessionPath,'^ASL4D_\d.nii+$',false,[],false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3.2 Iterate over runs
    for iRun = 1:(max(length(listRuns),1))
        xASL_imp_NII2BIDS_SubjectSessionRun(imPar, bidsPar, studyPar, [subjectLabel sessionLabel], nameSubject, listSessions{iSession}, inSessionPath, outSessionPath, listRuns, iRun);
    end

end


