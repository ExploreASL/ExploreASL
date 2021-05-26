function xASL_imp_NII2BIDS_SubjectSession(imPar, bidsPar, studyPar, iSubject, fSes, listSubjects, subjectLabel, kk)
%xASL_imp_NII2BIDS_SubjectSession NII2BIDS conversion for a single sessions.
%
% FORMAT: xASL_imp_NII2BIDS_SubjectSession(imPar, bidsPar, studyPar, iSubject, fSes, listSubjects, subjectLabel, kk)
% 
% INPUT:
%   imPar          - JSON file with structure with import parameter (STRUCT, REQUIRED)
%   bidsPar        - Output of xASL_imp_Config (STRUCT, REQUIRED)
%   studyPar       - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
%   iSubject       - Subject ID (INTEGER, REQUIRED)
%   fSes           - f sessions (CELL ARRAY, REQUIRED, e.g.: {'ASL_1'})
%   listSubjects   - list of subjects (CELL ARRAY, REQUIRED)
%   subjectLabel   - subject label (CHAR ARRAY, REQUIRED)
%   kk             - Session number (INTEGER, REQUIRED)
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
    if length(fSes)>1
        sessionLabel = ['ses-' fSes{kk}(5:end)];

        if ~exist(fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel),'dir')
            mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel));
            mkdir(fullfile(imPar.BidsRoot,['sub-' subjectLabel],sessionLabel,'asl'));
        end
        inSessionPath = fullfile(imPar.TempRoot,listSubjects{iSubject},fSes{kk});
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
        inSessionPath = fullfile(imPar.TempRoot,listSubjects{iSubject},fSes{kk});
        outSessionPath = fullfile(imPar.BidsRoot,['sub-' subjectLabel]);
    end

    % Check if there are multiple runs per session
    fRuns = xASL_adm_GetFileList(inSessionPath,'^ASL4D_\d.nii+$',false,[],false);
    nSes = length(fRuns);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3.2 Iterate over runs
    for mm = 1:(max(nSes,1))
        xASL_imp_NII2BIDS_SubjectSessionRun(imPar, bidsPar, studyPar, subjectLabel, sessionLabel, listSubjects, fSes, inSessionPath, outSessionPath, nSes, iSubject, kk, mm);
    end

end


