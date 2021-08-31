function xASL_imp_NII2BIDS_Run(imPar, bidsPar, studyPar, subjectSessionLabel, inSessionPath, outSessionPath, listRuns, iRun, nameSubjectSession)
%xASL_imp_NII2BIDS_Run NII2BIDS conversion for a single sessions, single run.
%
% FORMAT: xASL_imp_NII2BIDS_Run(bidsPar, studyPar, subjectSessionLabel, inSessionPath, outSessionPath, listRuns, iRun)
% 
% INPUT:
% imPar               - JSON file with structure with import parameter (STRUCT, REQUIRED)
% bidsPar             - Output of xASL_imp_Config (STRUCT, REQUIRED)
% studyPar            - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
% subjectSessionLabel - subject-session label (CHAR ARRAY, REQUIRED)
% inSessionPath       - input session path (CHAR ARRAY, PATH, REQUIRED)
% outSessionPath      - output session path (CHAR ARRAY, PATH, REQUIRED)
% listRuns            - list of runs (INTEGER, REQUIRED)
% iRun                - run number (INTEGER, REQUIRED)
% nameSubjectSession  - Name of subject & session (REQUIRED)
%
% OUTPUT:
% n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: NII2BIDS conversion for a single sessions, single runs.
%
% 1. Convert anat files
% 2. Convert perf files
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Anatomical files
    try
        xASL_imp_NII2BIDS_RunAnat(imPar, bidsPar, studyPar, subjectSessionLabel, outSessionPath, listRuns, iRun, nameSubjectSession);
    catch ME
        xASL_imp_NII2BIDS_Run_ISSUE_WARNING(ME, 'anatomical', subjectSessionLabel, iRun);
    end
    
    %% Perfusion files
    try
        xASL_imp_NII2BIDS_RunPerf(imPar, bidsPar, studyPar, subjectSessionLabel, inSessionPath, outSessionPath, listRuns, iRun);
    catch ME
        xASL_imp_NII2BIDS_Run_ISSUE_WARNING(ME, 'perfusion', subjectSessionLabel, iRun);
    end
    
end


%% ============================================================================================'

function xASL_imp_NII2BIDS_Run_ISSUE_WARNING(ME, Scantype, subjectSessionLabel, iRun)

    fprintf('\n\n\n%s\n', '============================================================================================');
    warning(['NII2BIDS went wrong for ' Scantype ' ' subjectSessionLabel '_run-' xASL_num2str(iRun)]);
    fprintf('\n%s\n', 'Error message:');
    fprintf('%s\n', ME.message);
    fprintf('\n%s\n', 'Continuing...');
    fprintf('%s\n\n\n\n', '============================================================================================');
    
end