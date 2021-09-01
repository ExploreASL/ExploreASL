function x = xASL_imp_NII2BIDS_Run(x, imPar, bidsPar, studyPar, subjectSessionLabel, inSessionPath, outSessionPath, listRuns, iRun, nameSubjectSession)
%xASL_imp_NII2BIDS_Run NII2BIDS conversion for a single sessions, single run.
%
% FORMAT: xASL_imp_NII2BIDS_Run(bidsPar, studyPar, subjectSessionLabel, inSessionPath, outSessionPath, listRuns, iRun)
% 
% INPUT:
% x                   - ExploreASL x structure (REQUIRED, STRUCT)
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
% x                   - ExploreASL x structure (STRUCT)
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
    catch loggingEntry
        [x] = xASL_qc_AddLoggingInfo(x, loggingEntry);
        xASL_imp_NII2BIDS_RunIssueWarning(loggingEntry, 'anatomical', subjectSessionLabel, iRun);
    end
    
    %% Perfusion files
    try
        xASL_imp_NII2BIDS_RunPerf(imPar, bidsPar, studyPar, subjectSessionLabel, inSessionPath, outSessionPath, listRuns, iRun);
    catch loggingEntry
        [x] = xASL_qc_AddLoggingInfo(x, loggingEntry);
        xASL_imp_NII2BIDS_RunIssueWarning(loggingEntry, 'perfusion', subjectSessionLabel, iRun);
    end
    
end


%% Issue failed run as a warning
function xASL_imp_NII2BIDS_RunIssueWarning(loggingEntry, Scantype, subjectSessionLabel, iRun)

    fprintf('\n==============================================================================================\n');
    fprintf(2,'NII2BIDS failed for %s image of %s_run-%s\n',Scantype,subjectSessionLabel,xASL_num2str(iRun));
    if size(loggingEntry.stack,1)>0
        fprintf(2,'Message: %s\n%s, line %d...\n',loggingEntry.message,loggingEntry.stack(1).name,loggingEntry.stack(1).line);
    end
    fprintf(2,'Continuing...\n');
    
end


