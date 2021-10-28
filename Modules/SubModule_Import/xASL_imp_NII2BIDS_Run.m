function x = xASL_imp_NII2BIDS_Run(x, imPar, bidsPar, studyPar, listRuns, nameSubjectSession, bidsLabel, iRun)
%xASL_imp_NII2BIDS_Run NII2BIDS conversion for a single run.
%
% FORMAT: x = xASL_imp_NII2BIDS_Run(x, imPar, bidsPar, studyPar, listRuns, nameSubjectSession, bidsLabel, iRun)
% 
% INPUT:
%   x                     - ExploreASL x structure (REQUIRED, STRUCT)
%   imPar                 - JSON file with structure with import parameter (STRUCT, REQUIRED)
%   bidsPar               - Output of xASL_imp_Config (STRUCT, REQUIRED)
%   studyPar              - JSON file with the BIDS parameters relevant for the whole study (STRUCT, REQUIRED)
%   listRuns              - list of runs (CELL ARRAY, REQUIRED, e.g.: {'ASL_1'})
%   nameSubjectSession    - name of the subject (CELL ARRAY, REQUIRED)
%   bidsLabel             - BIDS label (CHAR ARRAY, REQUIRED)
%   iRun                  - Run number (INTEGER, REQUIRED)
%
% OUTPUT:
%   x                     - ExploreASL x structure (STRUCT)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: NII2BIDS conversion for a single run.
% 
% 1. Make a subject directory with a correct session name
% 2. Convert structural and ASL runs
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     x = xASL_imp_NII2BIDS_Run(x, imPar, bidsPar, studyPar, listRuns, nameSubjectSession, bidsLabel, iRun);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Print the run that is being converted
    fprintf('\n======================================== CONVERT RUN =========================================\n');
    fprintf('Converting subject %s, session %s, run %s ', bidsLabel.subject, bidsLabel.visit, listRuns{iRun});

    %% 1. Make a subject directory with a correct session name
    if ~isempty(bidsLabel.visit)
        %if ~isempty(bidsLabel.visit)
        sessionLabel = ['ses-' bidsLabel.visit];
        %else
        %    sessionLabel = ['ses-' listRuns{iRun}(5:end)];
        %end
        
        % Create directory if it does not exist already
        sessionPerfusionDirectory = fullfile(imPar.BidsRoot,['sub-' bidsLabel.subject], sessionLabel,'perf');
        xASL_adm_CreateDir(sessionPerfusionDirectory);
        
        inSessionPath = fullfile(imPar.TempRoot, nameSubjectSession, listRuns{iRun});
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
        inSessionPath = fullfile(imPar.TempRoot, nameSubjectSession, listRuns{iRun});
        outSessionPath = fullfile(imPar.BidsRoot, ['sub-' bidsLabel.subject]);
    end


    %% 2. Convert structural and ASL runs
	subjectSessionLabel = [bidsLabel.subject sessionLabel];
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


