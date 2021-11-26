function x = xASL_imp_NII2BIDS_Subject(x, imPar, bidsPar, studyPar, nameSubjectSession)
%xASL_imp_NII2BIDS_Subject Run NII to ASL-BIDS for one individual subject.
%
% FORMAT: x = xASL_imp_NII2BIDS_Subject(x, imPar, bidsPar, studyPar, nameSubject)
% 
% INPUT:
%   x                   - ExploreASL x structure (REQUIRED, STRUCT)
%   imPar               - JSON file with structure with import parameter (REQUIRED, STRUCT)
%   bidsPar             - Output of xASL_imp_Config (REQUIRED, STRUCT)
%   studyPar            - JSON file with the BIDS parameters relevant for the whole study (REQUIRED, STRUCT)
%   nameSubjectSession  - name of the subject (REQUIRED, CELL STRUCT)
%
% OUTPUT:
%   x               - ExploreASL x structure (STRUCT)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run NII to ASL-BIDS for one individual subject.
%
% 1. Initialize
% 2. Process the anat & perfusion files
% - 1. Make a subject directory
% - 2. Iterate over sessions
% - 3. Iterate over runs
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     x = xASL_imp_NII2BIDS_Subject(x, imPar, bidsPar, studyPar, nameSubject);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% 1. Initialize
    bidsLabel = xASL_imp_CheckForAliasInVisit(imPar,nameSubjectSession);
    
    % Make a subject directory
    subjectDirectory = fullfile(imPar.BidsRoot,['sub-' bidsLabel.subject]);
    if xASL_exist(subjectDirectory,'dir')
        % Subject already exists and this is probably either a rerun or the pipeline crashed before
        fprintf(2,'Remove existing rawdata...\n');
        xASL_delete(subjectDirectory,true);
    end
    xASL_adm_CreateDir(subjectDirectory);
    
    
    %% 2. Process the anat & perfusion files
    % Subsequent code is based on having data per ASL scan, so we "fool" it
    % by renaming all runs into ASL_1 ASL_2 ASL_n and keeping the
    % unique runs only. Missing scans will issue a warning, not an error.
    listRuns = xASL_adm_GetFileList(fullfile(imPar.TempRoot,nameSubjectSession),'^(ASL|T1w|FLAIR).+$',false,[],true);
    listRuns = cellfun(@(y) y(end), listRuns, 'UniformOutput', false);
    listRuns = unique(listRuns);
    listRuns = cellfun(@(y) ['ASL_' y], listRuns, 'UniformOutput', false);
    
    % Go through all (ASL) runs
    for iRun = 1:length(listRuns)
        x = xASL_imp_NII2BIDS_Run(x, imPar, bidsPar, studyPar, listRuns, nameSubjectSession, bidsLabel, iRun);
    end
    
end



%% Check if there is a visit alias within the subject/visit name
function bidsLabel = xASL_imp_CheckForAliasInVisit(imPar,nameSubjectSession)

    % Get visitAliases from imPar
    if isfield(imPar,'tokenVisitAliases') && ~isempty(imPar.tokenVisitAliases) && size(imPar.tokenVisitAliases,2)>1
        visitAliases = imPar.tokenVisitAliases(:,2);
    else
        visitAliases = [];
    end
    
    % Separator subject/visit
    separator = '_';

    % Default
    subjectName = nameSubjectSession;
    visitName = '';

    % Iterate over aliases
    if ~isempty(visitAliases)
        for iAlias = 1:numel(visitAliases)
			if ~isempty(visitAliases{iAlias,1})
				if visitAliases{iAlias,1}(1) == separator
					checkExpression = regexp(nameSubjectSession, [visitAliases{iAlias,1} '$'], 'once');
				else
					checkExpression = regexp(nameSubjectSession, [separator visitAliases{iAlias,1} '$'], 'once');
				end
				if ~isempty(checkExpression) % nameSubject should end in the visit alias
					visitName = nameSubjectSession(checkExpression:end);
					subjectName = nameSubjectSession(1:checkExpression-1);
				end
			end
        end
    end
    
    bidsLabel.subject = xASL_adm_CorrectName(subjectName,2);
    bidsLabel.visit = xASL_adm_CorrectName(visitName,2);


end


