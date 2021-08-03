function xASL_imp_NII2BIDS_Subject(imPar, bidsPar, studyPar, nameSubjectSession)
%xASL_imp_NII2BIDS_Subject Run NII to ASL-BIDS for one individual subject.
%
% FORMAT: xASL_imp_NII2BIDS_Subject(imPar, bidsPar, studyPar, nameSubject)
% 
% INPUT:
%   imPar               - JSON file with structure with import parameter (REQUIRED, STRUCT)
%   bidsPar             - Output of xASL_imp_Config (REQUIRED, STRUCT)
%   studyPar            - JSON file with the BIDS parameters relevant for the whole study (REQUIRED, STRUCT)
%   nameSubjectSession  - name of the subject (REQUIRED, CELL STRUCT)
%
% OUTPUT:
%   n/a
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
% EXAMPLE:     xASL_imp_NII2BIDS_Subject(imPar, bidsPar, studyPar, nameSubject);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% 1. Initialize
    bidsLabel = xASL_imp_CheckForAliasInVisit(imPar,nameSubjectSession);
    
    % Make a subject directory
    xASL_adm_CreateDir(fullfile(imPar.BidsRoot,['sub-' bidsLabel.subject]));
    
    
    %% 2. Process the anat & perfusion files
    listSessions = xASL_adm_GetFileList(fullfile(imPar.TempRoot,nameSubjectSession),'^ASL.+$',false,[],true);
    
    % Go through all sessions
    for iSession = 1:length(listSessions)
        xASL_imp_NII2BIDS_Session(imPar, bidsPar, studyPar, listSessions, nameSubjectSession, bidsLabel, iSession);
    end
    
end



%% Check if there is a visit alias within the subject/visit name
function bidsLabel = xASL_imp_CheckForAliasInVisit(imPar,nameSubjectSession)

    % Get visitAliases from imPar
    visitAliases = imPar.tokenVisitAliases(:,2);
    
    % Separator subject/visit
    separator = '_';

    % Default
    subjectName = nameSubjectSession;
    visitName = '';

    % Iterate over aliases
    for iAlias = 1:numel(visitAliases)
        checkExpression = regexp(nameSubjectSession, [separator visitAliases{iAlias,1} '$'], 'once');
        if ~isempty(checkExpression) % nameSubject should end in the visit alias
            visitName = nameSubjectSession(checkExpression:end);
            subjectName = nameSubjectSession(1:checkExpression-1);
        end
    end
    
    bidsLabel.subject = xASL_adm_CorrectName(subjectName,2);
    bidsLabel.visit = xASL_adm_CorrectName(visitName,2);


end


