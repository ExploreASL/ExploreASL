function [x] = xASL_imp_DetermineStructureFromSourcedata(x)
%xASL_imp_DetermineStructureFromSourcedata Determine structure from sourcedata
%
% FORMAT: [x] = xASL_imp_DetermineStructureFromSourcedata(x)
%
% INPUT:
%   x                       - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   x.modules.import.tokens - Tokens for matching files/directories (REQUIRED)
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Determine structure from sourcedata.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2023 ExploreASL


    %% Read sourcedata
    x = xASL_imp_ReadSourceData(x);

	% Report missing tokenOrdering field
	if ~isfield(x.modules.import.imPar, 'tokenOrdering') || isempty(x.modules.import.imPar.tokenOrdering)
		error('tokenOrdering parameter not specified or empty');
	end
	
    %% VISITS
    if x.modules.import.imPar.tokenOrdering(2)==0
        % a zero means: no visits applicable
        x.modules.import.settings.bUseVisits = false;
        % vVisitIDs: each subject has a single visit
        x.modules.import.listsIDs.vVisitIDs = cellfun(@(y) '1', x.modules.import.listsIDs.vSubjectIDs, 'UniformOutput', false);
        x.modules.import.imPar.tokenVisitAliases = {'^1$', '_1'};
    else
        x.modules.import.settings.bUseVisits = true;
        % vVisitIDs: cell vector with extracted session IDs (for all subjects, sessions and scans)
        x.modules.import.listsIDs.vVisitIDs = x.modules.import.tokens(:, x.modules.import.imPar.tokenOrdering(2));
    end
    
    %% SESSIONS
    if x.modules.import.imPar.tokenOrdering(3)==0
        % a zero means: no sessions applicable
        x.modules.import.settings.bUseSessions = false;
        % vSessionIDs: each subject-visit has a single session
        x.modules.import.listsIDs.vSessionIDs = cellfun(@(y) '1', x.modules.import.listsIDs.vSubjectIDs, 'UniformOutput', false);
        x.modules.import.imPar.tokenSessionAliases = {'^1$', 'ASL_1'};
    else
        x.modules.import.settings.bUseSessions = true;
        % vSessionIDs: Cell vector with extracted session IDs (for all subjects and scans)
        x.modules.import.listsIDs.vSessionIDs = x.modules.import.tokens(:, x.modules.import.imPar.tokenOrdering(3));
    end
    
    %% SCANTYPES
    
    % vScanIDs: cell vector with extracted scan IDs (for all subjects, visits and sessions)
    x.modules.import.listsIDs.vScanIDs = x.modules.import.tokens(:, x.modules.import.imPar.tokenOrdering(4)); 
    
    % Convert the vectors to unique & sort sets by the output aliases
    x.modules.import.listsIDs.subjectIDs  = sort(unique(x.modules.import.listsIDs.vSubjectIDs));
    x.modules.import.nSubjects = length(x.modules.import.listsIDs.subjectIDs);
    
    % Determine x.SUBJECTS (required for xASL_Iteration)
    x.SUBJECTS = x.modules.import.listsIDs.subjectIDs;
    
    % Overview: get visits and sessions of each subject
    x.importOverview = struct;
    x = xASL_imp_AddSubjectOverview(x);
    
    % Print matching sourcedata
    fprintf('\n[\b========================================= SOURCEDATA =========================================]\b\n');
    
    % Print matching files
    if isfield(x.modules.import.imPar,'bVerbose') && x.modules.import.imPar.bVerbose
        fprintf('\nMatching files (#=%g):\n',length(x.modules.import.matches));
        for iMatch=1:size(x.modules.import.matches,1)
            fprintf('%s\n', x.modules.import.matches{iMatch,1});
        end
    end
    
    % Print overview
    xASL_imp_PrintOverview(x);


end


%% Print overview
function xASL_imp_PrintOverview(x)

    if isfield(x,'importOverview')
        % Get all main fields
        overviewFields = fieldnames(x.importOverview);
        for iField=1:numel(overviewFields)
            thisSubject = overviewFields{iField};
            if ~isempty(regexpi(thisSubject, 'subject_'))
                xASL_imp_PrintSubject(x, thisSubject);
            end
        end
    end
    
    fprintf('\n');

end
function xASL_imp_PrintSubject(x,thisSubject)

    % Print each individual subject
    if isfield(x.importOverview.(thisSubject),'name')
        subject = x.importOverview.(thisSubject).name;
    else
        subject = '';
    end
    fprintf('\nSubject: %s\n', subject);
    % Show subject content
    subjectLevelFields = fieldnames(x.importOverview.(thisSubject));
    for iSubjectField=1:numel(subjectLevelFields)
        thisSession = subjectLevelFields{iSubjectField};
        if ~isempty(regexpi(thisSession, 'visit_'))
            xASL_imp_PrintSession(x, thisSubject, thisSession);
        end
    end

end
function xASL_imp_PrintSession(x, thisSubject, thisSession)

    % Print each individual session (=visits in legacy terminology)
    if isfield(x.importOverview.(thisSubject).(thisSession),'name')
        session = x.importOverview.(thisSubject).(thisSession).name;
    else
        session = '';
    end
    fprintf('Session: %s\n', session);
    
    % Print runs (= session in legacy terminology)
    sessionLevelFields = fieldnames(x.importOverview.(thisSubject).(thisSession));
    for iSessionField=1:numel(sessionLevelFields)
        thisRun = sessionLevelFields{iSessionField};
        if ~isempty(regexpi(thisRun,'run_'))
            xASL_imp_PrintRun(x,thisSubject, thisSession, thisRun);
        end
    end

end
function xASL_imp_PrintRun(x,thisSubject,thisVisit,thisRun)

    % Print each individual run (= session in legacy terminology)
    if isfield(x.importOverview.(thisSubject).(thisVisit).(thisRun),'name')
        run = x.importOverview.(thisSubject).(thisVisit).(thisRun).name;
    else
        run = '';
    end
    fprintf('Run:     %s\n', run);

end



%% Add subjects to overview
function x = xASL_imp_AddSubjectOverview(x)
    
    for iSubject=1:numel(x.modules.import.listsIDs.subjectIDs)
        thisSubject = x.modules.import.listsIDs.subjectIDs{iSubject};
        x = xASL_imp_AddSubject(x, thisSubject, iSubject);
    end

end


%% Add single subject to overview
function x = xASL_imp_AddSubject(x,thisSubject,iSubject)

    % Add subject name
    subjectFieldName = ['subject_' num2str(iSubject,'%03.f')];
    x.importOverview.(subjectFieldName).name = thisSubject;

    % Get vSubjectIDs
    vSubjectIDs = strcmp(x.modules.import.listsIDs.vSubjectIDs,thisSubject);
    
    % Get visits of current subject
    currentVisitList = unique(x.modules.import.listsIDs.vVisitIDs(vSubjectIDs));
    
    % Determine visitIDs
    x.importOverview.(subjectFieldName).visitIDs = unique(x.modules.import.listsIDs.vVisitIDs(vSubjectIDs));
    
    for iVisit=1:numel(currentVisitList)
        thisVisit = currentVisitList{iVisit};
        x = xASL_imp_AddVisit(x,subjectFieldName,vSubjectIDs,thisVisit,iVisit);
    end
    

end



%% Add single visit to overview
function x = xASL_imp_AddVisit(x, sFieldName, vSubjectIDs, thisVisit, iVisit)


    %% Add basic visit fields

    % Determine visit field name
    vFieldName = ['visit_' num2str(iVisit, '%03.f')];
    
    % Get vVisitIDs
    vVisitIDs = strcmp(x.modules.import.listsIDs.vVisitIDs, thisVisit);
    vVisitIDs = vSubjectIDs & vVisitIDs;
    
    % Assign visit name and sessions
    x.importOverview.(sFieldName).(vFieldName).name = thisVisit;
    x.importOverview.(sFieldName).(vFieldName).sessions = unique(x.modules.import.listsIDs.vSessionIDs(vVisitIDs));
    
    % Get number of visits
    x.importOverview.(sFieldName).nVisits = length(x.importOverview.(sFieldName).visitIDs);
    
    % Sort by output
    if length(x.importOverview.(sFieldName).visitIDs)>1
        for iV=1:numel(x.importOverview.(sFieldName).visitIDs)
          	idVisit = cellfun(@(y) regexp(y, x.importOverview.(sFieldName).visitIDs{iV}), x.modules.import.imPar.tokenVisitAliases(:,1), 'UniformOutput', false);
			idVisit = find(cellfun(@(y) ~isempty(y), idVisit));
			if isempty(idVisit)
				error('Visit not identified');
			else
				IDrow(iV) = idVisit(1);
			end
        end
        %x.importOverview.(sFieldName).listsIDs.visitIDs = x.importOverview.(sFieldName).visitIDs(IDrow);
		x.importOverview.(sFieldName).listsIDs.visitIDs = x.modules.import.imPar.tokenVisitAliases(IDrow,1);
    end
    
    % Add additional fields
    x = xASL_imp_thisSubjectVisit(x,sFieldName,vVisitIDs,vFieldName);
    
    
    %% Sanity check for missing elements
    xASL_imp_DCM2NII_SanityChecks(x,x.importOverview.(sFieldName),x.importOverview.(sFieldName).(vFieldName));
    
    
    %% Preallocate space for (global) counts
    [x.importOverview.(sFieldName),x.importOverview.(sFieldName).(vFieldName)] = ...
        xASL_imp_PreallocateGlobalCounts(x.modules.import.nSubjects,x.importOverview.(sFieldName),x.importOverview.(sFieldName).(vFieldName));
    

end


%% Add fields of this subject/visit
function x = xASL_imp_thisSubjectVisit(x,sFieldName,vVisitIDs,vFieldName)


    % Visits (= sessions in BIDS)
    x = xASL_imp_AddVisitNames(x,sFieldName);
    
    
    % Sessions (= run in BIDS)
    x = xASL_imp_AddSessionNames(x,sFieldName,vFieldName,vVisitIDs);
    
    
    % Scans (= actual DICOM data)
    x = xASL_imp_AddScanNames(x,sFieldName,vFieldName,vVisitIDs);
    

    % Sessions (= runs in BIDS)
    x = xASL_imp_AddSessions(x,sFieldName,vFieldName);


end


%% Add visit names
function x = xASL_imp_AddVisitNames(x, sFieldName)

    if isempty(x.modules.import.imPar.visitNames)
        if isempty(x.importOverview.(sFieldName).visitIDs)
            x.modules.import.imPar.visitNames = cell(x.importOverview.(sFieldName).nVisits,1);
            for iVisit=1:x.importOverview.(sFieldName).nVisits
                x.modules.import.imPar.visitNames{iVisit} = sprintf('ASL_%g', iVisit);
            end
        else
            for iVisit=1:numel(x.importOverview.(sFieldName).visitIDs)
				% Find the name id of the visit according to its number
				idVisit = cellfun(@(y) regexp(y, x.importOverview.(sFieldName).visitIDs{iVisit}), x.modules.import.imPar.tokenVisitAliases(:,1), 'UniformOutput', false);
				idVisit = find(cellfun(@(y) ~isempty(y), idVisit));
				if isempty(idVisit)
					error('Visit not identified');
				else
					idVisit = idVisit(1);
				end

				% Resolve the new name of the visit
				nameVisit = x.modules.import.imPar.tokenVisitAliases{idVisit, 2};
				
				% Extract the number of the visit
				numberVisit = regexp(nameVisit, '\d+');
				if ~isempty(numberVisit)
					numberVisit = nameVisit(numberVisit:end);
				else
					numberVisit = sprintf('%g', idVisit);
				end
				x.modules.import.imPar.visitNames{iVisit} = sprintf('ASL_%s', numberVisit);
            end
        end
    end

end


%% Add session names
function x = xASL_imp_AddSessionNames(x,sFieldName,vFieldName,vVisitIDs)

    % Get number of session IDs, number of sessions, scan IDs, number of scans
    x.importOverview.(sFieldName).(vFieldName).sessionIDs  = sort(unique(x.modules.import.listsIDs.vSessionIDs(vVisitIDs)));

    if isempty(x.modules.import.imPar.sessionNames)
        % We can have human readble session names, but by default they are the same as the original tokens in the path
        if isempty(x.importOverview.(sFieldName).(vFieldName).sessionIDs)
            x.modules.import.imPar.sessionNames = cell(x.importOverview.(sFieldName).(vFieldName).nSessions,1);
            for iSession = 1:x.importOverview.(sFieldName).(vFieldName).nSessions
                x.modules.import.imPar.sessionNames{iSession} = sprintf('ASL_%g',iSession);
            end
        else
            for iSession = 1:numel(x.importOverview.(sFieldName).(vFieldName).sessionIDs)
                x.modules.import.imPar.sessionNames{iSession} = sprintf('ASL_%g',iSession);
            end
        end
    end

end


%% Add scan names
function x = xASL_imp_AddScanNames(x,sFieldName,vFieldName,vVisitIDs)

    % Add scan IDs and number of scans
    x.importOverview.(sFieldName).(vFieldName).scanIDs = sort(unique(lower(x.modules.import.listsIDs.vScanIDs(vVisitIDs))));
    x.importOverview.(sFieldName).(vFieldName).nScans = length(x.importOverview.(sFieldName).(vFieldName).scanIDs);
    x.importOverview.(sFieldName).(vFieldName).scanNames = x.importOverview.(sFieldName).(vFieldName).scanIDs;
    x.importOverview.(sFieldName).(vFieldName).scanSessionTable = horzcat(x.modules.import.listsIDs.vSessionIDs(vVisitIDs),lower(x.modules.import.listsIDs.vScanIDs(vVisitIDs)));

end


%% Add sessions
function x = xASL_imp_AddSessions(x,sFieldName,vFieldName)

    % Add runs fields
    if ~isfield(x.importOverview.(sFieldName).(vFieldName),'runs')
        x.importOverview.(sFieldName).(vFieldName).runs = {};
    end
    
    % If there are no tokenSessionAliases, we just say that each existing
    % and matched sessionID is a separate run...
    if ~isfield(x.modules.import.imPar,'tokenSessionAliases')
        numOfSessions = 0;
        
        % Iterate over all session names
        for iSession=1:numel(x.importOverview.(sFieldName).(vFieldName).sessions)
            thisSession = x.importOverview.(sFieldName).(vFieldName).sessions{iSession};
            x = xASL_imp_AddRun(x,sFieldName,vFieldName,thisSession,iSession);
            numOfSessions = numOfSessions+1;
        end
        
        % Assign x field
        x.importOverview.(sFieldName).(vFieldName).nSessions = numOfSessions;
    end
    
    % ... but if there is a tokenSessionAliases field, we can group the
    % sessions based on those tokens. We have to iterate over the tokens
    % and see if they were matched. If they were matched, we assign a name
    % and then we can assign all sessionIDs to their corresponding token.
    if isfield(x.modules.import.imPar,'tokenSessionAliases')
        numOfSessions = 0;
        tokenTable = x.modules.import.imPar.tokenSessionAliases;
        % Add empty last row for session (=run) names
        for iToken=1:size(tokenTable,1)
            tokenTable{iToken,3} = '';
        end
        
        % Iterate over session (=run) tokens
        for iToken=1:size(tokenTable,1)
            currentToken = tokenTable{iToken,1};
            % Iterate over session names
            for iSession=1:numel(x.importOverview.(sFieldName).(vFieldName).sessions)
                currentID = x.importOverview.(sFieldName).(vFieldName).sessions{iSession};
                if ~isempty(regexpi(currentID,currentToken))
                    % Check if session (=run) already exists
                    if isempty(tokenTable{iToken,3})
                        numOfSessions = numOfSessions+1;
                        tokenTable{iToken,3} = ['run_' num2str(numOfSessions,'%03.f')];
                    end
                end
            end
        end
        
        % Create runs
        iSession = 0;
        for iToken=1:size(tokenTable,1)
            if ~isempty(tokenTable{iToken,3})
                iSession=iSession+1;
                thisSession = tokenTable{iToken,2};
                x = xASL_imp_AddRun(x,sFieldName,vFieldName,thisSession,iSession,tokenTable{iToken,1});
            end
        end
        
        % Assign x field
        x.importOverview.(sFieldName).(vFieldName).nSessions = numOfSessions;
    end
    

end



%% Add runs to overview
function x = xASL_imp_AddRun(x,sFieldName,vFieldName,thisSession,iSession,thisRegExp)

    % Check if there was a regexp and token
    if nargin<6
        thisRegExp = '';
    end

    % Add field to overview
    vSessionName = ['run_' num2str(iSession,'%03.f')];
    % Make sure that the name starts with ASL_
    %if ~isempty(regexp(thisSession,'ASL_', 'once'))
        x.importOverview.(sFieldName).(vFieldName).(vSessionName).name = thisSession;
    %else
    %    x.importOverview.(sFieldName).(vFieldName).(vSessionName).name = ['ASL_' num2str(iSession)];
    %end
    x.importOverview.(sFieldName).(vFieldName).runs = vertcat(x.importOverview.(sFieldName).(vFieldName).runs,{thisSession});
    x.importOverview.(sFieldName).(vFieldName).(vSessionName).regexp = thisRegExp;
    
    % Add list of IDs
    if ~isempty(thisRegExp)
		x = xASL_imp_AddSessionIDListAndScansOfRun(x, sFieldName, vFieldName, vSessionName, thisRegExp);
	else
		x = xASL_imp_AddSessionIDListAndScansOfRun(x, sFieldName, vFieldName, vSessionName, thisSession);
    end


end

function x = xASL_imp_AddSessionIDListAndScansOfRun(x, sFieldName, vFieldName, vSessionName, thisRegExp)
% Iterate over scan IDs
x.importOverview.(sFieldName).(vFieldName).(vSessionName).scanIDs = {};
scanSessionTable = x.importOverview.(sFieldName).(vFieldName).scanSessionTable;
iSessionNew = 1;
for iID=1:size(scanSessionTable,1)
	thisToken = scanSessionTable{iID,1};
	thisScanID = scanSessionTable{iID,2};
	% Iterate over the sessions of this run
	for iSession=1:numel(x.importOverview.(sFieldName).(vFieldName).sessions)
		currentID = x.importOverview.(sFieldName).(vFieldName).sessions{iSession};
		if ~isempty(regexpi(currentID, thisRegExp)) && strcmp(currentID, thisToken)
			x.importOverview.(sFieldName).(vFieldName).(vSessionName).scanIDs{iSessionNew} = thisScanID;
			x.importOverview.(sFieldName).(vFieldName).(vSessionName).ids{iSessionNew} = currentID;
			iSessionNew = iSessionNew + 1;
			continue;
		end
	end
end

end
