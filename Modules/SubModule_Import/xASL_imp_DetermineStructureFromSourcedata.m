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
% Copyright 2015-2021 ExploreASL


    %% Read sourcedata
    x = xASL_imp_ReadSourceData(x);


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
        x.modules.import.listsIDs.vVisitIDs = x.modules.import.tokens(:,x.modules.import.imPar.tokenOrdering(2));
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
        x.modules.import.listsIDs.vSessionIDs = x.modules.import.tokens(:,x.modules.import.imPar.tokenOrdering(3));
    end
    
    %% SCANTYPES
    
    % vScanIDs: cell vector with extracted scan IDs (for all subjects, visits and sessions)
    x.modules.import.listsIDs.vScanIDs = x.modules.import.tokens(:,x.modules.import.imPar.tokenOrdering(4)); 
    
    % Convert the vectors to unique & sort sets by the output aliases
    x.modules.import.listsIDs.subjectIDs  = sort(unique(x.modules.import.listsIDs.vSubjectIDs));
    x.modules.import.nSubjects = length(x.modules.import.listsIDs.subjectIDs);
    
    % Determine x.SUBJECTS (required for xASL_Iteration)
    x.SUBJECTS = x.modules.import.listsIDs.subjectIDs;
    
    % Overview: get visits and sessions of each subject
    x.overview = struct;
    x = xASL_imp_AddSubjectOverview(x);
    
    % Print matching sourcedata
    fprintf('\n========================================= SOURCEDATA =========================================\n');
    
    % Print matching files
    if isfield(x.modules.import.imPar,'bVerbose') && x.modules.import.imPar.bVerbose
        fprintf('\nMatching files (#=%g):\n',length(x.modules.import.matches));
        for iMatch=1:size(x.modules.import.matches,1)
            fprintf('%s\n', x.modules.import.matches{iMatch,1});
        end
    end
    fprintf('\n');
    
    % Print overview
    xASL_imp_PrintOverview(x);


end


%% Print overview
function xASL_imp_PrintOverview(x)

    if isfield(x,'overview')
        % Get all main fields
        overviewFields = fieldnames(x.overview);
        for iField=1:numel(overviewFields)
            thisSubject = overviewFields{iField};
            if ~isempty(regexpi(thisSubject,'subject_'))
                xASL_imp_PrintSubject(x,thisSubject);
            end
        end
    end

end
function xASL_imp_PrintSubject(x,thisSubject)

    % Print each individual subject
    if isfield(x.overview.(thisSubject),'name')
        subject = x.overview.(thisSubject).name;
    else
        subject = '';
    end
    fprintf('Subject: %s\n', subject);
    % Show subject content
    subjectLevelFields = fieldnames(x.overview.(thisSubject));
    for iSubjectField=1:numel(subjectLevelFields)
        thisVisit = subjectLevelFields{iSubjectField};
        if ~isempty(regexpi(thisVisit,'visit_'))
            xASL_imp_PrintVisit(x,thisSubject,thisVisit);
        end
    end

end
function xASL_imp_PrintVisit(x,thisSubject,thisVisit)

    % Print each individual subject
    if isfield(x.overview.(thisSubject).(thisVisit),'name')
        visit = x.overview.(thisSubject).(thisVisit).name;
    else
        visit = '';
    end
    fprintf('Visit:   %s\n', visit);
    
    % Print sessions
    if isfield(x.overview.(thisSubject).(thisVisit),'sessions')
        for iSession=1:numel(x.overview.(thisSubject).(thisVisit).sessions)
            thisSession = x.overview.(thisSubject).(thisVisit).sessions{iSession};
            fprintf('Session: %s\n', thisSession);
        end
    end

end



%% Add subjects to overview
function x = xASL_imp_AddSubjectOverview(x)
    
    for iSubject=1:numel(x.modules.import.listsIDs.subjectIDs)
        thisSubject = x.modules.import.listsIDs.subjectIDs{iSubject};
        x = xASL_imp_AddSubject(x,thisSubject,iSubject);
    end

end


%% Add single subject to overview
function x = xASL_imp_AddSubject(x,thisSubject,iSubject)

    % Add subject name
    subjectFieldName = ['subject_' num2str(iSubject,'%03.f')];
    x.overview.(subjectFieldName).name = thisSubject;

    % Get vSubjectIDs
    vSubjectIDs = strcmp(x.modules.import.listsIDs.vSubjectIDs,thisSubject);
    
    % Get visits of current subject
    currentVisitList = unique(x.modules.import.listsIDs.vVisitIDs(vSubjectIDs));
    
    % Determine visitIDs
    x.overview.(subjectFieldName).visitIDs = unique(x.modules.import.listsIDs.vVisitIDs(vSubjectIDs));
    
    for iVisit=1:numel(currentVisitList)
        thisVisit = currentVisitList{iVisit};
        x = xASL_imp_AddVisit(x,subjectFieldName,vSubjectIDs,thisVisit,iVisit);
    end
    

end



%% Add single visit to overview
function x = xASL_imp_AddVisit(x,sFieldName,vSubjectIDs,thisVisit,iVisit)


    %% Add basic visit fields

    % Determine visit field name
    vFieldName = ['visit_' num2str(iVisit,'%03.f')];
    
    % Get vVisitIDs
    vVisitIDs = strcmp(x.modules.import.listsIDs.vVisitIDs,thisVisit);
    vVisitIDs = vSubjectIDs & vVisitIDs;
    
    % Assign visit name and sessions
    x.overview.(sFieldName).(vFieldName).name = thisVisit;
    x.overview.(sFieldName).(vFieldName).sessions = unique(x.modules.import.listsIDs.vSessionIDs(vVisitIDs));
    
    % Get number of visits
    x.overview.(sFieldName).nVisits = length(x.overview.(sFieldName).visitIDs);
    
    % Sort by output
    if length(x.overview.(sFieldName).visitIDs)>1
        for iV=1:numel(x.overview.(sFieldName).visitIDs)
            IDrow(iV) = find(cellfun(@(y) strcmp(y,x.overview.(sFieldName).visitIDs{iV}), x.modules.import.imPar.tokenVisitAliases(:,1)));
        end
        x.overview.(sFieldName).listsIDs.visitIDs = x.overview.(sFieldName).visitIDs(IDrow);
    end
    
    % Add additional fields
    x = xASL_imp_thisSubjectVisit(x,sFieldName,vVisitIDs,vFieldName);
    
    
    %% Sanity check for missing elements
    xASL_imp_DCM2NII_SanityChecks(x,x.overview.(sFieldName),x.overview.(sFieldName).(vFieldName));
    
    
    %% Preallocate space for (global) counts
    [x.overview.(sFieldName),x.overview.(sFieldName).(vFieldName)] = ...
        xASL_imp_PreallocateGlobalCounts(x.modules.import.nSubjects,x.overview.(sFieldName),x.overview.(sFieldName).(vFieldName));
    

end


%% Add fields of this subject/visit
function x = xASL_imp_thisSubjectVisit(x,sFieldName,vVisitIDs,vFieldName)


    %% SCAN NAMES
    
    % Get number of session IDs, number of sessions, scan IDs, number of scans
    x.overview.(sFieldName).(vFieldName).sessionIDs  = sort(unique(x.modules.import.listsIDs.vSessionIDs(vVisitIDs)));
    x.overview.(sFieldName).(vFieldName).nSessions = length(x.overview.(sFieldName).(vFieldName).sessionIDs);
    x.overview.(sFieldName).(vFieldName).scanIDs = sort(unique(lower(x.modules.import.listsIDs.vScanIDs(vVisitIDs))));
    x.overview.(sFieldName).(vFieldName).nScans = length(x.overview.(sFieldName).(vFieldName).scanIDs);
    
    
    %% VISIT NAMES (== session in BIDS)
    if isempty(x.modules.import.imPar.visitNames)
        if isempty(x.overview.(sFieldName).visitIDs)
            x.modules.import.imPar.visitNames = cell(x.overview.(sFieldName).nVisits,1);
            for kk=1:x.overview.(sFieldName).nVisits
                x.modules.import.imPar.visitNames{kk} = sprintf('ASL_%g', kk);
            end
        else
            x.modules.import.imPar.visitNames = x.overview.(sFieldName).visitIDs;
        end
    end
    
    %% SESSION NAMES (== run in BIDS)
    % optionally we can have human readble session names; by default they are the same as the original tokens in the path
    if isempty(x.modules.import.imPar.sessionNames)
        if isempty(x.overview.(sFieldName).(vFieldName).sessionIDs)
            x.modules.import.imPar.sessionNames = cell(x.overview.(sFieldName).(vFieldName).nSessions,1);
            for kk=1:x.overview.(sFieldName).(vFieldName).nSessions
                x.modules.import.imPar.sessionNames{kk}=sprintf('ASL_%g',kk);
            end
        else
            x.modules.import.imPar.sessionNames = x.overview.(sFieldName).(vFieldName).sessionIDs;
        end
    end
    
    
    %% SCAN NAMES
    x.overview.(sFieldName).(vFieldName).scanNames = x.overview.(sFieldName).(vFieldName).scanIDs;


end



