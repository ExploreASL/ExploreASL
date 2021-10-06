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
    x.modules.import.numOf.nSubjects = length(x.modules.import.listsIDs.subjectIDs);
    x.modules.import.listsIDs.visitIDs  = unique(x.modules.import.listsIDs.vVisitIDs);
    
    % Sort by output
    if length(x.modules.import.listsIDs.visitIDs)>1
        for iV=1:length(x.modules.import.listsIDs.visitIDs)
            IDrow(iV) = find(cellfun(@(y) strcmp(y,x.modules.import.listsIDs.visitIDs{iV}), x.modules.import.imPar.tokenVisitAliases(:,1)));
        end
        x.modules.import.listsIDs.visitIDs = x.modules.import.listsIDs.visitIDs(IDrow);
    end
    
    % Get number of visits, session IDs, number of sessions, scan IDs, number of scans
    x.modules.import.numOf.nVisits = length(x.modules.import.listsIDs.visitIDs);
    x.modules.import.listsIDs.sessionIDs  = sort(unique(x.modules.import.listsIDs.vSessionIDs));
    x.modules.import.numOf.nSessions = length(x.modules.import.listsIDs.sessionIDs);
    x.modules.import.listsIDs.scanIDs = sort(unique(lower(x.modules.import.listsIDs.vScanIDs)));
    x.modules.import.numOf.nScans = length(x.modules.import.listsIDs.scanIDs);
    
    %% VISIT NAMES (== session in BIDS)
    if isempty(x.modules.import.imPar.visitNames)
        if isempty(x.modules.import.listsIDs.visitIDs)
            x.modules.import.imPar.visitNames = cell(x.modules.import.numOf.nVisits,1);
            for kk=1:x.modules.import.numOf.nVisits
                x.modules.import.imPar.visitNames{kk} = sprintf('ASL_%g', kk);
            end
        else
            x.modules.import.imPar.visitNames = x.modules.import.listsIDs.visitIDs;
        end
    end
    
    %% SESSION NAMES (== run in BIDS)
    % optionally we can have human readble session names; by default they are the same as the original tokens in the path
    if isempty(x.modules.import.imPar.sessionNames)
        if isempty(x.modules.import.listsIDs.sessionIDs)
            x.modules.import.imPar.sessionNames = cell(x.modules.import.numOf.nSessions,1);
            for kk=1:x.modules.import.numOf.nSessions
                x.modules.import.imPar.sessionNames{kk}=sprintf('ASL_%g',kk);
            end
        else
            x.modules.import.imPar.sessionNames = x.modules.import.listsIDs.sessionIDs;
        end
    end
    
    %% SCAN NAMES
    x.modules.import.scanNames = x.modules.import.listsIDs.scanIDs;
    
    
    %% SUBJECTS
    x.SUBJECTS = x.modules.import.listsIDs.subjectIDs;
    
    
    % Sanity check for missing elements
    xASL_imp_DCM2NII_SanityChecks(x);
    
    % Preallocate space for (global) counts
    x = xASL_imp_PreallocateGlobalCounts(x);


end



