function [x] = xASL_init_DefineStudyStats(x)
%xASL_init_DefineStudyStats Define study subjects/parameters for this
% pipeline run.
%
% FORMAT: [x] = xASL_init_DefineStudyStats(x)
%
% INPUT:
%   x             - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x   - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This initialization wrapper initializes the statistical parameters for
%               this pipeline run. This function exists from the following parts:
%
% 9. Load & add statistical variables
% 10. Make x.S.SetsOptions horizontal if vertical by transposing
% 11. Create single site dummy, if there were no sites specified
% 14. Check what excluded from which TimePoints
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE: x = xASL_init_DefineStudyStats(x);
% __________________________________
% Copyright (c) 2015-2023 ExploreASL


% ------------------------------------------------------------------------------------------------



% ------------------------------------------------------------------------------------------------
%% 7) Add sessions as statistical variable

% Predefine SETS to avoid empty SETS & import predefined session settings as set settings
x.S.SetsName{1} = 'session';
x.S.SetsOptions = '';
x.S.SetsID = 0;
x.S.Sets1_2Sample(1) = 1; % sessions are always paired observations (e.g. belonging to same individual, looking at intra-individual changes)

% NOT SURE IF WE STILL USE THIS, WE MAY PHASE THIS OUT
% WE DO USE THIS FOR FEAST: CRUSHED VS NON-CRUSHED ASL
% Define SetsOptions for sessions
if ~isfield(x,'session')
    x.session = '';
end

if  isfield(x.session,'options')
    x.S.SetsOptions{1} = x.session.options; % with session options (e.g. morning, evening 2nd morning)
else
    for iS=1:x.dataset.nSessions
        x.S.SetsOptions{1}{iS} = ['ASL_' num2str(iS)];
    end
end        

% Create ID/data for session numbers
for iSubj=1:x.dataset.nSubjects
    for iSess=1:x.dataset.nSessions
        iSubjSess = (iSubj-1)*x.dataset.nSessions+iSess;

        x.S.SetsID(iSubjSess,1) = iSess;
    end
end

% ------------------------------------------------------------------------------------------------
%% 8) Add Longitudinal TimePoints (different T1 volumes) as covariate/set
% TimePoints should be indicated as different subjects, with a _1 _2 _3 _n suffix

[SubjectNList, TimePoint] = xASL_init_LongitudinalRegistration(x);

if isfield(x.S,'SetsID')
    iSetLong_TP = size(x.S.SetsID,2)+1;
else
    iSetLong_TP = 1;
end

x.S.SetsName{iSetLong_TP} = 'LongitudinalTimePoint';
x.S.SetsName{iSetLong_TP+1} = 'SubjectNList';
x.S.SetsID(:,iSetLong_TP) = TimePoint;
x.S.SetsID(:,iSetLong_TP+1) = SubjectNList;

for iT=1:max(TimePoint)
    x.S.SetsOptions{iSetLong_TP}{iT} = ['TimePoint_' num2str(iT)];
end

x.S.SetsOptions{iSetLong_TP+1} = {'SubjectNList'};
x.S.Sets1_2Sample(iSetLong_TP:iSetLong_TP+1) = 1; % Longitudinal TimePoints are always paired observations (e.g. belonging to same individual, looking at intra-individual changes)    
x.S.iSetLong_TP = iSetLong_TP;










% ================================================================================================
%% 9) Load & add statistical variables
% Keep some space here, to easily find this for debugging
x = xASL_init_LoadMetadata(x);







% ================================================================================================
%% 10) Make x.S.SetsOptions horizontal if vertical by transposing
if isfield(x.S,'SetsOptions')
    for iCell=1:length(x.S.SetsOptions)
        if size(x.S.SetsOptions{iCell},1)>size(x.S.SetsOptions{iCell},2)
            x.S.SetsOptions{iCell} = x.S.SetsOptions{iCell}';
        end
    end
else
    x.S.SetsOptions = {};
end


% ------------------------------------------------------------------------------------------------
%% 11) Create single site dummy, if there were no sites specified
if isfield(x.S, 'SetsName')
    % Search for a Site variable
    for iS=1:length(x.S.SetsName)
        if  strcmp(x.S.SetsName{iS},'Site')
            x.S.iSetSite=iS;
        end
    end
end
if ~isfield(x.S,'iSetSite')
    x.S.SetsName{end+1} = 'Site';
    x.S.SetsOptions{end+1} = {'SingleSite'};
    x.S.SetsID(:,end+1) = 1;
    x.S.Sets1_2Sample(end+1) = 2;
    x.S.iSetSite = length(x.S.SetsName);
end   


% ------------------------------------------------------------------------------------------------
%% 12) Check for time points sets
% Search for a variable specifying longitudinal time points
if isfield(x.S, 'SetsName')
    % Search for a variable specifying age at scanning
    iS = find(strcmp(x.S.SetsName,'Age'));
    if ~isempty(iS)
        x.S.iSetAge=iS;
    end
end




end
