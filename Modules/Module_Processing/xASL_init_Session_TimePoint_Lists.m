function [x] = xASL_init_Session_TimePoint_Lists(x)
%xASL_init_Session_TimePoint_Lists Generate sessions & visit lists and adding to x struct
%
% FORMAT: [x] = xASL_init_Session_TimePoint_Lists(x)
% 
% INPUT:
%   x          - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x          - ExploreASL x structure
%                         
%
%   Note that we want to run this part AFTER parallelization.
%
%               Note that ASL sessions are defined here as what BIDS calls "runs".
%
%               The "longitudinal_Registration functions here manage different
%               TimePoints, which is what BIDS calls "visits".
%               With different structural scans, from the same participant. This is
%               managed by subject name suffixes _1 _2 _n, and can be used for comparing
%               visits in the population module, or running SPM's longitudinal within-subject
%               registration.
%
% 1. Manage sessions
% 2. Manage TimePoint lists
% 3. Create list of baseline & follow-up subjects (i.e. after exclusion)
% 4. Check what excluded from which TimePoints
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [x] = xASL_init_Session_TimePoint_Lists(x);
% __________________________________
% Copyright (c) 2015-2024 ExploreASL

% ------------------------------------------------------------------------------------------------
%% 1. Manage sessions
fprintf('Automatically defining sessions...\n');

if isfield(x,'SESSIONS') && isstruct(x.SESSIONS)
    warning('Invalid x.SESSIONS structure, replacing this now');
    fprintf('%s\n', 'Check that the correct number of sessions were processed');
    x = rmfield(x,'SESSIONS');
end

if isfield(x,'SESSIONS') && ~isempty(x.SESSIONS)
    warning('Using predefined x.SESSIONS');
else
    x.SESSIONS = [];
    SessionPathList = xASL_adm_GetFileList(x.dir.xASLDerivatives, '^(ASL|func)_\d+$', 'FPListRec', [0 Inf],1);
    for iSess=1:length(SessionPathList)
        [~, x.SESSIONS{end+1}]  = fileparts(SessionPathList{iSess});
    end

    if isempty(x.SESSIONS)
        fprintf('%s\n', 'No sessions found, defaulting to a single ASL_1 session');
        x.SESSIONS{1} = 'ASL_1'; % default session
    else
        x.SESSIONS = xASL_adm_SortStringNumbers(unique(x.SESSIONS));
    end
end

x.dataset.nSessions = length(x.SESSIONS);
x.dataset.nSubjectsSessions = x.dataset.nSubjects .* x.dataset.nSessions;

% ------------------------------------------------------------------------------------------------
%% 2. Manage TimePoint lists
[~, TimePoint] = xASL_init_LongitudinalRegistration(x);

% TimePointTotalSubjects sorts subjects in different time point cells, for reporting purposes
% Note: 
% *TotalSubjects == all subjects found by the BIDS-matlab or subjectRegexp
% *Subjects      == all subjects after removing those defined in x.dataset.exclusion

for iT=unique(TimePoint)'
    x.dataset.TimePointTotalSubjects{iT} = cell(0);
end

for iSubj=1:x.dataset.nTotalSubjects
    iSess=1;
    iSubjSess = (iSubj-1)*x.dataset.nSessions + iSess;
    x.dataset.TimePointTotalSubjects{TimePoint(iSubjSess)}{end+1} = x.dataset.TotalSubjects{iSubj};
end

x.dataset.nTimePointsTotal = length(x.dataset.TimePointTotalSubjects);
for iT=1:x.dataset.nTimePointsTotal
    x.dataset.nTimePointTotalSubjects(iT) = length(x.dataset.TimePointTotalSubjects{iT});
end


% ------------------------------------------------------------------------------------------------
%% 3. Create list of baseline & follow-up subjects (i.e. after exclusion)

% Now we here create the same cells list of subjects per time point, but then after exclusion
% Note: 
% *TotalSubjects == all subjects found by the BIDS-matlab or subjectRegexp
% *Subjects      == all subjects after removing those defined in x.dataset.exclusion
for iCell=1:length(x.dataset.TimePointTotalSubjects)
    x.dataset.TimePointSubjects{iCell} = cell(0);
end
    
for iS=1:x.dataset.nSubjects
    iSession = 1; % append to accommodate sessions in SetsID
    iSubjSess = ((iS-1)*x.dataset.nSessions)+iSession;
    CurrentTimePoint = TimePoint(iS);
    if length(x.dataset.TimePointSubjects)<CurrentTimePoint % if this cell didn't exist yet
        x.dataset.TimePointSubjects{CurrentTimePoint} = ''; 
    end    
    x.dataset.TimePointSubjects{CurrentTimePoint}{end+1} = x.SUBJECTS{iS};
end
    
x.dataset.nTimePoints = length(x.dataset.TimePointSubjects);
for iT=1:x.dataset.nTimePoints
    if  length(x.dataset.TimePointSubjects)<iT
        % if an excluded later volume led to different composition
        % of TotalSubjects (i.e. before exclusion) & Subjects (i.e. after
        % exclusion)
        x.dataset.nTimePointSubjects(iT) = 0; 
    else
        x.dataset.nTimePointSubjects(iT) = length(x.dataset.TimePointSubjects{iT});
    end
end

% ------------------------------------------------------------------------------------------------
%% 4. Check what excluded from which TimePoints
for iT=1:x.dataset.nTimePoints
    x.dataset.TimePointExcluded{iT} = cell(0);
end

if x.dataset.nExcluded>0
    for iExcluded=1:x.dataset.nExcluded
        FoundE = false;
        iT = 0;
        while iT<x.dataset.nTimePoints && ~FoundE % excluded subject not found in previous TimePoint
            iT = iT+1;  % go to next TimePoint
            
            iS = find(strcmp(x.dataset.TimePointTotalSubjects{iT}, x.dataset.ExcludedSubjects{iExcluded}));
            if ~isempty(iS)
                x.dataset.TimePointExcluded{iT}{end+1} = x.dataset.TimePointTotalSubjects{iT}{iS};
                FoundE = true;
            end
        end
    end

    if ~FoundE
        warning('Did not find excluded subjects in any of the time points');
    end
end

for iT=1:x.dataset.nTimePoints
    x.dataset.nTimePointExcluded(iT) = length(x.dataset.TimePointExcluded{iT});
end


end