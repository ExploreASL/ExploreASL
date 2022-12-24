function [nSessions, bSessionsMissing, SESSIONS] = xASL_adm_GetPopulationSessions(x)
% xASL_adm_GetPopulationSessions(x) obtain number of Sessions by determining amount of input files present in the Population folder
% FORMAT: [nSessions, bSessionsMissing] = xASL_adm_GetPopulationSessions(x)
%
% INPUT:
%   x                   - struct containing statistical pipeline environment parameters (REQUIRED)
%
% OUTPUT:
%   nSessions           - Maximum amount of sessions present in Population folder
%   bSessionsMissing    - Boolean to show if no sessions can be found in the Population folder
%   SESSIONS            - x.SESSIONS
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function looks for the maximum amount of sessions that
% are present in selected processed files present in the Population folder.
%
%   1. Determine which files to look for in the Population folder
%   2. Obtain list of session files in the Population folder
%   3. Determine unique amount of session numbers present in list
%   4. Set nSessions as highest unique session number 
%   5. Check and provide warning of number of sesssions differs per subject
%
% EXAMPLE: n/a
% __________________________________
% Copyright (C) 2015-2021 ExploreASL

%% 1. Administration
if ~isfield(x, 'SESSIONS')
    SESSIONS = 'ASL_1';
else
    SESSIONS = x.SESSIONS;
end

currentSubjectRegExp = x.dataset.subjectRegexp;
if strcmp(currentSubjectRegExp(1), '^')
    currentSubjectRegExp = currentSubjectRegExp(2:end);
end
if strcmp(currentSubjectRegExp(end), '$')
    currentSubjectRegExp = currentSubjectRegExp(1:end-1);
end
% 2. Look for processed files from which to determine the amount of sessions present
SessionList = xASL_adm_GetFileList(x.D.PopDir,['^' x.S.InputDataStr '_' currentSubjectRegExp '_ASL_\d+(\.nii|\.nii\.gz)$'], 'FPList', [0 Inf]);

%% 3. Obtain nSessions
if isempty(SessionList) % If no files found, search for subject files instead of session files

    nSessions = 1;
    bSessionsMissing = 1;
    fprintf('%s\n','No session or subject files found');
    return;    
else % If files found, continue with defining sessions from SessionList
    
    for iSession = 1:numel(SessionList)
        [IndexStart, IndexEnd] = regexpi(SessionList{iSession}, 'ASL_\d+');
        NewList{iSession,1} = SessionList{iSession}(1,IndexStart:IndexEnd); % create a cell array containings characters of all session numbers present
    end
    
    UniqueSessions = unique(NewList); % determine unique session numbers
    % 4. define nSessions as highest unique session number
    nSessions = numel(UniqueSessions);
    bSessionsMissing = 0;

    for iSession = 1:nSessions
        iSessionName = char(UniqueSessions(iSession));
        iSessionNumber = iSessionName(5:end); % Session number
        CountSessionNumbers(iSession)= sum(~cellfun('isempty',strfind(NewList,iSessionNumber))); % Counts amount of individual session numbers
    end
    
    % 5. Check and provide warning of number of sessions differs per subject
    CompareSessions = ones(1,size(CountSessionNumbers,2)) .* x.dataset.nSubjects; % create an array to check differences in sessions per subject with maximum amount of sessions
    if ~isequal(CountSessionNumbers,CompareSessions) % Check if amount of sessions is similar for each subject and provide warning if not
        warning('Amount of Sessions differs between Subjects');
    end
end

if numel(SESSIONS)<nSessions
    for iSession=1:nSessions
        SESSIONS{iSession} = ['ASL_' num2str(iSession)];
    end
end


end