function [nSessions, bSessionsMissing] = xASL_adm_GetnSessions(x)
%xASL_adm_GetnSessions(x) obtain number of Sessions by determining amount of input files present in the Population folder
% FORMAT: [nSessions, bSessionsMissing] = xASL_adm_GetnSessions(x)
%
% INPUT:
%   x                   - struct containing statistical pipeline environment parameters (REQUIRED)
%
% OUTPUT:
%   nSessions           - Maximum amount of sessions present in Population folder
%   bSessionsMissing    - Variable to show if no sessions can be found in the Population folder
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function looks for the maximum amount of sessions that
% are present in the dataset.

% Determine which files to look in the Population folder
currentSubjectRegExp = x.subject_regexp;
if strcmp(currentSubjectRegExp(1), '^')
    currentSubjectRegExp = currentSubjectRegExp(2:end);
end
if strcmp(currentSubjectRegExp(end), '$')
    currentSubjectRegExp = currentSubjectRegExp(1:end-1);
end
% Look for processed files from which to determine the amount of sessions present
SessionList = xASL_adm_GetFileList(x.D.PopDir,['^' x.S.InputDataStr '_' currentSubjectRegExp '_ASL_\d*\.nii$'], 'FPList', [0 Inf]);

if isempty(SessionList) % If no files found, search for subjects instead of sessions
    if isempty(xASL_adm_GetFileList(x.D.PopDir,'^qCBF.*\.nii', 'FPList', [0 Inf]))
        fprintf('%s\n','No session or subject files found');
        return;
    end
    nSessions = 1;
    bSessionsMissing = 1;
else % If files found, continu with defining sessions from SessionList
    [IndexStart , ~] = cellfun(@(y) regexp(y,'ASL_\d+\.(nii|nii\.gz)'), SessionList);
    
    for n = 1:numel(SessionList)
        StrLoc = IndexStart(n);
        SessionsPopFolder = SessionList{n};
        NewList{n,1} = SessionsPopFolder(1,StrLoc:end); % create a cell containings strings of all session numbers present
    end
    
    UniqueSessions = unique(cell2mat(NewList),'rows'); % determine unique session numbers
    nSessions = size(UniqueSessions,1); % define nSessions as maximum amount of unique session numbers
    bSessionsMissing = 0;
    NewList = cellstr(NewList);
    for o = 1:nSessions
        CountSessionNumbers(o)= sum(~cellfun('isempty',strfind(NewList,num2str(o))));
    end
    
    CompareSessions = ones(1,size(CountSessionNumbers,2)) .* nSessions; % create an array to check differences in sessions per subject with maximum amount of sessions
    
    if ~isequal(CountSessionNumbers,CompareSessions) % Check if amount of sessions is similar for each subject and provide warning if not
        warning('Amount of Sessions per Subject different')
    end
end

end