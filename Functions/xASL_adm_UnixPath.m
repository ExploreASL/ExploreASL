function [PathIs] = xASL_adm_UnixPath(PathIs)
%xASL_adm_UnixPath Convert paths to Unix-compatible paths
%
% FORMAT: [PathIs] = xASL_adm_UnixPath(PathIs)
% 
% INPUT:
%   PathIs - string containing a single path to correct (REQUIRED)
%
% OUTPUT:
%   PathIs - corrected path (DEFAULT = uncorrected, same as input path.
%            Path is only changed when a Unix-filesystem is detected.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function performs the following steps to convert a path to a path that is compatible with the Unix-filesystem
%              as used in e.g. Linux/MacOS/Windows Subsystem for Linux (WSL).
% 1) Skip this function without Unix-filesystem
% 2) Remove leading spaces
% 3) Remove trailing spaces
% 4) Convert forward to backward slashes
% 5) Escape residual spaces
% 6) If WSL: add mounting prefix
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE xASL_adm_UnixPath('   \Users/User/Google Drive\ ');
% This should output '/Users/User/Google\ Drive/'
% __________________________________
% Copyright 2020 ExploreASL


%% ===================================================================================
%% 1) Skip this function without Unix-filesystem
% If we don't have a Unix-compatible filesystem, then skip this function
[statusWSL, result] = system('wsl ls');
if ~isunix && ~statusWSL==0
    return;
    % dont change path if we dont have unix installed, either as main OS or as
    % subsystem for Windows
end


%% ===================================================================================
%% 2) Remove leading spaces
NextN = 1;
SkipThis = false;
OrigLengthPathOut = length(PathIs);
while NextN<OrigLengthPathOut && SkipThis==false
    if strcmp(PathIs(1),' ')
        PathIs = PathIs(2:end);
        NextN = NextN+1;
    else
        SkipThis = true;
    end
end


%% ===================================================================================
%% 3) Remove trailing spaces
NextN = 1;
SkipThis = false;
OrigLengthPathOut = length(PathIs);
while NextN<OrigLengthPathOut && SkipThis==false
    if strcmp(PathIs(end),' ')
        PathIs = PathIs(1:end-1);
        NextN = NextN+1;
    else
        SkipThis = true;
    end
end


%% ===================================================================================
%% 4) Convert forward to backward slashes
PathIs = strrep(PathIs, '\','/'); % convert slashes
PathIs = strrep(PathIs, '/ ','\ '); % except for escaping spaces



%% ===================================================================================
%% 5) Escape residual spaces
StrInd = regexp(PathIs, '(?<!\\) '); % string indices with space that is not yet escaped
for iIndex=1:length(StrInd)
    PathIs = [PathIs(1:StrInd(iIndex)-1) '\ ' PathIs(StrInd(iIndex)+1:end)];
    % This added an escaping '\' so other indices should move 1 to the
    % right:
    StrInd(iIndex+1:end) = StrInd(iIndex+1:end)+1;
end


%% ===================================================================================
%% 6) If WSL: add mounting prefix
if ispc && statusWSL==0 && strcmp(PathIs(2), ':') % if we have Windows Subsystem for Linux
    if statusWSL==0 && strcmp(PathIs(2), ':')
        PathIs = ['/mnt/' lower(PathIs(1)) '/' PathIs(4:end)];
    end
end


end