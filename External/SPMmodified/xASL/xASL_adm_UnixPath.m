function [PathIs] = xASL_adm_UnixPath(PathIs, bWSL)
%xASL_adm_UnixPath Convert paths to Unix-compatible paths
%
% FORMAT: [PathIs] = xASL_adm_UnixPath(PathIs)
%
% INPUT:
%   PathIs - string containing a single path to correct (REQUIRED)
%   bWSL   - boolean, true for Windows Subsystem for Linux (WSL) functionality
%            (OPTIONAL, DEFAULT = false)
%
% OUTPUT:
%   PathIs - corrected path (DEFAULT = uncorrected, same as input path.
%            Path is only changed when a Unix-filesystem is detected.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function performs the following steps to convert a path to a path that is compatible with the Unix-filesystem
%              as used in e.g. Linux/MacOS.
%              It also has special support for Windows Subsystem for Linux (WSL),
%              though this should only be activated specifically for WSL calls.
%
% 1. Skip this function without Unix-filesystem
% 2. Trim whitespace
% 3. Selectively convert forward to backward slashes (ignore already escaped whitespace)
% 4. Escape characters and residual whitespaces (ignore already escaped whitespaces)
% 5. If WSL: add mounting prefix
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_adm_UnixPath('   \Users/User/Google Drive\My      Photos\Name(With)Brackets)  ');
% This should output '/Users/User/Google\ Drive/My\ \ \ \ \ \ Photos/Name\(With\)Brackets\)'
% __________________________________
% Copyright 2021 ExploreASL
        
    if nargin<2 || isempty(bWSL)
        bWSL = false;
    elseif ispc
        [statusWSL, ~] = system('wsl ls');
        if statusWSL==0
            bWSL = true;
        end
    end

    %% ===================================================================================
    %% 1) Skip this function without Unix-filesystem
    % If we don't have a Unix-compatible filesystem, then skip this function
    
    if ~isunix && ~bWSL
        return;
        % dont change path if we dont have unix installed, either as main OS or as
        % subsystem for Windows
    end
    
    %% ===================================================================================
    %% 2) Trim leading or trailing whitespace
    PathIs = strtrim(PathIs);
    
    %% ===================================================================================
    %% 3) Replace all \ which are not followed by a space with a /
    PathIs =regexprep(PathIs, '(\\)(?! )', '/');
    
    %% ===================================================================================
    %% 4) Escape characters in file name
    if ispc
        IllegalCharacters = '(?<!\\)( |\(|\)|[|]|{|}|*|;|+|=|,|<|>|!|~|@|#|%|^|&|*)';
    else
        IllegalCharacters = '(?<!\\)( |\(|\)|[|]|{|}|*|;|+|=|,|<|>|!|~|@|#|%|^|&|*|:)';
    end    
    
    PathIs = regexprep(PathIs, IllegalCharacters, '\\$1');
    
    %% ===================================================================================
    %% 5) If WSL: add mounting prefix
    if bWSL % if we have Windows Subsystem for Linux
        if strcmp(PathIs(2), ':')
            PathIs = ['/mnt/' lower(PathIs(1)) '/' PathIs(4:end)];
        end
    end
    
end