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
    % 2) Trim whitespace
    % 3) Selectively convert forward to backward slashes (ignore already escaped whitespace)
    % 4) Escape characters and residual whitespaces (ignore already escaped whitespaces)
    % 5) If WSL: add mounting prefix
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % EXAMPLE xASL_adm_UnixPath('   \Users/User/Google Drive\My      Photos\Name(With)Brackets)  ');
    % This should output '/Users/User/Google\ Drive/My\ \ \ \ \ \ Photos/Name\(With\)Brackets\)'
    % __________________________________
    % Copyright 2020 ExploreASL
        
    %% ===================================================================================
    %% 1) Skip this function without Unix-filesystem
    % If we don't have a Unix-compatible filesystem, then skip this function
    [statusWSL, ~] = system('wsl ls');
    if ~isunix && ~statusWSL==0
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
    PathIs = regexprep(PathIs, '(?<!\\)( |\(|\)|[|]|{|}|*|:|;|+|=|,|<|>|?|~|@|#|%|^|&|*)', '\\$1');
    
    %% ===================================================================================
    %% 5) If WSL: add mounting prefix
    if ispc && statusWSL==0 && strcmp(PathIs(2), ':') % if we have Windows Subsystem for Linux
        if statusWSL==0 && strcmp(PathIs(2), ':')
            PathIs = ['/mnt/' lower(PathIs(1)) '/' PathIs(4:end)];
        end
    end
end