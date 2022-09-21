function [PathIs] = xASL_adm_UnixPath(PathIs, bTryWSL)
%xASL_adm_UnixPath Convert paths to Unix-compatible paths
%
% FORMAT: [PathIs] = xASL_adm_UnixPath(PathIs[, bTryWSL])
%
% INPUT:
%   PathIs    - string containing a single path to correct (REQUIRED)
%   bTryWSL   - boolean, test for the presence of Windows Subsystem for Linux (WSL) functionality
%               and if present, provide a combined path to the mounted drive (OPTIONAL, DEFAULT = false)
%
% OUTPUT:
%   PathIs - corrected path (DEFAULT = uncorrected, same as input path).
%            Path is only changed when a Unix-filesystem is detected.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function performs the following steps to convert a path to a path that is compatible with the Unix-filesystem
%              as used in e.g. Linux/MacOS. It also has special support for Windows Subsystem for Linux (WSL),
%              though this should only be activated specifically for WSL calls.
%
%              Note that we want to use this function for most Unix calls
%              (to fix paths), but in the case of WSL only for some calls,
%              where Matlab in Windows calls Linux-code through WSL (e.g.
%              for FSL) - these have to be explicitly specified by the bTryWSL option.
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
% Copyright 2015-2022 ExploreASL
        
    %% ===================================================================================
	%% Input parameter administration
	
	% First initializing the optional input parameter
	if nargin < 2 || isempty(bTryWSL)
		bTryWSL = false;
	end
	
    bWSL = false;
    if bTryWSL && ispc
        % only in cases where bTryWSL is specifically called
        % and we have a PC system with WSL installed
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
