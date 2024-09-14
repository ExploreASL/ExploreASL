function xASL_bids_BIDS2xASL_CopyFile(pathOrig, pathDest, bOverwrite)
%xASL_bids_BIDS2xASL_CopyFile Copy files for BIDS to Legacy conversion.
%
% FORMAT:     xASL_bids_BIDS2xASL_CopyFile(pathOrig, pathDest, bOverwrite)
% 
% INPUT:      pathOrig - Origin path (STRING, REQUIRED)
%             pathDest - Destination path (STRING, REQUIRED)
%             bOverwrite  - Overwrite (BOOLEAN, REQUIRED)
%   
% OUTPUT:     n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Copy files for BIDS to Legacy conversion.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_bids_BIDS2xASL_CopyFile(pathOrig, pathDest, bOverwrite)
% __________________________________
% Copyright 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    % Create folder(s) if didnt exist
    xASL_adm_CreateDir(fileparts(pathDest));
    % check for existance & overwriting
    if ~xASL_exist(pathOrig)
        warning(['Couldnt find ' pathOrig, ' skipping']);
    elseif xASL_exist(pathDest) && ~bOverwrite
        warning([pathDest ' already existed, skipping']);
    else % if didnt already exist, or bOverwrite is true
        xASL_Copy(pathOrig, pathDest, 1);
    end
    
end