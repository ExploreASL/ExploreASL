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