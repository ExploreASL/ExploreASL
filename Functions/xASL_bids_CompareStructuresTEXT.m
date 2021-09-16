function [differences,identical,dn] = xASL_bids_CompareStructuresTEXT(differences,identical,bPrintReport,allFiles,iFile,dn,currentFileA,currentFileB)
%xASL_bids_CompareStructuresTEXT Compare TEXT files
%
% FORMAT: [differences,identical,dn] = xASL_bids_CompareStructuresTEXT(differences,identical,bPrintReport,allFiles,iFile,dn,currentFileA,currentFileB)
%
% INPUT:
%         differences    - Differences between datasets (CELL ARRAY, REQUIRED)
%         identical      - Datasets are identical (BOOLEAN, REQUIRED)
%         bPrintReport   - Print report to console (BOOLEAN, REQUIRED)
%         allFiles       - Array containing the file names (CELL ARRAY, REQUIRED)
%         iFile          - File number (INTEGER, REQUIRED)
%         dn             - Difference number for the row within differences (INTEGER, REQUIRED)
%         currentFileA   - Name of the current file in dataset A (STRING, REQUIRED) 
%         currentFileB   - Name of the current file in dataset B (STRING, REQUIRED) 
%
% OUTPUT:
%         differences    - Differences between datasets (CELL ARRAY, REQUIRED)
%         identical      - Datasets are identical (BOOLEAN, REQUIRED)
%         dn             - Difference number for the row within differences (INTEGER, REQUIRED)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      This script compares the content of two TEXT files for
%                   the BIDS flavor testing.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          ...
%
% REFERENCES:       ...
% __________________________________
% Copyright (c) 2015-2021 ExploreASL


    % Read files if they exist
    if (exist(currentFileA,'file') && exist(currentFileB,'file')) % xASL_exist somehow didn't work here (again)
        % Compare text files content directly
        currentFileTextA = fileread(currentFileA);
        currentFileTextB = fileread(currentFileB);
        if ~strcmp(currentFileTextA,currentFileTextB)
            [identical,differences,dn] = xASL_bids_CompareStoreDifference(bPrintReport,differences,dn,allFiles,iFile);
        end
    end

end





