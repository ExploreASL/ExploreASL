function [differences,identical,dn] = xASL_bids_CompareStructuresTSV(differences,identical,bPrintReport,allFiles,iFile,dn,currentFileA,currentFileB)
%xASL_bids_CompareStructuresTSV Compare TSV files
%
% FORMAT: [differences,identical,dn] = xASL_bids_CompareStructuresTSV(differences,identical,bPrintReport,allFiles,iFile,dn,currentFileA,currentFileB)
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
% DESCRIPTION:      This script compares the content of two TSV files for
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
    if xASL_exist(currentFileA,'file') && xASL_exist(currentFileB,'file')
        % Compare text files content directly
        currentTsvA = xASL_tsvRead(currentFileA);
        currentTsvB = xASL_tsvRead(currentFileB);
        if iscellstr(currentTsvA) && iscellstr(currentTsvB)
            if ~isempty(setdiff(currentTsvA,currentTsvB))
                [identical,differences,dn] = xASL_bids_CompareStoreDifference(bPrintReport,differences,dn,allFiles,iFile);
            end
        else
            % These cells contain both strings and other datatypes, we have
            % to do an iterative comparison if they have the same size
            if isequal(size(currentTsvA),size(currentTsvB))
                for iRow=1:size(currentTsvA,1)
                    for iColumn=1:size(currentTsvA,2)
                        if ~isequal(currentTsvA{iRow,iColumn},currentTsvB{iRow,iColumn})
                            [identical,differences,dn] = xASL_bids_CompareStoreDifference(bPrintReport,differences,dn,allFiles,iFile);
                        end
                    end
                end
            else
                [identical,differences,dn] = xASL_bids_CompareStoreDifference(bPrintReport,differences,dn,allFiles,iFile);
            end
        end
    end

end

