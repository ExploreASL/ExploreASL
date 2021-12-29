function [differences,identical,dn] = xASL_bids_CompareStructuresNIFTI(differences,identical,bPrintReport,detailedOutput,allFiles,iFile,dn,currentFileA,currentFileB,threshRmseNii)
%xASL_bids_CompareStructuresNIFTI Compare NIFTI files
%
% FORMAT: [differences,identical,dn] = xASL_bids_CompareStructuresNIFTI(differences,identical,bPrintReport,detailedOutput,allFiles,iFile,dn,currentFileA,currentFileB,threshRmseNii)
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
%         threshRmseNii  - Threshold for NIfTI comparison (DOUBLE, REQUIRED)
%
% OUTPUT:
%         differences    - Differences between datasets (CELL ARRAY, REQUIRED)
%         identical      - Datasets are identical (BOOLEAN, REQUIRED)
%         dn             - Difference number for the row within differences (INTEGER, REQUIRED)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      This script compares the content of two NIFTI files for
%                   the BIDS flavor testing.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          ...
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2021 ExploreASL

    % Read files if they exist
    if (exist(currentFileA,'file') && exist(currentFileB,'file'))
        [~,RMSE,~,~,dimCheck] = xASL_im_CompareNiftis(currentFileA,currentFileB,detailedOutput);
        
        % Check RMSE
        if (RMSE>threshRmseNii)
            if bPrintReport
                fprintf('File:      %s\n',allFiles{iFile});
                fprintf('           RMSE of NIFTIs above threshold.\n');
            end
            identical = false;
            differences{dn,1} = ['RMSE of NIFTIs above threshold: ', allFiles{iFile}, ' '];
            dn = dn+1;
        end
        
        % Check image dimensions
        if ~dimCheck
            if bPrintReport
                fprintf('File:      %s\n',allFiles{iFile});
                fprintf('           Matrix dimensions do not agree.\n');
            end
            identical = false;
            
            % Save difference
            differences{dn,1} = ['Matrix dimensions do not agree: ', allFiles{iFile}, ' '];
            dn = dn+1;
        end
                
    end

end


