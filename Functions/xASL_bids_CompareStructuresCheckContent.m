function [identical,differences] = xASL_bids_CompareStructuresCheckContent(filesDatasetA,filesDatasetB,pathDatasetA,pathDatasetB,identical,bPrintReport,detailedOutput,threshRmseNii)
 %xASL_bids_CompareStructuresCheckContent Check file contents
%
% FORMAT: [identical,differences] = xASL_bids_CompareStructuresCheckContent(filesDatasetA,filesDatasetB,pathDatasetA,pathDatasetB,identical,bPrintReport,detailedOutput,threshRmseNii)
%
% INPUT:
%         filesDatasetA  - List of files in dataset A (CELL ARRAY, REQUIRED)
%         filesDatasetB  - List of files in dataset B (CELL ARRAY, REQUIRED)
%         pathDatasetA   - Path to dataset A (STRING, REQUIRED)
%         pathDatasetB   - Path to dataset B (STRING, REQUIRED)
%         identical      - Datasets are identical (BOOLEAN, REQUIRED)
%         bPrintReport   - Print report to console (BOOLEAN, REQUIRED)
%         detailedOutput - Print an even more detailed report (BOOLEAN, REQUIRED)
%         threshRmseNii  - Threshold for NIfTI comparison (DOUBLE, REQUIRED)
%
% OUTPUT:
%         identical   - Datasets are identical (BOOLEAN)
%         differences - Differences between datasets (CELL ARRAY)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      This script iterates over the provided files (filesDatasetA,filesDatasetB).
%                   There are different comparisons implemented for JSON, TSV, TEXT, & NIFTI files.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:
%
%  filesDatasetA = {...};
%  filesDatasetB = {...};
%  pathDatasetA = '...';
%  pathDatasetB = '...';
%
%  [identical,differences] = ...
%                   xASL_bids_CompareStructuresCheckContent(...
%                   filesDatasetA,filesDatasetB,...
%                   pathDatasetA,pathDatasetB,...
%                   1,0,0,1e-3);
%
% REFERENCES:       ...
% __________________________________
% Copyright (c) 2015-2021 ExploreASL


    % All files
    allFiles = unique([filesDatasetA',filesDatasetB']');
    
    % Differences
    differences = cell(1,1);
    
    % Difference number
    dn = 1;
    
    % Iterate over list
    for iFile=1:length(allFiles)
        % Assign root directory of dataset A
        currentFileA = fullfile(pathDatasetA,allFiles{iFile});
        currentFileB = fullfile(pathDatasetB,allFiles{iFile});
        % Get extension
        [~,~,extension] = fileparts(allFiles{iFile});
        % Check extension
        if strcmp(extension,'.json')
			[differences,identical,dn] = ...
                xASL_bids_CompareStructuresJSON(differences,identical,bPrintReport,allFiles,iFile,dn,currentFileA,currentFileB);
        elseif strcmp(extension,'.tsv')
            [differences,identical,dn] = ...
                xASL_bids_CompareStructuresTSV(differences,identical,bPrintReport,allFiles,iFile,dn,currentFileA,currentFileB);
        elseif strcmp(extension,'.txt') || strcmp(extension,'.csv')
            [differences,identical,dn] = ...
                xASL_bids_CompareStructuresTEXT(differences,identical,bPrintReport,allFiles,iFile,dn,currentFileA,currentFileB);
        elseif strcmp(extension,'.nii') || ~isempty(strfind(allFiles{iFile},'.nii.gz'))
            [differences,identical,dn] = ...
                xASL_bids_CompareStructuresNIFTI(differences,identical,bPrintReport,detailedOutput,allFiles,iFile,dn,currentFileA,currentFileB,threshRmseNii);
        end
    end
    
    
    
end


