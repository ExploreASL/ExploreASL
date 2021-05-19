function [x] = xASL_imp_Import_UpdateStudyRoot(x, studyPath)
%xASL_imp_Import_UpdateStudyRoot Update x.opts.StudyRoot to dataset_description.json after NII2BIDS conversion
%
% FORMAT: [x] = xASL_imp_Import_UpdateStudyRoot(x, studyPath)
%
% INPUT:
%   x                   - ExploreASL x structure (STRUCT, REQUIRED)
%   studyPath           - Path to study root directory (CHAR ARRAY, PATH, REQUIRED)
%
% OUTPUT: 
%   x                   - ExploreASL x structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Update x.opts.StudyRoot to dataset_description.json after NII2BIDS conversion
% EXAMPLE:     n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Search for dataset_description.json within the rawdata subfolder
    foundFiles = xASL_adm_GetFileList(fullfile(studyPath,'rawdata'),'dataset_description.json');
    
    % Check if valid dataset_description.json exists within the rawdata folder
    if isempty(foundFiles)
        warning('No valid dataset_description.json found within the rawdata directory...');
    else
        x.dir.dataset_description = foundFiles{1};
    end

end


