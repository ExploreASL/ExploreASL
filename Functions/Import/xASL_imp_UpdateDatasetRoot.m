function x = xASL_imp_UpdateDatasetRoot(x)
%xASL_imp_UpdateDatasetRoot Update x.opts.DatasetRoot to dataset_description.json after NII2BIDS conversion
%
% FORMAT: x = xASL_imp_UpdateDatasetRoot(x)
%
% INPUT:
%   x     - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT: 
%   x     - ExploreASL x structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Update x.opts.DatasetRoot to dataset_description.json after NII2BIDS conversion
%
% EXAMPLE:     n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Search for dataset_description.json within the rawdata subfolder
    foundFiles = xASL_adm_GetFileList(fullfile(x.dir.RawData),'dataset_description.json');
    
    % Check if valid dataset_description.json exists within the rawdata folder
    if isempty(foundFiles)
        warning('No valid dataset_description.json found within the rawdata directory...');
        fprintf('%s\n', 'Ensure that your data is BIDS compatible by running the BIDS validator');
    else
        x.dir.dataset_description = foundFiles{1};
    end

end


