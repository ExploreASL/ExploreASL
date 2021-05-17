function [x] = ExploreASL_ImportMaster(x)
%ExploreASL_ImportMaster Multi-step import workflow for DCM2NII, NII2BIDS & BIDS2LEGACY.
%
% FORMAT: [x] = ExploreASL_ImportMaster(x)
%
% INPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Multi-step import workflow for DCM2NII, NII2BIDS & BIDS2LEGACY.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Import Workflow
    
    % We expect x.opts.DataParPath to be the sourceStructure.json file, so the root directory should be the study directory
	[x.dir.StudyRoot,~,~] = xASL_fileparts(x.opts.DataParPath);
	
    % Check if at least one of the three steps should be performed
    missingFields = false; % Fallback
    if sum(x.opts.ImportModules)>0
        % DICOM TO NII
        if x.opts.ImportModules(1)==1
            if ~isempty(x.dir.sourceStructure)
                xASL_module_Import(x.dir.StudyRoot, x.opts.DataParPath, [], [1 0 0], false, true, false, false, x);
            else
                missingFields = true;
            end
        end
        % NII TO BIDS
        if x.opts.ImportModules(2)==1
            if ~isempty(x.dir.sourceStructure)
                [x] = xASL_module_Import(x.dir.StudyRoot, x.opts.DataParPath, [], [0 1 0], false, true, false, false, x);
            else
                missingFields = true;
            end
        end
        % ANONYMIZE
        if x.opts.ImportModules(3)==1
            if ~isempty(x.dir.StudyRoot)
                [x] = xASL_module_Import(x.dir.StudyRoot, x.opts.DataParPath, [], [0 0 1], false, true, false, false, x);
            else
                missingFields = true;
            end
        end
        % BIDS TO LEGACY
        if x.opts.ImportModules(4)==1
            if ~isempty(x.dir.dataset_description)
                x = xASL_imp_BIDS2Legacy(x);
            else
                missingFields = true;
            end
        end
    end

    if missingFields
        warning('Import workflow is turned on, but at least one required JSON file is missing...');
    end
    
    % Reset the import parameter (for the second initialization including the loading of the dataset)
    x.opts.bImportData = 0;
    x.opts.ImportModules = [0 0 0 0];
    
end


