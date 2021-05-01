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
    
    % We expect x.DataParPath to be the sourceStructure.json file, so the root directory should be the study directory
    [x.StudyRoot,~,~] = xASL_fileparts(x.DataParPath);
    
    % Check if at least one of the three steps should be performed
    if sum(x.ImportModules)>0
        if ~isempty(x.DataParPath)
            if exist(x.DataParPath,'file')
                % DICOM TO NII
                if x.ImportModules(1)==1
                    xASL_module_Import(x.StudyRoot, x.DataParPath, [], [1 0 0], false, true, false, false, x);
                end
                % NII TO BIDS
                if x.ImportModules(2)==1
                    xASL_module_Import(x.StudyRoot, x.DataParPath, [], [0 1 0], false, true, false, false, x);
                end
                % ANONYMIZE
                if x.ImportModules(3)==1
                    xASL_module_Import(x.StudyRoot, x.DataParPath, [], [0 0 1], false, true, false, false, x);
                end
                % BIDS TO LEGACY
                if x.ImportModules(4)==1
                    x = xASL_Import_BIDS2LEGACY(x);
                end
            else
                warning('ImportModules was set to 1, but the sourceStructure file does not exist');
            end
        else
            warning('ImportModules was set to 1, but there is no DataParPath, Import will not be executed');
        end
    end
    
    % Reset the import parameter (for the second initialization including the loading of the dataset)
    x.bImportData = 0;
    x.ImportModules = [0 0 0 0];
    
end


