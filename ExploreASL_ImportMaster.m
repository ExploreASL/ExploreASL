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
    
    if x.opts.bImportData
        % We expect x.opts.DatasetRoot to be the study root directory, but if it is not defined, 
        % then the user probably used a path to a descriptive JSON file instead
        if isfield(x, 'dir') && isfield(x.dir, 'DatasetRoot') && isempty(x.dir.DatasetRoot)
            x.dir.DatasetRoot = xASL_fileparts(x.opts.DatasetRoot);
        end

        % Check if at least one of the three steps should be performed
        missingJSON = false; % Fallback
        if sum(x.opts.ImportModules)>0
            % DICOM TO NII
            if x.opts.ImportModules(1)==1
                if ~isempty(x.dir.sourceStructure)
                    try
                        xASL_module_Import(x.dir.DatasetRoot, x.dir.sourceStructure, x.dir.studyPar, [1 0 0], false, true, false, x);
                    catch loggingEntry
                        fprintf(2,'DICOM to NIfTI module failed...\n');
                        [x] = xASL_qc_AddLoggingInfo(x, loggingEntry);
                    end
                else
                    missingJSON = true;
                end
            end
            % NII TO BIDS
            if x.opts.ImportModules(2)==1
                if ~isempty(x.dir.sourceStructure)
                    try
                        [x] = xASL_module_Import(x.dir.DatasetRoot, x.dir.sourceStructure, [], [0 1 0], false, true, false, x);
                    catch loggingEntry
                        fprintf(2,'NIfTI to BIDS module failed:\n%s\n',loggingEntry.message);
                        [x] = xASL_qc_AddLoggingInfo(x, loggingEntry);
                    end
                else
                    missingJSON = true;
                end
            end
            % DEFACE
            if x.opts.ImportModules(3)==1
                if ~isempty(x.dir.DatasetRoot)
                    try
                        [x] = xASL_module_Import(x.dir.DatasetRoot, x.dir.sourceStructure, [], [0 0 1], false, true, false, x);
                    catch loggingEntry
                        fprintf(2,'Deface module failed...\n');
                        [x] = xASL_qc_AddLoggingInfo(x, loggingEntry);
                    end
                else
                    missingJSON = true;
                end
            end
            % BIDS TO LEGACY
            if x.opts.ImportModules(4)==1
                if ~isempty(x.dir.dataset_description)
                    try
                        x = xASL_imp_BIDS2Legacy(x);
                    catch loggingEntry
                        fprintf(2,'BIDS to legacy module failed...\n');
                        [x] = xASL_qc_AddLoggingInfo(x, loggingEntry);
                    end
                else
                    missingJSON = true;
                end
            end
        end

        % Print warning about missing JSON files
        if missingJSON
            warning('Import workflow is turned on, but at least one required JSON file is missing...');
        end

        % Reset the import parameters
        x.opts.bImportData = 0;
        x.opts.ImportModules = [0 0 0 0];
    end
    
end

