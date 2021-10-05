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
%
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Import Workflow

    % Run the initialization
    try
        x = ExploreASL_ImportMaster_Init(x);
    catch loggingEntry
        % Print user feedback if import crashed
        fprintf(2,'ExploreASL Import module failed...\n');
        % Check loggingEntry
        if size(loggingEntry.stack,1)>0
            fprintf(2,'%s\n%s, line %d...\n',loggingEntry.message,loggingEntry.stack(1).name,loggingEntry.stack(1).line);
        end
        [x] = xASL_qc_AddLoggingInfo(x, loggingEntry);
        % Turn off data loading and processing if import crashed
        x.opts.bLoadData = false;
        x.opts.bProcessData = false;
    end
    
    % Run the import submodules
    try
        [~, x] = xASL_Iteration(x,'xASL_module_Import');
    catch loggingEntry
        % Print user feedback if import crashed
        fprintf(2,'ExploreASL Import module failed...\n');
        % Check loggingEntry
        if size(loggingEntry.stack,1)>0
            fprintf(2,'%s\n%s, line %d...\n',loggingEntry.message,loggingEntry.stack(1).name,loggingEntry.stack(1).line);
        end
        [x] = xASL_qc_AddLoggingInfo(x, loggingEntry);
        % Turn off data loading and processing if import crashed
        x.opts.bLoadData = false;
        x.opts.bProcessData = false;
    end
    
    % Reset the import parameters
    x.opts.bImportData = 0;
    x.opts.ImportModules = [0 0 0 0];

    
end



%% Init function for ExploreASL_ImportMaster
function x = ExploreASL_ImportMaster_Init(x)


    % Initialize x struct
    x = xASL_init_SubStructs(x);
    
    % We expect x.opts.DatasetRoot to be the study root directory, but if it is not defined,
    % then the user probably used a path to a descriptive JSON file instead
    if isfield(x, 'dir') && isfield(x.dir, 'DatasetRoot') && isempty(x.dir.DatasetRoot)
        x.dir.DatasetRoot = xASL_fileparts(x.opts.DatasetRoot);
    end
    
    % We are running the import here, but if BIDS to Legacy will not run,
    % we can not load the data afterwards!
    if ~x.opts.ImportModules(4)
        x.opts.bLoadData = false;
    end
    
    % Currently fixed import settings
    x.modules.import.settings.bCopySingleDicoms = false;
    x.modules.import.settings.bUseDCMTK = true;
    x.modules.import.settings.bCheckPermissions = false;
    
    % Basic parameter checks
    x = xASL_imp_BasicParameterChecks(x);

    % Initialize the import setup
    if x.opts.ImportModules(1) || x.opts.ImportModules(2) || x.opts.ImportModules(3)
        % Load the sourceStructure.json and initialize the corresponding struct
        x.modules.import.imPar = xASL_imp_Initialize(x.dir.DatasetRoot, x.dir.sourceStructure);
    else
        x.modules.import.imPar = NaN;
    end

    % Determine subject/session/run structure from sourcedata or temp data
    x = xASL_imp_DetermineSubjectStructure(x);


end




