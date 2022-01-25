function [x] = ExploreASL_Import(x)
%ExploreASL_Import Multi-step import workflow for DCM2NII & NII2BIDS.
%
% FORMAT: [x] = ExploreASL_Import(x)
%
% INPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Multi-step import workflow for DCM2NII & NII2BIDS.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
%
% __________________________________
% Copyright (c) 2015-2022 ExploreASL


    %% Import Initialization
    try
        % Before we run the subject-wise xASL_module_Import we need to initialize some x structure fields. For DCM2NII & NII2BIDS 
        % we try to read the sourceStructure.json and studyPar.json as well as the general file structure. We determine the general
        % subject/visit/session structure and store everything required for import/processing in x.
        x = xASL_imp_ImportInitialization(x);
    catch loggingEntry
        % Print user feedback if import crashed
        fprintf(2,'ExploreASL Import initialization failed...\n');
        % Check loggingEntry
        if size(loggingEntry.stack,1)>0
            fprintf(2,'%s\n%s, line %d...\n',loggingEntry.message,loggingEntry.stack(1).name,loggingEntry.stack(1).line);
        end
        [x] = xASL_qc_AddLoggingInfo(x, loggingEntry);
        % Turn off data loading and processing if import crashed
        x.opts.bLoadData = false;
        x.opts.bProcessData = false;
    end
    
    
    %% Import workflow
    try
        % Here we run the subject-wise ExploreASL xASL_module_Import
        [~, x] = xASL_init_Iteration(x,'xASL_module_Import');
        % Turn off the import
        x.opts.bImportData = false;
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

    
end




