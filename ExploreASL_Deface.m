function [x] = ExploreASL_Deface(x)
%ExploreASL_Deface Run the defacing.
%
% FORMAT: [x] = ExploreASL_Deface(x)
%
% INPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x      - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Deface master script.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
%
% __________________________________
% Copyright (c) 2015-2022 ExploreASL


    %% Deface Initialization
    try
        % Before we run the subject-wise xASL_module_Deface we need to initialize some x structure fields.
        % We determine the general subject/visit/session structure and store everything required for import/defacing/processing in x.
        x = xASL_init_Import(x);
    catch loggingEntry
        % Print user feedback if deface crashed
        fprintf(2,'ExploreASL Deface initialization failed...\n');
        % Check loggingEntry
        if size(loggingEntry.stack,1)>0
            fprintf(2,'%s\n%s, line %d...\n',loggingEntry.message,loggingEntry.stack(1).name,loggingEntry.stack(1).line);
        end
        [x] = xASL_qc_AddLoggingInfo(x, loggingEntry);
        % Turn off data loading and processing if import crashed
        x.opts.bLoadData = false;
        x.opts.bProcessData = false;
    end
    
    
    %% Deface workflow
    try
        % Here we run the subject-wise ExploreASL xASL_module_Deface.
        [~, x] = xASL_init_Iteration(x,'xASL_module_Deface');
    catch loggingEntry
        % Print user feedback if deface crashed
        fprintf(2,'ExploreASL Deface module failed...\n');
        % Check loggingEntry
        if size(loggingEntry.stack,1)>0
            fprintf(2,'%s\n%s, line %d...\n',loggingEntry.message,loggingEntry.stack(1).name,loggingEntry.stack(1).line);
        end
        [x] = xASL_qc_AddLoggingInfo(x, loggingEntry);
        % Turn off data loading and processing if import crashed
        x.opts.bLoadData = false;
        x.opts.bProcessData = false;
    end
    
    % Reset the deface parameters
    x.opts.bDefaceData = 0;

    
end




