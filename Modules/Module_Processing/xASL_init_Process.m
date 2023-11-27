function x = xASL_init_Process(x)
%xASL_init_Process Initialization before ExploreASL_Process
%
% FORMAT: x = xASL_init_Process(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Initialization before ExploreASL_Process.
% 1. 
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright (c) 2015-2022 ExploreASL


    %% Initialization

    %% 1. Print the hyperlink
	if ~isdeployed && usejava('desktop') % true if the Matlab GUI is loaded, false when in CLI with or without Java VM
        disp('<a href="https://exploreasl.github.io/Documentation/latest/Tutorials-Processing; ">Click here for the ExploreASL processing tutorial</a>');
        disp('<a href="https://exploreasl.github.io/Documentation/latest/ProcessingParameters; ">Click here for the ExploreASL processing settings overview</a>');
	else % text only
        fprintf('Examples of processing-parameter settings are at: https://exploreasl.github.io/Documentation/latest/Tutorials-Processing\n');
        fprintf('A full explanation of processing parameters is @: https://exploreasl.github.io/Documentation/latest/ProcessingParameters\n');
	end
    

    %% Data loading

	% Go to ExploreASL folder
    cd(x.opts.MyPath);


    %% Initialize x struct
    x = xASL_init_SubStructs(x);

    
    % by default we read rawdata, if this folder doesn't exist, we read the derivatives
    x.opts.bReadRawdata = true;
    if ~xASL_exist(x.dir.RawData, 'dir')
        x.opts.bReadRawdata = false;
        fprintf('%s\n', ['rawdata folder missing: ' x.dir.RawData]);
        
        if ~xASL_exist(x.dir.xASLDerivatives)
            fprintf('%s\n', ['derivatives ExploreASL folder missing: ' x.dir.xASLDerivatives]);

            error('Both rawdata and derivatives folders missing, cannot continue processing');

            x.opts.bLoadData = false;
            x.opts.bLoadableData = false;
        end
    end
    
    x.opts.bLoadableData = true;

    % Create logging directory if it does not exist already
    xASL_adm_CreateDir(fullfile(x.dir.xASLDerivatives, 'log'));



    %% 1. Check basic directories
    if isempty(x.dir.DatasetRoot)
	    error('x.dir.DatasetRoot is a required parameter.');
    end
    
    % Verify that the rawdata subfolder exists
    if ~exist(fullfile(x.dir.DatasetRoot,'rawdata'), 'dir')
        warning(['Rawdata folder missing: ' fullfile(x.dir.DatasetRoot,'rawdata')]);
        return;
    end
    
    % Check derivatives of ExploreASL
    if ~isfield(x.dir,'xASLDerivatives')
        error('Missing xASL derivatives field...');
    end

    %% 2. Create the derivatives directory
    if exist(x.dir.xASLDerivatives, 'dir')
        fprintf('%s\n', [x.dir.xASLDerivatives ' already exists']);
        fprintf('%s\n', 'Note that all pre-existing derivative subject folders will be overwritten,');
        fprintf('%s\n', 'unless BIDS2Legacy lock files already exist for a subject');
    else
        xASL_adm_CreateDir(x.dir.xASLDerivatives);
    end
    
    
    %% 3. Load BIDS configuration for file renaming
    x.modules.bids2legacy.bidsPar = xASL_bids_Config;
    
        
    % This is the line needed by xASL_init_Iteration for BIDS2Legacy
    x.D.ROOT = x.dir.DatasetRoot;
    



end