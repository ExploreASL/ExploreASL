function [x] = xASL_init_DataLoading(x)
%xASL_init_DataLoading Load dataset by adding relevant fields to xASL x struct
%
% FORMAT: [x] = xASL_init_DataLoading(x)
% 
% INPUT:
%   x          - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x          - ExploreASL x structure
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Here, we load all the ExploreASL derivatives data into the Matlab x structure, 
% such that we can use it for loading the data/subjects.
%
% 1. Print the hyperlink
% 2. Go to ExploreASL folder
% 3. Initialize x struct
% 4. Which data to read
% 5. Check basic directories
% 6. Create the derivatives directory
% 7. Load BIDS configuration for file renaming
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright (c) 2015-2023 ExploreASL


    %% Initialization

    %% 1. Print the hyperlink
	if ~isdeployed && usejava('desktop') % true if the Matlab GUI is loaded, false when in CLI with or without Java VM
        disp('<a href="https://exploreasl.github.io/Documentation/latest/Tutorials-Processing; ">Click here for the ExploreASL processing tutorial</a>');
        disp('<a href="https://exploreasl.github.io/Documentation/latest/ProcessingParameters; ">Click here for the ExploreASL processing settings overview</a>');
	else % text only
        fprintf('Examples of processing-parameter settings are at: https://exploreasl.github.io/Documentation/latest/Tutorials-Processing\n');
        fprintf('A full explanation of processing parameters is @: https://exploreasl.github.io/Documentation/latest/ProcessingParameters\n');
	end
    

    %% 2. Go to ExploreASL folder
    cd(x.opts.MyPath);


    %% 3. Initialize x struct
    x = xASL_init_SubStructs(x);

    
    %% 4. Which data to read
    if ~isfield(x.opts, 'bReadRawdata')
        x.opts.bReadRawdata = true;
        %  by default we read rawdata
    end
    
    if x.opts.bReadRawdata && ~xASL_exist(x.dir.RawData, 'dir')
        warning(['Rawdata folder missing: ' x.dir.RawData]);
        error('rawdata did not exist, if you want to only process /derivatives/ExploreASL, set x.opts.bReadRawdata to false in dataPar.json');
    elseif ~x.opts.bReadRawdata && ~xASL_exist(x.dir.xASLDerivatives, 'dir')
        warning(['ExploreASL derivatives folder missing: ' x.dir.xASLDerivatives]);
        error('/derivatives/ExploreASL did not exist, if you want to load BIDS data from /rawdata, set x.opts.bReadRawdata to true in dataPar.json');
    end

    % Create logging directory if it does not exist already
    xASL_adm_CreateDir(fullfile(x.dir.xASLDerivatives, 'log'));



    %% 5. Check basic directories
    if isempty(x.dir.DatasetRoot)
	    error('x.dir.DatasetRoot is a required parameter');
    end
    
    % Check derivatives of ExploreASL
    if ~isfield(x.dir,'xASLDerivatives')
        error('Missing xASL derivatives field...');
    end

    %% 6. Create the derivatives directory
    if exist(x.dir.xASLDerivatives, 'dir')
        fprintf('%s\n', [x.dir.xASLDerivatives ' already exists']);
        fprintf('%s\n', 'Note that all pre-existing derivative subject folders will be overwritten,');
        fprintf('%s\n', 'unless BIDS2Legacy lock files already exist for a subject');
    else
        xASL_adm_CreateDir(x.dir.xASLDerivatives);
    end
    
    
    %% 7. Load BIDS configuration for file renaming
    x.modules.bids2legacy.bidsPar = xASL_bids_Config;

    
end