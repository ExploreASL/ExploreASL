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
% 5. Create directories
% 6. Load BIDS configuration for file renaming
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright (c) 2015-2023 ExploreASL



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
    if ~isfield(x.opts, 'bReadRawdata') || isempty(x.opts.bReadRawdata) || ~islogical(x.opts.bReadRawdata)
        % If the developer has not set this parameter, we try to detect whether we should read from /rawdata or /derivatives/ExploreASL 
        % Usually, this parameter will not exist
        if xASL_exist(x.dir.RawData, 'dir')
            x.opts.bReadRawdata = true; %  by default we read rawdata
            fprintf('%s\n', '/rawdata folder detected, trying to load subjects from BIDS');
        elseif xASL_exist(x.dir.xASLDerivatives, 'dir')
            x.opts.bReadRawdata = false;
            fprintf('%s\n', '/derivatives/ExploreASL folder detected, trying to load subjects from previous ExploreASL processing');
        else
            error('No rawdata or derivatives found to load subjects from, first convert dicoms 2 BIDS');
        end
    elseif x.opts.bReadRawdata && ~xASL_exist(x.dir.RawData, 'dir')
        % Optionally, developers have set this parameter to force reading subjects from /rawdata
        warning(['rawdata folder missing: ' x.dir.RawData]);
        error('rawdata did not exist, if you want to read subjects from /derivatives/ExploreASL, set x.opts.bReadRawdata to false in dataPar.json');
    elseif ~x.opts.bReadRawdata && ~xASL_exist(x.dir.xASLDerivatives, 'dir')
        % Optionally, developers have set this parameter to force reading subjects from /derivatives/ExploreASL
        warning(['ExploreASL derivatives folder missing: ' x.dir.xASLDerivatives]);
        error('/derivatives/ExploreASL did not exist, if you want to load BIDS subjects from /rawdata, set x.opts.bReadRawdata to true in dataPar.json');
    end


    %% 5. Create directories
    xASL_adm_CreateDir(fullfile(x.dir.xASLDerivatives, 'log')); % Create logging directory if it does not exist
    xASL_adm_CreateDir(x.dir.xASLDerivatives); % Create ExploreASL derivatives directory if it does not exist
    
    
    %% 6. Load BIDS configuration for file renaming
    x.modules.bids2legacy.bidsPar = xASL_bids_Config;

    
end