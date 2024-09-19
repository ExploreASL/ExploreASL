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
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________




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
    if ~isfield(x.opts, 'bReadRawdata') || isempty(x.opts.bReadRawdata)
        % If the developer has not set this parameter, we try to detect whether we should read from /rawdata or /derivatives/ExploreASL 
        % Usually, this parameter will not exist

        if xASL_exist(x.dir.RawData, 'dir')
            % we assume we start with BIDS data
            
            if xASL_exist(x.dir.xASLDerivatives, 'dir') && isequal(x.opts.bProcess(:), [0; 0; 1])
                 % when we also have pre-existing derivatives
                % when only running the population module, we skip BIDS2Legacy
			    x.opts.bReadRawdata = false; %
			    fprintf('%s\n', '/rawdata and /derivatives/ExploreASL folders detected, and running Population module only --> We read only derivatives and ignore rawdata');
            else
                % By default we read rawdata
                x.opts.bReadRawdata = true;
                fprintf('%s\n', '/rawdata folder detected, trying to load subjects from BIDS');
            end
        elseif xASL_exist(x.dir.xASLDerivatives, 'dir')
			% Only derivatives exists, we skip rawdata and read only derivatives
            x.opts.bReadRawdata = false;
            fprintf('%s\n', '/derivatives/ExploreASL folder detected, trying to load subjects from previous ExploreASL processing');
        else
            error('No rawdata or derivatives found to load subjects from, first convert DICOMs 2 BIDS');
        end
    elseif x.opts.bReadRawdata && ~xASL_exist(x.dir.RawData, 'dir')
        % Optionally, developers have set this parameter to force reading subjects from /rawdata
        warning(['rawdata folder missing: ' x.dir.RawData]);
        error('rawdata did not exist, if you want to read subjects from /derivatives/ExploreASL, set x.opts.bReadRawdata to false in dataPar.json');
    elseif x.opts.bReadRawdata && isequal(x.opts.bProcess(:), [0; 0; 1])
        % We force reading rawdata, but we run population module only
        warning('Rawdata reading set to true, but running a Population model only. Any newly converted data will not be processed before the population module');
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