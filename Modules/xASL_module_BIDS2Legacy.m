function [result, x] = xASL_module_BIDS2Legacy(x, bOverwrite, bVerbose)
%xASL_module_BIDS2Legacy Convert BIDS rawdata to ExploreASL legacy format
%
% FORMAT: [x] = xASL_module_BIDS2Legacy(x [,bOverwrite]);
%
% INPUT:
%   x             - ExploreASL x structure, containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%   bOverwrite    - Boolean indicating whether to overwrite existing files (OPTIONAL, DEFAULT = true)
%   bVerbose   - boolean, true for verbose output (OPTIONAL, DEFAULT = false)
%
% OUTPUT:
%   x             - ExploreASL x structure, containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%  result         - Boolean indicating whether module succeeded
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function converts BIDS rawdata (in pathStudy/rawdata/) 
% to xASL legacy derivative format (e.g. pathStudy/derivatives/ExploreASL/)
%
% Can be updated step-by-step when ExploreASL's derivative structure moves to BIDS
% NB: ask how Visits/session layer is defined in bids-matlab (should be
% separate layer within subjects, but now isn't?)
%
% This function performs the following steps:
%
% 0. Initialization

% 6. Finalize and unlock mutex for this module
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        x = xASL_module_BIDS2Legacy(x);
%
% __________________________________
% Copyright (c) 2015-2023 ExploreASL


    %% 0. Initialization
    if nargin<3 || isempty(bVerbose)
        bVerbose = 0; % default to not have verbose output for each subject (the mutex set up will already do this)
    end
    if nargin<2 || isempty(bOverwrite)
        bOverwrite = 1; % Default is to overwrite existing files
    end

    % Make sure that logging is still active
    if isfield(x.dir,'diaryFile')
        diary(x.dir.diaryFile); % Start diary
    end
    
    x = xASL_init_InitializeMutex(x, 'BIDS2Legacy'); % Start Mutex
    result = true; % Default for result
    xASL_adm_BreakString('BIDS to ExploreASL LEGACY CONVERSION'); % Print feedback
	

    if ~x.mutex.HasState('010_BIDS2LEGACY')







        x = xASL_wrp_BIDS2Legacy(x, bOverwrite, bVerbose); % Run BIDS2Legacy conversion








        x.mutex.AddState('010_BIDS2LEGACY');
    elseif x.mutex.HasState('010_BIDS2LEGACY')
        fprintf('%s\n', ['BIDS2Legacy already done, skipping ' x.SUBJECT '.    ']);
    end

    
    %% 6. Finalize and unlock mutex for this module
    x.mutex.AddState('999_ready'); % Add ready state
    x.mutex.Unlock(); % Unlock mutex
        
end