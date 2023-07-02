function [result, x] = xASL_module_BIDS2Legacy(x)
%xASL_module_BIDS2Legacy BIDS2LEGACY conversion script which calls xASL_wrp_BIDS2Legacy.
%
% FORMAT: [x] = xASL_module_BIDS2Legacy(x);
%
% INPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    BIDS to Legacy conversion script which calls xASL_module_BIDS2Legacy.
%
% 0. Input check
% 1. Run legacy conversion
% 2. Finalize and unlock mutex for this module
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
%
% __________________________________
% Copyright (c) 2015-2023 ExploreASL


    %% 0. Initialization
    
    % Make sure that logging is still active
    if isfield(x.dir,'diaryFile')
        diary(x.dir.diaryFile);
    end
    
    x = xASL_init_InitializeMutex(x, 'BIDS2Legacy'); % Start Mutex
    result = true; % Default for result
    xASL_adm_BreakString('BIDS to ExploreASL LEGACY CONVERSION'); % Print feedback
	

    %% 1. Run legacy conversion
    if ~x.mutex.HasState('010_BIDS2LEGACY')
        x = xASL_wrp_BIDS2Legacy(x, 1);
        x.mutex.AddState('010_BIDS2LEGACY');
    elseif x.mutex.HasState('010_BIDS2LEGACY')
        fprintf('%s\n', ['BIDS2Legacy already done, skipping ' x.SUBJECT '.    ']);
    end

    
    %% 2. Finalize and unlock mutex for this module
    x.mutex.AddState('999_ready');
    x.mutex.Unlock();
        
end