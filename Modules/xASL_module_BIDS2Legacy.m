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
% 1. Input check
% 2. Run legacy conversion
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
%
% __________________________________
% Copyright (c) 2015-2023 ExploreASL


    %% 1. Input check
    
    % Make sure that logging is still active
    if isfield(x.dir,'diaryFile')
        diary(x.dir.diaryFile);
    end
    
    % Start Mutex
    x = xASL_init_InitializeMutex(x, 'BIDS2Legacy');
    
    % Define lock states
    StateName{1} = '010_BIDS2LEGACY';
    
    % Default for result
    result = true;
    
    % Print feedback
    xASL_adm_BreakString('BIDS to ExploreASL LEGACY CONVERSION');
	
    %% 2. Run legacy conversion
    iState = 1;
    if ~x.mutex.HasState(StateName{1})
        x = xASL_wrp_BIDS2Legacy(x, 1);
        x.mutex.AddState(StateName{iState});
    elseif x.mutex.HasState(StateName{1})
        fprintf('%s\n', ['BIDS2Legacy already done, skipping ' x.SUBJECT '...']);
    end
    
    %% 3. Clean-Up
    x.mutex.AddState('999_ready');
    x.mutex.Unlock();
        
end
