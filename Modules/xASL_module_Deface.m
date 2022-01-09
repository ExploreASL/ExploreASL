function [result, x] = xASL_module_Deface(x)
%xASL_module_Deface Run defacing.
%
% FORMAT: [result, x] = xASL_module_Deface(x)
% 
% INPUT:
%   x          - ExploreASL x structure (REQUIRED, STRUCT)
%
% OUTPUT:
%   x          - ExploreASL x structure
%   result     - True for successful run of this module, false for insuccessful run
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Run defacing.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     [~, x] = xASL_init_Iteration(x,'xASL_module_Deface');
% __________________________________
% Copyright 2015-2022 ExploreASL


    %% 1. Input check
    
    % Make sure that logging is still active
    if isfield(x.dir,'diaryFile')
        diary(x.dir.diaryFile);
    end
    
    % Start Mutex
    x = xASL_init_InitializeMutex(x, 'Deface');
    
    % Define lock states
    StateName{1} = '010_DEFACE';
    StateName{2} = '020_CleanUp';
    
    % Default for result
    result = true;
    
    % Run Defacing
    iState = 1;
    if ~x.mutex.HasState(StateName{1})
        xASL_wrp_Deface(x);
        x.mutex.AddState(StateName{iState});
    elseif x.mutex.HasState(StateName{1})
        fprintf('Defacing was run before...   \n');
    end
    
    %% Clean-Up
    
    x.mutex.AddState('999_ready');
    x.mutex.Unlock();
    
    
end


