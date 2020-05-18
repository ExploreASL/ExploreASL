function [x] = xASL_init_InitializeMutex(x, ModuleName)
%xASL_init_InitializeMutex Initialize mutex system for ExploreASL modules
%
% FORMAT: [x] = xASL_init_InitializeMutex(x, ModuleName)
%
% INPUT:
%   x           - struct containing pipeline environment parameters (REQUIRED)
%   ModuleName  - name of the module (options are 'Structural', 'ASL' 'Population')
%                 expert options include 'DARTEL', "LongReg'
%                 module name are prefixed by 'xASL_module_'
%                 (REQUIRED)
% OUTPUT:
%   x           - struct containing pipeline environment parameters
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function initializes the mutex/lock system of
% ExploreASL for a module. Mutex (for mutual exclusion) is a
% synchronization mechanism for enforcing limits of access to data (here a
% module for a single scan) to allow parallelization. It also allows
% stopping and continuing of ExploreASL. This function runs the following
% steps:
% 1) Lock folder management
% 2) Initialize mutex object
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [x] = xASL_init_InitializeMutex(x, 'xASL_module_Structural')
% __________________________________
% Copyright 2015-2020 ExploreASL



%% --------------------------------------------------------
%% 1) Lock folder management
x.ModuleName  = ModuleName;
x.result = false;
if x.RERUN
    [status, message] = rmdir(x.LockDir,'s');
    if status~=1 && exist(x.LockDir,'dir') % backwards compatibility
        fprintf(2,['ERROR in module_' x.ModuleName ': could not remove lock folder:\n%s\n'], message);
        return; % exit with an error
    end
    x.bOverwrite = true; % re-running makes no sense if you're not overwriting existing files...
end

if xASL_adm_CreateDir(x.LockDir,3)<0
    fprintf(2,['ERROR in module_' x.ModuleName ': could not create lock folder: \n%s\n'],x.LockDir);
    return; % exit with an error
end

if ~exist(x.LockDir,'dir')
    fprintf(['ERROR in module_' x.ModuleName ': wrong LockDir assigned\n']);
end

%% --------------------------------------------------------
%% 2) Initialize mutex object
x.mutex = xASL_GoNoGo(x.LockDir);
if ~x.mutex.Lock(x.MUTEXID)
    fprintf('Also, check that there is no filesystem permission issue\n');
	fprintf(2,['ERROR in module_' x.ModuleName ': mutex is locked: %s in %s\n'], x.MUTEXID, x.LockDir);
	fprintf('This means that this module is currently being parallel processed by another Matlab instance/worker\n');
	fprintf('If this is not the case, the locked folder needs to be removed before proceeding\n');
	fprintf('Also, check that there is no filesystem permission issue\n');
    fprintf('Otherwise this error can be ignored\n');
	return;
end


end

