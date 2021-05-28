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
%
% 1. Lock folder management
% 2. Initialize mutex object
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [x] = xASL_init_InitializeMutex(x, 'xASL_module_Structural')
% __________________________________
% Copyright 2015-2020 ExploreASL

% Check inputs
if ~isfield(x.settings,'RERUN') || ~isfield(x.settings,'MUTEXID') || ~isfield(x.dir,'LockDir')
    warning(['Seemingly you are using xASL_module_' ModuleName ' without initialized dbSettings, consider running xASL_Iteration instead...']);
end

%% --------------------------------------------------------
%% 1) Lock folder management
x.ModuleName  = ModuleName;
x.result = false;
if isfield(x.settings,'RERUN')
    if x.settings.RERUN
        [status, message] = rmdir(x.dir.LockDir,'s');
        if status~=1 && exist(x.dir.LockDir,'dir') % backwards compatibility
            fprintf(2,['ERROR in module_' x.ModuleName ': could not remove lock folder:\n%s\n'], message);
            return; % exit with an error
        end
        x.settings.bOverwrite = true; % re-running makes no sense if you're not overwriting existing files...
    end
else
    error('RERUN field of x structure not defined...');
end

if isfield(x.dir,'LockDir')
    if xASL_adm_CreateDir(x.dir.LockDir,3)<0
        fprintf(2,['ERROR in module_' x.ModuleName ': could not create lock folder: \n%s\n'],x.dir.LockDir);
        return; % exit with an error
    end

    if ~exist(x.dir.LockDir,'dir')
        fprintf(['ERROR in module_' x.ModuleName ': wrong LockDir assigned\n']);
    end
else
    error('LockDir field of x structure not defined...');
end

%% --------------------------------------------------------
%% 2) Initialize mutex object
x.mutex = xASL_GoNoGo(x.dir.LockDir);
if ~x.mutex.Lock(x.settings.MUTEXID)
    fprintf('Also, check that there is no filesystem permission issue\n');
	fprintf(2,['ERROR in module_' x.ModuleName ': mutex is locked: %s in %s\n'], x.settings.MUTEXID, x.dir.LockDir);
	fprintf('This means that this module is currently being parallel processed by another Matlab instance/worker\n');
	fprintf('If this is not the case, the locked folder needs to be removed before proceeding\n');
	fprintf('Also, check that there is no filesystem permission issue\n');
    fprintf('Otherwise this error can be ignored\n');
	return;
end


end

