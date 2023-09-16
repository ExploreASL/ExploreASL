function [x, bLocked] = xASL_init_InitializeMutex(x, ModuleName)
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

%% 3) Check if this module is locked by another process
bLocked = 0;

if ~x.mutex.Lock(x.settings.MUTEXID)
	fprintf(2,['ERROR in module_' x.ModuleName ': mutex is locked: %s in %s\n'], x.settings.MUTEXID, x.dir.LockDir);
	fprintf('\n');
    fprintf('This means that this module is currently being parallel processed by another Matlab instance/worker\n');
	fprintf('If this is not the case, the locked folder needs to be removed before proceeding\n');
	fprintf('Also, check that there is no filesystem permission issue\n');
    fprintf('Otherwise this error can be ignored\n');
    fprintf('\n');

    bLocked = 1;
else
    
    %% 4) Check if any module for this subject is locked
    if ~strcmpi(x.ModuleName, 'population')
        % the population module is subject-independent
        if ~isfield(x, 'mutex') || ~isprop(x.mutex, 'Root')
            warning('mutex field missing, something odd going wrong');
        else
            [~, subjectName] = fileparts(x.mutex.Root);
            lockDir = fullfile(x.dir.xASLDerivatives, 'lock');
            moduleDirs = xASL_adm_GetFileList(lockDir, '.*', 'FPList', [], 1);

            for iModule=1:length(moduleDirs)
                if ~strcmpi(moduleDirs{iModule}, 'population') || ~strcmpi(moduleDirs{iModule}, 'import')
                    thisModuleDir = fullfile(moduleDirs{iModule}, subjectName);

                    otherLockedFolders = xASL_adm_GetFileList(thisModuleDir, '^locked$', 'FPListRec', [], 1);
                    if ~isempty(otherLockedFolders)

                        fprintf(['ERROR in module_' x.ModuleName ', there is another module locked for the same subject:\n']);
                        for iOther=1:length(otherLockedFolders)
                            fprintf('%s\n', otherLockedFolders{iOther});
                        end
                        fprintf('\n');
	                    fprintf('This means that this module is currently being parallel processed by another Matlab instance/worker\n');
	                    fprintf('If this is not the case, the locked folder needs to be removed before proceeding\n');
	                    fprintf('Also, check that there is no filesystem permission issue\n');
                        fprintf('Otherwise this error can be ignored\n');
                        fprintf('\n');
        
                        bLocked = 1;
                    end % ~isempty(otherLockedFolders)
                end % ~strcmpi
            end % iModule=1:length(moduleDirs)

        end

    end

end


end