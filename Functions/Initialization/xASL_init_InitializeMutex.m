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
% Copyright 2015-2023 ExploreASL

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
% xASL_GoNoGo initializes the x.mutex object and states, and it deals with
% locking of the current module only
x.mutex = xASL_GoNoGo(x.dir.LockDir);

%% 3) Check if this module is locked by another process
if ~x.mutex.Lock(x.settings.MUTEXID)
	fprintf(2,['ERROR in module_' x.ModuleName ': mutex is locked: %s in %s\n'], x.settings.MUTEXID, x.dir.LockDir);
	fprintf('\n');
    fprintf('This means that this module is currently being parallel processed by another Matlab instance/worker\n');
	fprintf('If this is not the case, the locked folder needs to be removed before proceeding\n');
	fprintf('Also, check that there is no filesystem permission issue\n');
    fprintf('Otherwise this error can be ignored\n');
    fprintf('\n');

    x.mutex.bAnyModuleLocked = true;
else
    
    %% 4) Check if any module for this subject is locked
    % Here we check if any other module for the same subject is locked.
    % So whereas xASL_GoNoGo in step 2 above only checks the current
    % module, here we check if ANY module for this subject is locked.

    if ~strcmpi(x.ModuleName, 'population')

        % If the current module is the population module, we don't need to
        % perform this check because the population module is cannot be
        % iterated over subjects, it is always run once.
        % Plus, at the start of the population module itself, we already
        % check if nWorkers>1, in which case we skip the population module
        % (== the population module should not run in parallel).
        if ~isfield(x, 'mutex') || ~isprop(x.mutex, 'Root')
            warning('mutex field missing, something odd going wrong');
        else
            [~, subjectName] = fileparts(x.mutex.Root);
            lockDir = fullfile(x.dir.xASLDerivatives, 'lock');
						
			% Go through all module that have lock dirs already and are thus potentially locked
			checkModuleList = xASL_adm_GetFileList(lockDir, '.*', 'FPList', [], 1);
            for iCheckModule=1:length(checkModuleList)
                [~, checkModuleName] = fileparts(checkModuleList{iCheckModule});

				% Verify that it is really a directory with an ASL module lock files
                if strcmp(checkModuleName(1:12), 'xASL_module_')
                    checkModuleName = checkModuleName(13:end); % Extract the name of the module
                    if strcmpi(checkModuleName, x.ModuleName)
                        % this is the current module, so we can continue if this is locked
						
						% Note that we are checking here only for cross-module activity. We do not care if there is an activity within the same module, because 
						% a parallel activity within the same module only matters if exactly the same visit/session/subject would have been processed and this behavior
						% is protected by the module locking itself.
						
                    elseif strcmpi(checkModuleName, 'population') || strcmpi(checkModuleName, 'import')
                        % for population or import modules, we can continue if other modules are locked
					else
						% We are in a BIDS2Legacy, structural, or ASL module and we are checking for activity in one of those other modules
						
						% By default, no locks are detected
						bLockedFolders = false;
						
						if ~strcmp(checkModuleName, 'BIDS2Legacy')
							% We are checking for activity in one of the processing modules
							
							if strcmp(x.ModuleName, 'BIDS2Legacy')
								% If we are now in the BIDS2Legacy module, we don't have a visit suffix, so we have to search through all visits to find locks
								% In the future, BIDS2Legacy will also have visits and this option will be removed
								checkModulePath = xASL_adm_GetFileList(checkModuleList{iCheckModule}, ['^' subjectName '_\d+$'], 'FPListRec', [], 1);
							else
								% If we are in a processing module, we only check for activity in the same visit
								checkModulePath = xASL_adm_GetFileList(checkModuleList{iCheckModule}, ['^' subjectName '$'], 'FPListRec', [], 1);
							end
							
							% Check all the necessary visits for the presence of locks
							for iVisit = 1:length(checkModulePath)
								if strcmpi(checkModuleName, 'structural') || strcmpi(checkModuleName, 'LongReg') || strcmpi(checkModuleName, 'DARTEL') || strcmpi(checkModuleName, 'Analyze')
									% For structural module, we only check a single session
									if ~isempty(xASL_adm_GetFileList(fullfile(checkModulePath{iVisit}, ['xASL_module_' checkModuleName]), '^locked$', 'FPListRec', [], 1))
										bLockedFolders = true;
									end
								elseif strcmpi(checkModuleName, 'asl') || strcmpi(checkModuleName, 'func') || strcmpi(checkModuleName, 'dwi') 
									% For other modules, we check for multiple sessions
									checkSessionPath = xASL_adm_GetFileList(checkModulePath{iVisit}, ['xASL_module_' checkModuleName '_' checkModuleName '_\d+$'], 'FPListRec', [], 1);
									for iSession = 1:length(checkSessionPath)
										if ~isempty(xASL_adm_GetFileList(checkSessionPath{iVisit}, '^locked$', 'FPListRec', [], 1))
											bLockedFolders = true;
										end
									end
								else
									error(['Unknown module ' checkModuleName]);
								end
								
								
							end
									
						elseif ~strcmp(x.ModuleName, 'BIDS2Legacy') && strcmp(checkModuleName, 'BIDS2Legacy')
							% In case we are in a processing module (not in BIDS2Legacy) and we are checking for locked BIDS2Legacy, we have to remove the suffix and check for all cases for this subject
							% In the future, we will update the BIDS2Legacy to lock visit wise, we will not remove this lock and we will only check for a specific visit
							ivisitSuffix = regexp(subjectName, '_\d+$');
							if  ~isempty(ivisitSuffix)
								subjectName = subjectName(1:ivisitSuffix-1);
							end
							checkModulePath = fullfile(checkModuleList{iCheckModule}, subjectName, ['xASL_module_' checkModuleName]);
							bLockedFolders = ~isempty(xASL_adm_GetFileList(checkModulePath, '^locked$', 'FPListRec', [], 1));
						end
                        
                        if bLockedFolders
    
                            fprintf(['ERROR in module_' x.ModuleName ', there is another module locked for the same subject: ' checkModuleName '\n']);
                            fprintf('\n');
	                        fprintf('This means that this module is currently being parallel processed by another Matlab instance/worker\n');
	                        fprintf('If this is not the case, the locked folder needs to be removed before proceeding\n');
	                        fprintf('Also, check that there is no filesystem permission issue\n');
                            fprintf('Otherwise this error can be ignored\n');
                            fprintf('\n');
            
                            x.mutex.bAnyModuleLocked = true;
                        end
                    end
                end
            end

        end

    end

end


end