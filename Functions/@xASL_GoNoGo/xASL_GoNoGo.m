% Class that implements a process locking mechanism.
%
% DESCRIPTION: % xASL_GoNoGo is a class which implements a process locking mechanism to prevent
% (parts of) scripts to run twice. Simple files are used as flags/semaphores.
% However, since matlab doesn't support any form of file locking, a subfolder
% is used as semaphore. A xASL_GoNoGo lock is identified by a unique ID. The ID should
% be a valid directory name for the filesystem that is used to store the folder
% structure. The root folder must be specified before locking can occur.
% After a lock is set, one or more status files can be created to set one or 
% more independend states.
% EXAMPLE: 
% x.S = xASL_GoNoGo('C:\Temp\')
% x.S.Lock('MySemaphore')
% x.S.State = 'busy';
% x.S.State = 'ready';
% x.S.State
% x.S.State = { 'busy', 'run1' }
% x.S.HasState('run1')
% x.S.HasState({ 'busy', 'run1', 'ready' })
% x.S.Unlock
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL
%

classdef xASL_GoNoGo < handle
    properties
        Root = []; % root directory for storing lock files
    end
    
    properties (SetAccess = private)
        ID = [];
        Locked = false;
    end
    
    properties (GetAccess = private, SetAccess = private)
    end
   
    properties (Dependent = true)
        State; % assigned state, stored as empty file using <State> as filename
    end
    
    methods
        function obj = xASL_GoNoGo(folder)
            if nargin == 1 && ischar(folder) && exist(folder,'dir')
                % if the only input argument is a directory, try to load it.
                obj.Root = folder;
            else
                % otherwise do nothing; creates empty xASL_GoNoGo object
                disp('Empty xASL_GoNoGo object created. Use obj.path = folder to set the folder.')
            end
        end % xASL_GoNoGo
        
        function delete(obj)
            obj.Unlock;
        end
        
        function p = get.Root(obj)
            p = obj.Root;
        end
       
        function set.Root(obj, folder)
            if obj.IsLocked
                % ERROR
            else
                if ~exist(folder,'dir')
                    mkdir(folder);
                end
                obj.Root = folder;
            end
        end
        
        function id = get.ID(obj)
            id = obj.ID;
        end

        function state = get.State(obj)
            state = {};
            if obj.Locked && isdir(obj.Root) && ~isempty(obj.ID)
                E = dir(fullfile(obj.Root, obj.ID, '*.status'));
                if ~isempty(E)
                    state = { E.name };
                    state = regexprep(state,'\.status$','');
                end
            else
                % WARNING: not locked yet
            end
        end
        
        function set.State(obj, newState)
            % Add state(s) if not already set, and remove previously set state(s)
            % newState can be either a string or cell array.
            if obj.Locked && isdir(obj.Root) && ~isempty(obj.ID)
                oldState = obj.State;
                if ~iscell(newState)
                    newState = { newState };
                end
                addStates = setdiff(newState,oldState);
                for ii=1:length(addStates)
                    fid = fopen(fullfile(obj.Root, obj.ID, [addStates{ii} '.status']), 'w');
                    fclose(fid);
                end
                delStates = setdiff(oldState,newState);
                for ii=1:length(delStates)
                    delete(fullfile(obj.Root, obj.ID, [delStates{ii} '.status']));
                end
            else
                % WARNING: not locked yet
            end
        end

        function AddState(obj, someState)
            % Add state(s) to set
            if obj.Locked && isdir(obj.Root) && ~isempty(obj.ID)
                oldState = obj.State;
                if ~ismember(someState,oldState)
                    obj.State = union(oldState, someState);
                end
            else
                % WARNING: not locked yet
            end
        end
        
        function DelState(obj, someState)
            % Add state(s) to set
            if obj.Locked && isdir(obj.Root) && ~isempty(obj.ID)
                oldState = obj.State;
                if ismember(someState,oldState)
                    obj.State = setdiff(oldState, someState);
                end
            else
                % WARNING: not locked yet
            end
        end
        
        function b = HasState(obj, someState)
            % Check if state(s) are set
            % Returns a logical vector if someState is a cell array.
            if obj.Locked && isdir(obj.Root) && ~isempty(obj.ID)
                oldState = obj.State;
                if ~iscell(someState)
                    someState = { someState };
                end
                b = ismember(someState,oldState);
            else
                % WARNING: not locked yet
            end
        end
        
        function b = Lock(obj, ID)
            b = false;
            if isempty(obj.Root)
                % ERROR
            else
                if ~isempty(obj.ID)
                    obj.Unlock()
                end
                [status, msg] = mkdir(obj.Root, ID);
                if status==1
%                     if isempty(msg)
%                         % a new folder was created, so we are free to go...
%                     else
%                         % folder already exists, this could become tricky...
%                     end
                    [status, msg] = mkdir(fullfile(obj.Root, ID), 'locked');
                    if status==1
                        if isempty(msg)
                            % a new folder was created, so we are free to go...
                            obj.ID = ID;
                            obj.Locked = true;
                            b = true;
                        else
                            % the lock folder already exists, refuse...
                        end
                    else
                        % ERROR
                        % could not create lock directory
                    end
                else
                    % ERROR
                    % could not create ID directory
                end
            end
        end
     
        function b = Unlock(obj)
            b       = false;
            
            if obj.Locked && isdir(obj.Root) && ~isempty(obj.ID)
                LockDir        = fullfile(obj.Root, obj.ID, 'locked');
                if      exist(LockDir,'dir')
                        [status, msg] = rmdir(LockDir,'s');
                elseif  exist(LockDir)>0
                        status  = 0;
                else    % if there is no LockDir, assume it has been removed before
                        status  = 1;
                end
                
                if  status==1
                    b = true;
                else
                    error('xASL_GoNoGo:unlockFailed', 'xASL_GoNoGo:Unlock failed: %s', msg);
                end
            end
            % even on errors...
            obj.ID = [];
            obj.Locked = false; 
        end
        
        function b = get.Locked(obj)
            b = obj.Locked;
        end
        
        function b = IsLocked(obj)
            % for internal use in properties (set.Root)
            b = obj.Locked;
        end
    end
    
    methods (Access = 'private') % Access by class members only
        
    end
    
end % classdef xASL_GoNoGo
