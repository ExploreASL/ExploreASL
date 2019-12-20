function x = xASL_init_GenericMutexModules(x, ModuleName)
%xASL_init_GenericMutexModules Runs the generic loopDB stuff that Paul
% created to start a module with


%% Mutex
x.ModuleName  = ModuleName;
x.result = false;
if  x.RERUN
    [status, message] = rmdir(x.LOCKDIR,'s');
    if status~=1 && exist(x.LOCKDIR,'dir') % backwards compatibility
        fprintf(2,['ERROR in module_' x.ModuleName ': could not remove touch folder:\n%s\n'],message);
        return % exit with an error
    end
    x.OVERWRITE = true; % re-running makes no sense if you're not overwriting existing files...
end

if  xASL_adm_CreateDir(x.LOCKDIR,3)<0
    fprintf(2,['ERROR in module_' x.ModuleName ': could not create touch folder: \n%s\n'],x.LOCKDIR);
    return % exit with an error
end

%% 0.3 create a 'lock' folder to keep track of which analysis step ran (or still runs)
% Create a mutex folder to assure that everything will run only once
x.mutex = xASL_GoNoGo(x.LOCKDIR);
if ~x.mutex.Lock(x.MUTEXID)
    fprintf(2,['ERROR in module_' x.ModuleName ': mutex is locked: %s in %s\n'], x.MUTEXID, x.LOCKDIR);
    return;
end


end

