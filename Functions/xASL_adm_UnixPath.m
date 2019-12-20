function [PathOut] = xASL_adm_UnixPath(PathIn)
%xASL_adm_UnixPath Converts PC path to Linux path, also taking into
%account Windows Subsystem for Linux (WSL)
% Can be useful for inline use, avoiding to have to declare paths multiple
% times for different OSes
% WARNING: when using in non-Linux OS (e.g. Windows), this only works
% when this path is used in a Linux subsystem (e.g. WSL)

PathOut = strrep(PathIn, '\','/'); % first convert slashes

if ispc % only change when we find WSL
    [status, result] = system('wsl ls');
        if status==0 && strcmp(PathOut(2),':')
            PathOut = ['/mnt/' lower(PathOut(1)) '/' PathOut(4:end)];
        end
end



end

