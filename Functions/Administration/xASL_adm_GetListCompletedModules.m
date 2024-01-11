function [listStructDone, listASLDone] = xASL_adm_GetListCompletedModules(x)
%xASL_adm_GetListCompletedModules This function obtains a list of completed modules
%
% FORMAT:       [listStructDone, listASLDone] = xASL_adm_GetListCompletedModules(x)
% 
% INPUT:        x               - x structure of ExploreASL
%
% OUTPUT:       listStructDone  - list of subjects for which the structural module is completed
%               listASLDone     - list of sessions/runs for which the ASL module is completed
%
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  For both structural and ASL modules, we get a list of subject/sessions of which these modules are fully completed.
%               This can be useful for subsequent group processing in e.g., the ExploreASL population module.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      [listStructDone, listASLDone] = xASL_adm_GetListCompletedModules(x);
% __________________________________
% Copyright 2015-2024 ExploreASL

%% First, we create a list of subjects and sessions of which the structural and ASL modules, respectively, are completed

dirLockRoot = fullfile(x.dir.xASLDerivatives, 'lock');
dirLockStruct = fullfile(dirLockRoot, 'xASL_module_Structural');
dirLockASL = fullfile(dirLockRoot, 'xASL_module_ASL');

listStructSubjects = xASL_adm_GetFileList(dirLockStruct, '.*', 'List', [], 1);
listStructDone = {};

for iList = 1:length(listStructSubjects)
    lockFile999 = fullfile(dirLockStruct, listStructSubjects{iList}, 'xASL_module_Structural', '999_ready.status');
    if exist(lockFile999, 'file')
        listStructDone{end+1} = listStructSubjects{iList};
    end
end

listASLSubjects = xASL_adm_GetFileList(dirLockASL, '.*', 'List', [], 1);
listASLDone = {};
for iList = 1:length(listASLSubjects)
    sessions = xASL_adm_GetFileList(fullfile(dirLockASL, listASLSubjects{iList}), '^xASL_module_ASL_.*', 'List', [], 1);
    for iSession=1:length(sessions)
        lockFile999 = fullfile(dirLockASL, listASLSubjects{iList}, sessions{iSession}, '999_ready.status');
        if exist(lockFile999, 'file')
            listASLDone{end+1} = [listASLSubjects{iList} '_' sessions{iSession}(length('xASL_module_ASL_')+1:end)];
        end
    end
end


end