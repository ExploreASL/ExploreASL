function [x] = xASL_init_RemoveLockDirs(x)
%xASL_init_RemoveLockDirs Remove 'lock-dir' if present from aborted previous run, for current subjects only
%
% FORMAT: [x] = xASL_init_RemoveLockDirs(x)
% 
% INPUT:
%   x        - ExploreASL x struct (STRUCT, REQUIRED)
%
% OUTPUT:
%   x        - ExploreASL x struct 
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Remove 'lock-dir' if present from aborted previous run, for current subjects only.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        [x] = xASL_init_RemoveLockDirs(x);
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% LockDir within 2 directories (e.g. T1, FLAIR or ASL)
    LockDir = fullfile(x.D.ROOT, 'lock');

    if exist(LockDir, 'dir')
        % fprintf('%s\n','Searching for locked previous ExploreASL image processing');
        LockDirFound = 0;
        LockDir = xASL_adm_FindByRegExp(fullfile(x.D.ROOT, 'lock'), {'(ASL|Structural|LongReg_T1)', x.subject_regexp, '.*module.*','^(locked)$'}, 'Match', 'Directories');
        if ~isempty(LockDir)
            warning('Locked folders were found, consider removing them before proceeding');
        end

        % LockDir within 2 directories (e.g. DARTEL)
        LockDir = xASL_adm_FindByRegExp(fullfile(x.D.ROOT, 'lock'), {'(Population|DARTEL_T1)', '.*module.*','^(locked)$'}, 'Match','Directories');
        if ~isempty(LockDir)
            warning('Locked folders were found, consider removing them before proceeding');
        end

        if LockDirFound==0
            % fprintf('%s\n', 'No locked folders found from previous ExploreASL image processing');
        end
    end

end


