function status = xASL_adm_CreateDir(varargin)
% Recursively creates missing subdirectories of the given path
%
% FORMAT: status = xASL_adm_CreateDir(strPath)                  create all missing subdirs
%         status = xASL_adm_CreateDir(strPath, strBranch)       create strBranch (if missing) under existing strPath
%         status = xASL_adm_CreateDir(strPath, nMaxNewDirs)     impose limit on the number of new directories
%
% INPUT:
%   strPath	    - base path to start with or to be created
%   strBranch   - the directories to be created at strPath/strBranch
%   nMaxNewDirs - limit on the number of directories to be created
% OUTPUT:
%   status      - number of create subdirectories
%                 0 if nothing was created 
%                 -1 if creation failed.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Recursively creates missing directories at the given path or for given subdirectories, with an option
%              to limit the number of newly created directories.
% EXAMPLE: status = xASL_adm_CreateDir('c:\User\path\dir\');
%          status = xASL_adm_CreateDir('c:\User\path\file\dir',2);
%          status = xASL_adm_CreateDir('c:\User\path\','file\dir');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2019 ExploreASL

   % Admin
    if nargin == 1
		% We give single path
        strPath = varargin{1};
        nMaxNewDirs = Inf;
    elseif nargin == 2
        if ischar(varargin{2})
			% A path and subpath is given
			% Checks for the existence of the base directory
			if ~exist(varargin{1},'dir')
				error(['xASL_adm_CreateDir: Base directory does not exist: ', varargin{1}]);
			end
			
			% how many subdirs are added (ignore leading/trailing filesep)
            strBranch = varargin{2};
            cSubDirs = regexp(strBranch,'[\\/]+','split');
			nMaxNewDirs = length(cSubDirs) -  sum(cellfun('isempty',cSubDirs));
			
			% Merge the paths
            strPath = fullfile(varargin{1},strBranch);
			
        elseif isnumeric(varargin{2})
			% A path and limit is given
            strPath = varargin{1};
            nMaxNewDirs = varargin{2};
        else
            error('xASL_adm_CreateDir: Second parameter should be char or numeric.');
        end
    else
        error('xASL_adm_CreateDir: Requires one or two input parameters.');
    end
    
    ROOTdir = fileparts(fileparts(fileparts(fileparts(fileparts(fileparts(fileparts(strPath)))))));
    if ~exist(ROOTdir,'dir') && ~isempty(ROOTdir)
        warning('This is an invalid directory path');
        status = 0;
        ret = 0;
        return;
    
    elseif exist(strPath,'dir')
		% Already exists
        ret = 0;
    elseif nMaxNewDirs>0
        strParent = fileparts(strPath);
        if ~isempty(strParent) && ~exist(strParent,'dir')
            ret = xASL_adm_CreateDir(strParent,nMaxNewDirs-1);
        else
            ret = 0;
        end
        if ret>=0
             % disp(strPath)
            if ~isempty(strPath) && mkdir(strPath)==1
                ret = ret + 1;
            else
                ret = -1;
            end
        end
    else
        ret = -1; % not allowed to create so many dir levels
    end
    if nargout>=1
        status = ret;
    end
end

