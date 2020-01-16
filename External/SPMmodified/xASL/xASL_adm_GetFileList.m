function filepaths = xASL_adm_GetFileList(strDirectory, strRegEx, mode, nRequired, bGetDirNames)
% List the list of files (or directories) in a directory, using regular expressions
%
% FORMAT: filepaths = xASL_adm_GetFileList(strDirectory[, strRegEx, mode, nRequired, bGetDirNames])
%
% INPUT:
%   strDirectory   Directory to search (REQUIRED)
%                  Optionally use wildcards to be more specific (i.e. 'C:\data\*.img')
%   strRegEx       Regular expression that will be used to filter the returned list (OPTIONAL)
%                  Examples: 
%                            '^.+$'       - matches all files (DEFAULT)
%                            '^\d+$'      - matches names that contains only digits
%                            '^.+\.nii$'  - matches names that ends with .nii
%                            '^\d{3}$'    - matches names that exist of exactly 3 digits
%                            '^\d{3}_.+$' - matches names that start with 3 digits, followed by an underscore and other characters
%                            '^PP\d+$'    - matches names that start with PP, followed by a number
%   mode           true to return full paths with path+name (OPTIONAL, DEFAULT = 'FPList')
%                  false to return only the name (= 'List')
%                  Or directly assign the spm_select compatible mode 'List' (only filenames), 'FPList' (fullpath),
%                  or 'FPListRec' (full path and recursively).
%   nRequired      Set to a value or range [a b] to define a minimum/maximum number of required files. (OPTIONAL)
%                  An exception will be thrown if this requirement is not met. Use [n Inf] to accept at least n.
%                  Default is to accept any count. You can also provide a scalar [a] requiring to accept exactly 'a' matches
%   bGetDirNames   Set to true to get directory names. (OPTIONAL)
%                  Set to false to get filenames (DEFAULT).
%
% OUTPUT:
%   filepaths      returned result: a cell array with matching names
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: List files or directories from a given path. And optionally uses regular expressions to filter the result
%              with option to set a minimal requirement on the number of results.
% EXAMPLE: 
%   filepaths = xASL_adm_GetFileList('D:\data','^\d$')                    get files with fullpath
%   filepaths = xASL_adm_GetFileList('D:\data','^.+\.nii$',false)         get filenames only
%   filepaths = xASL_adm_GetFileList('D:\data','^\d{5}_.+$',true,[1 Inf]) get files with fullpath and require at least one file
%   filepaths = xASL_adm_GetFileList('D:\data','^\d{5}_.+$',true,[0 10])  ..., require at maximum ten files
%   filepaths = xASL_adm_GetFileList('D:\data','^\d{5}_.+$',true,[5])     ..., require exactly 5 files
%   filepaths = xASL_adm_GetFileList('D:\data','^\d{5}_.+$',[],[],true)   get any number of directories with full path
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2019 ExploreASL

    % Checks the input arguments
	if nargin<2 || isempty(strRegEx)
		strRegEx = '^.+$'; % list all files by default
	end
	
    if nargin>=3 && islogical(mode)
		% If mode is true, the return the path+name
		if mode
			mode='FPList';
		else
			% Otherwise only name
			mode='List';
		end
    elseif nargin<3 || isempty(mode)
        mode='FPList'; % default mode is to return full paths
    end

    if nargin<4 || isempty(nRequired)
        nRequired=[0 Inf]; % default mode is to accept any file count
    elseif nargin>=4 && isnumeric(nRequired)
        n=numel(nRequired);
        if n<1 || n>2
            error('xASL_adm_GetFileList:invalidArgument','range should be a vector with one or two values')
        end
        if min(nRequired)<0
            error('xASL_adm_GetFileList:invalidArgument','range should only have values >=0')
        end
    else
        error('xASL_adm_GetFileList:invalidArgument','range should be a numerical vector with one or two values')
    end
    
    % Manage searching for \.nii (including \.gz)
    [StartInd, EndInd] = regexp(strRegEx, '\\.nii(?!.*\.gz)');
    if ~isempty(StartInd) && ~isempty(EndInd)
        strRegEx = [strRegEx(1:StartInd-1) '(\.nii|\.nii\.gz)' strRegEx(EndInd+1:end)];
    end
    
    if ispc
        % spm_select will malfunction if you use / on windows ...
        strDirectory = strrep(strDirectory,'/',filesep); 
    end
    
    if  nargin>=5 && bGetDirNames
        filepaths = spm_select(mode, strDirectory, 'dir', strRegEx);
    else
        filepaths = spm_select(mode, strDirectory, strRegEx);
    end        
    
    % return the list as cell array of strings because SPM cfg_files structure requires this.
    if ~isempty(filepaths)
        filepaths = cellstr(filepaths);
    end
    
    if isempty(nRequired)
        % check if any file was found
        if isempty(filepaths)
            warning('xASL_adm_GetFileList:nofilefound','Didn''t find any files: [%s] in [%s]',strRegEx,strDirectory)
        end
    else
        n=length(filepaths);
        mincount=min(nRequired);
        maxcount=max(nRequired);
        if n<mincount
            error('xASL_adm_GetFileList:missingfiles','Expected at least %i files, but found %i\nUsing pattern [%s] in [%s]',mincount,n,strRegEx,strDirectory);
        elseif n>maxcount
            error('xASL_adm_GetFileList:toomanyfiles','Expected at most %i files, but found %i\nUsing pattern [%s] in [%s]',maxcount,n,strRegEx,strDirectory);
        end
    end
end
