function filepaths = xASL_adm_DeleteFileList(strDirectory, strRegEx, bRecurse, nRequired)
% Deletes files in the specified directory that match the regular expression.
%
% FORMAT: filepaths = xASL_adm_DeleteFileList(strDirectory, strRegEx[, bRecurse, nRequired])
%                     xASL_adm_DeleteFileList(strDirectory, strRegEx[, bRecurse, nRequired])
%
% INPUT:
%   strDirectory   Directory to search
%                  Optionally use wildcards to be more specific (i.e. 'C:\data\*.img')
%   strRegEx       Regular expression that will be used to filter the returned list.
%                  Examples: 
%                            '^.+$'       - matches all files (DEFAULT)
%                            '^\d+$'      - matches names that contains only digits
%                            '^.+\.nii$'  - matches names that ends with .nii
%                            '^\d{3}$'    - matches names that exist of exactly 3 digits
%                            '^\d{3}_.+$' - matches names that start with 3 digits, followed by an underscore and other characters
%                            '^PP\d+$'    - matches names that start with PP, followed by a number
%   bRecurse       When FALSE then deletes file normaly (DEFAULT)
%                  When TRUE  then deletes file recursively
%   nRequired      Set to a value or range [a b] to define a minimum/maximum number of required files. 
%                  An exception will be thrown if this requirement is not met. Use [n Inf] to accept at least n.
%                  Default is to accept any count. You can also provide a scalar [a] requiring to accept exactly 'a' matches
% OUTPUT:
%   filepaths - list of the deleted files
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Delete the files that match regular expression STRREGEXP in the given directory STRDIRECTORY.
%              Deletes recursively if specified in BRECURSE. Deletes all files unless the number is specified
%              by NREQUIRED, if the number is not met, then does not delete anything and throws an error.
% EXAMPLE: filepaths = xASL_adm_DeleteFileList('c:\User\path', '^\d$'); Non-recursively deletes all matching files
%          filepaths = xASL_adm_DeleteFileList('c:\User\path', '^\d$',true,[0 Inf]); Recursively deletes all matching files
%                      xASL_adm_DeleteFileList('c:\User\path', '^\d$',true,[5 10]); Recursively deletes matching files only if between 5 and 10 matching files are found
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2019 ExploreASL

    % Check the input arguments
	if nargin<2
		error('xASL_adm_DeleteFileList:Requires at least two input arguments.');
	end

	% Non recursive by default
    mode='FPList';
	
	% If the BRECURSE is set, then select recursive listing -> leading to recursive delete
	if nargin>=3 && ~isempty(bRecurse)
        if bRecurse
            mode='FPListRec';
        end
	end

	% allow zero or more files to be deleted by default
    if nargin<4 || isempty(nRequired)
        nRequired = [0 Inf]; 
    end
    
    % return the list as cell array of strings because SPM cfg_files structure requires this.
    filepaths = xASL_adm_GetFileList(strDirectory, strRegEx, mode, nRequired, false);
	
    % check if any file was found and then delete them
    if ~isempty(filepaths)
        % first disable warnings for double files (in case of symbolic links)
        warning('off', 'MATLAB:DELETE:FileNotFound');
        
        try
            delete(filepaths{:});
            filepaths = xASL_adm_GetFileList(strDirectory, strRegEx, mode, nRequired, false);
        catch ME
            % try deleting residual paths file-by-file
            % this can help in cases of symbolic links
            filepaths = xASL_adm_GetFileList(strDirectory, strRegEx, mode, nRequired, false);
            for iPath=1:numel(filepaths)
                xASL_delete(filepaths{iPath});
            end
            % then check if really everything was deleted, otherwise throw the
            % warning
            filepaths = xASL_adm_GetFileList(strDirectory, strRegEx, mode, nRequired, false);
            if ~isempty(filepaths)
                fprintf('%s\n', ME);
            end
        end
        
        if ~isempty(filepaths)
            warning('Something went wrong deleting files');
        end
        
        warning('on', 'MATLAB:DELETE:FileNotFound');
    end

end