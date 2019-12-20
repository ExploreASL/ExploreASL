function [result, files] = xASL_adm_CheckFileCount(path, expr, mincount, failifmissing)
% Checks if the number of files matching the expression is more than mincount.
%
% FORMAT: 
%         [result, files] = xASL_adm_CheckFileCount(path, expr[, mincount, failifmissing])
%         [result]        = xASL_adm_CheckFileCount(...)
%
% INPUT:
%   path 	     - path to be checked
%   expr         - is an spm_select compatible regular expression.
%   mincount     - minimum number of files expected (default 1)
%   filifmissing - if true, then throw error when the number of files is less than mincount (default TRUE)
% OUTPUT:
%   result       - returns TRUE if the number of files is higher or equal to mincount
%   files        - xASL_adm_GetFileList structure of the files corresponding to the expression
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Checks the given PATH for files corresponding to the SPM_SELECT regular expression EXPR. Returns if the number of files is equal to or higher
%              than MINCOUNT. If FAILIFMISSING is true and not enough files, then throw and error. If everything goes ok and second output argument is specified
%              then return also the list of files.
% EXAMPLE: [result, files] = xASL_adm_CheckFileCount('/user/test','*.nii',2,0);
%          [result]        = xASL_adm_CheckFileCount('/user/test','*.nii',3);
%          [result, files] = xASL_adm_CheckFileCount('/user/test','*.nii');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL

    if nargin<2
		error('xASL_adm_CheckFileCount:Need at least 2 input arguments');
	end
	
	% Set the defaults
    if nargin<3 || isempty(mincount)
        mincount = 1;
	end
	
    if nargin<4 || isempty(failifmissing)
        failifmissing = true;
	end
    
	% Load files and check if there is enough of them
    filepaths = xASL_adm_GetFileList( path, expr, 'List' );
    n = length(filepaths);
    bOK = n>=mincount ;
    
	% Report error if not enough of files
    if ~failifmissing
        result = bOK;
    elseif ~bOK
        error('xASL_adm_CheckFileCount:missingfiles','Expected %d files, but found %d',mincount,n);
	end

	% Output the files if the second output argument is given
    if nargout>=2
        files = filepaths;
    end
end
