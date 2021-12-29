function filepaths = xASL_adm_ZipFileList(strDirectory, strRegExp, bRecurse, bUseGzip, nRequired, bDelete)
% Zip files in the specified directory that match the regular expression.
%
% FORMAT: filepaths = xASL_adm_ZipFileList(strDirectory, strRegExp[, bRecurse, bUseGzip, nRequired])
%                     xASL_adm_ZipFileList(strDirectory, strRegExp[, bRecurse, bUseGzip, nRequired])
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
%   bUseGzip       When TRUE then uses GZIP (DEFAULT)
%                  When FALSE uses ZIP command
%   nRequired      Set to a value or range [a b] to define a minimum/maximum number of required files. 
%                  An exception will be thrown if this requirement is not met. Use [n Inf] to accept at least n.
%                  Default is to accept any count. You can also provide a scalar [a] requiring to accept exactly 'a' matches
% OUTPUT:
%   filepathsOut - list of the zipped files
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Zip the files that match regular expression STRREGEXP in the given directory STRDIRECTORY.
%              Zips recursively if specified in BRECURSE. Zips all files unless the number is specified
%              by NREQUIRED, if the number is not met, then does not zip anything and throws an error.
% EXAMPLE: filepaths = xASL_adm_ZipFileList('c:\User\path', '^\d$'); Non-recursively zips all matching files
%          filepaths = xASL_adm_ZipFileList('c:\User\path', '^\d$',true,[],[0 Inf]); Recursively zips all matching files using GZIP
%                      xASL_adm_ZipFileList('c:\User\path', '^\d$',true,false,[5 10]); Recursively zips matching files only if between 5 and 10 matching files are found using ZIP
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2019 ExploreASL

    % Check input arguments
    if nargin<3 || isempty(bRecurse)
        bRecurse = false;
    end
	if nargin<4 || isempty(bUseGzip)
		bUseGzip = true;
	end
	if nargin<5 || isempty(nRequired)
		nRequired=[0 Inf]; % default mode is to accept any file count
    end
    if nargin<6 || isempty(bDelete)
        bDelete = false; % by default don't delete the zipped source file(s)
    end

	% Set recursive or non-recursive delete
    if bRecurse
        mode='FPListRec'; 
    else
        mode='FPList'; 
    end
    
    % Load the list of files using the xASL_adm_GetFileList function
    filepaths = xASL_adm_GetFileList(strDirectory, strRegExp, mode, nRequired);
	
	% Zips the files from the list
    if bUseGzip
        for iFile=1:length(filepaths)
            filepathsOut(iFile) = gzip(filepaths{iFile});
            if bDelete
                delete(filepaths{iFile});
            end
        end
    else
        for iFile=1:length(filepaths)
            zipfile = [filepaths{iFile} '.zip'];
            zip(zipfile,filepaths{iFile},strDirectory); 
            if bDelete
                delete(filepaths{iFile});
            end
            filepathsOut{iFile} = zipfile;
        end
    end
end
