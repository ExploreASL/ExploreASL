function  RES = xASL_adm_GetFsList(strDirectory, strRegEx, bGetDirNames, bExcludeHidden, bIgnoreCase, nRequired)
% List the content of a directory, using regular expressions, and sorts the output.
%
% FORMAT: RES = xASL_adm_GetFsList([strDirectory, strRegEx, bGetDirNames, bExcludeHidden, bIgnoreCase, nRequired])
%
% INPUT:
%   strDirectory   Directory to search
%                  (use '' or [] for current dir)
%                  Optionally use wildcards to be more specific (i.e. 'C:\data\*.img')
%   strRegEx       Regular expression that will be used to filter the returned list.
%                  Note that the returned strings (RES) only contain the characters that match the expression.
%                  Examples: 
%                            '^.+$'       - matches all files (DEFAULT)
%                            '^\d+$'      - matches names that contains only digits
%                            '^.+\.nii$'  - matches names that ends with .nii
%                            '^\d{3}$'    - matches names that exist of exactly 3 digits
%                            '^\d{3}_.+$' - matches names that start with 3 digits, followed by an underscore and other characters
%                            '^PP\d+$'    - matches names that start with PP, followed by a number
%   bGetDirNames   Set to true to get directory names.
%                  Set to false to get filenames (DEFAULT).
%   bExcludeHidden Set to true to skip hidden objects. (DEFAULT)
%                  Set to false to include hidden objects.
%   bIgnoreCase    Set to true to perform case insensitive matches (DEFAULT on windows) 
%                  Set to false to perform case sensitive matches (DEFAULT on Linux and Mac)
%   nRequired      Set to a value or range [a b] to define a minimum/maximum number of required files. 
%                  An exception will be thrown if this requirement is not met. Use [n Inf] to accept at least n.
%                  Default is to accept any count. You can also provide a scalar [a] requiring to accept exactly 'a' matches
%
% OUTPUT:
%   RES            returned result: a cell array with matching names
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: List files or directories from a given path. And optionally uses regular expressions to filter the result
%              with options to exclude hidden files, ignore case, and set a minimal requirement on the number of results.
%              Sorts the results at the end.
% EXAMPLE: 
%   RES = xASL_adm_GetFsList()                                    all directories from the current directory, exclude hidden
%   RES = xASL_adm_GetFsList('D:\data','^\d$')                    only files, exclude hidden
%   RES = xASL_adm_GetFsList('D:\data','^.+\.nii$',false)         only files, exclude hidden
%   RES = xASL_adm_GetFsList([],'^.+\.nii$',false)                only files from the current directory, exclude hidden
%   RES = xASL_adm_GetFsList('D:\data','^\d{5}_.+$',true)         only directories, exclude hidden
%   RES = xASL_adm_GetFsList('D:\data','^\d{5}_.+$',true,false)   only directories, include hidden
%   RES = xASL_adm_GetFsList('D:\data','^\d{5}_.+$',true,false,false) ..., case sensitive
%   RES = xASL_adm_GetFsList('D:\data','^\d{5}_.+$',true,false,false,[1 Inf]) ..., require at least one file
%   RES = xASL_adm_GetFsList('D:\data','^\d{5}_.+$',true,false,false,[0 10]) ..., require at maximum ten files
%   RES = xASL_adm_GetFsList('D:\data','^\d{5}_.+$',true,false,false,[5]) ..., require exactly 5 files
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2020 ExploreASL

    % Check input arguments
    if nargin<1 || isempty(strDirectory)
        % use current directory if not specified
        strDirectory=pwd;
    end
    if nargin<2 || isempty(strRegEx)
        % look for all names by default
        strRegEx='^.+$';
    end
    if nargin<3 || isempty(bGetDirNames)
        % look for files by default
        bGetDirNames = false;
    end
    if nargin<4 || isempty(bExcludeHidden)
        % ignore hidden folders by default
        bExcludeHidden = true;
	end

    if nargin<5 || isempty(bIgnoreCase)
		% By default ignore case on Windows, and case sensitive on Linux and Mac
		bIgnoreCase = ispc();
	end
    
	if nargin<6 || isempty(nRequired)
		nRequired=[]; % default mode is to accept any file count
	elseif nargin>=6 && isnumeric(nRequired)
		n=numel(nRequired);
		if n<1 || n>2
			error('xASL_adm_GetFsList:invalidArgument','range should be a vector with one or two values')
		end
		if min(nRequired)<0
			error('xASL_adm_GetFsList:invalidArgument','range should only have values >=0')
		end
	else
		error('xASL_adm_GetFsList:invalidArgument','range should be a numerical vector with one or two values')
	end
	
	% Set the ignore case string
	if bIgnoreCase
		casearg = 'ignorecase';
	else
		casearg = 'matchcase';
	end
    
	if ispc
        % spm_select will malfunction if you use / on windows ...
        strDirectory = strrep(strDirectory,'/',filesep); 
    end
	
    %% get the contents of all files and directories that match D
    if isempty(regexp(strDirectory,'[\*]', 'once', casearg)) % ? is not supported as wildcard here
        if isdir(strDirectory)
            % D is a directory: 
            % just get all items, but make sure to add filesep to be able to use fileparts below
            strDirectory = fullfile(strDirectory,filesep);
        else
            % weird: not looking in a directory?
            warning('xASL_adm_GetFsList:Directory doesn''t exist: %s',strDirectory)
        end
    end
    E = dir(strDirectory); 
    strDirectory = fileparts(strDirectory); % only need folder name itself in script below
    
    %% collect all names from struct E in cell array N
    RES = {E.name};
    
    %% only keep directories OR filenames
    B = [E.isdir];
    if bGetDirNames
        % only keep the directories
        B = ~B;
    end
    
    %% remove hidden objects
    if bExcludeHidden
        % on all systems we check if files start with a dot
        M = regexp(RES, '^\.', 'match', 'once', casearg);  
        L = strcmp(M, '');
        B = B | ~L; % combine the logical vectors
        if ispc
            % on windows, there is also a specific file attribute to mark a file as hidden:
            for iFile=1:length(RES)
                if ~B(iFile)
                    [stat, attr]=fileattrib(fullfile(strDirectory,RES{iFile}));
                    if stat && attr.hidden
                        B(iFile) = 1;
                    end
                end
            end
        end
    end
    
    %% remove unwanted objects
    if sum(B)>0
        RES(B) = []; % only keep files OR directories
    end
    
    % get only those that match the regular expression
    RES = regexp(RES, strRegEx, 'match', 'once', casearg);  
    
    %% remove (non matching) empty cellï¿½s from result
    RES(strcmp(RES, ''))= []; 
    
    %% sort the result
    RES = sort(RES); 
    
    %% check count
    if isempty(nRequired)
        % check if any file was found
        if isempty(RES)
            warning('xASL_adm_GetFsList::nofilefound','Didn''t find any files: [%s] in [%s]',strRegEx,strDirectory)
        end
    else
        n=length(RES);
        mincount=min(nRequired);
        maxcount=max(nRequired);
        if n<mincount
            error('xASL_adm_GetFsList::missingfiles','Expected at least %i files, but found %i\nUsing pattern [%s] in [%s]',mincount,n,strRegEx,strDirectory);
        elseif n>maxcount
            error('xASL_adm_GetFsList::toomanyfiles','Expected at most %i files, but found %i\nUsing pattern [%s] in [%s]',maxcount,n,strRegEx,strDirectory);
        end
    end
end
