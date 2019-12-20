function unpackedFiles = xASL_adm_UnzipOrCopy(srcDir, wildCard, destDir, bOverwrite)
% Unpacks (or copy if unpacked) one or more files matching the regular expresion
%
% FORMAT: unpackedFiles = xASL_adm_UnzipOrCopy(srcDir, wildCard, destDir [, bOverwrite])
%
% INPUT:
%   srcDir 	      - source path
%   wildCard      - wild card to identify the files
%   destDir       - destination directory
%   bOverwrite    - (TRUE/FALSE) overwrite on destination (default FALSE)
% OUTPUT:
%   unpackedFiles - the list of the unpacked files
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This is a simple wrapper function to (g)unzip one or more files to the specified destination
%              directory. Existing files or directories will not be overwritten, unless forced with bOverwrite.
%              A regular file-copy will be used if the source files don't have gz or zip filename extensions.

% EXAMPLE: unpackedFiles = xASL_adm_UnzipOrCopy('c:\User\path\', '*.nii','c:\User\path2\',0);
%          unpackedFiles = xASL_adm_UnzipOrCopy('c:\User\path\', 'file*.*','c:\User\path2\',0);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2019 ExploreASL

    % Admin
	if nargin<3 || isempty(srcDir) || isempty(destDir)
		error('xASL_adm_UnzipOrCopy:missingArgs','Need at least 3 arguments');
	end
	
	if nargin<4 || isempty(bOverwrite)
		bOverwrite = false;
	end
	
	% Search for all files
    unpackedFiles = [];
    entries = dir(fullfile(srcDir,wildCard));
    for ii=1:length(entries)
		if ~(entries(ii).isdir)
			[ ~, name, ext ] = fileparts(entries(ii).name);
			srcpath = fullfile(srcDir, entries(ii).name);
			X = {};
			if strcmpi(ext,'.gz')
				destpath = fullfile(destDir, name );
				if bOverwrite || ~exist(destpath,'file')
					X = gunzip(srcpath, destDir);
				end
			elseif strcmpi(ext,'.zip')
				destpath = fullfile(destDir, name );
				if bOverwrite || ~exist(destpath,'file')
					X = unzip(srcpath, destDir);
				end
			else
				destpath = fullfile(destDir, entries(ii).name);
				if bOverwrite || ~exist(destpath,'file')
					copyfile(srcpath, destDir);
					X = destpath;
				end
			end
			if ~isempty(X)
				if iscell(X)
					assert(length(X)==1,'xASL_adm_UnzipOrCopy: Expected exactly one file');
					X = X{1};
				end
				unpackedFiles{end+1} = X;
			end
		end
    end
end
