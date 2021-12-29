function filepaths = xASL_adm_DeleteFilePair(path, varargin)
% Deletes a file and corresponding sibblings with same basename and specified extension.
%
% FORMAT: filepaths = xASL_adm_DeleteFilePair(path, ext1[, ext2 [, ext3 ...]])
%                     xASL_adm_DeleteFilePair(path, ext1[, ext2 [, ext3 ...]])
%
% INPUT:
%   path  	  - path and filename to the file to be deleted
%   ext1      - an extension of a sibling file to be deleted. Either given as '.ext' or 'ext'
%   ext2,...  - other extensions for deletion
% OUTPUT:
%   filepaths - list of the deleted files
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Delete the file given in PATH, and also deletes files with the same name, but with extension
%              given in EXT1, and potentially also EXT2, EXT3... 
% EXAMPLE: filepaths = xASL_adm_DeleteFilePair('c:\User\path\file.nii', 'mat');
%          filepaths = xASL_adm_DeleteFilePair('c:\User\path\file.nii', 'mat','jpg');
%                      xASL_adm_DeleteFilePair('c:\User\path\file.nii', 'mat');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2019 ExploreASL

    % The list of deleted paths is empty
    filepaths = {};
	
	% Deletes the main file
	if xASL_exist(path,'file')
		xASL_delete(path);
		filepaths{1} = path;
	end
       
	% Extracts the filename and extension from the main file
	[thepath, base, ~] = xASL_fileparts(path);
	
	% Goes through all the given extensions
	for iExt=1:length(varargin)
		ext = varargin{iExt};
		% Adds fullstop before the extension if it is not there
		if ext(1)~='.'
			ext = ['.' ext];
		end
		
		% If the file exists, then delete it and add to the list of deleted files
		path = fullfile(thepath, [ base, ext ]);
		if xASL_exist(path,'file')
			xASL_delete(path);
			filepaths{end+1} = path;
		end
	end
end
