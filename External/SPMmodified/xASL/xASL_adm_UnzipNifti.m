function pathOut = xASL_adm_UnzipNifti(pathIn, varargin)
% Takes the pathIn file, unzips it and deletes the original ZIP file.
%
% FORMAT: pathOut = xASL_adm_UnzipNifti(pathIn [,bOverwrite])
%
% INPUT:
%   pathIn     - path and filename to the file to be unzipped (.NII or .NII.GZ)
%   bOverwrite - indicates if we should overwrite the .NII file (OPTIONAL, default = FALSE)
% OUTPUT:
%   pathOut    - path to the unzipped file
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Takes the input file, unzips if needed, delete the zipped file and return the path to the unzipped file.
%              If the input is already unzipped, then does nothing, but returns the original filename - so it
%              can be run just to be sure a file is unzipped without much overhead.
%              Returns error if more than one file is in the archive, if the filename does not exist, is a directory etc.
%              If there's a NII and NII.GZ already existing, then return error, or just overwrite in case overwrite is set to 1
%
% EXAMPLE: pathOut = xASL_adm_UnzipNifti('test.nii',1)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright (C) 2015-2020 ExploreASL

% Check for the optional parameter overwrite
if nargin < 2
	bOverwrite = 0;
else
	bOverwrite = varargin{1};
end

if nargin > 2
	error('Maximum number of parameters is 2.');
end

% Check if the extension pathIn extension is .nii, .gz - otherwise gives an error
[pathstr, name0, ext0] = fileparts(pathIn);

switch (ext0)
	case '.nii'
		isGZ = 0;
	case '.gz'
		isGZ = 1;
	otherwise
		error(['Handles only .nii and .nii.gz files: ' pathIn]);
end

% Gets the correct paths for the NII and NII.GZ files
if isGZ
	pathNII = fullfile(pathstr,name0);
	nameNII = name0;
	pathGZ  = pathIn;
else
	pathNII = pathIn;
	nameNII = [name0 ext0];
	pathGZ  = [pathIn '.gz'];
end

% Checks that none of those is a directory
if exist(pathNII,'dir')
	error(['Nifti file is a directory: ' pathNII]);
end

if exist(pathGZ,'dir')
	error(['GZ file is a directory: ' pathGZ]);
end

% Checks for the existence of the .nii and .nii.gz files
existGZ = 0;
existNII = 0;
if exist(pathNII,'file')
	existNII = 1;
end

if exist(pathGZ,'file')
	existGZ = 1;
end

% If both .nii and .nii.gz exist, then report an error unless the overwrite option is set to 1
if existNII && existGZ
	% We try the last option that both unzipped and zipped are the same
	% We call the ZipFileNameHandling function with pathNII - if both are the same, then
	% it deletes the GZ
	xASL_adm_ZipFileNameHandling(pathNII);
	
	% Update the file existing booleans
	if ~exist(pathNII,'file')
		existNII = 0;
	end
	if ~exist(pathGZ,'file')
		existGZ = 0;
	end
	if existNII && existGZ
		if bOverwrite
			% Overwrite == 1, delete NII and continue unzipping GZ
			delete(pathNII);
			existNII = 0;			
		else
			% If both files still exist, then they were different and we report an error
			error(['Both .nii and .nii.gz exist and bOverwrite is set to 0: ' pathIn]);
		end
	end
	% Otherwise if the files were the same and either NII or NII.GZ was removed - to be addressed later
end

% Already unzipped, return the name of the .nii file and do nothing else
if existNII && ~existGZ
	pathOut = pathNII;
	return;
end

% None of the files exist
if ~existNII && ~existGZ
	warning(['Cannot locate .nii[.gz] files: ' pathIn ', skipping...']);
	return;
end

%% When it comes here - unzip the GZ file and delete it
bDelete = true;
try
    if ~ispc
        [result1, result2] = system(['gunzip -f ' pathGZ]); % use system unzipping, is faster & doesn't need JVM
        extracted{1} = pathNII;
        if result1~=0
            warning('Couldnt unzip NIfTI with CLI');
            fprintf('%s\n', result2);
        end
    else % pc, use matlab's JVM
        extracted = gunzip(pathGZ);
    end

catch ME % if unzipping didnt work, try reading it without unzipping
    try 
        nifti(pathGZ);
        % if this didnt crash, proceed
        xASL_SysMove(pathGZ, pathNII, true);
        warning('Couldnt unzip file, but seemed to be a normal NIfTI file, renamed instead');
        extracted{1} = pathNII;
        bDelete = false;
    catch
        warning('Couldnt unzip file, illegal file?');
        error('%s\n', ME.message);
    end
end


% A single file only is extracted
if length(extracted)==1 && bDelete
	pathOut = extracted{1};
	[~, name1, ext1] = fileparts(pathOut);
	% The file name has to match the expected .nii.gz file needs to be extracted as .nii file.
	% If not then delete it and throw an error
	if ~strcmp([name1 ext1],nameNII)
		delete(extracted{1});
		error(['Name after extraction does not match: ' pathIn]);
	end
	
	% If everything was successful, then delete the unpacked .nii.gz file
	if exist(pathGZ, 'file')
        delete(pathGZ);
    end
elseif bDelete
	% Multiple files are extracted - delete extracted files and report an error
	for ii = 1:length(extracted)
		delete(extracted{ii});
	end
	error(['Archive contained more than 1 file: ' pathIn]);
end


end