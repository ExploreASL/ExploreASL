function pathOut = xASL_adm_GzipNifti(pathIn, bOverwrite)
%xASL_adm_GzipNifti Take the pathIn file, gzip it and delete the original Nifti file
%
% FORMAT: pathOut = xASL_adm_GzipNifti(pathIn [,bOverwrite])
%
% INPUT:
%   pathIn     - path and filename to the file to be zipped (.NII or .NII.GZ)
%   bOverwrite - indicates if we should overwrite the .NII file (OPTIONAL, default = FALSE)
% OUTPUT:
%   pathOut    - path to the zipped file
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Take the input file, zips it, overwriting any existing zipped file and return the path of the zipped file.
%
% EXAMPLE: pathOut = xASL_adm_GzipNifti('test.nii',1)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright (C) 2015-2020 ExploreASL


%% ------------------------------------------------------
%% 0. Admin

% Check for the optional parameter overwrite
if nargin < 2 || isempty(bOverwrite)
    bOverwrite = 1;
end

pathOut = '';

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

% Get the correct paths for the NII and NII.GZ files
if isGZ
	pathNII = fullfile(pathstr,name0);
	pathGZ  = pathIn;
else
	pathNII = pathIn;
	pathGZ  = [pathIn '.gz'];
end

% Checks that none of those is a directory
if exist(pathNII, 'dir')
	error(['Nifti file is a directory: ' pathNII]);
end

if exist(pathGZ, 'dir')
	error(['GZ file is a directory: ' pathGZ]);
end


%% 1. Deal with overwrite
if ~exist(pathNII, 'file') && exist(pathGZ, 'file')
    xASL_adm_UnzipNifti(pathGZ);
end

pathNII = xASL_adm_ZipFileNameHandling(pathNII);
if exist(pathNII, 'file') && exist(pathGZ, 'file')
    if bOverwrite
        delete(pathGZ);
    else
        error('pathGZ already existed, not overwriting');
    end
end


%% ------------------------------------------------------
%% 2. Zip the NIfTI
if ~ispc
    [result1, result2] = system(['gzip -1 -f ' pathNII]); % use system unzipping, is faster & doesn't need JVM

    if result1~=0
        warning('Couldnt gzip NIfTI with CLI');
        fprintf('%s\n', result2);
    else
        pathOut = pathGZ;
    end
else % pc, use matlab's JVM
    gzip(pathNII);
    delete(pathNII);
    pathOut = pathGZ;
end


end
