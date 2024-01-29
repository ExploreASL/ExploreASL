function [NiftiObject, pathOut, json] = xASL_io_ReadNifti(pathIn)
%xASL_io_ReadNifti Wrapper around SPM nifti function to read Nifti file
%
% FORMAT: [NiftiObject, pathOut, json] = xASL_io_ReadNifti(pathIn)
%
% INPUT:
%   pathIn - The path to the image. String or a single cell with a string. (REQUIRED)
%
% OUTPUT:
%   NiftiObject - Nifti structure of the loaded file.
%   pathOut     - Modified path to the image in the form that is valid for the current (possibly unzipped) file.
%   json        - Loaded JSON sidecar, in ExploreASL legacy format
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Read Nifti file given by the path. Return the NII object. And also return the actual path to the loaded 
%              Nifti if by any reason the name changed during the function runtime (e.g. unzipping).
%
% EXAMPLE:  
%   NiftiObject = xASL_io_ReadNifti('/home/tmp/CBF.nii'); % for only loading the nifti
%   [NiftiObject, ~, json] = xASL_io_ReadNifti('/home/tmp/CBF.nii'); % for also loading the json sidecar
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright ExploreASL 2015-2024 ExploreASL

% Admin
% Set defaults
NiftiObject = [];
pathOut = [];
json = [];

[fPath, fFile] = xASL_fileparts(pathIn);
pathUnzipped = fullfile(fPath, [fFile '.nii']);
pathZipped = fullfile(fPath, [fFile '.nii.gz']);
pathJson = fullfile(fPath, [fFile '.json']);

if iscell(pathIn)
	if length(pathIn)==1
		pathIn = pathIn{1};
	else
		error('NiftiPath is a cell with more than one element, not a character string.');
	end
end

% remove SPM image extension
IndexComma  = strfind(pathIn,',');
if ~isempty(IndexComma)
	ImExt       = pathIn(IndexComma:end);
	pathIn   = pathIn(1:IndexComma-1);
else
	ImExt       = '';
end
    
if ~xASL_exist(pathIn,'file')
	error([pathIn ' does not exist as .nii[.gz]']);
end

%% Checking whether .nii & .nii.gz are equal (if both exist)
xASL_adm_UnzipNifti(pathUnzipped);

pathIn = pathUnzipped;
    
xASL_adm_ManageMoCoMat(pathIn); % manage .mat orientation sidecars

%% Load NIfTI
% Check first if NIfTI is valid
try
    NiftiObject = nifti(pathIn); % load NIfTI
catch ME
    try % try unzipping
        gunzip(pathIn);

        warning('NIfTI seemed to be zipped, have unzipped this');
        try % if works, rename file
            nifti(pathUnzipped);
            xASL_SysMove(pathUnzipped, pathIn, true);                
        catch % otherwise throw error
            warning(['Couldnt read ' pathIn]);
            if exist(pathUnzipped, 'file') && exist(pathIn, 'file')
                delete(pathUnzipped);
            end
            error(ME);
        end
    end
end

pathOut = [pathIn ImExt]; % put it back
    
%% Load JSON sidecar if it exists

if exist(pathJson)
    json = xASL_bids_parms2BIDS([], xASL_io_ReadJson(pathJson), 0);
end


end