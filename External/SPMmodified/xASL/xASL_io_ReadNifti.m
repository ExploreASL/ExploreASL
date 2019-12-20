function [NiftiObject, pathOut] = xASL_io_ReadNifti(pathIn)
% Wrapper around SPM nifti function to read Nifti file.
%
% FORMAT: [NiftiObject, pathIn] = xASL_io_ReadNifti(pathIn)
%
% INPUT:
%   pathIn - The path to the image. String or a single cell with a string. (REQUIRED)
%
% OUTPUT:
%   NiftiObject - Nifti structure of the loaded file.
%   pathOut   - Modified path to the image in the form that is valid for the current (possibly unzipped) file.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Read Nifti file given by the path. Return the NII object. And also return the actual path to the loaded 
%              Nifti if by any reason the name changed during the function runtime (e.g. unzipping).
%
% EXAMPLE: 
%     [NiftiObject, pathOut] = xASL_io_ReadNifti('/home/tmp/CBF.nii')
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright ?? 2015-2019 ExploreASL

% Admin
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
    
[SrcPath, SrcFile] = xASL_fileparts(pathIn);
src_fileNII = fullfile(SrcPath, [SrcFile '.nii']);

%% Checking whether .nii & .nii.gz are equal (if both exist)
xASL_adm_UnzipNifti(src_fileNII);
    
pathIn = src_fileNII;
    
xASL_adm_ManageMoCoMat(pathIn); % manage .mat orientation sidecars

%% Load NIfTI
% Check first if NIfTI is valid
try
    NiftiObject = nifti(pathIn); % load NIfTI
catch ME
    try % try unzipping
        gunzip(pathIn);

        [Fpath, Ffile] = fileparts(pathIn);
        TempUnzippedFile = fullfile(Fpath, Ffile);
        warning('NIfTI seemed to be zipped, have unzipped this');
        try % if works, rename file
            nifti(TempUnzippedFile);
            xASL_SysMove(TempUnzippedFile, pathIn, true);                
        catch % otherwise throw error
            warning(['Couldnt read ' pathIn]);
            if exist(TempUnzippedFile,'file') && exist(pathIn,'file')
                delete(TempUnzippedFile);
            end
            error(ME);
        end
    end
end

pathOut = [pathIn ImExt]; % put it back
    
end
