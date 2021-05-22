function xASL_im_FlipNifti(pathInput, flipAxis, bOverwrite)
%xASL_im_FlipNifti Flip the NifTI image matrix
%
% FORMAT: xASL_im_FlipNifti(pathInput[, flipAxis, bOverwrite])
% 
% INPUT:
%   pathInput   - path to NifTI file
%   flipAxis    - dimension of the image matrix that will be flipped, or axis
%                 that will be flipped over, can be:
%                 1: x-axis, left-right in MNI
%                 2: y-axis, anterior-posterior in MNI
%                 3: z-axis, inferior-posterior in MNI
%                 (OPTIONAL, DEFAULT = 1)
%   bOverwrite  - boolean specifying if an existing output NIfTI
%                 
%
% OUTPUT: n/a
% OUTPUT FILES  - input path but with 'flipped' suffix
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function allows correcting an inappropriate flip in the image matrix. It
%               will not change the orientation matrix in the header but the image
%               itself. So any NifTI program will not be aware of this flip!
%
%               This function runs the following steps:
%               1. Manage if we overwrite the new NIfTI
%               2. Manage if we zip the new NIfTI
%               3. Load image from NIfTI
%               4. Flip image
%               5. Save image to NIfTI
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_im_FlipNifti('CICERO_Nolan/analysis/C3T-S007_1/RPI_4V/4V.nii');
% __________________________________
% Copyright 2021 ExploreASL


%% Admin

if nargin<3 || isempty(bOverwrite)
    bOverwrite = false;
end
if nargin<2 || isempty(flipAxis)
    flipAxis = 1; % x-axis
end
if nargin<1 || isempty(pathInput)
    error('Input path missing');
end

[Fpath, Ffile, Fext] = xASL_fileparts(pathInput);
pathOutput = fullfile(Fpath, [Ffile '_flipped' Fext]);

%% 1. Manage if we overwrite the new NIfTI
if xASL_exist(pathOutput, 'file') && bOverwrite
    xASL_delete(pathOutput);
elseif xASL_exist(pathOutput, 'file') && ~bOverwrite
    warning('Output file already existed, skipping:');
    fprintf('%s\n', pathOutput);
    fprintf('%s\n', 'Hint, either delete this file manually or set bOverwrite to true');
end
    

%% 2. Manage if we zip the new NIfTI
if strcmp(Fext, '.nii')
    bGZip = 0;
elseif strcmp(Fext, '.nii.gz')
    bGZip = 1;
else
    warning('Unknown extension!');
end

%% 3. Load image from NIfTI
imageIn = xASL_io_Nifti2Im(pathInput);

%% 4. Flip image
imageOut = xASL_im_Flip(imageIn, flipAxis);

%% 5. Save image to NIfTI
xASL_io_SaveNifti(pathInput, pathOutput, imageOut, [], bGZip);


end