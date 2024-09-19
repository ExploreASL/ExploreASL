function xASL_im_Move3DVolume(pathNifti, pathResult, dim, volNum, numVoxel)
%xASL_im_MoveImage Move zero padded image in spatial direction.
%
% FORMAT:       xASL_im_Move3DVolume(pathNifti, pathResult, dim, volNum, numVoxel)
% 
% INPUT:        pathNifti (CHAR ARRAY, REQUIRED)
%               pathResult (CHAR ARRAY, REQUIRED)
%               dim (CHAR, REQUIRED)
%               volNum (INTEGER, REQUIRED)
%               numVoxel (INTEGER, REQUIRED)
%
% OUTPUT:       NewList  - vertical list of things
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Move zero padded image in spatial direction. Was written to help with motion correction testing.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      n/a
% __________________________________
% Copyright 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


    %% Checks
    if nargin<5
        error('Missing input arguments...');
    end

    %% Load image
    imageNifti = xASL_io_Nifti2Im(pathNifti);
    originalVolume = imageNifti(:,:,:,volNum);

    %% Create new volume
    newVolume = zeros(size(originalVolume));

    % Get image dimensions
    dimensions = size(originalVolume);
    dimX = dimensions(1);
    dimY = dimensions(2);
    dimZ = dimensions(3);

    %% Move 3D volume
    switch dim
        case 'x'
            newVolume(1+numVoxel:dimX,:,:) = originalVolume(1:dimX-numVoxel,:,:);
        case 'y'
            newVolume(:,1+numVoxel:dimY,:) = originalVolume(:,1:dimY-numVoxel,:);
        case 'z'
            newVolume(:,:,1+numVoxel:dimZ) = originalVolume(:,:,1:dimZ-numVoxel);
        otherwise
            error('Unknown dimension...');
    end

    %% Store volume
    imageNifti(:,:,:,volNum) = newVolume;
    xASL_io_SaveNifti(pathNifti, pathResult, imageNifti);
    

end

