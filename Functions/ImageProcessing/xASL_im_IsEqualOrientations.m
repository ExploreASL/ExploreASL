function [isEqualdim, isEqualMat0, isEqualMat] = xASL_im_IsEqualOrientations(pathNifti1, pathNifti2)
%xASL_im_IsEqualOrientations Check if NIfTI images have equal orientations
%
% FORMAT: [isEqualdim, isEqualMat0, isEqualMat] = xASL_im_IsEqualOrientations(pathNifti1, pathNifti2)
% 
% INPUT:
%
%   Definitions:
%
%   pathNifti1 - Path to the first NIfTI (REQUIRED, STRING)
%   pathNifti2 - Path to the second NIfTI (REQUIRED, STRING)
%
% OUTPUT:
%
%   isEqualdim  - Boolean for equal image matrix dimensions (true = equal, false = not equal)
%   isEqualMat0 - Boolean for equal orientation mat0 (the original image orientation, i.e. before realignments)
%   isEqualMat  - Boolean for equal orientation mat (the current image orientation, i.e. after realignments)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:          This function checks is NIfTI images have equal orientations, which can be useful e.g., when derivatives and original
%                       images need to stay aligned (e.g., FLAIR and WMH_SEGM)
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE: [isEqualdim, isEqualMat0, isEqualMat] = xASL_im_IsEqualOrientations('pathToMyNifti', 'pathToMyOtherNifti');
%
% __________________________________
% Copyright (c) 2015-2025 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


%% 1. Load NIfTIs
objectNifti1 = xASL_io_ReadNifti(pathNifti1);
objectNifti2 = xASL_io_ReadNifti(pathNifti2);


%% 2. Compare image dimensions
isEqualdim = isequal(objectNifti1.hdr.dim, objectNifti2.hdr.dim);


%% 3. Compare orientations
isEqualMat0 = isequal(objectNifti1.mat0, objectNifti2.mat0);
isEqualMat = isequal(objectNifti1.mat, objectNifti2.mat);


end