function [structOut] = xASL_qc_ComputeNiftiOrientation(PathNIfTI, structIn)
%xASL_qc_ComputeNiftiOrientation Adds basic volume information about a Nifti file to a structure
%
% FORMAT:       [structOut] = xASL_qc_ComputeNiftiOrientation(PathNIfTI[, structIn])
% 
% INPUT:        PathNIfTI - Path to the file (REQUIRED)
%               structIn  - Adds the calculated results to this structure and keeps all fields intact (OPTIONAL, DEFAULT = empty struct)
%
% OUTPUT:       The following fields are calculated and added to the output structure
%                   Matrix - A vector with the dimensions of the image 
%                   VoxelSize_mm - A vector with voxel size of the image
%                   RigidBody2Anat_mm - A Net Displacement Vector (RMS) (mm) between the original and new image orientation
%                   If PathNIfTI does exist, then 'n/a', 'n/a', and NaN are returned in the structure
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It loads the input Nifti, finds its dimension, voxel size and a net vector distance from its 
%              original position before registration. Adds all these information into an output structure structOut
%              while copying all from structIn and keeping it intact.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      structOut = xASL_qc_ComputeNiftiOrientation('test.nii')
% __________________________________
% Copyright 2015-2021 ExploreASL

% Initializes with an empty structure if missing
if nargin < 2 || isempty(structIn)
	structOut = struct();
else
	structOut = structIn;
end

    if ~xASL_exist(PathNIfTI,'file')
		% File doesn't exist - n/a values are returned
        warning(['Missing file: ' PathNIfTI]);
        structOut.Matrix = 'n/a';
        structOut.VoxelSize_mm = 'n/a';
        structOut.RigidBody2Anat_mm = NaN;
	else
		% Reads the file
        nii = xASL_io_ReadNifti(PathNIfTI);
		
		% Extracts the file dimensions
        structOut.Matrix = num2str(size(nii.dat));

		% Extracts the voxel size
        SPMmat = spm_imatrix(nii.mat);
        SPMmat0 = spm_imatrix(nii.mat0);
        VoxelSize_mm = abs(SPMmat(7:9));
        structOut.VoxelSize_mm = [num2str(VoxelSize_mm(1)) ' ' num2str(VoxelSize_mm(2)) ' ' num2str(VoxelSize_mm(3))]; % this avoids too much whitespace

        % Now calculate Net Displacement Vector (RMS) from ASL to T1w image (mm)
        RigidBody2anat = SPMmat(1:6) - SPMmat0(1:6);
        part_translation = sum(RigidBody2anat(1:3).^2);
        mean_radius = 50; % typical distance center head to cerebral cortex (Power et al., NeuroImage 2012)
        rx = RigidBody2anat(4); ry = RigidBody2anat(5); rz = RigidBody2anat(6);

        part_rotation = 0.2*mean_radius^2* ((cos(rx)-1).^2 + (sin(rx)).^2 + (cos(ry)-1).^2 + (sin(ry)).^2 + (cos(rz)-1).^2 + (sin(rz)).^2);
        structOut.RigidBody2Anat_mm = xASL_round(sqrt(part_translation + part_rotation),3);
    end
end
