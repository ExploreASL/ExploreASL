function [Struct] = xASL_qc_ComputeNiftiOrientation(x, PathNIfTI, Struct)
%xASL_qc_ComputeNiftiOrientation ...
%
% FORMAT:       [Struct] = xASL_qc_ComputeNiftiOrientation(x, PathNIfTI, Struct)
% 
% INPUT:        x         - ...
%               PathNIfTI - ...
%               Struct    - ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2021 ExploreASL

    if ~xASL_exist(PathNIfTI,'file')
        warning(['Missing file: ' PathNIfTI]);
        Struct.Matrix = 'n/a';
        Struct.VoxelSize_mm = 'n/a';
        Struct.RigidBody2Anat_mm = NaN;
    else
        nii = xASL_io_ReadNifti(PathNIfTI);
        Struct.Matrix = num2str(size(nii.dat));

        SPMmat = spm_imatrix(nii.mat);
        SPMmat0 = spm_imatrix(nii.mat0);
        VoxelSize_mm = abs(SPMmat(7:9));
        Struct.VoxelSize_mm = [num2str(VoxelSize_mm(1)) ' ' num2str(VoxelSize_mm(2)) ' ' num2str(VoxelSize_mm(3))]; % this avoids too much whitespace

        % Now calculate Net Displacement Vector (RMS) from ASL to T1w image (mm)
        RigidBody2anat = SPMmat(1:6) - SPMmat0(1:6);
        part_translation = sum(RigidBody2anat(1:3).^2);
        mean_radius = 50; % typical distance center head to cerebral cortex (Power et al., NeuroImage 2012)
        rx = RigidBody2anat(4); ry = RigidBody2anat(5); rz = RigidBody2anat(6);

        part_rotation = 0.2*mean_radius^2* ((cos(rx)-1).^2 + (sin(rx)).^2 + (cos(ry)-1).^2 + (sin(ry)).^2 + (cos(rz)-1).^2 + (sin(rz)).^2);
        Struct.RigidBody2Anat_mm = xASL_round(sqrt(part_translation + part_rotation),3);
    end

end

