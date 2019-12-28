function xASL_im_EqualizeDimsOrder(PathSrc, PathRef)
%xASL_im_EqualizeDimsOrder Sets the order of dims of SrcPath to RefPath
%
% FORMAT: xASL_im_EqualizeDimsOrder(PathSrc, PathRef)
%
% INPUT:
% PathSrc  - path to source NIfTI file (REQUIRED)
% PathRef  - path to reference NIfTI file (REQUIRED)
%
% OUTPUT: n/a
%         But WARNING: SrcPath will be overwritten
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function ensures that the orientation dimension order
% is the same for SrcPath as for RefPath. It takes the order of RefPath and
% applies this to SrcPath. This is useful when wanting to resample one
% scantype (e.g. T1w segmentation) to the space of another scantype (e.g.
% ASL), without having to apply any realignment (which would be the case with SPMs reslice function).
% A good example is the MPRAGE T1w that is usually
% acquired sagittal, whereas the ASL image is acquired transversal.
%
% WARNING: this code assumes that the voxelsize is larger than rotations/shearing
%
% This function performs the following steps:
%
% 1) Obtain X Y Z order of reference NIfTI orientation matrix
% 2) Obtain X Y Z order of source NIfTI orientation matrix
% 3) Throw warning when the voxelsize is not significantly larger than other shearing/rotations
% 4) Obtain relative X Y Z orientation matrix order
% 5) Apply relative orientation matrix order
% 6) Apply shiftdim of image matrix and concomitant rotations
% 7) Save the new NIfTI
%
% EXAMPLE: xASL_im_EqualizeDimsOrder('/MyStudy/MySubject/ASL_1/sPVgm.nii', '/MyStudy/MySubject/ASL_1/ASL4D.nii');
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


niiSrc = xASL_io_ReadNifti(PathSrc);
niiRef = xASL_io_ReadNifti(PathRef);

%% 1) Obtain X Y Z order of reference NIfTI orientation matrix
RefMat = zeros(3,3);
for iDim=1:3
    RefMat(iDim,1:3) = abs(niiRef.mat(iDim,1:3)./norm(niiRef.mat(iDim,1:3)));
    RefInd(iDim) = find(RefMat(iDim,:)==max(RefMat(iDim,:)));
end

%% 2) Obtain X Y Z order of source NIfTI orientation matrix
SrcMat = zeros(3,3);
for iDim=1:3
    SrcMat(iDim,1:3) = abs(niiSrc.mat(iDim,1:3)./norm(niiSrc.mat(iDim,1:3)));
    SrcInd(iDim) = find(SrcMat(iDim,:)==max(SrcMat(iDim,:)));
end

%% 3) Throw warning when the voxelsize is not significantly larger than other shearing/rotations
for iDim=1:3
    VoxelSizeIs = RefMat(iDim,RefInd(iDim));
    RestIs = RefMat(iDim,find(RefMat(iDim,:)~=VoxelSizeIs));
    if VoxelSizeIs<=abs(1.5*sum(RestIs))
        warning('Finding dim order may have gone wrong, voxelsize is not significantly larger than shearing/rotations');
    end
end

%% 4) Obtain relative X Y Z orientation matrix order
niiSrc = xASL_io_ReadNifti(PathSrc);
SrcMat = niiSrc.mat;
for iDim=1:3
    MoveInd(iDim) = find(SrcInd==RefInd(iDim));
end
ShiftDimN=find(MoveInd==1)-1;

%% 5) Apply relative orientation matrix order
SrcMat(1:3,:) = SrcMat(MoveInd,:);

%% 6) Apply shiftdim of image matrix and concomitant rotations
if ShiftDimN~=0
    fprintf(['Shifting NIfTI with ' num2str(ShiftDimN) ' dimensions, and rotating with ' num2str(90*ShiftDimN) ' degrees\n']);
    IM = shiftdim(niiSrc.dat(:,:,:), ShiftDimN);
    IM = xASL_im_rotate(IM, 90*ShiftDimN);

    %% 7) Save the new NIfTI
    clear niiSrc;
    xASL_io_SaveNifti(PathSrc, PathSrc, IM, [], false, SrcMat);
end

end