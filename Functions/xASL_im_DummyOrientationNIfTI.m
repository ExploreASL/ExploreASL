function xASL_im_DummyOrientationNIfTI(PathSrc, PathRef, PathDummyOut, bApplyRotationMinor, bApplyRotation90, bApplyZoom, bApplyTranslation)
%xASL_im_DummyOrientationNIfTI Create a dummy NIfTI with desired
%orientation matrix
% FORMAT: xASL_im_DummyOrientationNIfTI(PathSrc, PathRef, PathDummyOut[, bApplyRotationMinor, bApplyRotation90, bApplyZoom, bApplyTranslation])
%
% INPUT:    PathSrc     - path to source NIfTI that we want to transform
%           PathRef     - path to reference NIfTI containing orientation matrix that we want to (partly) transform to
%           PathDummyOut- path to dummy NIfTI with desired new orientation matrix
%                         applied but without an image matrix (empty). This dummy NIfTI
%                         can be used as reference for xASL_spm_reslice to put (the source) NIfTI in this
%                         desired orientation
%           bApplyRotationMinor - boolean if we want to apply subtle
%                                 rotations/shearing (OPTIONAL, DEFAULT=true)
%           bApplyRotation90    - boolean if we want to apply 90 degree
%                                 rotations or dim order shifts (OPTIONAL, DEFAULT=true)
%           bApplyZoom          - boolean if we want to zoom/stretch (OPTIONAL, DEFAULT=true)
%           bApplyTranslation   - boolean if we want to apply the translation (OPTIONAL, DEFAULT=true)
%   
% OUTPUT: n/a
% OUTPUT NIfTI: PathDummyOut - NIfTI as above
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function creates a dummy image as reference for xASL_spm_reslice,
%               allowing to only apply specific parts of the transformation between the
%               two images. E.g. only the rotation, or only the zooming.
%               This can be useful to correct for any erroneous rotations from registration,
%               or to put two images in the same space without applying their
%               realignment. This function performs the following steps:
%               1) Load orientations & calculate transformation
%               2) Calculate the desired transformation
%               3) Calculate new orientation matrix
%               4) Calculate the new image size
%               5) Save the dummy NIfTI
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% THEORY:   Any rotation will always swap the other dims (X rotation swaps Y/Z, Y
%           rotation swaps X/Z etc.) because they are perpendicular (haaks)
%           Dims X Y Z care for LR, AP and IS translation
%
%           X-rotation will rotate the transverse slice (LR <-> AP)
%           swapping Y (coronal) & Z (saggital)
%           Y-rotation will rotate the coronal slice (IS <-> LR) slice,
%           swapping X (transversal) and Z (sagittal)
%           Z-rotation will rotate the sagittal slice (AP <-> IS)
%           swapping X (transversal) and Y (sagittal)
%
%           So, there are no differences in dim orders, only rotations.
%           E.g., MPRAGE is acquired in sagittal slices, and ASL/fMRI/BOLD in
%           transversal slices. This is an Y rotation (you look into the coronal
%           plane, rotate this, which will swap the sagittal slices into transversal)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: See: https://nipy.org/nibabel/coordinate_systems.html,
% And also the explanation inside xASL_im_DecomposeAffineTransformation.m
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: for only applying zooming & dim shifts: xASL_im_DummyOrientationNIfTI('../c1T1.nii', 'ASL4D.nii', OrientationDummy.nii, 0, 1, 1, 0]);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL




%% 0) Admin: by default, apply all orientations
if nargin<4 || isempty(bApplyRotationMinor)
    bApplyRotationMinor = true;
end
if nargin<5 || isempty(bApplyRotation90)
    bApplyRotation90 = true;
end
if nargin<6 || isempty(bApplyZoom)
    bApplyZoom = true;
end
if nargin<7 || isempty(bApplyTranslation)
    bApplyTranslation = true;
end
    
    

%% 1) Load orientations & calculate transformation
niiSrc = xASL_io_ReadNifti(PathSrc);
niiRef = xASL_io_ReadNifti(PathRef);

Mtransformation = niiRef.mat/niiSrc.mat; % obtain transformation between 2 orientations


%% 2) Calculate the desired transformation
% by replacing the disabled ones by the identity matrix eye(4,4)
[M, P] = xASL_im_DecomposeAffineTransformation(Mtransformation);

if ~bApplyTranslation
    M.Translation = eye(4,4);
end
if ~bApplyRotation90
    M.rotations90 = eye(4,4);
end
if ~bApplyRotationMinor
    M.rotationsMinor = eye(4,4);
end
if ~bApplyZoom
    M.Zoom = eye(4,4);
end

TransformMat = M.Translation*M.rotationsMinor*M.rotations90*M.Zoom; % Transformation matrix = Magnification x Rotations x Translation



%% 3) Calculate new orientation matrix
% Here we apply the net transformation matrix as calculated above, to the
% old orientation matrix
NewMat = TransformMat*niiSrc.mat;

%% 4) Calculate the new image size
% Here we create a dummy image, because we only need to know the image
% matrix size, and xASL_spm_reslice will do the rest when resampling to the
% dummy NIfTI

ZoomApply = [norm(TransformMat(1:3,1)) norm(TransformMat(1:3,2)) norm(TransformMat(1:3,3))];
NewImageSize = round(size(niiSrc.dat(:,:,:))./ZoomApply);
NewImageSize([1 2 3]) = NewImageSize([3 2 1]);

DummyImage = zeros(NewImageSize);

% ->>>>>>>>>> HERE WE STILL NEED TO ROTATE THE DUMMY IMAGE MATRIX INTO THE
% CORRECT IMAGE MATRIX SIZE
% I PUT HERE CODE THAT DOESNT WORK, BUT IT SHOULD BE SOMETHING SIMILAR
% X rotations
if P.rotations90(1)~=0
    DummyImage = xASL_im_rotate(DummyImage,P.rotations90(1));
end
% Y rotations
if P.rotations90(2)~=0
    DummyImage = shiftdim(DummyImage, 1); % First go to the Y plane by shifting dim once
    xASL_im_rotate(DummyImage, P.rotations90(2)); % perform the rotations
    DummyImage = shiftdim(DummyImage,2); % Bring back to the original plane
end
 % Z rotations
if P.rotations90(3)~=0
    DummyImage = shiftdim(DummyImage, 2); % First go to the Z plane by shifting dim twice
    xASL_im_rotate(DummyImage, P.rotations90(3)); % perform the rotations
    DummyImage = shiftdim(DummyImage, 1); % Bring back to the original plane
end


%% 5) Save the dummy NIfTI
xASL_io_SaveNifti(PathSrc, PathDummyOut, DummyImage, [], 0, NewMat);


end