function [M, P] = xASL_im_DecomposeAffineTransformation(Mtransformation)
%xASL_im_DecomposeAffineTransformation Decompose transformation matrix
% FORMAT: [M, P] = xASL_im_DecomposeAffineTransformation(Mtransformation)
%
% INPUT:  Mtransformation - SPM orientation matrix containing the transformation to be decomposed (REQUIRED)
%   
% OUTPUT: M = decomposed orientation matrix:
%         M.Full = original input transformation matrix
%         M.Translation = matrix for X Y Z translations
%         M.Zoom = matrix for X Y Z zoom (& flip)
%         M.rotations90 = matrix for 90 degree rotations
%         M.rotationsMinor = matrix for residual rotations & shearing
%
%         P = decomposed orientation parameters, with the same fieldnames
%         as described for M above
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function splits a transformation matrix into individual
%               components, which can be useful to guide the SPM reslicing.
%               The components are the same as in spm_(i)matrix.m, except for
%               the shearing: these are included in the rotations, and 
%               the 90 degree rotations, these are separated.
% 
%               Reason for the separation of the 90 degree rotations, is
%               that these indicate if orientations (transversal, coronal &
%               sagittal) have been switched in the NIfTI.
%
%               This can be useful to correct for any erroneous 90degree rotations from registration,
%               or to put two images in the same orientation order or voxelsize without applying
%               their subtle realignment (e.g. for manipulating registration references)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% THEORY 90 degree rotations:
%               Any rotation will always swap the other dims (X rotation swaps Y/Z, Y
%               rotation swaps X/Z etc.) because they are perpendicular (haaks)
% 
%               Dims X Y Z care for LR, AP and IS translation.
%             - X-rotation will rotate the transverse slice (LR <-> AP)
%               swapping Y (coronal) & Z (saggital)
%             - Y-rotation will rotate the coronal slice (IS <-> LR) slice,
%               swapping X (transversal) and Z (sagittal)
%             - Z-rotation will rotate the sagittal slice (AP <-> IS)
%               swapping X (transversal) and Y (sagittal)
%
%               E.g., MPRAGE is acquired in sagittal slices, and ASL/fMRI/BOLD in
%               transversal slices. This is an Y rotation (you look into the coronal
%               plane, rotate this, which will swap the sagittal slices into transversal)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function performs the following steps:
%
%               0. Start with an output P and M struct
%               1. Obtain translations
%               2. Obtain zoom
%               3. Obtain 90degree rotations
%               4. Obtain subtle rotations & shearing
%               5. Check the rounding errors of the decomposition
%
% EXAMPLE: [M, P] = xASL_im_DecomposeAffineTransformation(Mtransformation);
% See also: https://nipy.org/nibabel/coordinate_systems.html
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: T1w MPRAGE (acquired sagitally) has niiSrc.mat:
% (note the 1x1x0.56 mm resolution, & dim order shift)
%   -0.0079    0.0175    0.5599  -92.1031
%   -0.9444    0.2544   -0.0070   82.4730
%    0.2545    0.9443   -0.0085 -146.3968
%         0         0         0    1.0000
%   ASL CBF.nii has niiRef.mat:
% (note the 3x3x7 mm resolution)
%    0.0861   -0.7175    6.7090  -28.7328
%   -0.0179    2.8620    1.6206 -148.1755
%    3.0277    0.0375   -0.1895 -120.8767
%         0         0         0    1.0000
% Transformation matrix to go from T1w to ASL:
% MTransformation = niiRef.mat/niiSrc.mat = 
%    11.9642   -0.4249   -0.8667  981.3742
%     2.9457    0.7426    2.7763  468.3356
%    -0.3627   -2.9745    0.8477  215.1320
%          0         0         0    1.0000
% Where M.Zoom:
% (7/0.56=12.3mm, 3/1 = 3mm)
%    12.3269         0         0         0
%          0    3.0916         0         0
%          0         0    3.0224         0
%          0         0         0    1.0000
% and M.rotations90:
% (this is a single negative rotation in the first plane,
%  which moves the sagittal plane into transversal plane:
%     1.0000         0         0         0
%          0    cos(1)    sin(1)         0
%          0   -sin(1)    cos(1)         0
%          0         0         0    1.0000
%
%     1.0000         0         0         0
%          0    0.5403    0.8415         0
%          0   -0.8415    0.5403         0
%          0         0         0    1.0000
% and M.rotationsMinor:
%     0.9706   -0.3155   -0.0393         0
%     0.2390    0.9027    0.2942         0
%    -0.0294   -0.2838    0.9611         0
%          0         0         0    1.0000
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


%% 0) Start with an output P and M struct
M.Full = Mtransformation;
P.Full = spm_imatrix(M.Full); % translate into parameters

%% 1) Obtain translations
M.Translation = eye(4,4);
M.Translation(1:3,4) = M.Full(1:3,4);
P.Translation = M.Translation(1:3,4);

%% 2) Obtain zoom
PZoom = zeros(1,12);
PZoom(7:9) = P.Full(7:9);
M.Zoom = spm_matrix(PZoom);
P.Zoom = diag(M.Zoom(1:3,1:3));

%% 3) Obtain 90degree rotations
Rad90 = 90/(180/pi); % this is a ninety degree rotation in radians
Rotations90 = round(P.Full(4:6)./Rad90); % These are the round 90 degrees rotations
Protations90 = zeros(1,12);
Protations90(7:9) = 1;
Protations90(4:6) = Rotations90; % This contains the 90 degree rotations only
M.rotations90 = spm_matrix(Protations90);
P.rotations90 = Protations90(4:6);

%% 4) Obtain subtle rotations & shearing
M.rotationsMinor = eye(4,4); % identity matrix
M.rotationsMinor(1:3,1:3) = M.Full(1:3,1:3)/M.Zoom(1:3,1:3)/M.rotations90(1:3,1:3);

%% 5) Check the rounding errors of the decomposition
DiffMatrix = M.Full - M.Translation*M.rotationsMinor*M.rotations90*M.Zoom;
SumRoundingErrors = xASL_stat_SumNan(xASL_stat_SumNan(DiffMatrix./M.Full));

fprintf('%s\n',['Orientation matrix decomposed with a sum of ' xASL_num2str(SumRoundingErrors) ' rounding error']);


end

