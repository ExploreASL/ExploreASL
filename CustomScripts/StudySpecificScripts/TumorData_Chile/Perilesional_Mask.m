%% Trial extrapolation as perilesional

SaveFile    = 'C:\Backup\ASL\Carolina\Trial_Lesion_nii\dartel\rLesion_T1_1_001_Lizama.nii';

LesionIM    = xASL_io_Nifti2Im(SaveFile);
LesionIM    = LesionIM(:,:,:,1);
T1File      = 'C:\Backup\ASL\Carolina\Trial_Lesion_nii\dartel\rT1_001_Lizama.nii';
pGMfile     = 'C:\Backup\ASL\Carolina\Trial_Lesion_nii\dartel\rc1T1_001_Lizama.nii';
pWMfile     = 'C:\Backup\ASL\Carolina\Trial_Lesion_nii\dartel\rc2T1_001_Lizama.nii';

LesionIM    = RobustMap2Mask(LesionIM);

[DistanceMap, x, y, z] = xASL_distanceTransform(LesionIM);
DistanceMap = 1.5 .* DistanceMap; % convert voxels to 1.5 mm
PeriMask2   = DistanceMap>0 & DistanceMap<=15;

PeriMask    = new_mask-LesionIM;

% BrainMasking
BrainMask   = (xASL_io_Nifti2Im(pGMfile)+xASL_io_Nifti2Im(pWMfile))>0.1;
PeriMask    = PeriMask.*BrainMask;
AddMask     = LesionIM+PeriMask;
ContraMask  = xASL_im_Flip(AddMask,1);

OutIm1      = visual_registration_MNI_1_5_INPUT(xASL_io_Nifti2Im(T1File),   LesionIM, 0.75,0.35,'red');
OutIm2      = visual_registration_MNI_1_5_INPUT(xASL_io_Nifti2Im(T1File),   PeriMask, 0.75,0.5,'green');
OutIm3      = visual_registration_MNI_1_5_INPUT(xASL_io_Nifti2Im(T1File), ContraMask, 0.75,0.75,'blue');

% Check masks, where to combine the images
OutIm1Mask  = OutIm1(:,:,1)~=OutIm1(:,:,2) | OutIm1(:,:,2)~=OutIm1(:,:,3)  | OutIm1(:,:,1)~=OutIm1(:,:,3);
OutIm1Mask  = repmat(OutIm1Mask,[1 1 3]);

OutIm2Mask  = OutIm2(:,:,1)~=OutIm2(:,:,2) | OutIm2(:,:,2)~=OutIm2(:,:,3)  | OutIm2(:,:,1)~=OutIm2(:,:,3);
OutIm2Mask  = repmat(OutIm2Mask,[1 1 3]);

OutIm3Mask  = OutIm3(:,:,1)~=OutIm3(:,:,2) | OutIm3(:,:,2)~=OutIm3(:,:,3)  | OutIm3(:,:,1)~=OutIm3(:,:,3);
OutIm3Mask  = repmat(OutIm3Mask,[1 1 3]);

% Combine OutIm1 2 & 3
OutIm4                  = zeros(size(OutIm1));
OutIm4(OutIm1==OutIm2)  = OutIm1(OutIm1==OutIm2);
OutIm4(OutIm1Mask)      = OutIm1(OutIm1Mask);
OutIm4(OutIm2Mask)      = OutIm2(OutIm2Mask);
OutIm4(OutIm3Mask)      = OutIm3(OutIm3Mask);


figure(1);imshow(OutIm4)

% Save segmentations
LesionIM(:,:,:,2)       = PeriMask;
LesionIM(:,:,:,3)       = ContraMask;
LesionIM                = logical(LesionIM);
xASL_io_SaveNifti(SaveFile,SaveFile,LesionIM,8,1);

%% Where red is the lesion segmentation, green perilesional, blue contralesional
