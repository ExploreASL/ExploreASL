% This created the overview of M0 division artifacts, how to affect the qCBF images,
% with and without the techniques used in ExploreASL

x = ExploreASL_Master('',0);

x.P.PopDir   = 'C:\BackupWork\ASL\EPAD_xASL_Paper\analysis_EPAD\Population';

qCBF_new        = xASL_io_Nifti2Im(fullfile(x.P.PopDir,'qCBF_040EPAD00040_ASL_1.nii'));
PWI             = xASL_io_Nifti2Im(fullfile(x.P.PopDir,'PWI_040EPAD00040_ASL_1.nii'));
M0              = xASL_io_Nifti2Im(fullfile(x.P.PopDir,'noSmooth_M0_040EPAD00040_ASL_1.nii'));
M0_bf           = xASL_io_Nifti2Im(fullfile(x.P.PopDir,'M0_040EPAD00040_ASL_1.nii'));


%%%%% % RESCALING HERE
ScaleF1 = 3.5;
ScaleF2 = 1.3;
PWI = PWI.*ScaleF1;
%%%%%
Slice1 = 55;
Slice2 = 70;
Slice3 = 35;


MaskM0          = M0>0.3*max(M0(:));
% MaskWide        = dilate_erode_full(MaskM0,'dilate',dilate_erode_sphere(5));
MeanM0          = mean(M0(MaskM0));

MaskCBF         = ((PWI./M0).*MaskM0).*MeanM0.*(83/55);
NonMaskCBF      = ((PWI./M0).*1     ).*MeanM0.*(83/55);

cA  = 10;
cB  = 109;
cC  = 14;
cD  = 136;

M0              = xASL_im_ndnanfilter(M0,'gauss',[3.76 3.76 3.76]);

IM_NonMask      = imrotate(NonMaskCBF(cA:cB,cC:cD,Slice1),90);
IM_Mask         = imrotate(MaskCBF(cA:cB,cC:cD,Slice1),90);
IM_new          = imrotate(qCBF_new(cA:cB,cC:cD,Slice1),90);

% IM_NonMask2      = imrotate(NonMaskCBF(cA:cB,cC:cD,Slice2),90);
% IM_Mask2         = imrotate(MaskCBF(cA:cB,cC:cD,Slice2),90);
% IM_new2          = imrotate(qCBF_new(cA:cB,cC:cD,Slice2),90);
%
% IM_NonMask3      = imrotate(NonMaskCBF(cA:cB,cC:cD,Slice3),90);
% IM_Mask3         = imrotate(MaskCBF(cA:cB,cC:cD,Slice3),90);
% IM_new3          = imrotate(qCBF_new(cA:cB,cC:cD,Slice3),90);

PWI             = imrotate(PWI(cA:cB,cC:cD,Slice1),90);
M0_bf           = imrotate(M0_bf(cA:cB,cC:cD,Slice1),90) ./ MeanM0.*60;
M0              = imrotate(M0(cA:cB,cC:cD,Slice1),90)    ./ MeanM0.*60;
MaskM0          = imrotate(MaskM0(cA:cB,cC:cD,Slice1),90) .*M0;

IM_NonMask = IM_NonMask./ScaleF2;
IM_Mask = IM_Mask./ScaleF2;

% FullSliceIM      = [[IM_NonMask3 IM_NonMask IM_NonMask2];[IM_Mask3 IM_Mask IM_Mask2];[IM_new3 IM_new IM_new2]];
FullSliceIM      = [[PWI M0 IM_NonMask];[PWI MaskM0 IM_Mask];[PWI M0_bf IM_new]];
FullSliceIM2      = [IM_NonMask; IM_Mask; IM_new];
figure(1);imshow(FullSliceIM,[0 100],'InitialMagnification',250,'border','tight')
figure(2);imshow(FullSliceIM2,[0 100],'InitialMagnification',250,'border','tight','colormap',x.S.jet256)
