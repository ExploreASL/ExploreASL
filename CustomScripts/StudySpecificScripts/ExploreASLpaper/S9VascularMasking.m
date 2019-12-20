% This script assumes that we are running within debugging
% xASL_wrp_CreateAnalysisMask.m within 040EPAD00040

NegativeMaskMNI = xASL_im_MaskNegativeVascularSignal(x, 2); % standard space
PositiveMaskMNI = xASL_im_MaskPeakVascularSignal(x.P.Pop_Path_qCBF, x.P.Pop_Path_M0);

pGM = xASL_io_Nifti2Im(x.P.Pop_Path_rc1T1);
pWM = xASL_io_Nifti2Im(x.P.Pop_Path_rc2T1);
pCSF = xASL_io_Nifti2Im(x.P.Pop_Path_rc3T1);

NonGMmask = pWM>0.9 | pCSF>0.9;
WB = (pGM+pWM+pCSF)>0.5; % for visualization only

IM1 = xASL_io_Nifti2Im(x.P.Pop_Path_qCBF);
IM1(IM1>100) = 100;
IM1(IM1<0) = 0;

NegativeMaskMNI = xASL_im_DilateErodeFull(NegativeMaskMNI,'dilate',xASL_im_DilateErodeSphere(1));
PositiveMaskMNI = xASL_im_DilateErodeFull(PositiveMaskMNI,'dilate',xASL_im_DilateErodeSphere(1));
NegativeMaskMNI(NonGMmask) = 0;
PositiveMaskMNI(NonGMmask) = 0;
NegativeMaskMNI(~WB) = 0;
PositiveMaskMNI(~WB) = 0;

x.S.TraSlices = x.S.slicesLarge([6 8 10]);
x.S.ConcatSliceDims = 0;
x.S.Square = 0;

Image1 = xASL_im_CreateVisualFig(x, IM1);
Image2 = xASL_im_CreateVisualFig(x, {IM1 NegativeMaskMNI});
Image3 = xASL_im_CreateVisualFig(x, {IM1 PositiveMaskMNI});

figure(1);imshow([Image1 Image2 Image3],'InitialMagnification',250)

%% Check visualization native space
% pGM = xASL_io_Nifti2Im(x.P.Path_rc1T1);
% pWM = xASL_io_Nifti2Im(x.P.Path_rc2T1);
% pCSF = xASL_io_Nifti2Im(x.P.Path_rc3T1);
%
% NonGMmask = pWM>0.8 | pCSF>0.8;
% WB = (pGM+pWM+pCSF)>0.5; % for visualization only
%
% IM1 = xASL_io_Nifti2Im(x.P.Path_CBF);
% IM1(IM1>100) = 100;
% IM1(IM1<0) = 0;
%
% % NegativeMaskNative = xASL_im_DilateErodeFull(NegativeMaskNative,'dilate',xASL_im_DilateErodeSphere(1));
%
% NegativeMaskNative(NonGMmask) = 0;
% PositiveMaskNative(NonGMmask) = 0;
% NegativeMaskNative(~WB) = 0;
% PositiveMaskNative(~WB) = 0;
%
% dip_image([IM1 IM1+NegativeMaskNative.*60 IM1+PositiveMaskNative.*60])
