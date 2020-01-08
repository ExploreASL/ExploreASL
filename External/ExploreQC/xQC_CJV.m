function [ CJV_Ratio ] = xQC_CJV(T1im, GMim, WMim)
% xASL_qc_CJV calculates the Coefficient of Joint Variance 
%
% FORMAT: [ CJV_Ratio ] = xASL_qc_CJV(GMim, WMim)
%
% INPUT:
%   T1Im         - Main structural image (Can be a path or an image matrix, REQUIRED)
%   GMim  	     - Segmentation Image of the Gray Matter (can be either path to NIfTI or image matrix, REQUIRED)
%   WMim  	     - Segmentation Image of the White Matter (can be either path to NIfTI or image matrix, REQUIRED)
%  
%
% OUTPUT:
%   CJV          - Coefficient of Joint Variance
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Description : This function calculates the Coefficient of Joint Variance within a T1 image, higher values are related to the presence of heavy head motion and
% large INU artefacts
%
% 
% EXAMPLE: [CJV] = xASL_qc_CJV('/examplepath/subjectdir/c1T1.nii', '/examplepath/subjectdir/c2T1.nii')
%
%

% Load Images
T1im = xASL_io_Nifti2Im(T1im);
GMim = xASL_io_Nifti2Im(GMim);
WMim = xASL_io_Nifti2Im(WMim);

GMmask = GMim>0.5;
WMmask = WMim>0.5;

GM = T1im(GMmask);
WM = T1im(WMmask);

% Compute mean intensity
meanGM = xASL_stat_MeanNan(GM(:));
meanWM = xASL_stat_MeanNan(WM(:));

% Comute Standard Deviation
sdGM = xASL_stat_StdNan(GM(:));
sdWM = xASL_stat_StdNan(WM(:));

%CJV 
CJV_Ratio = (sdWM+sdGM)/(meanWM-meanGM);
end