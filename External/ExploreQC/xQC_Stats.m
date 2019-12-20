function [Stats] = xQC_Stats(InputPath, WM2MAX, C1Path, C2Path, C3Path)
% xASL_qc_Stats computes several summary statistics 
%
% FORMAT: [Stats] = xASL_qc_Stats(InputPath, C1Path, C2Path, C3Path, WM2MAX)
%
% INPUT:
%    InputPath                              - Input image on which stats are computed can either be an image or a path (REQUIRED)
%    WM2MAX                                 - Set to 1 extract the WM2MAX coefficient (usually on  T1 images) (REQUIRED)
%    C1Path                                 - GM segmentation of the main image, either an image or a path (REQUIRED)
%    C2Path                                 - WM segmentation of the main image, either an image or a path (REQUIRED)
%    C3Path                                 - CSF segmentation of the main image, either an image or a path (REQUIRED)
%    
%
% OUTPUT:
%   Stats.(SegmentationName).mean           - Mean intensity in the specified segmentation of the input immage
%   Stats.(SegmentationName).Median         - Median intensity in the specified segmentation of the input immage 
%   Stats.(SegmentationName).SD             - Standard Deviation of intensity in the specified segmentation of the input immage
%   Stats.(SegmentationName).Value5         - 5% percentile of intensity distribution in the specified segmentation of the input immage
%   Stats.(SegmentationName).Value95        - 95% percentile of intensity distribution in the segmentation of the input immage
%   Stats.(SegmentationName).Skewnss        - Skewness of intensity distribution in the specified segmentation of the input immage
%   Stats.(SegmentationName).Kurtss         - Kurtosis of intensity  distribution in the specified segmentation of the input immage
%   Stats.WM2MAX                            - WM2MAX coefficient, i.e. the median intensity within the WM mask over the 95% percintile of the full intensity distribution  
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Description : This function computes several summary statistics on
% theintensity distribution of the input image. For each of the
% segmentation masks and for the full brain mask. if WM2MAX is set to 1 it
% computes this qc parameter (see above), usally for Structural imaging
%
% 
% EXAMPLE: 
%           [Stats] = xASL_qc_Stats(/pathToInput/bold.nii, 0)
%           [Stats] = xASL_qc_Stats(/pathToInput/T1.nii, 1, PathToGM/c1T1.nii,PathToGM/c2T1.nii, PathToGM/c3T1.nii)
%
%


%Load Input image if is a path 
InputImage = xASL_io_Nifti2Im(InputPath);
%Subject Directory
SubjectDir = xASL_fileparts(InputPath);

% Prepare mask values extraction
Masks = {};
MaskNames = {};
   
% reslice the masks into the image space
xASL_spm_reslice(InputPath, C1Path, [], [], [], fullfile(SubjectDir, 'GM_tmp_reslicedmask'), 0);
xASL_spm_reslice(InputPath, C2Path, [], [], [], fullfile(SubjectDir, 'WM_tmp_reslicedmask'), 0);
xASL_spm_reslice(InputPath, C3Path, [], [], [], fullfile(SubjectDir, 'CSF_tmp_reslicedmask'), 0);

% Load  the temporary resliced image
GMmask = xASL_io_Nifti2Im(fullfile(SubjectDir, 'GM_tmp_reslicedmask.nii'));
WMmask = xASL_io_Nifti2Im (fullfile(SubjectDir, 'WM_tmp_reslicedmask.nii'));
CSFmask = xASL_io_Nifti2Im(fullfile(SubjectDir, 'CSF_tmp_reslicedmask.nii'));

% Threshold Mask
GMmask = GMmask>0.5;
WMmask = WMmask>0.5;
CSFmask = CSFmask>0.5;

% Delete temporary files
xASL_delete(fullfile(SubjectDir, 'GM_tmp_reslicedmask.nii'))
xASL_delete(fullfile(SubjectDir, 'WM_tmp_reslicedmask.nii'))
xASL_delete(fullfile(SubjectDir, 'CSF_tmp_reslicedmask.nii'))

Masks = {GMmask , WMmask  , CSFmask};
MaskNames=  {'GM' , 'WM', 'CSF' };


   
% Create full brain mask 
AllBrainMask  = GMmask+WMmask+CSFmask;
AllBrainMask = AllBrainMask>0.5;
Masks{end+1} = AllBrainMask;
MaskNames{end+1} = 'AllBrain';

%Here it should make the 4d image into a 3D image but IT DOESNT WORK
if ndims(InputImage)>4
    InputImage = xASL_im_IM2Column(InputImage, AllBrainMask);
    

end     
    
    
% Extract summary statistics 
Stats = struct();

% Loop across Segmentations
for iMask = 1:length(Masks)
    
    % Select Mask 
    Segmentation = Masks{iMask};
    SegName = MaskNames{iMask};
    
    % Extract value from input and sort it
    SegmentedImage = (InputImage(Segmentation));
    SortedIm = sort(SegmentedImage(:));
    
    % Stats
    Stats.(SegName).Mean = xASL_stat_MeanNan(SortedIm);
    Stats.(SegName).Median = xASL_stat_MedianNan(SortedIm); 
    Stats.(SegName).SD = xASL_stat_StdNan(SortedIm);
    Stats.(SegName).Value5 = SortedIm(round(0.05*length(SortedIm)));
    Stats.(SegName).Value95 = SortedIm(round(0.95*length(SortedIm)));
    Stats.(SegName).Skewnss = (xASL_stat_SumNan((SegmentedImage-xASL_stat_MeanNan(SortedIm)).^3)./length(SortedIm)) ./ (xASL_stat_VarNan(SortedIm).^1.5);
    Stats.(SegName).Kurtss = (xASL_stat_SumNan((SegmentedImage-xASL_stat_MeanNan(SortedIm)).^4)./length(SortedIm)) ./ (xASL_stat_VarNan(SortedIm).^2);
end

% if WM2MAX is one, compute 
if WM2MAX == 1
    Stats.WM2MAX = Stats.WM.Median/Stats.AllBrain.Value95;
end 

end 









