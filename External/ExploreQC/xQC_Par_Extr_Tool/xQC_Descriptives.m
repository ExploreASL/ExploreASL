
function [Descriptives] = xQC_Descriptives(InputPath, bWM2MAX, C1Path, C2Path, C3Path)
% xASL_qc_Stats computes several summary statistics 
%
% FORMAT: [Stats] = xASL_qc_Stats(InputPath, C1Path, C2Path, C3Path, WM2MAX)
%
% INPUT:
%    InputPath                              - Input image on which stats are computed can either be an image or a path (REQUIRED)
%    bWM2MAX                                - Boolean, set to true to extract the WM2MAX coefficient (usually on  T1 images) (OPTIONAL, DEFAULT=false)
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
%   Stats.(SegmentationName).Skewness       - Skewness of intensity distribution in the specified segmentation of the input immage
%   Stats.(SegmentationName).Kurtosis       - Kurtosis of intensity  distribution in the specified segmentation of the input immage
%   Stats.WM2MAX                            - WM2MAX coefficient, i.e. the median intensity within the WM mask over the 95% percintile of the full intensity distribution  
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Description: This function computes several summary statistics on
%              the intensity distribution of the input image. For each of the
%              segmentation masks and for the full brain mask. if WM2MAX is set to true it
%              computes this qc parameter (see above), usally for Structural imaging
%
% 
% EXAMPLE: 
% Stats = xASL_qc_Stats('/pathToInput/bold.nii', 0);
% Stats = xASL_qc_Stats('/pathToInput/T1.nii', 1, 'PathToGM/c1T1.nii', 'PathToGM/c2T1.nii', 'PathToGM/c3T1.nii');
 


%% Admin
if nargin<2 || isempty(bWM2MAX)
    bWM2MAX = false;
end

% Load Input image if is a path 
InputImage = xASL_io_Nifti2Im(InputPath);
% Subject Directory
SubjectDir = xASL_fileparts(InputPath);

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

% Prepare mask values extraction
Masks = {GMmask, WMmask, CSFmask};
MaskNames=  {'GM', 'WM', 'CSF'};
   
% Create full brain mask 
WholeBrainMask  = GMmask+WMmask+CSFmask;
WholeBrainMask = WholeBrainMask>0.5;
Masks{end+1} = WholeBrainMask;
MaskNames{end+1} = 'WholeBrain';

    
% Extract summary statistics 
Descriptives = struct();

% Loop across Segmentations
for iMask = 1:length(Masks)
    
    
    if ndims(InputImage)==4
        SegmentedImage = bsxfun(@times,InputImage,Masks{iMask}); 
        SortedIm = nonzeros(sort(SegmentedImage(:)));

    else 
        SegmentedImage = InputImage(Masks{iMask});
        SortedIm = sort(SegmentedImage(:)); 
        
        % Skewness and kurtosis just on 3D Images
        Descriptives.(MaskNames{iMask}).Skewness = (xASL_stat_SumNan((SegmentedImage-xASL_stat_MeanNan(SortedIm)).^3)./length(SortedIm)) ./ (xASL_stat_VarNan(SortedIm).^1.5);
        Descriptives.(MaskNames{iMask}).Kurtosis = (xASL_stat_SumNan((SegmentedImage-xASL_stat_MeanNan(SortedIm)).^4)./length(SortedIm)) ./ (xASL_stat_VarNan(SortedIm).^2);
   
    end 
    % Stats
    Descriptives.(MaskNames{iMask}).Mean = xASL_stat_MeanNan(SortedIm);
    Descriptives.(MaskNames{iMask}).Median = xASL_stat_MedianNan(SortedIm); 
    Descriptives.(MaskNames{iMask}).SD = xASL_stat_StdNan(SortedIm);
    Descriptives.(MaskNames{iMask}).Value5 = SortedIm(round(0.05*length(SortedIm)));
    Descriptives.(MaskNames{iMask}).Value95 = SortedIm(round(0.95*length(SortedIm)));
    
end 
% if bWM2MAX is true, compute 
if bWM2MAX
    Descriptives.WM2MAX = Descriptives.WM.Median/Descriptives.WholeBrain.Value95;
end 


end