function [MeanAI_PreTopUp_Perc, MeanAI_PostTopUp_Perc] = xASL_im_VisualQC_TopUp(PathPopB0, PathPopUnwarped, x, iSubject, CheckDir)
%xASL_im_VisualQC_TopUp Compute & visualize the effect of TopUp
%
% FORMAT: [MeanAI_PreTopUp_Perc, MeanAI_PostTopUp_Perc] = xASL_im_VisualQC_TopUp(PathPopB0, PathPopUnwarped, x, iSubject, CheckDir)
%
% INPUT:
%   PathPopB0 	    - path to NIfTI containing NormPE & RevPE images resampled to standard space (called B0 by FSL TopUp) (REQUIRED)
%   PathPopUnwarped - path to NIfTI containing NormPE & RevPE images after TopUp, resampled to standard space (called Unwarped by FSL TopUp) (REQUIRED)
%   x               - structure containing fields with all information required to run this submodule (REQUIRED)
%   iSubject        - index of current subject (REQUIRED)
%   CheckDir        - path to folder in Population dir where the Figure will be saved (REQUIRED)
%
% OUTPUT:
%   MeanAI_PreTopUp_Perc  - average relative voxel-wise asymmetry index (AI) between NormPE & RevPE before TopUp (%) (OPTIONAL)
%   MeanAI_PostTopUp_Perc - average relative voxel-wise asymmetry index (AI) between NormPE & RevPE after TopUp (%) (OPTIONAL)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function creates a Figure that showes the effect of TopUp
%              with 6 images with axial slices: the NormPE, RevPE and
%              their difference image in colorscale, and this before (upper
%              row) & after (lower row) TopUp.
%
%
% EXAMPLE: xASL_im_VisualQC_TopUp('analysis/Sub-001/dwi/B0.nii', 'analysis/Sub-001/dwi/Unwarped.nii', x, 10, 'analysis/Population/TopUpCheckDir');
%
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


%% Define images
Image{1} = xASL_io_Nifti2Im(PathPopB0);
Image{2} = xASL_io_Nifti2Im(PathPopUnwarped);

%% Reinforce transversal slices
if isfield(x.S, 'CorSlices')
    x.S = rmfield(x.S, 'CorSlices');
end
if isfield(x.S, 'SagSlices')
    x.S = rmfield(x.S, 'SagSlices');
end
x.S.TraSlices = x.S.slicesLarge;

%% Create the figure parts
IM1 = xASL_im_CreateVisualFig(x, Image{1}(:,:,:,1));
IM2 = xASL_im_CreateVisualFig(x, Image{1}(:,:,:,2));
IM3 = xASL_im_CreateVisualFig(x, Image{2}(:,:,:,1));
IM4 = xASL_im_CreateVisualFig(x, Image{2}(:,:,:,2));

%% Create the asymmetry figure parts
for ii=1:2
    Image{ii}(Image{ii}<0) = 0;
    DeltaIm{ii} = abs(Image{ii}(:,:,:,1)-Image{ii}(:,:,:,2));
    MeanIm{ii} = xASL_stat_MeanNan(Image{ii},4);
    AsymIndexIm{ii} = DeltaIm{ii} ./ MeanIm{ii};
end

%% Masking
MaskInt = sort(MeanIm{1}(isfinite(MeanIm{1})));
MaskInt = MaskInt(round(0.75*length(MaskInt)));
BrainMask = MeanIm{1}>MaskInt;

%% Calculate mean asymmetry indices
for ii=1:2
    MeanAsym(ii) = 100*xASL_stat_MeanNan(AsymIndexIm{ii}(BrainMask));
end

MeanAI_PreTopUp_Perc = MeanAsym(1);
MeanAI_PostTopUp_Perc = MeanAsym(2);

fprintf('%s\n','TopUp reduced the mean voxel-wise asymmetry index between NormPE & RevPE');
fprintf('%s\n',['within the BrainMask from ' xASL_num2str(MeanAI_PreTopUp_Perc) '% to ' xASL_num2str(MeanAI_PostTopUp_Perc) '%']);

%% Image clipping
% Clip both images manually now, at the same top
ClipValue = 0.95;
SortInt = sort(DeltaIm{1}(isfinite(DeltaIm{1})));
SortInt = SortInt(round(ClipValue*length(SortInt)));

for ii=1:2
    DeltaIm{ii}(DeltaIm{ii}>SortInt) = SortInt; % clip
    RMSimS{ii} = xASL_im_CreateVisualFig(x, DeltaIm{ii}, [], [], [], x.S.jet256, false);
end

%% Save the Figure
ComposIm = [IM1,IM2,RMSimS{1};IM3,IM4,RMSimS{2}];

OutputFile = fullfile(CheckDir, ['TopUp_QC_' x.SUBJECTS{iSubject} '.jpg']);
xASL_delete(OutputFile);
xASL_imwrite(ComposIm, OutputFile);


end
