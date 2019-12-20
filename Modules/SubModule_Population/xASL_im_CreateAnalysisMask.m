function [x] = xASL_im_CreateAnalysisMask(x, Threshold)
%xASL_im_CreateAnalysisMask Creates population analysis mask
%
% FORMAT: [x] = xASL_im_CreateAnalysisMask(x, Threshold)
%
% INPUT:
%   x           - struct containing pipeline environment parameters, useful
%                 when only initializing ExploreASL/debugging (REQUIRED)
%   Threshold   - percentage of population that needs to have data within a
%                 voxel (OPTIONAL, DEFAULT=0.95)
%
% OUTPUT:
%   x                       - same as input
%   x.S.VBAmask             - final analysis mask, to be used later in the pipeline
%   x.S.MaskVascular        - vascular mask, to be used later in the pipeline
%   x.S.MaskSusceptibility  - susceptibility mask, to be used later in the pipeline
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function takes the mean population-based probability
%              maps of masks, thresholds and combines them:
%              A) Creation GM, WM & WholeBrain masks by p>0.5
%              B) Create, combine & save vascular, susceptibity & FoV
%                 masks:
%                 - MaskVascular
%                 - MaskSusceptibility = MaskSusceptibility & MaskFoV
%              C) Create & save VBA mask
%                 - MaskAnalysis = MaskVascular & MaskSusceptibility
%                 - x.S.VBAmask = MaskAnalysis & GMmask
%              D) Visualization: Creates a figure with columns being
%                 following maps/masks overlaid over mean population T1w:
%                 1) FoV probability 0-50% missing voxels
%                 2) Vascular 0-7.5% missing voxels
%                 3) Susceptibility 0-50% missing voxels
%                 4) Analysis mask
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_im_CreateAnalysisMask(x);
% __________________________________
% Copyright 2015-2019 ExploreASL


%% Admin
if nargin<2 || isempty(Threshold)
    Threshold = 0.95; % default threshold
end

% Define paths to create
VBA_maskPath = fullfile(x.D.PopDir,'VBA_mask_final.nii'); % later move to statsdir
MaskVascularPath = fullfile(x.S.StatsDir,'MaskVascular.nii');
MaskSusceptibilityPath = fullfile(x.S.StatsDir,'MaskSusceptibility.nii');

% Define pre-existing paths
PathFoV = fullfile(x.D.TemplatesStudyDir, 'FoV_bs-mean_Unmasked.nii');
PathVascularMask = fullfile(x.D.TemplatesStudyDir, 'MaskVascular_bs-mean.nii');
PathSusceptibilityMask = fullfile(x.D.TemplatesStudyDir, 'MaskSusceptibility_bs-mean.nii');
PathpGM = fullfile(x.D.TemplatesStudyDir, 'pGM_bs-mean.nii');
PathpWM = fullfile(x.D.TemplatesStudyDir, 'pWM_bs-mean.nii');
PathpCSF = fullfile(x.D.TemplatesStudyDir, 'pCSF_bs-mean.nii');
PathT1 = fullfile(x.D.TemplatesStudyDir, 'T1_bs-mean_Unmasked.nii');

if x.nSubjectsSessions<16
    fprintf('%s\n', ['Too few subjects (' num2str(x.nSubjectsSessions) ') to create population-based analysis mask']);
    x.S.MaskSusceptibility = xASL_im_IM2Column(ones(121,145,121),x.WBmask);
    x.S.VBAmask = x.S.MaskSusceptibility;
    return;
elseif ~xASL_exist(PathSusceptibilityMask, 'file') || ~xASL_exist(PathpGM, 'file')
    warning('Missing susceptibility maps, skipping...');
    return;
end

%% Creation GM, WM & WholeBrain masks by p>0.5
GMmask = xASL_io_Nifti2Im(PathpGM)>0.5;
WMmask = xASL_io_Nifti2Im(PathpWM)>0.5;

if xASL_exist(PathpCSF,'file')
    pCSF = xASL_io_Nifti2Im(PathpCSF)>0.5;
    WholeBrain = (GMmask+WMmask+pCSF)>0.5;
else
    WholeBrain = (GMmask+WMmask)>0.5;
end

GMmask = GMmask>0.5;
WMmask = WMmask>0.5;

%% B) Create, combine & save vascula, susceptibity & FoV masks
MaskVascular = xASL_io_Nifti2Im(PathVascularMask)>=Threshold;
MaskSusceptibility = xASL_io_Nifti2Im(PathSusceptibilityMask);


DoSusceptibility = true;
if ~isfield(x,'Sequence')
    error('x.Sequence parameter missing!');
elseif strcmp(lower(x.Sequence),'2d_epi')
    Path_Template = fullfile(x.D.MapsDir,'Templates','Susceptibility_pSignal_2D_EPI.nii');
elseif strcmp(lower(x.Sequence),'3d_grase')
    Path_Template = fullfile(x.D.MapsDir,'Templates','Susceptibility_pSignal_3D_GRASE.nii');
else
    % don't mask susceptibility artifacts
    DoSusceptibility = false;
end

MaskFoV = xASL_io_Nifti2Im(PathFoV)>=Threshold;
if DoSusceptibility
    TemplateMask = xASL_io_Nifti2Im(Path_Template)~=1; % MaskSusceptibility<0.85 &
    TemplateMask(:,:,1:10) = 0;
    TemplateMask = TemplateMask & MaskSusceptibility<0.8;
    ThrValues = MaskSusceptibility(TemplateMask);
    ThresholdSuscept = Threshold.*max(ThrValues(:));

    % Remove some extremes using 95% percentile,
    % as here we aim to obtain the 95% threshold within the probability map
    % (i.e. excluding the manually created mask to isolate the skull
    % regions of air cavities)
    MaskSusceptibility = MaskSusceptibility > ThresholdSuscept;
    
    % Combine susceptibility & FoV
    MaskSusceptibility = MaskSusceptibility & MaskFoV;
else
    MaskSusceptibility = MaskFoV; % if we don't have susceptibility artifacts
end


% Save them
xASL_io_SaveNifti(PathVascularMask, MaskVascularPath, MaskVascular, 8, true);
xASL_io_SaveNifti(PathSusceptibilityMask, MaskSusceptibilityPath, MaskSusceptibility, 8, true);
% this is used in stats:
x.S.MaskSusceptibility = xASL_im_IM2Column(MaskSusceptibility,x.WBmask);

%% C) Create & save VBA mask
MaskAnalysis = MaskSusceptibility;
x.S.VBAmask = MaskAnalysis & GMmask; % this should be an image matrix, not a column
xASL_io_SaveNifti(PathSusceptibilityMask, VBA_maskPath, x.S.VBAmask, 8, true);

%% D) Visualization
% 0) Admin
IntScale = [0.75 0.65];
ColorMap = {x.S.gray x.S.jet256};
bClip = [];
bWhite = false;
x.S.TraSlices = [53 35]; % 35]; % Z slice
x.S.CorSlices = 77; %
x.S.SagSlices = 63; % left right
x.S.ConcatSliceDims = 0;
x.S.bCrop = -15;
MasksAre = {ones(121,145,121) WholeBrain};


% 1) FoV probability
FoVim = 1-xASL_io_Nifti2Im(PathFoV);
FoVim(FoVim>0.5) = 0.5; % show 0-50% missing voxels in colors
ImOut{1} = xASL_im_CreateVisualFig(x, {PathT1 FoVim}, [], IntScale, [], ColorMap, bClip, [], bWhite);

% % 2) Vascular
% ImVascular = 1-xASL_io_Nifti2Im(PathVascularMask);
% ImVascular(ImVascular>0.075) = 0.075; % show
% ImOut{2} = xASL_im_CreateVisualFig(x, {PathT1 ImVascular}, [], IntScale, [], ColorMap, bClip, MasksAre, bWhite);

% 2) Susceptibility
ImSusceptibility = 1-xASL_io_Nifti2Im(PathSusceptibilityMask);
if DoSusceptibility
    % ImSusceptibility(ImSusceptibility>0.5) = 0.5; % ThresholdSuscept
    ImSusceptibility(ImSusceptibility<(1-ThresholdSuscept)) = 1-ThresholdSuscept;
    ImSusceptibility = ImSusceptibility-(1-ThresholdSuscept);
    ImSusceptibility(~TemplateMask) = 0;
end

ImOut{2} = xASL_im_CreateVisualFig(x, {PathT1 ImSusceptibility}, [], IntScale, [], ColorMap, bClip, MasksAre, bWhite);

% 3) Resulting analysis mask
ImOut{3} = xASL_im_CreateVisualFig(x, {PathT1 MaskAnalysis}, [], IntScale, [], ColorMap, bClip, MasksAre, bWhite);

% 5) Save image
xASL_imwrite([ImOut{1},ImOut{2},ImOut{3}], fullfile(x.S.StatsDir,'AnalysisMaskCreation.jpg'),[], false);

%     figure(1);  imshow([ImOut{1},ImOut{2},ImOut{3}],'border','tight','InitialMagnification',400)

%% Reset visualization settings
x.S.TraSlices = x.S.slicesLarge;
x.S.CorSlices = [];
x.S.SagSlices = [];

% %% Create the colorbars
% MaxValue = [0.5 0.75];
% ColorMap = ColorMap{2};
% for iBar=1:length(MaxValue)
%     if isfinite(MaxValue(iBar))
%         DummyIm = repmat([0:0.01:1].*MaxValue(iBar),[101 1]);
%         figure(iBar); imshow(DummyIm,[],'colormap',ColorMap,'InitialMagnification',400);
%         colorbar;
%     end
% end


%     Unused T1 mask
%     T1mask = xASL_io_Nifti2Im(PathT1);
%     SortInt = sort(T1mask(:));
%     Thr = SortInt(round(0.5*length(SortInt)));
%     T1mask = T1mask>Thr;


end
