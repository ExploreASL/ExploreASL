function xASL_im_CreatePseudoCBF(x, spatialCoV, bPVC)
%xASL_im_CreatePseudoCBF Create a pseudoCBF reference image for CBF-based
%ASL->T1w registration
%
% FORMAT: xASL_im_CreatePseudoCBF(x, spatialCoV[,bPVC])
%
% INPUT:
%   x          - structure containing fields with all information required to run this submodule (REQUIRED)
%   spatialCoV - estimated spatialCoV of CBF image, that determines the mix
%                of mean CBF, ATT biasfield and vascular artifacts
%                (REQUIRED). When this parameter is set to 0 or lower, this
%                function will skip creating the pseudoCBF NIfTI
%   bPVC       - Does not set fixed values of GM and WM CBF, but rather calculate the regional PVEC and set 
%                the values of the PseudoTissue based on that (OPTIONAL, DEFAULT = FALSE)
% OUTPUT: n/a
% OUTPUT FILES: in the ASL Session Folder:
%              - PseudoCBF.nii: final pseudoCBF image used for registration
%              - PseudoTissue.nii: pseudoCBF image created from pGM+pWM
%                from the subject only. Is blended in the pseudoCBF image
%              - Mean_CBF_Template.nii: native space copy of the CBF
%                template
%              - ATT_Biasfield.nii: same as the CBF template but with a
%                biasfield for transit times, simulating a CBF image of a
%                subject with high ATT and high spatial CoV
%              - VascularArtifact_Template.nii: average template of
%                vascular artifacts, simulating the vascular artifacts
%                observed in a subject with high spatial CoV, high ATT,
%                where the label still resides in the proximal arteries.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% RATIONALE:   This function creates a pseudoCBF image as reference image
%              for CBF-based registration. This allows to have a reliable
%              CBF-based registration even when there is less tissue CBF contrast
%              and more vascular contrast. The alternative is to use the
%              (mean) control image for registration, but this only
%              outperforms the CBF-based registration for images with 
%              high effective spatial resolution (2D EPI), and some ASL
%              sequence do not provide the Control images. Smooth ASL
%              sequences still have sufficient CBF contrast, but the
%              control image can be nearly flat, especially with background suppression and 
%              incomplete T1 relaxation. See Mutsaerts, Petr et al, JMRI 2018
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function creates a pseudo-CBF image from mean CBF template,
%              arterial transit time (ATT) bias field & vascular artifacts, weighted through spatial CoV
%              The first part of this code puts templates in the native space and
%              creates a pseudoCBF image from a downsampled pGM & pWM tissue (PseudoTissue). The latter
%              is used for registration but also as reference for the template
%              registration, to speed this up.
%              The second part of this code computes a pseudoCBF image based on the
%              pseudoTissue & the CBF templates of CBF, ATT biasfield and vascular peaks, based on spatial CoV.
%
%              This submodule performs the following steps:
%              1) Create the pseudoTissue CBF reference image, if it doesnt exist already
%              2) Create the native space copies of ASL templates, if they dont exist already
%              3) Spatial CoV input argument check
%              4) Load native space copies of templates
%              5) Create pseudoTissue from segmentations, mix this with the mean CBF template depending on spatial CoV
%              6) Create pseudoCBF reference image used for CBF-based registration
%              7) Scale mean_PWI_Clipped source image to the same range as PseudoCBF
%
% EXAMPLE: xASL_im_CreatePseudoCBF(x, 0.4);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL

if nargin < 3 || isempty(bPVC)
	bPVC = false;
end

%% ----------------------------------------------------------------------------------------
% 1) Create the pseudoTissue CBF reference image, if it doesnt exist already
if ~xASL_exist(x.Path_PseudoTissue,'file') || bPVC
	% Obtain the ASL VoxelSize (although the effective spatial resolution can
    % be lower). For isotropy and simplicity we take the highest resolution of
    % the 3 dimensions (smallest VoxelSize)
    nii = xASL_io_ReadNifti(x.P.Path_ASL4D);
    NewVoxelSize = repmat(min(nii.hdr.pixdim(2:4)),[1,3]);
    % Smooth the high resolution pGM & pWM images with the ASL voxelsize, and    
    
	% Pre-smooth the high resolution pGM & pWM images to match the ASL voxelsize, and
    % save them as temporary resampled copy prefixed with "r"
    xASL_spm_smooth(x.P.Path_c1T1, NewVoxelSize, x.P.Path_rc1T1);
    xASL_spm_smooth(x.P.Path_c2T1, NewVoxelSize, x.P.Path_rc2T1);
    % Downsample these temporary images to the ASL VoxelSize
    xASL_im_Upsample(x.P.Path_rc1T1, x.P.Path_rc1T1, NewVoxelSize);
    xASL_im_Upsample(x.P.Path_rc2T1, x.P.Path_rc2T1, NewVoxelSize);
	
    % Calculate the pseudoCBF image from the downsampled tissue posteriors,
    % called PseudoTissue
    PseudoTissue = xASL_io_Nifti2Im(x.P.Path_rc1T1).*80 + xASL_io_Nifti2Im(x.P.Path_rc2T1).*26.7;
    % Save this PseudoTissue NIfTI
    xASL_io_SaveNifti(x.P.Path_rc1T1, x.Path_PseudoTissue, PseudoTissue, [], false);
	if ~bPVC
		xASL_delete(x.P.Path_rc1T1);
		xASL_delete(x.P.Path_rc2T1);
	end
end

% Prepare the PWI image in the same space as the PseudoTissue
if exist(x.P.Path_mean_PWI_Clipped_sn_mat,'file')
	xASL_spm_reslice(x.Path_PseudoTissue, x.P.Path_mean_PWI_Clipped_ORI, x.P.Path_mean_PWI_Clipped_sn_mat, 0, x.Quality, x.P.Path_rPWI);
else
	xASL_spm_reslice(x.Path_PseudoTissue, x.P.Path_mean_PWI_Clipped_ORI, [], 0, x.Quality, x.P.Path_rPWI);
end

% Calculate the pseudoTissue again using PVC
if bPVC
	% Load the PV tissue priors
	imGM = xASL_io_Nifti2Im(x.P.Path_rc1T1);
	imWM = xASL_io_Nifti2Im(x.P.Path_rc2T1);
	xASL_delete(x.P.Path_rc1T1);
	xASL_delete(x.P.Path_rc2T1);
	imPWI = xASL_io_Nifti2Im(x.P.Path_rPWI);
	
	% Run PV-correction with smooth kernel
	imGM(imGM<0) = 0;
	imGM(imGM>1) = 1;
	imWM(imWM<0) = 0;
	imWM(imWM>1) = 1;
	imPV = [];
	imPV(:,:,:,1) = imGM;
	imPV(:,:,:,2) = imWM;
	[imPVEC,~,imResidual] = xASL_im_PVCkernel(imPWI,imPV,[5 5 5],'asllani');
	
	% Remove from mask all voxels with too high residuals after PVEc
	imMask = (imGM+imWM)>0.1;
	imErr = imResidual(imMask);
	meanErr = mean(imErr);
	stdErr = std(imErr);
	imMask = (imResidual < (meanErr + 2*stdErr));
	imMask = imMask>0;
	
	% Create a pseudoCBF image
	PseudoTissue = imGM.*imPVEC(:,:,:,1) + imWM.*imPVEC(:,:,:,2);
	PseudoTissue(PseudoTissue<0) = 0;
	PseudoTissue(isnan(PseudoTissue)) = 0;
	xASL_io_SaveNifti(x.P.Path_rc1T1, x.Path_PseudoTissue, PseudoTissue, [], false);
end

%% ----------------------------------------------------------------------------------------
%% 2) Create the native space copies of ASL templates, if they dont exist already
if ~xASL_exist(x.Mean_Native,'file') || ~xASL_exist(x.Mask_Native,'file') || ~xASL_exist(x.Vasc_Native,'file')
    % We assume the same structural reference for all ASL sessions, so only have to warp these templates to native space once
    % Trilinear interpolation is fine for smooth template
    xASL_spm_deformations(x,{x.Mean_MNI;x.Bias_MNI;x.Vasc_MNI;x.Mask_MNI},{x.Mean_Native;x.Bias_Native;x.Vasc_Native;x.Mask_Native}, 1, x.Path_PseudoTissue, [], x.P.Path_y_ASL);
end



%% ----------------------------------------------------------------------------------------
%% 3) Spatial CoV input argument check
% We discontinue here for spatial CoV = 0 or smaller
if numel(spatialCoV)~=1 % spatial CoV should be a single value
    error('Wrong definition of spatial CoV, this should be a single number');
elseif spatialCoV<eps
    % when the spatial CoV argument is 0 or lower, we skip creating the
    % pseudoCBF image
    return;
end



%% ----------------------------------------------------------------------------------------
%% 4) Load native space copies of templates
Mean_IM = xASL_io_Nifti2Im(x.Mean_Native);
Bias_IM = xASL_io_Nifti2Im(x.Bias_Native);
Bias_IM(isnan(Bias_IM)) = 0;
PWIim = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped);


%% ----------------------------------------------------------------------------------------
%% 5) Create pseudoTissue from segmentations, mix this with the mean CBF template depending on spatial CoV
PseudoTissue = xASL_io_Nifti2Im(x.Path_PseudoTissue);
% Here we mix the template-based (Mean_IM) and the tissue-based (PseudoProper)
% For spatial CoV 0.4 we can do 60% tissue-based, for spatial CoV 0 we can
% do 100% tissue-based.
PseudoMixFactor = min(spatialCoV,1);
Mean_IM = Mean_IM.*PseudoMixFactor+PseudoTissue.*(1-PseudoMixFactor);
Mean_IM(isnan(Mean_IM)) = 0;

%% ----------------------------------------------------------------------------------------
%% 6) Create pseudoCBF reference image used for CBF-based registration
% PseudoCBF equation
PseudoCBFim = Mean_IM - (Bias_IM.*5.*spatialCoV.^2);
% clip at 0
PseudoCBFim = max(PseudoCBFim, 0);
% add residual, for if no signal is left
PseudoCBFim = PseudoCBFim+Mean_IM./25;

if ~strcmp(x.Sequence,'3D_spiral')
    % With 3D spiral, we nearly see no vascular artifacts because of low
    % effective spatial resolution
    Vasc_IM = xASL_io_Nifti2Im(x.Vasc_Native);
	Vasc_IM(isnan(Vasc_IM)) = 0;
    PseudoCBFim = PseudoCBFim + Vasc_IM./5.*spatialCoV.^2;
end

% Save the reference NIfTI
xASL_io_SaveNifti(x.Mean_Native, x.P.Path_PseudoCBF, PseudoCBFim, [], 0);


%% ----------------------------------------------------------------------------------------
%% 7) Scale mean_PWI_Clipped source image to the same range as PseudoCBF
pseudoIM = xASL_io_Nifti2Im(x.P.Path_PseudoCBF);

% Use 5 and 95% percentiles instead of min and max
sortPWI = sort(PWIim(:),'ascend');
minPWI = sortPWI(ceil(0.05*length(sortPWI)));
maxPWI = sortPWI(ceil(0.95*length(sortPWI)));
sortPseu = sort(pseudoIM(:),'ascend');
minPseu = sortPseu(ceil(0.05*length(sortPseu)));
maxPseu = sortPseu(ceil(0.95*length(sortPseu)));
PWIim = PWIim-minPWI+minPseu;
% Use the whole range (max-min) not only the max for the scaling
RatioMax = (maxPseu-minPseu)/(maxPWI-minPWI);
PWIim = PWIim.*RatioMax;

xASL_io_SaveNifti(x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped, PWIim, [], 0);

% Calculate the proper scaling of rPWI to PseudoTissue
imMask = pseudoIM > 10;
imMask(:,:,[1:2,(end-1):(end)]) = 0;
PWIIM = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped_ORI);
rPWIIM = xASL_io_Nifti2Im(x.P.Path_rPWI);
rPWIIM(rPWIIM<0) = 0;
rPWIIMsort = sort(rPWIIM(imMask(:)),'ascend');
rPWIIMmax = rPWIIMsort(floor(0.97*length(rPWIIMsort)));
PWIIM(PWIIM<0) = 0;
PWIIM(PWIIM>rPWIIMmax) = rPWIIMmax;
rPWIIM(rPWIIM>rPWIIMmax) = rPWIIMmax;rPWIIM(isnan(rPWIIM)) = 0;
X = rPWIIM(imMask);
X = [ones(length(X),1),X];

Y = pseudoIM(imMask);
sol = pinv(X)*Y;
rPWIIM = rPWIIM*sol(2)+sol(1);
PWIIM = PWIIM*sol(2)+sol(1);
					
xASL_io_SaveNifti(x.P.Path_rPWI, x.P.Path_rPWI, rPWIIM, [], 0);
xASL_io_SaveNifti(x.P.Path_mean_PWI_Clipped_ORI, x.P.Path_PWI, PWIIM, [], 0);

%% Print what we did
fprintf('%s\n',['PseudoCBF.nii created for spatial CoV=' num2str(100*spatialCoV,3) '% & rescaled mean_PWI_Clipped.nii to this']);


end

%% This example shows the different blending of mean CBF images
% ii=1
% PseudoMixFactor = ii/10;
% MeanImIs{ii} = Mean_IM.*PseudoMixFactor+PseudoTissue.*(1-PseudoMixFactor);
% NewIm = MeanImIs{ii};
% 
% for ii=2:10
%     PseudoMixFactor = ii/10;
%     MeanImIs{ii} = Mean_IM.*PseudoMixFactor+PseudoTissue.*(1-PseudoMixFactor);
%     NewIm = [NewIm MeanImIs{ii}];
% end
% dip_image(NewIm)
