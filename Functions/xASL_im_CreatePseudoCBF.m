function xASL_im_CreatePseudoCBF(x, spatialCoV, bPVC)
%xASL_im_CreatePseudoCBF Create a pseudoCBF reference image for CBF-based
%ASL->T1w registration
%
% FORMAT: xASL_im_CreatePseudoCBF(x, spatialCoV[, bPVC])
%
% INPUT:
%   x             - structure containing fields with all information required to run this submodule (STRUCT, REQUIRED)
%   spatialCoV    - estimated spatialCoV of CBF image, that determines the mix
%                   of mean CBF, ATT biasfield and vascular artifacts
%                   (REQUIRED). When this parameter is set to 0 or lower, this
%                   function will skip creating the pseudoCBF NIfTI
%   bPVC          - boolean for performing regional partial volume correction and to set
%                   the values of the PseudoTissue based on its results, rather than a global scaling.
%                   Improved scaling saved to x.P.Path_mean_PWI_Clipped_DCT
%                   to be used, especially for non-NMI cost functions (OPTIONAL, DEFAULT = FALSE)
% OUTPUT: n/a
% OUTPUT FILES: in the ASL native space folder:
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
%
%              1. Create the pseudoTissue CBF reference image, if it doesnt exist already
%              2. Create the native space copies of ASL templates, if they dont exist already
%              3. Spatial CoV input argument check
%              4. Load native space copies of templates
%              5. Create pseudoTissue from segmentations, mix this with the mean CBF template depending on spatial CoV
%              6. Create pseudoCBF reference image used for CBF-based registration
%              7. Scale mean_PWI_Clipped source image to the same range as PseudoCBF
%
% EXAMPLE: xASL_im_CreatePseudoCBF(x, 0.4);
% __________________________________
% Copyright (C) 2015-2021 ExploreASL

if nargin < 3 || isempty(bPVC)
	bPVC = false;
end

%% ----------------------------------------------------------------------------------------
% 1) Create the pseudoTissue CBF reference image, if it doesnt exist already
if ~xASL_exist(x.D.Path_PseudoTissue,'file') || bPVC
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
	xASL_io_SaveNifti(x.P.Path_rc1T1, x.D.Path_PseudoTissue, PseudoTissue, [], false);

	% The rc1T1 can be deleted, unless you do PVC and it is used later in the code - it is then also deleted later
	if ~bPVC
		xASL_delete(x.P.Path_rc1T1);
		xASL_delete(x.P.Path_rc2T1);
	end

end

if bPVC
	% Prepare the PWI image in the same space as the PseudoTissue
	if exist(x.P.Path_mean_PWI_Clipped_sn_mat,'file')
		xASL_spm_reslice(x.D.Path_PseudoTissue, x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped_sn_mat, 0, x.Quality, x.P.Path_mean_PWI_Clipped_DCT);
	else
		xASL_spm_reslice(x.D.Path_PseudoTissue, x.P.Path_mean_PWI_Clipped, [], 0, x.Quality, x.P.Path_mean_PWI_Clipped_DCT);
	end
end

% Calculate the pseudoTissue again using PVC - that means it does the contrast matching locally and takes into account PV
if bPVC
	% Load the PV tissue priors
	imGM = xASL_io_Nifti2Im(x.P.Path_rc1T1);
	imWM = xASL_io_Nifti2Im(x.P.Path_rc2T1);
	imPWI = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped_DCT);

	% Run PV-correction with smooth kernel
	imGM(imGM<0) = 0;
	imGM(imGM>1) = 1;
	imWM(imWM<0) = 0;
	imWM(imWM>1) = 1;
	imPV = [];
	imPV(:,:,:,1) = imGM;
	imPV(:,:,:,2) = imWM;
	[imPVEC,~,~] = xASL_im_PVCkernel(imPWI,imPV,[5 5 5],'asllani');

	% Create a pseudoCBF image that has the local contrast variation based on the estimated PVEC
	PseudoTissue = imGM.*imPVEC(:,:,:,1) + imWM.*imPVEC(:,:,:,2);
	PseudoTissue(PseudoTissue<0) = 0;
	xASL_io_SaveNifti(x.P.Path_rc1T1, x.D.Path_PseudoTissue, PseudoTissue, [], false);

	% Delete the temporary files
	xASL_delete(x.P.Path_rc1T1);
	xASL_delete(x.P.Path_rc2T1);
end

%% ----------------------------------------------------------------------------------------
%% 2) Create the native space copies of ASL templates, if they dont exist already
if ~xASL_exist(x.D.Mean_Native,'file') || ~xASL_exist(x.D.Mask_Native,'file') || ~xASL_exist(x.D.Vasc_Native,'file')
    % We assume the same structural reference for all ASL sessions, so only have to warp these templates to native space once
    % Trilinear interpolation is fine for smooth template
    xASL_spm_deformations(x,{x.D.Mean_MNI;x.D.Bias_MNI;x.D.Vasc_MNI;x.D.Mask_MNI},{x.D.Mean_Native;x.D.Bias_Native;x.D.Vasc_Native;x.D.Mask_Native}, 1, x.D.Path_PseudoTissue, [], x.P.Path_y_ASL);
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
Mean_IM = xASL_io_Nifti2Im(x.D.Mean_Native);
Bias_IM = xASL_io_Nifti2Im(x.D.Bias_Native);
PWIim = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped);

%% ----------------------------------------------------------------------------------------
%% 5) Create pseudoTissue from segmentations, mix this with the mean CBF template depending on spatial CoV
PseudoTissue = xASL_io_Nifti2Im(x.D.Path_PseudoTissue);

% Here we mix the template-based (Mean_IM) and the tissue-based (PseudoProper)
% For spatial CoV 0.4 we can do 60% tissue-based, for spatial CoV 0 we can
% do 100% tissue-based.
PseudoMixFactor = min(spatialCoV,1);
Mean_IM = Mean_IM.*PseudoMixFactor+PseudoTissue.*(1-PseudoMixFactor);
%% ----------------------------------------------------------------------------------------
%% 6) Create pseudoCBF reference image used for CBF-based registration
% PseudoCBF equation
PseudoCBFim = Mean_IM - (Bias_IM.*5.*spatialCoV.^2);
% clip at 0
PseudoCBFim = max(PseudoCBFim, 0);
% add residual, for if no signal is left
PseudoCBFim = PseudoCBFim+Mean_IM./25;

if ~strcmpi(x.Sequence,'3D_spiral')
	% With 3D spiral, we nearly see no vascular artifacts because of low
	% effective spatial resolution
	Vasc_IM = xASL_io_Nifti2Im(x.D.Vasc_Native);
	PseudoCBFim = PseudoCBFim + Vasc_IM./5.*spatialCoV.^2;
end

% Save the reference NIfTI
xASL_io_SaveNifti(x.D.Mean_Native, x.P.Path_PseudoCBF, PseudoCBFim, [], 0);


%% ----------------------------------------------------------------------------------------
%% 7) Scale mean_PWI_Clipped source image to the same range as PseudoCBF for the rigid and affine registration
% rPWI (x.P.Path_mean_PWI_Clipped_DCT) is used for scaling and PWI.nii is used for DCT
if ~bPVC
	pseudoIM = xASL_io_Nifti2Im(x.P.Path_PseudoCBF);
	pseudoIM(isnan(pseudoIM)) = 0;
	PWIim(isnan(PWIim)) = 0;

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

	% Save the image for the rigid registration
	xASL_io_SaveNifti(x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped, PWIim, [], 0);
else
	% Calculate the proper scaling of rPWI (x.P.Path_mean_PWI_Clipped_DCT) to PseudoTissue, but leave the other mean_PWI_Clipped unaffected
	% Masks the perfused brain
	pseudoIM = xASL_io_Nifti2Im(x.P.Path_PseudoCBF);
	pseudoIM(isnan(pseudoIM)) = 0;
	imMask = pseudoIM > 10;
	imMask(:,:,[1:2,(end-1):(end)]) = 0;

	% Load only the perfusion values in a relevant range and mask out the air
    % "rPWI" here is x.P.Path_mean_PWI_Clipped_DCT
	PWIIM = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped);
	rPWIIM = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped_DCT);
	PWIIM(isnan(PWIIM)) = 0;
	rPWIIM(isnan(rPWIIM)) = 0;
	rPWIIM(rPWIIM<0) = 0;
	rPWIIMsort = sort(rPWIIM(imMask(:)),'ascend');
	rPWIIMmax = rPWIIMsort(floor(0.97*length(rPWIIMsort)));
	PWIIM(PWIIM<0) = 0;
	PWIIM(PWIIM>rPWIIMmax) = rPWIIMmax;
	rPWIIM(rPWIIM>rPWIIMmax) = rPWIIMmax;rPWIIM(isnan(rPWIIM)) = 0;
	X = rPWIIM(imMask);

	% Linear regression for least squares fit of the contrasts
	X = [ones(length(X),1),X];
	Y = pseudoIM(imMask);
	sol = pinv(X)*Y;

	% Modify and save the images
	rPWIIM = rPWIIM*sol(2)+sol(1);
	PWIIM = PWIIM*sol(2)+sol(1);
	xASL_io_SaveNifti(x.P.Path_mean_PWI_Clipped_DCT, x.P.Path_mean_PWI_Clipped_DCT, rPWIIM, [], 0);
	xASL_io_SaveNifti(x.P.Path_mean_PWI_Clipped, x.P.Path_mean_PWI_Clipped, PWIIM, [], 0);
end

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
