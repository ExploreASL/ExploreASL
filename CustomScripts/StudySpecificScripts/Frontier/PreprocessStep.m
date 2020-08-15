%% Initialize ExploreASL
% First import data with ImportStep
% Then run the normal xASL analysis
% Then run the preprocess step to do the additional resampling and co-registration - something to be later integrated directly in xASL
%ExploreASL_Master('/pet/projekte/asl/data/FRONTIER/analysis/DataParFrontier.json');
ExploreASL_Master('',0);
%% Initialize the paths
rawDir    = '/pet/projekte/asl/data/FRONTIER';
PETresol  = [5.5 5.5 5.5];
patientNameList = xASL_adm_GetFileList(fullfile(rawDir,'analysis'), '^P\d{2}$', 'List', [], 1);
%patientNameList = xASL_adm_GetFileList(fullfile(rawDir,'analysis'), '^P05$', 'List', [], 1);
x.Quality = 1;
%% Downsample GM and WM and ROIs to the ASL space
for iL = 1:length(patientNameList)
	% Pre-smooth GM and WM
	%vol1 = spm_vol(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'));
	%vol2 = spm_vol(fullfile(rawDir,'analysis',patientNameList{iL},'T1.nii'));

	xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'c1T1.nii'));
	xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'c2T1.nii'));

	xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sc1T1.nii'),...
		[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_GM.nii'), 1);
	xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sc2T1.nii'),...
		[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_WM.nii'), 1);

	% Delete the pre-smoothed images
	xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'sc1T1.nii'));
	xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'sc2T1.nii'));

	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_T1_1.nii'))
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_T1_1.nii'));
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sLesion_T1_1.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Lesion_T1.nii'), 1);
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'sLesion_T1_1.nii'));
	end

	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_FLAIR_1.nii'))
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_FLAIR_1.nii'));
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sLesion_FLAIR_1.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Lesion_FLAIR.nii'), 1);
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'sLesion_FLAIR_1.nii'));
	end
end
%% Run the alignment of PET-CBF to ASL-CBF
for iL = 1:length(patientNameList)
	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'))
		xASL_im_CenterOfMass(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'));
		matlabbatch = [];
		matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii')};
		%matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_GM.nii')};
		matlabbatch{1}.spm.spatial.coreg.estwrite.source = {fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii')};
		matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
		matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
		matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
		matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
		matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
		matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
		matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
		matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
		matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
		spm_jobman('run',matlabbatch);
		xASL_Move(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','rCBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_PET.nii'),1);
	end
end
%% DSC - align all
for iL = 1:length(patientNameList)
	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'))
		xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'));
		xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','K2.nii'));
		xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','MTT.nii'));
		xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBV.nii'));
		xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBV_correct.nii'));
		xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC4D.nii'));
		xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','tMIP.nii'));
		xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','Tmax.nii'));
		xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','TTP.nii'));
		matlabbatch = [];
		%matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii')};
		matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fullfile(rawDir,'analysis',patientNameList{iL},'c1T1.nii')};
		matlabbatch{1}.spm.spatial.coreg.estimate.source = {fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii')};
		matlabbatch{1}.spm.spatial.coreg.estimate.other = {fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','K2.nii')
														   fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','MTT.nii')
														   fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBV.nii')
														   fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBV_correct.nii')
														   fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC4D.nii,1')
														   fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','tMIP.nii')
														   fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','Tmax.nii')
														   fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','TTP.nii')};
		matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
		matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
		matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
		matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
		spm_jobman('run',matlabbatch);
	end
end
%% DSC smooth and reslice
for iL = 1:length(patientNameList)
	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'))
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'));
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','srBF.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_DSC_rBF.nii'), 1)
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','srBF.nii'));
	end
end
%% Downsample GM and WM and ROIs to the DSC space
for iL = 1:length(patientNameList)
	if xASL_exist (fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'))
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'c1T1.nii'));
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'c2T1.nii'));
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'c3T1.nii'));

		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sc1T1.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_GM.nii'), 1);
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sc2T1.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_WM.nii'), 1);
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sc3T1.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_CSF.nii'), 1);

		% Delete the pre-smoothed images
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'sc1T1.nii'));
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'sc2T1.nii'));
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'sc3T1.nii'));

		if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_T1_1.nii'))
			xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_T1_1.nii'));
			xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sLesion_T1_1.nii'),...
				[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Lesion_T1.nii'), 1);
			xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'sLesion_T1_1.nii'));
		end

		if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_FLAIR_1.nii'))
			xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_FLAIR_1.nii'));
			xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sLesion_FLAIR_1.nii'),...
				[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Lesion_FLAIR.nii'), 1);
			xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'sLesion_FLAIR_1.nii'));
		end
	end
end
%% Run further alignment of ASL and T1w
for iL = 1:length(patientNameList)
	% Read CBF and GM/WM for intensity scaling
	imCBF = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'));
	imGM = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_GM.nii'));
	imWM = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_WM.nii'));

	% Otherwise align to a pseudoCBF
	% Run PV-correction with smooth kernel
	imGM(imGM<0) = 0;
	imGM(imGM>1) = 1;
	imWM(imWM<0) = 0;
	imWM(imWM>1) = 1;
	imPV = [];
	imPV(:,:,:,1) = imGM;
	imPV(:,:,:,2) = imWM;
	[imPVEC,~,imResidual] = xASL_im_PVCkernel(imCBF,imPV,[5 5 5],'asllani');

	% Remove from mask all voxels with too high residuals after PVEc
	imMask = (imGM+imWM)>0.1;
	imErr = imResidual(imMask);
	meanErr = mean(imErr);
	stdErr = std(imErr);
	imMask = (imResidual < (meanErr + 2*stdErr));
	imMask = imMask>0;

	% Create a pseudoCBF image
    imPVEC(imPVEC<-20) = -20;
    imPVEC(imPVEC>200) = 200;
    imPseudoCBF = imGM.*imPVEC(:,:,:,1) + imWM.*imPVEC(:,:,:,2);
    imPseudoCBF(imPseudoCBF<0) = 0;
    imPseudoCBF(isnan(imPseudoCBF)) = 0;

	% Combine the Lesion ROIs to exclude them
	imLesion = zeros(size(imMask));
	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Lesion_T1.nii'))
		imLesion = imLesion + (xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Lesion_T1.nii'))>0);
	end
	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Lesion_FLAIR.nii'))
		imLesion = imLesion + (xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Lesion_FLAIR.nii'))>0);
	end
	imLesion = 1-(imLesion>0);
	imLesion((imGM+imWM)<0.01) = 0;
	imLesion(isnan(imLesion)) = 0;
	xASL_io_SaveNifti(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_GM.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_PseudoMask.nii'), imLesion, [], 0);

	% Load the M0 image - also use for exclusions to mask the skull
	% Exclude the first and last three slices.
	imMask(:,:,[1:2,(end-1):(end)]) = 0;
	xASL_io_SaveNifti(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_GM.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Mask.nii'), imMask, [], 0);

	% Find the ideal scaling and apply it
	imMask = (imMask.*((imGM+imWM)>0.3)) > 0;
	X = imPseudoCBF(imMask);
	X = [ones(length(X),1),X];
	Y = imCBF(imMask);
	sol = pinv(X)*Y;
	imPseudoCBF = imPseudoCBF*sol(2)+sol(1);
	xASL_io_SaveNifti(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_GM.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Pseudo.nii'), imPseudoCBF, [], 0);
	xASL_io_SaveNifti(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_GM.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Affine.nii'), imCBF, [], 0);

	% Estimate the affine registration
	flags.smosrc = 1;
	flags.smoref = 1;
	flags.regtype = 'subj';
	flags.cutoff = 25;
	flags.nits = 16;
	flags.reg = 1;
	paramsAffine = spm_normalise(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Pseudo.nii'),...
		fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Affine.nii'),...
		fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Affine_src_sn.mat'),...
		fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_PseudoMask.nii'),...
		fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Mask.nii'),flags);

	% Apply the affine registration
	matlabbatch = [];
	matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname = {fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Affine_src_sn.mat')};
	matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = [NaN NaN NaN];
	matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = [NaN NaN NaN
		NaN NaN NaN];
	matlabbatch{1}.spm.util.defs.comp{2}.id.space = {fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_GM.nii')};
	matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {
		fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Affine.nii')
		};
	matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
	matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
	spm_jobman('run',matlabbatch);

	% Delete temporary files
	xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Affine_src_sn.mat'));
	xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Mask.nii'));
	xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_PseudoMask.nii'));

	xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Pseudo.nii'));
	xASL_Move(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','wCBF_Affine.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Deform.nii'),1);
	xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Affine.nii'));
end

%% Run further alignment of DSC and T1w
for iL = 1:length(patientNameList)
	if xASL_exist (fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'))
		% Read CBF and GM/WM for intensity scaling
		imCBF = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'));
		imGM = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_GM.nii'));
		imWM = xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_WM.nii'));

		% Otherwise align to a pseudoCBF
		% Run PV-correction with smooth kernel
		imGM(imGM<0) = 0;
		imGM(imGM>1) = 1;
		imWM(imWM<0) = 0;
		imWM(imWM>1) = 1;

		imPV = [];
		maskPV = (imGM+imWM)>0.2;
		imPV(:,:,:,1) = imGM;
		imPV(:,:,:,2) = imWM;
		[imPVEC,~,imResidual] = xASL_im_PVCkernel(imCBF.*maskPV,imPV.*repmat(maskPV,[1 1 1 2]),[5 5 5],'asllani');

		% Remove from mask all voxels with too high residuals after PVEc
		imMask = (imGM+imWM)>0.1;
		imErr = imResidual(imMask);
		meanErr = mean(imErr);
		stdErr = std(imErr);
		imMask = (imResidual < (meanErr + 2*stdErr));
		imMask = imMask>0;

		% Create a pseudoCBF image
		imPVEC(imPVEC<-20) = -20;
		imPVEC(imPVEC>200) = 200;
		imPseudoCBF = imPV(:,:,:,1).*imPVEC(:,:,:,1) + imPV(:,:,:,2).*imPVEC(:,:,:,2);
		imPseudoCBF(imPseudoCBF<0) = 0;
		imPseudoCBF(isnan(imPseudoCBF)) = 0;

		% Combine the Lesion ROIs to exclude them
		imLesion = zeros(size(imMask));
		if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Lesion_T1.nii'))
			imLesion = imLesion + (xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Lesion_T1.nii'))>0);
		end
		if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Lesion_FLAIR.nii'))
			imLesion = imLesion + (xASL_io_Nifti2Im(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Lesion_FLAIR.nii'))>0);
		end
		imLesion = 1-(imLesion>0);
		imLesion((imGM+imWM)<0.9) = 0;
		imLesion(isnan(imLesion)) = 0;
		xASL_io_SaveNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_GM.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_PseudoMask.nii'), imLesion, [], 0);

		% Load the M0 image - also use for exclusions to mask the skull
		% Exclude the first and last three slices.
		imMask(:,:,[1:4,(end-3):(end)]) = 0;
		imMask = (imMask.*(imCBF>0))>0;
		xASL_io_SaveNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_GM.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Mask.nii'), imMask, [], 0);

		% Find the ideal scaling and apply it
		imMask = (imMask.*((imGM+imWM)>0.8)) > 0;
		X = imPseudoCBF(imMask);
		X = [ones(length(X),1),X];
		Y = imCBF(imMask);
		sol = pinv(X)*Y;
		imPseudoCBF = imPseudoCBF*sol(2)+sol(1);% + imCSF*valCSF;
		xASL_io_SaveNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_GM.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Pseudo.nii'), imPseudoCBF, [], 0);
		xASL_io_SaveNifti(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_GM.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Affine.nii'), imCBF, [], 0);

		% Estimate the affine registration
		flags.smosrc = 1;
		flags.smoref = 1;
		flags.regtype = 'subj';
		flags.cutoff = 25;
		flags.nits = 16;
		flags.reg = 1;
		paramsAffine = spm_normalise(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Pseudo.nii'),...
			fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Affine.nii'),...
			fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Affine_src_sn.mat'),...
			fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_PseudoMask.nii'),...
			fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Mask.nii'),flags);

		% Apply the affine registration
		matlabbatch = [];
		matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname = {fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Affine_src_sn.mat')};
		matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = [NaN NaN NaN];
		matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = [NaN NaN NaN
			NaN NaN NaN];
		matlabbatch{1}.spm.util.defs.comp{2}.id.space = {fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_GM.nii')};
		matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {
			fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Affine.nii')
			};
		matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
		matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
		matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
		matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
		spm_jobman('run',matlabbatch);

		% Delete temporary files
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Affine_src_sn.mat'));
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Mask.nii'));
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_PseudoMask.nii'));

		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Pseudo.nii'));
		xASL_Move(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','wDSC_Affine.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Deform.nii'),1);
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Deform_rBF.nii'));
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Affine.nii'));
	end
end

%% Correctly downsample all to the ASL final resolution
for iL = 1:length(patientNameList)
	% Pre-smooth GM and WM
	xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'c1T1.nii'),[],PETresol,[]);
	xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'c2T1.nii'),[],PETresol,[]);
	xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sc1T1.nii'),...
		[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','Final_GM.nii'), 1);
	xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sc2T1.nii'),...
		[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','Final_WM.nii'), 1);
	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'))
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'c1T1.nii'),[],PETresol,[]);
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'c2T1.nii'),[],PETresol,[]);
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sc1T1.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','Final_GM.nii'), 1);
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'sc2T1.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','Final_WM.nii'), 1);
	end
	xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'sc1T1.nii'));
	xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'sc2T1.nii'));

	% Smooth the PET images
	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'))
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','Final_PET.nii'), 1);
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','Final_PET.nii'), 1);
	end

	% Smooth the CBF images
	xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),[],PETresol,[]);
	xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','sCBF.nii'),...
		[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','Final_CBF.nii'), 1);
	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'))
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),[],PETresol,[]);
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','sCBF.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','Final_CBF.nii'), 1);
	end
	xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','sCBF.nii'));

	% Smooth the CBF images
	xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Deform.nii'),[],PETresol,[]);
	xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','sCBF_Deform.nii'),...
		[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','Final_CBF_Deform.nii'), 1);
	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'))
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF_Deform.nii'),[],PETresol,[]);
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','sCBF_Deform.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','Final_CBF_Deform.nii'), 1);
	end
	xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','sCBF_Deform.nii'));

	if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'))
		% Smooth the DSC images
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'),[],PETresol,[]);
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','srBF.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','Final_DSC.nii'), 1);
		if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'))
			xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','rBF.nii'),[],PETresol,[]);
			xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','srBF.nii'),...
				[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','Final_DSC.nii'), 1);
		end
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','srBF.nii'));

		% Smooth the DSC images
		xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Deform.nii'),[],PETresol,[]);
		xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','sDSC_Deform.nii'),...
			[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','Final_DSC_Deform.nii'), 1);
		if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'))
			xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','DSC_Deform.nii'),[],PETresol,[]);
			xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','sDSC_Deform.nii'),...
				[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','Final_DSC_Deform.nii'), 1);
		end
		xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1','sDSC_Deform.nii'));
	end
end

%% Transform the multiple-ROI back to T1w. And then downsample to ASL and PET space

% Standard to T1w space
for iL = 1:length(patientNameList)
	for iN = 1:2
		if iN == 1
			fnm = 'Lesion_T1';
		else
			fnm = 'Lesion_FLAIR';
		end
		if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},[fnm '_1.nii']))
			% Unzip all the needed files
			xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'y_T1.nii'));
			xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'T1.nii'));
			xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},[fnm '_1.nii']));
			xASL_adm_UnzipNifti(fullfile(rawDir,'analysis','Population',['r' fnm '_1_' patientNameList{iL} '.nii']));
			matlabbatch = [];
			matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(rawDir,'analysis',patientNameList{iL},'y_T1.nii')};
			matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {fullfile(rawDir,'analysis',patientNameList{iL},'T1.nii')};
			matlabbatch{1}.spm.util.defs.comp{2}.id.space = {fullfile(rawDir,'analysis',patientNameList{iL},[fnm '_1.nii'])};
			matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(rawDir,'analysis','Population',['r' fnm '_1_' patientNameList{iL} '.nii'])};
			matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(rawDir,'analysis',patientNameList{iL})};
			matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 1;
			matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
			matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
			matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
			spm_jobman('run',matlabbatch);

			xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},['wr' fnm '_1_' patientNameList{iL} '.nii']),[],PETresol,[]);
			xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},['swr' fnm '_1_' patientNameList{iL} '.nii']),...
				[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'ASL_1',['Final_' fnm '.nii']), 1);
			if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'))
				xASL_im_PreSmooth(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),fullfile(rawDir,'analysis',patientNameList{iL},['wr' fnm '_1_' patientNameList{iL} '.nii']),[],PETresol,[]);
				xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'), fullfile(rawDir,'analysis',patientNameList{iL},['swr' fnm '_1_' patientNameList{iL} '.nii']),...
					[],[],x.Quality, fullfile(rawDir,'analysis',patientNameList{iL},'PET_1',['Final_' fnm '.nii']), 1);
			end
			xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},['wr' fnm '_1_' patientNameList{iL} '.nii']));
			xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},['swr' fnm '_1_' patientNameList{iL} '.nii']));
		end
	end
end

%% Transform the everything and realign in the T1 space

% Standard to T1w space
for iL = 1:length(patientNameList)
	fnm = 'Lesion_T1';
	if ~xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},[fnm '_1.nii']))
		fnm = 'Lesion_FLAIR';
	end
	%if xASL_exist(fullfile(rawDir,'analysis',patientNameList{iL},[fnm '_1.nii']))
	% Unzip all the needed files
	xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'y_T1.nii'));
	xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'T1.nii'));
	xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},[fnm '_1.nii']));
	xASL_adm_UnzipNifti(fullfile(rawDir,'analysis','Population',['r' fnm '_1_' patientNameList{iL} '.nii']));
	
	matlabbatch = [];
	matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(rawDir,'analysis',patientNameList{iL},'y_T1.nii')};
	matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {fullfile(rawDir,'analysis',patientNameList{iL},'T1.nii')};
	matlabbatch{1}.spm.util.defs.comp{2}.id.space = {fullfile(rawDir,'analysis',patientNameList{iL},[fnm '_1.nii'])};
	matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(rawDir,'analysis','Population',['r' fnm '_1_' patientNameList{iL} '.nii'])};
	matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(rawDir,'analysis',patientNameList{iL})};
	matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
	matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
	spm_jobman('run',matlabbatch);

	xASL_Move(fullfile(rawDir,'analysis',patientNameList{iL},['wr' fnm '_1_' patientNameList{iL} '.nii']),fullfile(rawDir,'analysis',patientNameList{iL},['Final_' fnm '_1_.nii']),1);
	
	xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},['Final_' fnm '_1_.nii']),fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','Final_CBF.nii'),[],[],1,fullfile(rawDir,'analysis',patientNameList{iL},'Final_CBF.nii'),[]);
	xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},['Final_' fnm '_1_.nii']),fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','Final_DSC.nii'),[],[],1,fullfile(rawDir,'analysis',patientNameList{iL},'Final_DSC.nii'),[]);
	xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},['Final_' fnm '_1_.nii']),fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','Final_PET.nii'),[],[],1,fullfile(rawDir,'analysis',patientNameList{iL},'Final_PET.nii'),[]);
	xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},['Final_' fnm '_1_.nii']),fullfile(rawDir,'analysis',patientNameList{iL},'T1.nii'),[],[],1,fullfile(rawDir,'analysis',patientNameList{iL},'Final_T1.nii'),[]);
	xASL_spm_reslice(fullfile(rawDir,'analysis',patientNameList{iL},['Final_' fnm '_1_.nii']),fullfile(rawDir,'analysis',patientNameList{iL},'T1_ORI.nii'),[],[],1,fullfile(rawDir,'analysis',patientNameList{iL},'Final_T1_ORI.nii'),[]);
		
end
