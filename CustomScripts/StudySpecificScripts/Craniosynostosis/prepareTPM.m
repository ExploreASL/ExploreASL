pathTPM = '/pet/projekte/asl/data/Craniosynostosis-test/UNCInfant012Atlas_20140325';

nameTPM{1} = 'infant-1yr';
nameTPM{2} = 'infant-2yr';
nameTPM{3} = 'infant-neo';

%% Prepare a 6-component template that can be used for SPM12 segmentation
% TPM 1GM,2WM,3CSF,4cavity,5skull,6air
for ii=1:length(nameTPM)
	imAll = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '.nii']));
	imWithSkull = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '-withSkull.nii']));
	imWithCerebellum = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '-withCerebellum.nii']));
	
	imWM = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '-seg-wm.nii']));
	imGM = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '-seg-gm.nii']));
	imCSF = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '-seg-csf.nii']));
	
	imWM = imWM/250;
	imGM = imGM/150;
	imCSF = imCSF/10;
	
	imAir = (imWithSkull < 25).*(imWM==0).*(imCSF==0).*(imGM==0);
	imCavitySkull = (imWithCerebellum < 25).*(imAir==0).*(imWM==0).*(imCSF==0).*(imGM==0);
	imCavity = (imWithSkull<120) .* imCavitySkull;
	imSkull = (imWithSkull>=120) .* imCavitySkull;
	
	imTPM = zeros([size(imAll), 6]);
	imTPM(:,:,:,1) = imGM;
	imTPM(:,:,:,2) = imWM;
	imTPM(:,:,:,3) = imCSF;
	imTPM(:,:,:,4) = imCavity;
	imTPM(:,:,:,5) = imSkull;
	imTPM(:,:,:,6) = imAir;
	
	xASL_io_SaveNifti(fullfile(pathTPM,[nameTPM{ii} '.nii']),fullfile(pathTPM,[nameTPM{ii} '-TPM.nii']),imTPM,16,0);
end

%% Smooth to prepare the DARTEL template
for ii=1:length(nameTPM)
	% Presmooth
	xASL_im_PreSmooth(fullfile(pathTPM,'Template_6_IXI555_MNI152.nii'),fullfile(pathTPM,[nameTPM{ii} '-seg-wm.nii']));
	xASL_im_PreSmooth(fullfile(pathTPM,'Template_6_IXI555_MNI152.nii'),fullfile(pathTPM,[nameTPM{ii} '-seg-gm.nii']));
	
	% And resample to 1.5mm DARTEL space
	matlabbatch = [];
	matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathTPM,'Template_6_IXI555_MNI152.nii')};
	matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM,['s' nameTPM{ii} '-seg-gm.nii']),
		                                               fullfile(pathTPM,['s' nameTPM{ii} '-seg-wm.nii'])};
	matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
	matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
	spm_jobman('run',matlabbatch);

	xASL_delete(fullfile(pathTPM,['s' nameTPM{ii} '-seg-gm.nii']));
	xASL_delete(fullfile(pathTPM,['s' nameTPM{ii} '-seg-wm.nii']));

	% Load GM and WM and save as the step 6 of DARTEL
	imGM = xASL_io_Nifti2Im(fullfile(pathTPM,['ws' nameTPM{ii} '-seg-gm.nii']));
	imWM = xASL_io_Nifti2Im(fullfile(pathTPM,['ws' nameTPM{ii} '-seg-wm.nii']));
	
	xASL_delete(fullfile(pathTPM,['ws' nameTPM{ii} '-seg-gm.nii']));
	xASL_delete(fullfile(pathTPM,['ws' nameTPM{ii} '-seg-wm.nii']));
	
	imGM = imGM./150;
	imGM(:,:,:,2) = imWM./250;
	xASL_io_SaveNifti(fullfile(pathTPM,'Template_6_IXI555_MNI152.nii'),fullfile(pathTPM,[ 'Template_6_' nameTPM{ii} '_DARTEL.nii']),imGM,[],0);
	
	% Create the DARTEL steps by smoothing

	StartSmooth     = 1;
	EndSmooth       = 6;
	StepsN          = 4;

	xASL_Copy(fullfile(pathTPM,[ 'Template_6_' nameTPM{ii} '_DARTEL.nii']),fullfile(pathTPM,[ 'Template_5_' nameTPM{ii} '_DARTEL.nii']));
	for jj=1:5
		fwhm(jj,:) = repmat(EndSmooth-(jj-1)*((EndSmooth-StartSmooth)/StepsN),[1 3]);

		matlabbatch = [];
		INPUTim{1,1} = fullfile(pathTPM,[ 'Template_5_' nameTPM{ii} '_DARTEL.nii,1']);
		INPUTim{2,1} = fullfile(pathTPM,[ 'Template_5_' nameTPM{ii} '_DARTEL.nii,2']);

		matlabbatch{1}.spm.spatial.smooth.data      = INPUTim;
		matlabbatch{1}.spm.spatial.smooth.fwhm      = fwhm(jj,:);
		matlabbatch{1}.spm.spatial.smooth.dtype     = 0;
		matlabbatch{1}.spm.spatial.smooth.im        = 0;
		matlabbatch{1}.spm.spatial.smooth.prefix    = 's';

		spm_jobman('run',matlabbatch);
	
		xASL_Move( fullfile(pathTPM,['sTemplate_5_' nameTPM{ii} '_DARTEL.nii']) , fullfile(pathTPM,['Template_' num2str(jj-1) '_' nameTPM{ii} '_DARTEL.nii']), true);
	end
	
	StartSmooth     = 1;
	EndSmooth       = 4;
	StepsN          = 3;

	xASL_Copy(fullfile(pathTPM,[ 'Template_6_' nameTPM{ii} '_DARTEL.nii']),fullfile(pathTPM,[ 'Template_4_' nameTPM{ii} '_CAT.nii']));
	for jj=1:4
		fwhm(jj,:) = repmat(EndSmooth-(jj-1)*((EndSmooth-StartSmooth)/StepsN),[1 3]);

		matlabbatch = [];
		INPUTim{1,1} = fullfile(pathTPM,[ 'Template_4_' nameTPM{ii} '_CAT.nii,1']);
		INPUTim{2,1} = fullfile(pathTPM,[ 'Template_4_' nameTPM{ii} '_CAT.nii,2']);

		matlabbatch{1}.spm.spatial.smooth.data      = INPUTim;
		matlabbatch{1}.spm.spatial.smooth.fwhm      = fwhm(jj,:);
		matlabbatch{1}.spm.spatial.smooth.dtype     = 0;
		matlabbatch{1}.spm.spatial.smooth.im        = 0;
		matlabbatch{1}.spm.spatial.smooth.prefix    = 's';

		spm_jobman('run',matlabbatch);
	
		xASL_Move( fullfile(pathTPM,['sTemplate_4_' nameTPM{ii} '_CAT.nii']) , fullfile(pathTPM,['Template_' num2str(jj-1) '_' nameTPM{ii} '_CAT.nii']), true);
	end
	
end