pathTPM = '/pet/projekte/asl/data/Craniosynostosis-test/UNCInfant012Atlas_20140325';

nameTPM{1} = 'infant-1yr';
nameTPM{2} = 'infant-2yr';
nameTPM{3} = 'infant-neo';

pathMaps = '/home/janpetr/code/ExploreASL/External/SPMmodified/MapsAdded';
pathSPM = '/home/janpetr/code/ExploreASL/External/SPMmodified';
pathTemplates = '/home/janpetr/code/ExploreASL/Maps/Templates';
pathAtlases = '/home/janpetr/code/ExploreASL/External/AtlasesNonCommercial';

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

%% Adapt the whole xASL templates including vascular territories/masks/etc to be transformed to pediatric sizes

for ii=1:length(nameTPM)
	mkdir(fullfile(pathMaps,nameTPM{ii}));
	mkdir(fullfile(pathMaps,nameTPM{ii},'VascularTerritories'));
	mkdir(fullfile(pathTemplates,nameTPM{ii}));
	mkdir(fullfile(pathAtlases,nameTPM{ii}));
	
	% First copy from the source to xASL directory
	xASL_Copy(fullfile(pathTPM,[nameTPM{ii} '.nii']),fullfile(pathMaps,nameTPM{ii},[nameTPM{ii} '.nii']));
	
	% Then align the pediatric map to the MNI of an adult
	matlabbatch = [];
	matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {fullfile(pathMaps,nameTPM{ii},[nameTPM{ii} '.nii,1'])};
	matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {fullfile(pathMaps,nameTPM{ii},[nameTPM{ii} '.nii,1'])};
	matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
	matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
	matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {fullfile(pathSPM,'tpm','TPM.nii')};
	matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
	matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
	matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
	matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
	matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
		78 76 85];
	matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
	matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
	matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
	spm_jobman('run',matlabbatch);
	
	% Cleanup
	xASL_delete(fullfile(pathMaps,nameTPM{ii},['w' nameTPM{ii} '.nii']));
	
	% Copy all the adult maps and transform to pediatric size
	listMap = {'rc1T1' 'rc2T1' 'rT1' 'rc1T1_ASL_res' 'rc2T1_ASL_res' 'rgrey' 'rbrainmask'};
	matlabbatch = [];
	matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(pathMaps,nameTPM{ii},['y_' nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {fullfile(pathMaps,nameTPM{ii},[nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{2}.id.space = {fullfile(pathMaps,'rc1T1.nii')};
	for jj=1:length(listMap)
		%xASL_Copy(fullfile(pathMaps,[listMap{jj} '.nii']), fullfile(pathMaps,nameTPM{ii},[listMap{jj} '.nii']),true);
		xASL_adm_UnzipOrCopy(pathMaps,[listMap{jj} '.nii'], fullfile(pathMaps,nameTPM{ii}),true);
		matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{jj,1} = fullfile(pathMaps,nameTPM{ii},[listMap{jj} '.nii']);
	end
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(pathMaps,nameTPM{ii})};
	matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
	matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
	matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
	spm_jobman('run',matlabbatch);
	
	% Cleanup
	for jj=1:length(listMap)
		xASL_Move(fullfile(pathMaps,nameTPM{ii},['w' listMap{jj} '.nii']),fullfile(pathMaps,nameTPM{ii},[listMap{jj} '.nii']),true);
	end

	% Copy all the adult masks and transform to pediatric size
	listMask = {'brainmask' 'TotalGM' 'TotalWM' 'WholeBrain'  'brainmask_supratentorial' 'ParenchymNarrow' 'LeftRight' 'DeepWM' 'CentralWM_QC' 'GhostSignalRatio' 'MNI_structural'};
	matlabbatch = [];
	matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(pathMaps,nameTPM{ii},['y_' nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {fullfile(pathMaps,nameTPM{ii},[nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{2}.id.space = {fullfile(pathMaps,'rc1T1.nii')};
	for jj=1:length(listMask)
		if jj == 1
			xASL_Copy(fullfile(pathMaps,[listMask{jj} '.nii']), fullfile(pathMaps,nameTPM{ii},[listMask{jj} '.nii']),true);
		else
			xASL_adm_UnzipOrCopy(pathMaps,[listMask{jj} '.nii.gz'], fullfile(pathMaps,nameTPM{ii}),true);
		end
		matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{jj,1} = fullfile(pathMaps,nameTPM{ii},[listMask{jj} '.nii']);
	end
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(pathMaps,nameTPM{ii})};
	matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
	matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
	matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
	spm_jobman('run',matlabbatch);
	
	% Cleanup
	for jj=1:length(listMask)
		xASL_Move(fullfile(pathMaps,nameTPM{ii},['w' listMask{jj} '.nii']),fullfile(pathMaps,nameTPM{ii},[listMask{jj} '.nii']),true);
	end

	% Save as .mat and zip
	for jj=1:length(listMask)
		IM = xASL_io_Nifti2Im(fullfile(pathMaps,nameTPM{ii},[listMask{jj} '.nii']));
		save(fullfile(pathMaps,nameTPM{ii},[listMask{jj} '.nii.mat']),'IM');
		gzip(fullfile(pathMaps,nameTPM{ii},[listMask{jj} '.nii']));
		delete(fullfile(pathMaps,nameTPM{ii},[listMask{jj} '.nii']));
	end	
	
	% Copy normal files
	listFile = {'TotalGM.tsv' 'TotalWM.tsv' 'WholeBrain.tsv' 'MNI_structural.tsv' 'LeftRight.tsv' 'LabelColors.mat' 'GhostSignalRatio.tsv' 'DeepWM.tsv' 'Identity_sn.mat' 'Identity_Deformation_y_T1.nii.gz'...
		        fullfile('VascularTerritories','TatuICA_PCA.tsv'),fullfile('VascularTerritories','LabelingTerritories.tsv'),fullfile('VascularTerritories','CortVascTerritoriesTatu.tsv'),fullfile('VascularTerritories','ATTbasedFlowTerritories.tsv')};
	for jj=1:length(listFile)
		xASL_Copy(fullfile(pathMaps,listFile{jj}), fullfile(pathMaps,nameTPM{ii},listFile{jj}),true);
	end
	
	% Copy all the vascular masks and transform to pediatric size
	listMask = {'ATTbasedFlowTerritories' 'CortVascTerritoriesTatu' 'LabelingTerritories' 'TatuICA_PCA'};
	matlabbatch = [];
	matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(pathMaps,nameTPM{ii},['y_' nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {fullfile(pathMaps,nameTPM{ii},[nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{2}.id.space = {fullfile(pathMaps,'rc1T1.nii')};
	for jj=1:length(listMask)
		% xASL_Copy(fullfile(pathMaps,'VascularTerritories',[listMask{jj} '.nii']), fullfile(pathMaps,nameTPM{ii},'VascularTerritories',[listMask{jj} '.nii']),true);
		xASL_adm_UnzipOrCopy(fullfile(pathMaps,'VascularTerritories'),[listMask{jj} '.nii.gz'], fullfile(pathMaps,nameTPM{ii},'VascularTerritories'),true);
		matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{jj,1} = fullfile(pathMaps,nameTPM{ii},'VascularTerritories',[listMask{jj} '.nii']);
	end
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(pathMaps,nameTPM{ii},'VascularTerritories')};
	matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
	matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
	matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
	spm_jobman('run',matlabbatch);
	
	% Cleanup
	for jj=1:length(listMask)
		xASL_Move(fullfile(pathMaps,nameTPM{ii},'VascularTerritories',['w' listMask{jj} '.nii']),fullfile(pathMaps,nameTPM{ii},'VascularTerritories',[listMask{jj} '.nii']),true);
	end
	
	% Save as .mat and zip
	for jj=1:length(listMask)
		IM = xASL_io_Nifti2Im(fullfile(pathMaps,nameTPM{ii},'VascularTerritories',[listMask{jj} '.nii']));
		save(fullfile(pathMaps,nameTPM{ii},'VascularTerritories',[listMask{jj} '.nii.mat']),'IM');
		gzip(fullfile(pathMaps,nameTPM{ii},'VascularTerritories',[listMask{jj} '.nii']));
		delete(fullfile(pathMaps,nameTPM{ii},'VascularTerritories',[listMask{jj} '.nii']));
	end

	% Copy all the templates and transform to pediatric size
	listTemp = {'ATT_BiasField' 'GE_3Dspiral_Product_CBF' 'MaxVesselTemplate' 'Philips_2DEPI_Bsup_CBF' 'Philips_2DEPI_noBsup_CBF' 'Philips_2DEPI_noBsup_Control'...
		        'Siemens_2DEPI_PCASL_noBsup_CBF' 'Siemens_2DEPI_PCASL_noBsup_Control' 'Siemens_3DGRASE_PASL_CBF' 'Siemens_3DGRASE_PCASL_Control_BiasfieldCorr_MoodStudy'...
				'Susceptibility_pSignal_2D_EPI' 'Susceptibility_pSignal_3D_GRASE'};
	matlabbatch = [];
	matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(pathMaps,nameTPM{ii},['y_' nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {fullfile(pathMaps,nameTPM{ii},[nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{2}.id.space = {fullfile(pathMaps,'rc1T1.nii')};
	for jj=1:length(listTemp)
		xASL_adm_UnzipOrCopy(fullfile(pathTemplates),[listTemp{jj} '.nii.gz'], fullfile(pathTemplates,nameTPM{ii}),true);
		matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{jj,1} = fullfile(pathTemplates,nameTPM{ii},[listTemp{jj} '.nii']);
	end
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(pathTemplates,nameTPM{ii})};
	matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
	matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
	matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
	spm_jobman('run',matlabbatch);
	
	% Cleanup
	for jj=1:length(listTemp)
		xASL_Move(fullfile(pathTemplates,nameTPM{ii},['w' listTemp{jj} '.nii']),fullfile(pathTemplates,nameTPM{ii},[listTemp{jj} '.nii']),true);
	end
	
	for jj=1:length(listTemp)
		if sum(jj==[4,6,11,12])
			IM = xASL_io_Nifti2Im(fullfile(pathTemplates,nameTPM{ii},[listTemp{jj} '.nii']));
			save(fullfile(pathTemplates,nameTPM{ii},[listTemp{jj} '.nii.mat']),'IM');
		end
		gzip(fullfile(pathTemplates,nameTPM{ii},[listTemp{jj} '.nii']));
		delete(fullfile(pathTemplates,nameTPM{ii},[listTemp{jj} '.nii']));
	end
	
	% Copy all the template masks and transform to pediatric size
	listTemp = {'Philips_2DEPI_Bsup_QC_mask' 'Siemens_3DGRASE_PASL_QC_mask'};
	matlabbatch = [];
	matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(pathMaps,nameTPM{ii},['y_' nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {fullfile(pathMaps,nameTPM{ii},[nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{2}.id.space = {fullfile(pathMaps,'rc1T1.nii')};
	for jj=1:length(listTemp)
		xASL_adm_UnzipOrCopy(fullfile(pathTemplates),[listTemp{jj} '.nii.gz'], fullfile(pathTemplates,nameTPM{ii}),true);
		matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{jj,1} = fullfile(pathTemplates,nameTPM{ii},[listTemp{jj} '.nii']);
	end
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(pathTemplates,nameTPM{ii})};
	matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
	matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
	matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
	spm_jobman('run',matlabbatch);
	
	% Cleanup
	for jj=1:length(listTemp)
		xASL_Move(fullfile(pathTemplates,nameTPM{ii},['w' listTemp{jj} '.nii']),fullfile(pathTemplates,nameTPM{ii},[listTemp{jj} '.nii']),true);
	end
	
	% Save as .mat and zip
	for jj=1:length(listTemp)
		gzip(fullfile(pathTemplates,nameTPM{ii},[listTemp{jj} '.nii']));
		delete(fullfile(pathTemplates,nameTPM{ii},[listTemp{jj} '.nii']));
	end	
		
	% Copy all the Atlases and transform to pediatric size
	listMask = {'Hammers' 'HOcort_CONN' 'HOsub_CONN' 'Thalamus'};
	matlabbatch = [];
	matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(pathMaps,nameTPM{ii},['y_' nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {fullfile(pathMaps,nameTPM{ii},[nameTPM{ii} '.nii'])};
	matlabbatch{1}.spm.util.defs.comp{2}.id.space = {fullfile(pathMaps,'rc1T1.nii')};
	for jj=1:length(listMask)
		xASL_adm_UnzipOrCopy(pathAtlases,[listMask{jj} '.nii.gz'], fullfile(pathAtlases,nameTPM{ii}),true);
		matlabbatch{1}.spm.util.defs.out{1}.pull.fnames{jj,1} = fullfile(pathAtlases,nameTPM{ii},[listMask{jj} '.nii']);
	end
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(pathAtlases,nameTPM{ii})};
	matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
	matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
	matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
	spm_jobman('run',matlabbatch);
	
	% Cleanup
	for jj=1:length(listMask)
		xASL_Move(fullfile(pathAtlases,nameTPM{ii},['w' listMask{jj} '.nii']),fullfile(pathAtlases,nameTPM{ii},[listMask{jj} '.nii']),true);
	end
	
	% Save as .mat and zip
	for jj=1:length(listMask)
		IM = xASL_io_Nifti2Im(fullfile(pathAtlases,nameTPM{ii},[listMask{jj} '.nii']));
		save(fullfile(pathAtlases,nameTPM{ii},[listMask{jj} '.nii.mat']),'IM');
		gzip(fullfile(pathAtlases,nameTPM{ii},[listMask{jj} '.nii']));
		delete(fullfile(pathAtlases,nameTPM{ii},[listMask{jj} '.nii']));
	end
		
	% Copy normal files
	listFile = {'Hammers.tsv' 'Hammers.txt' 'HOcort_CONN.tsv' 'Thalamus.tsv	'};
	for jj=1:length(listFile)
		xASL_Copy(fullfile(pathAtlases,listFile{jj}), fullfile(pathAtlases,nameTPM{ii},listFile{jj}),true);
	end

	
end
	