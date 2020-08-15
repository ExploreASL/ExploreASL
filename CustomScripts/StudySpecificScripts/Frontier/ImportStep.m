%% You need to also load the P05, ASL4D.nii and change the (1,4) field in the matrix to 100 and save again - to prealign T1 and ASL - otherwise the registration subsequently fails
% This seems to work now with the new version of ExploreASL and no changes are necessary
%% Initialize ExploreASL

ExploreASL_Master('',0);
%% Initialize the paths
rawDir    = '/pet/projekte/asl/data/FRONTIER';
rawDirDSC = fullfile(rawDir,'rawDSC');
rawDirPET = fullfile(rawDir,'rawPET');
rawDirStruct = fullfile(rawDir,'rawStruct');
rawDirROI = fullfile(rawDir,'rawROIs');

%% Import the raw ASL and DSC data
imPar = ExploreASL_ImportConfig(rawDir);
ExploreASL_Import(imPar, 0,1,0);

jsonParms = spm_jsonread(fullfile(rawDir,'analysis','P05','ASL_1','ASL4D.json'));
jsonParms.MRScaleSlope = jsonParms.MRScaleSlope*2;
spm_jsonwrite(fullfile(rawDir,'analysis','P05','ASL_1','ASL4D.json'),jsonParms);

% Now we have them!
% Remove DSC of subjects 2 and 6 since these were not quantified
%xASL_delete(fullfile(rawDir,'analysis','P02','DSC_1','DSC4D.nii'));
%xASL_delete(fullfile(rawDir,'analysis','P02','DSC_1','DSC4D.json'));
%xASL_delete(fullfile(rawDir,'analysis','P02','DSC_1'));
%xASL_delete(fullfile(rawDir,'analysis','P06','DSC_1','DSC4D.nii'));
%xASL_delete(fullfile(rawDir,'analysis','P06','DSC_1','DSC4D.json'));
%xASL_delete(fullfile(rawDir,'analysis','P06','DSC_1'));

%% Copy the quantified DSC to the
patientNameList = xASL_adm_GetFileList(fullfile(rawDir,'analysis'), '^P\d{2}$', 'List', [], 1);

for iL = 1:length(patientNameList)
	fileNameList = xASL_adm_GetFileList(fullfile(rawDirDSC), ['^' patientNameList{iL} '.+$'], 'List', [], 0);
	for iF = 1:length(fileNameList)
		xASL_Copy(fullfile(rawDirDSC,fileNameList{iF}),fullfile(rawDir,'analysis',patientNameList{iL},'DSC_1',fileNameList{iF}(4:end)),1,0);
	end
end

%% Import raw PET
patientNameList = xASL_adm_GetFileList(fullfile(rawDir,'analysis'), '^P\d{2}$', 'List', [], 1);
for iL = 1:length(patientNameList)
	fileNameList = xASL_adm_GetFileList(fullfile(rawDirPET), ['^' patientNameList{iL} '.*\.v$'], 'List', [], 0);
	for iF = 1:length(fileNameList)
		[im mhead shead] = readEcat(fullfile(rawDirPET,fileNameList{iF}));
		im = im*100;
		imNew = zeros(size(im,2),size(im,1),size(im,3));
		for iS = 1:size(im,3)
			imNew(:,:,end-iS+1) = im(end:-1:1,end:-1:1,iS)';
		end
		if ~exist(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1'),'dir')
			mkdir(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1'));
		end
		xASL_io_CreateNifti( fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'), imNew, [shead.y_pixel_size, shead.x_pixel_size, shead.z_pixel_size]*10,32,0);
	end
end

%% Import raw PET that was in NII
patientNameList = xASL_adm_GetFileList(fullfile(rawDir,'analysis'), '^P\d{2}$', 'List', [], 1);
for iL = 1:length(patientNameList)
	fileNameList = xASL_adm_GetFileList(fullfile(rawDirPET), ['^' patientNameList{iL} '.*\.nii$'], 'List', [], 0);
	for iF = 1:length(fileNameList)
		im = xASL_io_Nifti2Im(fullfile(rawDirPET,fileNameList{iF}));
		im = im*100;
		if ~exist(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1'),'dir')
			mkdir(fullfile(rawDir,'analysis',patientNameList{iL},'PET_1'));
		end
		xASL_io_SaveNifti(fullfile(rawDirPET,fileNameList{iF}),fullfile(rawDir,'analysis',patientNameList{iL},'PET_1','CBF.nii'),im,32,0,[]); 
	end
end

%% Import FLAIR and T1

patientNameList = xASL_adm_GetFileList(fullfile(rawDir,'analysis'), '^P\d{2}$', 'List', [], 1);
for iL = 1:length(patientNameList)
	% Copy the T1w
	fileNameList = xASL_adm_GetFileList(fullfile(rawDirStruct), ['^' patientNameList{iL} 'T1G.+$'], 'List', [], 0);
	for iF = 1:length(fileNameList)
		xASL_Copy(fullfile(rawDirStruct,fileNameList{iF}),fullfile(rawDir,'analysis',patientNameList{iL},'T1.nii.gz'),1,0);
	end

	% Copy the FLAIR
	fileNameList = xASL_adm_GetFileList(fullfile(rawDirStruct), ['^' patientNameList{iL} 'FLR.+$'], 'List', [], 0);
	for iF = 1:length(fileNameList)
		xASL_Copy(fullfile(rawDirStruct,fileNameList{iF}),fullfile(rawDir,'analysis',patientNameList{iL},'FLAIR.nii.gz'),1,0);
	end
end
%% Import ROIs

patientNameList = xASL_adm_GetFileList(fullfile(rawDir,'analysis'), '^P\d{2}$', 'List', [], 1);
for iL = 1:length(patientNameList)
	switch(patientNameList{iL})
		% T1w ROI and T1w are perfectly aligned
		case {'P02','P04','P07','P08'}
			imROI = xASL_io_Nifti2Im(fullfile(rawDirROI,[patientNameList{iL} 'T1G.nii.gz']));
			imROI = imROI > 65000;
			imROI = shiftdim(imROI,2);
			imROI = imROI(:,end:-1:1,:);
			xASL_io_SaveNifti(fullfile(rawDir,'analysis',patientNameList{iL},'T1.nii.gz'), fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_T1_1.nii'), imROI, 32, 0);
		% FLAIR ROI and FLAIR are perfectly aligned
		case {'P03','P05','P06'}
			imROI = xASL_io_Nifti2Im(fullfile(rawDirROI,[patientNameList{iL} 'FLR.nii.gz']));
			imROI = imROI > 65000;
			imROI = shiftdim(imROI,2);
			imROI = imROI(:,end:-1:1,:);
			xASL_io_SaveNifti(fullfile(rawDir,'analysis',patientNameList{iL},'FLAIR.nii.gz'), fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_FLAIR_1.nii'), imROI, 32, 0);
	end

	% Realign the image with ROI with the resulting image

	switch(patientNameList{iL})
		% T1w ROI drawn on T1w aligned to T1w
		case 'P01'
			% Read the ROI file
			imT1w = xASL_io_Nifti2Im(fullfile(rawDirROI,[patientNameList{iL} 'T1G.nii.gz']));
			% Save the ROI
			imROI = imT1w > 65000;
			xASL_io_SaveNifti(fullfile(rawDirROI,[patientNameList{iL} 'T1G.nii.gz']), fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_T1_1.nii'), imROI, 32, 0);
			% Save the ROI-reference with ROI set to 0
			imT1w(imT1w > 65000) = 0;
			xASL_io_SaveNifti(fullfile(rawDirROI,[patientNameList{iL} 'T1G.nii.gz']), fullfile(rawDir,'analysis',patientNameList{iL},'Ref_T1_1.nii'), imT1w, 32, 0);

			xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'T1.nii'));
			% Align the T1w and the ROI-reference T1w and apply the same transformation to the lesion
			matlabbatch = [];
			matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(rawDir,'analysis',patientNameList{iL},'T1.nii,1')};
			matlabbatch{1}.spm.spatial.coreg.estwrite.source = {fullfile(rawDir,'analysis',patientNameList{iL},'Ref_T1_1.nii,1')};
			matlabbatch{1}.spm.spatial.coreg.estwrite.other = {fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_T1_1.nii,1')};
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
			spm_jobman('run',matlabbatch);

			% Delete the temporary files and rename the Lesion file
			xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'Ref_T1_1.nii'));
			xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'rRef_T1_1.nii'));
			xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_T1_1.nii'));
			xASL_Move(fullfile(rawDir,'analysis',patientNameList{iL},'rLesion_T1_1.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_T1_1.nii'));
	end

	switch(patientNameList{iL})
		% FLAIR ROI drawn on T1w aligned to FLAIR
		case {'P01','P02','P04','P07','P08'}
			% Read the ROI file
			imFLAIR = xASL_io_Nifti2Im(fullfile(rawDirROI,[patientNameList{iL} 'FLR.nii.gz']));
			% Save the ROI
			imROI = imFLAIR > 65000;
			xASL_io_SaveNifti(fullfile(rawDirROI,[patientNameList{iL} 'FLR.nii.gz']), fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_FLAIR_1.nii'), imROI, 32, 0);
			% Save the ROI-reference with ROI set to 0
			imFLAIR(imFLAIR > 65000) = 0;
			xASL_io_SaveNifti(fullfile(rawDirROI,[patientNameList{iL} 'FLR.nii.gz']), fullfile(rawDir,'analysis',patientNameList{iL},'Ref_FLAIR_1.nii'), imFLAIR, 32, 0);

			% Unzip the FLAIR
			xASL_adm_UnzipNifti(fullfile(rawDir,'analysis',patientNameList{iL},'FLAIR.nii'));

			% Align the T1w and the ROI-reference T1w and apply the same transformation to the lesion
			matlabbatch = [];
			matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {fullfile(rawDir,'analysis',patientNameList{iL},'FLAIR.nii,1')};
			matlabbatch{1}.spm.spatial.coreg.estwrite.source = {fullfile(rawDir,'analysis',patientNameList{iL},'Ref_FLAIR_1.nii,1')};
			matlabbatch{1}.spm.spatial.coreg.estwrite.other = {fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_FLAIR_1.nii,1')};
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
			matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 0;
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
			matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
			spm_jobman('run',matlabbatch);

			% Delete the temporary files and rename the Lesion file
			xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'Ref_FLAIR_1.nii'));
			xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'rRef_FLAIR_1.nii'));
			xASL_delete(fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_FLAIR_1.nii'));
			xASL_Move(fullfile(rawDir,'analysis',patientNameList{iL},'rLesion_FLAIR_1.nii'),fullfile(rawDir,'analysis',patientNameList{iL},'Lesion_FLAIR_1.nii'));
	end
end

%% Fix the T1contrast for P08
cd /pet/projekte/asl/data/FRONTIER/analysis/P08/
xASL_spm_reslice('T1.nii','Lesion_FLAIR_1.nii',[],[],1,'LF.nii',[]);
LT = xASL_io_Nifti2Im('Lesion_T1_1.nii');
LF = xASL_io_Nifti2Im('LF.nii');
T1 = xASL_io_Nifti2Im('T1.nii');
T1(LF>0.5) = 740;
T1(LT>0.5) = 740;
xASL_io_SaveNifti('T1.nii','T1new.nii',T1,[],1,[]);
xASL_delete('LF.nii');
xASL_Move('T1.nii','T1old.nii');
xASL_Move('/pet/projekte/asl/data/FRONTIER/analysis/P08/T1new.nii','/pet/projekte/asl/data/FRONTIER/analysis/P08/T1.nii');