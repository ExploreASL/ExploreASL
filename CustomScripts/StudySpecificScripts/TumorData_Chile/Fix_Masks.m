%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes the finished calculated data a recalculate several masks, mainly to remove the edemas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPath = '/pet/projekte/asl/data/Chile/analysis';

pathNewQC = fullfile(dataPath,'Population/ExtraQC');
pathNewFiles = fullfile(dataPath,'Population/ExtraFiles');

doInitialize = 0;
doPrepareTemplates = 0;
doFixMasks = 0;
doRunPVEC = 0;
doCalculateResults = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the dirs etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doInitialize
	cd ~/code/ExploreASLREVAMP/ExploreASL/
	x = ExploreASL_Initialize('/pet/projekte/asl/data/Chile/DATA_PAR.m');
	cd ~/code/ExploreASLREVAMP/ExploreASL/CustomScripts/TumorData_CarolinaMendez/
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the WM from the mean image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doPrepareTemplates
	if ~exist(pathNewQC,'dir')
		mkdir(pathNewQC);
	end
	if ~exist(pathNewFiles,'dir')
		mkdir(pathNewFiles);
	end

	% Load the mean WM mask
	fileWM = fullfile(dataPath,'Population/Templates','pWM_bs-mean.nii');
	WMim = xASL_io_Nifti2Im(fileWM);

	% Threshold and erode the WM mask to have a region that should be WM only - i.e. free of edemas
	for ii = 1:2
		WMimNew = WMim>0.8;
		[WMimNew,~,~,~] = xASL_im_DistanceTransform(1-WMimNew);

		if ii == 1
			WMimNew = (WMimNew>1).*WMim;
		else
			WMimNew = (WMimNew>2).*WMim;
		end

		xASL_io_SaveNifti(fileWM,fullfile(pathNewFiles,['pWMtoRemove' num2str(ii) '.nii']),WMimNew);

		% Show the new WM overlayed over the mean T1w
		xASL_im_CreateVisualFig(x, {fullfile(dataPath,'Population/Templates','T1_bs-mean.nii'),fullfile(pathNewFiles,['pWMtoRemove' num2str(ii) '.nii'])},pathNewQC,[0.8 1],'rem',[]);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix masks for all subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doFixMasks
	% Load the WM to remove template
	imWMremove1 = xASL_io_Nifti2Im(fullfile(pathNewFiles,'pWMtoRemove1.nii'));
	imWMremove2 = xASL_io_Nifti2Im(fullfile(pathNewFiles,'pWMtoRemove2.nii'));

	% Go through all images
	for ii = 1:x.nSubjects %25 sub,[46,51,60,61]vasc

		% Load GM and WM files (in original resolution)
		imGM = xASL_io_Nifti2Im(fullfile(dataPath,'Population',['rc1T1_' x.SubjectNameList{ii} '.nii']));
		imWM = xASL_io_Nifti2Im(fullfile(dataPath,'Population',['rc2T1_' x.SubjectNameList{ii} '.nii']));

		imCBFun = xASL_io_Nifti2Im(fullfile(dataPath,'Population',['qCBF_untreated_' x.SubjectNameList{ii} '_ASL_1.nii']));
		imCBF = xASL_io_Nifti2Im(fullfile(dataPath,'Population',['qCBF_' x.SubjectNameList{ii} '_ASL_1.nii']));

		% Load lesion files %1 lesion, 2 peri, 3lesion+peri, 4contra+lesion+peri, 5ipsi-lesion-peri, 6contra
		if xASL_exist(fullfile(dataPath,'Population',['rLesion_T1_1_' x.SubjectNameList{ii} '.nii']))
			isLesion = 1;
		else
			isLesion = 0;
		end

		% Load the 2nd lesion
		if xASL_exist(fullfile(dataPath,'Population',['rLesion_T1_2_' x.SubjectNameList{ii} '.nii']))
			isEdema = 1;
		else
			isEdema = 0;
		end

		if isLesion
			imLesion = xASL_io_Nifti2Im(fullfile(dataPath,'Population',['rLesion_T1_1_' x.SubjectNameList{ii} '.nii']));

			imGMORI = double(xASL_io_Nifti2Im(fullfile(dataPath,'Population',['rc1T1_ORI_' x.SubjectNameList{ii} '.nii'])))/255;
			imWMORI = double(xASL_io_Nifti2Im(fullfile(dataPath,'Population',['rc2T1_ORI_' x.SubjectNameList{ii} '.nii'])))/255;

			if isEdema
				imEdema = xASL_io_Nifti2Im(fullfile(dataPath,'Population',['rLesion_T1_2_' x.SubjectNameList{ii} '.nii']));
			end

			% Define the ipsilateral and contralateral masks
			imMaskIpsi = zeros(size(imGMORI));
			imMaskContra = zeros(size(imGMORI));
			if sum(sum(sum(imLesion(62:end,:,:,1))))>sum(sum(sum(imLesion(1:60,:,:,1))))
				imMaskIpsi(62:end,:,:) = 1;
				imMaskContra(1:60,:,:) = 1;
			else
				imMaskIpsi(1:60,:,:) = 1;
				imMaskContra(62:end,:,:) = 1;
			end

			imLesion(:,:,:,2) = imLesion(:,:,:,2).*imMaskIpsi;
			imLesion(:,:,:,3) = imLesion(:,:,:,3).*imMaskIpsi;
			imLesion(:,:,:,4) = imLesion(:,:,:,4).*imMaskContra;

			% First checks the difference between GM and GM-ORI (outside of Lesion1) as this is mostly an edema
			maskEdema = (imGM-imGMORI);
			maskEdema(maskEdema<0.01) = 0; % to compensate for small rounding errors, and for bigger imGM
			maskEdema = maskEdema.*(1-imLesion(:,:,:,1));

			imGM = imGMORI.*(1-imLesion(:,:,:,1));
			imWM = imWMORI;

			% If the edema is manually segmented, then add it to the mask
			if isEdema
				imGM(imEdema(:,:,:,1)>0.1) = 0;
				imWM(imEdema(:,:,:,1)>0.1) = 1;
				maskEdema = maskEdema + (imEdema(:,:,:,1)>0.1);
			end
			maskEdema(maskEdema>1) = 1;

			% Look for regions that should belong to WM, but are in GM
			% And outside of lesion
 			imWMedemaIpsi   = imLesion(:,:,:,5).*imGM.*imWMremove2;
 			imWMedemaContra = imLesion(:,:,:,6).*imGM.*imWMremove2;
 			countLesVox = [sum(sum(sum(imWMedemaIpsi))), sum(sum(sum(imWMedemaContra)))];
 			disp(['Subject ' x.SubjectNameList{ii} ' GM in WM on L/R: ' num2str(countLesVox)])

			% We remove those from GM and put them do edema mask
			% To remove, use the less restricted WM mask
			imWMedemaIpsiRem   = imMaskIpsi.*imGM.*imWMremove1;
			imWMedemaContraRem = imMaskContra.*imGM.*imWMremove1;

			% If the count exceeds 950, then use the less restrictive template to flip GM/WM
			if countLesVox(1) > 2950
				imWM = imWM + imWMedemaIpsiRem;
				imGM = imGM - imWMedemaIpsiRem;
			end

			if countLesVox(2) > 2950
				imWM = imWM + imWMedemaContraRem;
				imGM = imGM - imWMedemaContraRem;
			end

			% Look for regions that should belong to WM, but are neither GM or WM
			% And outside of edema and lesions
			imGMWMinv = (1-imGM-imWM);
			imGMWMinv(imGMWMinv<0) = 0;
			imGMWMedemaIpsi   = imLesion(:,:,:,5).*imWMremove1.*imGMWMinv;
 			imGMWMedemaContra = imLesion(:,:,:,6).*imWMremove1.*imGMWMinv;
			countLesVox = [sum(sum(sum(imGMWMedemaIpsi))), sum(sum(sum(imGMWMedemaContra)))];
 			disp(['Subject ' x.SubjectNameList{ii} ' non GM/WM in WM on L/R: ' num2str(countLesVox)])

			% We can safely add these to edema and nothing else, as they are already not in GM or WM
			imWM = imWM + imGMWMedemaIpsi + imGMWMedemaContra;
			maskEdema = maskEdema + (imGMWMedemaIpsi + imGMWMedemaContra)>0.5;
			% If there are no lesions, then skip the exclusion step
			imLesionNoSuper = imLesion;
		end

		% Downsample the image to the correct ASL resolution
		% imGM and imWM are now in 1.5x1.5x1.5, we smooth them to 2.75x2.75x5
		% Get the resolution to voxels by dividing by 1.5, then -1 to get the smoothing kernel, then get the sigma from FWHM
		sigmaSmoothing = sqrt(([2.75,2.75,5]).^2 - 1)/2.3551/1.5;

		isSuper = 0;

		% Remove the superhigh-signal in the qCBF untreated from the GM and WM and lesion masks as well
		if (~isempty(strfind(x.SubjectNameList{ii},'0002'))) ||...
				(~isempty(strfind(x.SubjectNameList{ii},'0005'))) ||...
				(~isempty(strfind(x.SubjectNameList{ii},'0017'))) ||...
				(~isempty(strfind(x.SubjectNameList{ii},'0018'))) ||...
				(~isempty(strfind(x.SubjectNameList{ii},'0019'))) ||...
				(~isempty(strfind(x.SubjectNameList{ii},'0024'))) ||...
				(~isempty(strfind(x.SubjectNameList{ii},'0049'))) ||...
				(~isempty(strfind(x.SubjectNameList{ii},'0051'))) ||...
				(~isempty(strfind(x.SubjectNameList{ii},'0037'))) ||...
				(~isempty(strfind(x.SubjectNameList{ii},'0043'))) ||...
				(~isempty(strfind(x.SubjectNameList{ii},'0054'))) ||...
				(~isempty(strfind(x.SubjectNameList{ii},'0055')))
			maskSuper = double(imCBFun > 250);
			%maskSuper = xASL_im_DilateErodeSeparable(maskSuper,'dilate',[1 1 1 1 1 1 1],[1 1 1 1 1 1 1],1);
			maskSuper = spm_dilate(maskSuper,ones(3,3,3));
			%maskSuper = spm_dilate(maskSuper,ones(3,3,3));

			maskSuper = maskSuper > 0;
			imWM(maskSuper) = 0;
			imGM(maskSuper) = 0;

			if isLesion
				imLesion(repmat(maskSuper,[1 1 1 6])) = 0;
				maskEdema(maskSuper) = 0;
			end

			[imSuperSmo,~,~,~] = xASL_im_Smooth3D(sigmaSmoothing,double(maskSuper),{'gaussian','gaussian','gaussian'});

			imCBF(imSuperSmo>0.1) = 0;
			isSuper = 1;
		end

		% Downsampling of the rest

		[imGMSmo,~,~,~] = xASL_im_Smooth3D(sigmaSmoothing,double(imGM),{'gaussian','gaussian','gaussian'});
		[imWMSmo,~,~,~] = xASL_im_Smooth3D(sigmaSmoothing,double(imWM),{'gaussian','gaussian','gaussian'});
		if isLesion
			[imEdemaSmo,~,~,~] = xASL_im_Smooth3D(sigmaSmoothing,double(maskEdema),{'gaussian','gaussian','gaussian'});
		end

		% And save
		xASL_io_SaveNifti(fullfile(dataPath,'Population',['rc1T1_' x.SubjectNameList{ii} '.nii']),...
			fullfile(pathNewFiles,['rc1T1_Smo_' x.SubjectNameList{ii} '.nii']),imGMSmo);
		xASL_io_SaveNifti(fullfile(dataPath,'Population',['rc1T1_' x.SubjectNameList{ii} '.nii']),...
			fullfile(pathNewFiles,['rc2T1_Smo_' x.SubjectNameList{ii} '.nii']),imWMSmo);
		if isLesion
			xASL_io_SaveNifti(fullfile(dataPath,'Population',['rLesion_T1_1_' x.SubjectNameList{ii} '.nii']),...
				fullfile(pathNewFiles,['rLesion_T1_1_' x.SubjectNameList{ii} '.nii']),imLesion);
			xASL_io_SaveNifti(fullfile(dataPath,'Population',['rLesion_T1_1_' x.SubjectNameList{ii} '.nii']),...
				fullfile(pathNewFiles,['rLesion_T1_1_NoSuper' x.SubjectNameList{ii} '.nii']),imLesionNoSuper);
			xASL_io_SaveNifti(fullfile(dataPath,'Population',['rLesion_T1_1_' x.SubjectNameList{ii} '.nii']),...
				fullfile(pathNewFiles,['rEdema_' x.SubjectNameList{ii} '.nii']),imEdemaSmo);
			xASL_io_SaveNifti(fullfile(dataPath,'Population',['rLesion_T1_1_' x.SubjectNameList{ii} '.nii']),...
				fullfile(pathNewFiles,['rAddEdema1_' x.SubjectNameList{ii} '.nii']),imWMedemaIpsi+imWMedemaContra);
			xASL_io_SaveNifti(fullfile(dataPath,'Population',['rLesion_T1_1_' x.SubjectNameList{ii} '.nii']),...
				fullfile(pathNewFiles,['rAddEdema2_' x.SubjectNameList{ii} '.nii']),imGMWMedemaIpsi+imGMWMedemaContra);
		end
		if isSuper
			xASL_io_SaveNifti(fullfile(dataPath,'Population',['rLesion_T1_1_' x.SubjectNameList{ii} '.nii']),...
				fullfile(pathNewFiles,['rSuper_' x.SubjectNameList{ii} '.nii']),imSuperSmo);
		end
		xASL_io_SaveNifti(fullfile(dataPath,'Population',['qCBF_' x.SubjectNameList{ii} '_ASL_1.nii']),...
			fullfile(pathNewFiles,['qCBF_' x.SubjectNameList{ii} '_ASL_1.nii']),imCBF);

	    % Show QC before and after
		if isLesion
			ImOut1 = xASL_im_CreateVisualFig(x, {fullfile(dataPath,'Population',['rT1_' x.SubjectNameList{ii} '.nii']) fullfile(pathNewFiles,['rLesion_T1_1_' x.SubjectNameList{ii} '.nii'])},[],[0.8 1],[],[]);
			%ImOut1 = xASL_im_CreateVisualFig(x, {fullfile(dataPath,'Population',['rT1_' x.SubjectNameList{ii} '.nii']) fullfile(pathNewFiles,['rEdema_' x.SubjectNameList{ii} '.nii'])},[],[0.8 1],[],[]);
			%ImOut1 = xASL_im_CreateVisualFig(x, {fullfile(dataPath,'Population',['rT1_' x.SubjectNameList{ii} '.nii']) fullfile(pathNewFiles,'pWMtoRemove1.nii')},[],[0.8 1],[],[]);
			ImOut2 = xASL_im_CreateVisualFig(x, {fullfile(pathNewFiles,['rc1T1_Smo_' x.SubjectNameList{ii} '.nii']),fullfile(pathNewFiles,['rc2T1_Smo_' x.SubjectNameList{ii} '.nii'])},[],[0.8 1],[],[]);
		%	if isLesion
		%		ImOut2 = xASL_im_CreateVisualFig(x, {fullfile(pathNewFiles,['rc1T1_Smo_' x.SubjectNameList{ii} '.nii']),fullfile(pathNewFiles,['rLesion_T1_1_' x.SubjectNameList{ii} '.nii'])},[],[0.8 1],[],[]);
		%	else
		%		ImOut2 = xASL_im_CreateVisualFig(x, {fullfile(pathNewFiles,['rc1T1_Smo_' x.SubjectNameList{ii} '.nii'])},[],[0.8 1],[],[]);
			%end
			%ImOut3 = xASL_im_CreateVisualFig(x, {fullfile(dataPath,'Population',['rT1_' x.SubjectNameList{ii} '.nii']) fullfile(pathNewFiles,['rAddEdema2_' x.SubjectNameList{ii} '.nii'])},[],[0.8 1],[],[]);
			ImOut3 = xASL_im_CreateVisualFig(x, {fullfile(pathNewFiles,['qCBF_' x.SubjectNameList{ii} '_ASL_1.nii'])},[],[0.8 1],[],[]);
			IM          = [ImOut1,ImOut2,ImOut3];

			xASL_imwrite((IM+eps)./max(IM(:)), fullfile( pathNewQC ,['EdemaCheck_' x.SubjectNameList{ii} '.jpg']));
		end

		if isSuper
			IM = xASL_im_CreateVisualFig(x, {fullfile(dataPath,'Population',['rT1_' x.SubjectNameList{ii} '.nii']) fullfile(pathNewFiles,['rSuper_' x.SubjectNameList{ii} '.nii'])},[],[0.8 1],[],[]);
			xASL_imwrite((IM+eps)./max(IM(:)), fullfile( pathNewQC ,['VascularArtifact_' x.SubjectNameList{ii} '.jpg']));
		end
	end
end

if doRunPVEC
	for ii = 1:x.nSubjects
		imGM = xASL_io_Nifti2Im(fullfile(pathNewFiles,['rc1T1_Smo_' x.SubjectNameList{ii} '.nii']));
		imWM = xASL_io_Nifti2Im(fullfile(pathNewFiles,['rc2T1_Smo_' x.SubjectNameList{ii} '.nii']));
		imFOV = xASL_io_Nifti2Im(fullfile(dataPath,'Population',['FoV_' x.SubjectNameList{ii} '_ASL_1.nii']));
		imCBF = xASL_io_Nifti2Im(fullfile(pathNewFiles,['qCBF_' x.SubjectNameList{ii} '_ASL_1.nii']));

		imExcludeMask = ((imFOV<0.9) + ((imGM+imWM)<0.05))>0.1;
		imGM(imExcludeMask) = 0;
		imWM(imExcludeMask) = 0;
		imCBF(imExcludeMask) = 0;

		imPV = imGM;
		imPV(:,:,:,2) = imWM;

		[imPVEC,imCBFrec,imResidual,FWHM] = xASL_im_PVCbspline(imCBF,imPV,[21 23 21]);

		xASL_io_SaveNifti(fullfile(dataPath,'Population',['rc1T1_' x.SubjectNameList{ii} '.nii']),...
			fullfile(pathNewFiles,['rGMPVEC_' x.SubjectNameList{ii} '.nii']),imPVEC(:,:,:,1).*(1-imExcludeMask));
		xASL_io_SaveNifti(fullfile(dataPath,'Population',['rc1T1_' x.SubjectNameList{ii} '.nii']),...
			fullfile(pathNewFiles,['rWMPVEC_' x.SubjectNameList{ii} '.nii']),imPVEC(:,:,:,2).*(1-imExcludeMask));
	end
end

if doCalculateResults
	imFOVmean = xASL_io_Nifti2Im(fullfile(dataPath,'Population','Templates','FoV_bs-mean.nii'));
	imFOVmean = imFOVmean>0.9;

	% Load the HO atlas
	imHO = xASL_io_Nifti2Im(fullfile(x.D.AtlasDir,'HOcort_CONN.nii'));
	maskBroca = ((imHO == 5) + (imHO == 6)) > 0;

	imMNI = xASL_io_Nifti2Im(fullfile(x.D.AtlasDir,'MNI_structural.nii'));

	% Set the left and right masks
	maskLeft = zeros(size(imHO));
	maskRight = zeros(size(imHO));
	maskLeft(1:60,:,:) = 1;
	maskRight(62:end,:,:) = 1;

	for ii = 1:x.nSubjects
		imGM = xASL_io_Nifti2Im(fullfile(pathNewFiles,['rc1T1_Smo_' x.SubjectNameList{ii} '.nii']));
		imWM = xASL_io_Nifti2Im(fullfile(pathNewFiles,['rc2T1_Smo_' x.SubjectNameList{ii} '.nii']));
		imGMPVEC = xASL_io_Nifti2Im(fullfile(pathNewFiles,['rGMPVEC_' x.SubjectNameList{ii} '.nii']));
		imWMPVEC = xASL_io_Nifti2Im(fullfile(pathNewFiles,['rWMPVEC_' x.SubjectNameList{ii} '.nii']));
		imFOV = xASL_io_Nifti2Im(fullfile(dataPath,'Population',['FoV_' x.SubjectNameList{ii} '_ASL_1.nii']));
		imCBF = xASL_io_Nifti2Im(fullfile(pathNewFiles,['qCBF_' x.SubjectNameList{ii} '_ASL_1.nii']));
		imFOV = imFOV>0;

		if xASL_exist(fullfile(dataPath,'Population',['rLesion_T1_1_' x.SubjectNameList{ii} '.nii']))
			isLesion = 1;
			imLesion = xASL_io_Nifti2Im(fullfile(pathNewFiles,['rLesion_T1_1_' x.SubjectNameList{ii} '.nii']));
			imLesionOrig = xASL_io_Nifti2Im(fullfile(pathNewFiles,['rLesion_T1_1_NoSuper' x.SubjectNameList{ii} '.nii']));
			imEdema = xASL_io_Nifti2Im(fullfile(pathNewFiles,['rEdema_' x.SubjectNameList{ii} '.nii']));
		else
			isLesion = 0;
		end

		if xASL_exist(fullfile(pathNewFiles,['rSuper_' x.SubjectNameList{ii} '.nii']))
			isSuper = 1;
			imSuper = xASL_io_Nifti2Im(fullfile(pathNewFiles,['rSuper_' x.SubjectNameList{ii} '.nii']));
		else
			isSuper = 0;
		end

		switch(4)
			case {1}
				if isLesion
					fprintf([x.SubjectNameList{ii} ',']);

					% Mean/90%quantile within tumor
					dat = imCBF((imFOV.*imLesion(:,:,:,1))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
					fprintf([num2str(xASL_stat_QuantileNan(dat,0.9,1),'%.2f'),','])
					fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])
					if isSuper
						dat = imSuper((imFOV.*imLesionOrig(:,:,:,1))>0);
						dat(dat<=0) = NaN;
						fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])
					else
						fprintf([num2str(0,'%.2f'),','])
					end

					% Mean/90%quantile perilesion
					dat = imCBF((imFOV.*imLesion(:,:,:,2))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
					fprintf([num2str(xASL_stat_QuantileNan(dat,0.9,1),'%.2f'),','])
					fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])
					if isSuper
						dat = imSuper((imFOV.*imLesionOrig(:,:,:,2))>0);
						dat(dat<=0) = NaN;
						fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])
					else
						fprintf([num2str(0,'%.2f'),','])
					end

					% Mean/90%quantile perilesion GM
					dat = imCBF((imFOV.*imLesion(:,:,:,2).*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
					fprintf([num2str(xASL_stat_QuantileNan(dat,0.9,1),'%.2f'),','])
					fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

					% Mean/90%quantile Contralateral
					dat = imCBF((imFOV.*imLesion(:,:,:,4))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
					fprintf([num2str(xASL_stat_QuantileNan(dat,0.9,1),'%.2f'),','])
					fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

					% Mean/90%quantile Contralateral GM
					dat = imCBF((imFOV.*imLesion(:,:,:,4).*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
					fprintf([num2str(xASL_stat_QuantileNan(dat,0.9,1),'%.2f'),','])
					fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

					% Mean Ipsi hemi
					dat = imCBF((imFOV.*imLesion(:,:,:,5))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

					% Mean Ipsi hemi GM
					dat = imCBF((imFOV.*imLesion(:,:,:,5).*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

					% Mean Ipsi hemi GM-PVEC
					dat = imGMPVEC((imFOV.*imLesion(:,:,:,5).*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

					% Mean Contra hemi
					dat = imCBF((imFOV.*imLesion(:,:,:,6))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

					% Mean Contra hemi GM
					dat = imCBF((imFOV.*imLesion(:,:,:,6).*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

					% Mean Contra hemi GM-PVEC
					dat = imGMPVEC((imFOV.*imLesion(:,:,:,6).*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

					% Mean Ipsi hemi on global FOV
					dat = imCBF((imFOVmean.*imFOV.*imLesion(:,:,:,5))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

					% Mean Ipsi hemi GM on global FOV
					dat = imCBF((imFOVmean.*imFOV.*imLesion(:,:,:,5).*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

					% Mean Ipsi hemi GM-PVEC on global FOV
					dat = imGMPVEC((imFOVmean.*imFOV.*imLesion(:,:,:,5).*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

					% Mean Contra hemi on global FOV
					dat = imCBF((imFOVmean.*imFOV.*imLesion(:,:,:,6))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

					% Mean Contra hemi GM on global FOV
					dat = imCBF((imFOVmean.*imFOV.*imLesion(:,:,:,6).*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

					% Mean Contra hemi GM -PVEC on global FOV
					dat = imGMPVEC((imFOVmean.*imFOV.*imLesion(:,:,:,6).*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),'\n'])

				end
			case {2}
				% Subject name
				fprintf([x.SubjectNameList{ii} ',']);

				if isLesion
					maskLesionExcl = (1-imLesion(:,:,:,3))>0;
				else
					maskLesionExcl = (ones(size(imFOV)))>0;
				end

				% Mean Left hemi GM
				dat = imCBF((imFOV.*maskLeft.*maskLesionExcl.*(imGM>0.7))>0);
				dat(dat<=0) = NaN;
				fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

				% Mean Right hemi GM
				dat = imCBF((imFOV.*maskRight.*maskLesionExcl.*(imGM>0.7))>0);
				dat(dat<=0) = NaN;
				fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

				% Mean Left hemi GM global FOV
				dat = imCBF((imFOV.*imFOVmean.*maskLeft.*maskLesionExcl.*(imGM>0.7))>0);
				dat(dat<=0) = NaN;
				fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

				% Mean Right hemi GM global FOV
				dat = imCBF((imFOV.*imFOVmean.*maskRight.*maskLesionExcl.*(imGM>0.7))>0);
				dat(dat<=0) = NaN;
				fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])

				% Broca + individual FOV + L/R +-GM ++ size
				% Broca + individual FOV + L
				dat = imCBF((imFOV.*maskLeft.*maskBroca.*maskLesionExcl)>0);
				dat(dat<=0) = NaN;
				fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
				fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

				% Broca + individual FOV + R
				dat = imCBF((imFOV.*maskRight.*maskBroca.*maskLesionExcl)>0);
				dat(dat<=0) = NaN;
				fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
				fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

				% Broca + individual FOV + L + GM
				dat = imCBF((imFOV.*maskLeft.*maskBroca.*(imGM>0.7).*maskLesionExcl)>0);
				dat(dat<=0) = NaN;
				fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
				fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

				% Broca + individual FOV + R + GM
				dat = imCBF((imFOV.*maskLeft.*maskBroca.*(imGM>0.7).*maskLesionExcl)>0);
				dat(dat<=0) = NaN;
				fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
				fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),'\n'])
			case {3}
				% Subject name
				fprintf([x.SubjectNameList{ii} ',']);

				if isLesion
					maskLesionExcl = (1-imLesionOrig(:,:,:,1))>0;
				else
					maskLesionExcl = (ones(size(imFOV)))>0;
				end

				if isSuper
					maskLesionExcl = maskLesionExcl + (imSuper>0.5);
					maskLesionExcl = (maskLesionExcl>0);
				end

				for jj=1:11
					switch (jj)
						case {1} % Local FoV
							maskLoc = (imGM+imWM)>0.7;
						case {2} % Global Fov
							maskLoc = imFOVmean.*((imGM+imWM)>0.7);
						case {3} % Cerebellum
							maskLoc = (imMNI == 2);
						case {4} % 1. IFG (opercularis, triangularis, orbitalis)
							maskLoc = ((imHO == 5) + (imHO == 6))>0;
						case {5} % 2. Juxtapositional Lobule Cortex (formerly Supplementary Motor Cortex)
							maskLoc = ((imHO == 26))>0;
						case {6} % 3. Insula
							maskLoc = ((imHO == 2))>0;
						case {7} %4. STG (superior temporar gyrus) - AD+PD
							maskLoc = ((imHO == 9) + (imHO == 10))>0;
						case {8} % 5. MTG Middle Temporal gyrus - AD+PD+TO
							maskLoc = ((imHO == 11) + (imHO == 12) + (imHO==13))>0;
						case {9} %6. Supramarginal gyrus - AD+PD
							maskLoc = ((imHO == 19) + (imHO == 20))>0;
						case {10} % 7. Angular gyrus
							maskLoc = ((imHO == 21))>0;
						case {11} % 8. all
							maskLoc = (((imHO==2)+(imHO==5)+(imHO==6)+(imHO==26)+(imHO==2)+(imHO==9)+(imHO==10)+(imHO==11)+(imHO==12)+(imHO==13)+(imHO==19)+(imHO==20)+(imHO==21))>0);
					end
					dat = imCBF((maskLoc.*imFOV.*maskLeft.*maskLesionExcl)>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
					fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

					dat = imCBF((maskLoc.*imFOV.*maskRight.*maskLesionExcl)>0);
					dat(dat<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
					fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

					dat = imCBF((maskLoc.*imFOV.*maskLeft.*maskLesionExcl.*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					dat2 = imGMPVEC((maskLoc.*imFOV.*maskLeft.*maskLesionExcl.*(imGM>0.7))>0);
					dat2(dat2<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
					if jj<=2
						fprintf([num2str(xASL_stat_MeanNan(dat2),'%.2f'),','])
					end
					fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

					dat = imCBF((maskLoc.*imFOV.*maskRight.*maskLesionExcl.*(imGM>0.7))>0);
					dat(dat<=0) = NaN;
					dat2 = imGMPVEC((maskLoc.*imFOV.*maskRight.*maskLesionExcl.*(imGM>0.7))>0);
					dat2(dat2<=0) = NaN;
					fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
					if jj<=2
						fprintf([num2str(xASL_stat_MeanNan(dat2),'%.2f'),','])
					end
					fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])



				end

			case {4}
				if isLesion
					% Subject name
					fprintf([x.SubjectNameList{ii} ',']);

					if (sum(sum(sum(imLesion(:,:,:,1).*maskLeft))) > sum(sum(sum(imLesion(:,:,:,1).*maskRight))))
						maskIpsi = maskLeft;
						maskContra = maskRight;
					else
						maskIpsi = maskRight;
						maskContra = maskLeft;
					end

					maskLesionBilateral = imLesionOrig(:,:,:,1) + imLesionOrig(end:-1:1,:,:,1);
					maskLesionBilateral(maskLesionBilateral>0) = 1;
					maskLesionBilateralExcl = 1-(imLesionOrig(:,:,:,1) + imLesionOrig(end:-1:1,:,:,1));
					maskLesionBilateralExcl(maskLesionBilateralExcl<1) = 0; % in case of overlap of the two masks, we would get negative numbers
					maskPeriBilateralExcl   = 1-(imLesionOrig(:,:,:,3) + imLesionOrig(end:-1:1,:,:,3));
					maskPeriBilateralExcl(maskPeriBilateralExcl<1) = 0;

					% The mask of the language areas
					maskLang = (((imHO==2)+(imHO==5)+(imHO==6)+(imHO==26)+(imHO==2)+(imHO==9)+(imHO==10)+(imHO==11)+(imHO==12)+(imHO==13)+(imHO==19)+(imHO==20)+(imHO==21))>0);

					for jj=1:6
						switch (jj)
							case {1} % Local FoV, hemisphere, without lesion+peri
								maskLoc = ((imGM+imWM)>0.7).*maskPeriBilateralExcl;
							case {2} % Global Fov, hemisphere, without lesion+peri
								maskLoc = imFOVmean.*((imGM+imWM)>0.7).*maskPeriBilateralExcl;
							case {3}
								maskLoc = ((imGM+imWM+maskLesionBilateral)>0.7).*maskLang;
							case {4}
								maskLoc = ((imGM+imWM)>0.7).*maskLang.*maskLesionBilateralExcl;
							case {5}
								maskLoc = ((imGM+imWM)>0.7).*maskLang.*maskPeriBilateralExcl;
							case {6}
								maskLoc = maskLang;

						end
						if jj<6
							dat = imCBF((maskLoc.*imFOV.*maskIpsi)>0);
							dat(dat<=0) = NaN;
							fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
							fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

							dat = imCBF((maskLoc.*imFOV.*maskContra)>0);
							dat(dat<=0) = NaN;
							fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
							fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

							if jj~=3 % When calculating whole language area also within lesion then do not calculate GM>0.7
								dat = imCBF((maskLoc.*imFOV.*maskIpsi.*(imGM>0.7))>0);
								dat(dat<=0) = NaN;
								fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
								if jj<=2
									dat2 = imGMPVEC((maskLoc.*imFOV.*maskIpsi.*(imGM>0.7))>0);
									dat2(dat2<=0) = NaN;
									fprintf([num2str(xASL_stat_MeanNan(dat2),'%.2f'),','])
								end
								fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])

								dat = imCBF((maskLoc.*imFOV.*maskContra.*(imGM>0.7))>0);
								dat(dat<=0) = NaN;
								fprintf([num2str(xASL_stat_MeanNan(dat),'%.2f'),','])
								if jj<=2
									dat2 = imGMPVEC((maskLoc.*imFOV.*maskContra.*(imGM>0.7))>0);
									dat2(dat2<=0) = NaN;
									fprintf([num2str(xASL_stat_MeanNan(dat2),'%.2f'),','])
								end
								fprintf([num2str(xASL_stat_SumNan(dat>0)*1.5*1.5*1.5/1000,'%.2f'),','])
							end
						else
							dat = 0;
							if isSuper
								dat = imSuper((maskLoc.*imFOV.*maskIpsi)>0);
								dat = sum(sum(sum(dat)));
							end
							fprintf([num2str(dat,'%.2f'),','])
							dat = 0;
							if isSuper
								dat = imSuper((maskLoc.*imFOV.*maskContra)>0);
								dat = sum(sum(sum(dat)));
							end
							fprintf([num2str(dat,'%.2f'),','])
						end
					end
				end


		end
		fprintf('\n');
	end
end

% 1 local FoV hemi without lesion + peri, normal, GM, PVEC
% 2 global FoV hemi without lesion + peri, normal, GM, PVEC

% 3 language excluding nothing - everything

% 4 language excluding lesion - whole+GM
% 5 language excluding lesion+peri - whole+GM



% Hemi-GM-PVEC 32, 43,47 - asymmetry above 10%
