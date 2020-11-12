%% Load the DICOMs to the BIDS format and fill in all the missing info to the JSON
ExploreASL_Initialize([],false);

% FILLIN
baseDir = '/pet/projekte/asl/data/BIDS/';
basePath = [baseDir 'BIDS']; % The raw dicoms should be here in /StudyName/raw/...
parPath  = [baseDir 'BIDSpar']; % Put it your data_par.m file to load your parameters with name StudyName.m or StudyName.json (it converts all to JSON on the first run)
outputPath = [baseDir 'BIDS']; % The directory to save the DICOM2NII import
finalPath = [baseDir 'BIDSfinal']; % Takes files in NIFTI+JSON from outputPath and saves the complete BIDS format to finalPath
anonymPath = [baseDir 'BIDSanonymized']; % Takes files in NIFTI+JSON from outputPath and saves the complete BIDS format to finalPath

%% Load the list of the directories
fList = xASL_adm_GetFileList(basePath,[],false,[],true);

%% Go through all studies and import them
% This simply runs the ExploreASL_Import
for ii = 1:length(fList)
	% Import the whole session to JSON and NIFTI
	importStr{ii}.dirName = fList{ii};
	%ExploreASL_Import(imPar,false, true, false, true, false);
	switch(fList{ii})
		case 'Siemens_PCASL_3DGRASE_vascular'
			ExploreASL_ImportBIDS(fullfile(basePath,fList{ii}), [],[], [1 0 0], false, true, false, false);
			system(['rm ' basePath '/Siemens_PCASL_3DGRASE_vascular/analysis/Sub1/ASL_1/ASL4D_0170*']);
			system(['rm ' basePath '/Siemens_PCASL_3DGRASE_vascular/analysis/Sub1/ASL_1/ASL4D_2460*']);
			nii_files = xASL_adm_GetFileList([basePath '/Siemens_PCASL_3DGRASE_vascular/analysis/Sub1/ASL_1/'],'^*.nii$','FPList',[],false);
			nii_files = xASL_bids_MergeNifti(nii_files, 'ASL');
			ExploreASL_ImportBIDS(fullfile(basePath,fList{ii}), [],[], [0 1 0], false, true, false, false);
			
		case 'Philips_PCASL_3DGRASE_R5.4_TopUp'
			ExploreASL_ImportBIDS(fullfile(basePath,fList{ii}), [],[], [1 0 0], false, true, false, false);
			system(['mv ' basePath '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0_1.nii ' basePath '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0.nii']);
			system(['mv ' basePath '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0_1.json ' basePath '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0.json']);
			system(['mv ' basePath '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0_1_2.nii ' basePath '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0PERev.nii']);
			system(['mv ' basePath '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0_1_2.json ' basePath '/Philips_PCASL_3DGRASE_R5.4_TopUp/analysis/Sub1/ASL_1/M0PERev.json']);
			ExploreASL_ImportBIDS(fullfile(basePath,fList{ii}), [],[], [0 1 0], false, true, false, false);
			
		case 'Siemens_PCASL_volunteer'
			ExploreASL_ImportBIDS(fullfile(basePath,fList{ii}), [],[], [1 0 0], false, true, false, false);
			if xASL_exist([basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_NS.nii'])
				system(['rm ' basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_NS.json']);
				system(['mv ' basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_SS.json ' basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D.json']);
				imNS = xASL_io_Nifti2Im([basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_NS.nii']);
				imSS = xASL_io_Nifti2Im([basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_SS.nii']);
				imNS(:,:,:,2) = imSS;
				xASL_io_SaveNifti([basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_NS.nii'],[basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D.nii'],imNS/10,[],1,[]);
				system(['rm ' basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_NS.nii']);
				system(['rm ' basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/ASL4D_SS.nii']);
				system(['mv ' basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/M0_2.json ' basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/M0PERev.json']);
				system(['mv ' basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/M0_2.nii ' basePath '/Siemens_PCASL_volunteer/analysis/Sub1/ASL_1/M0PERev.nii']);
			end
			ExploreASL_ImportBIDS(fullfile(basePath,fList{ii}), [],[], [0 1 0], false, true, false, false);
			
		case 'Siemens_PASL_multiTI'
			ExploreASL_ImportBIDS(fullfile(basePath,fList{ii}), [],[], [1 0 0], false, true, false, false);
			if xASL_exist([basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_300.nii'])
				mTIvec = [300,600,900,1200,1500,1800,2100,2400,2700,3000];
				for jj = 1:length(mTIvec)
					if jj>1
						system(['rm ' basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(jj)) '.json']);
						system(['rm ' basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_SS_' num2str(mTIvec(jj)) '.json']);
						imNSSS(:,:,:,2*(jj-1)+1) = xASL_io_Nifti2Im([basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(jj)) '.nii']);
						imNSSS(:,:,:,2*(jj-1)+2) = xASL_io_Nifti2Im([basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_SS_' num2str(mTIvec(jj)) '.nii']);
					else
						system(['mv ' basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(jj)) '.json ' basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D.json']);
						system(['rm ' basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_SS_' num2str(mTIvec(jj)) '.json']);
						imNSSS = xASL_io_Nifti2Im([basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(jj)) '.nii']);
						imNSSS(:,:,:,2) = xASL_io_Nifti2Im([basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_SS_' num2str(mTIvec(jj)) '.nii']);
					end
				end
				xASL_io_SaveNifti([basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(1)) '.nii'],[basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D.nii'],imNSSS/10,[],1,[]);
				for jj = 1:length(mTIvec)
					system(['rm ' basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(jj)) '.nii']);
					system(['rm ' basePath '/Siemens_PASL_multiTI/analysis/Sub1/ASL_1/ASL4D_SS_' num2str(mTIvec(jj)) '.nii']);
				end
			end
			ExploreASL_ImportBIDS(fullfile(basePath,fList{ii}), [],[], [0 1 0], false, true, false, false);
			
		otherwise
			ExploreASL_ImportBIDS(fullfile(basePath,fList{ii}), [],[], [1 1 0], false, true, false, false);
	end
end











%% Specify the missing study parameters
% This is the most important part - it takes the imported JSONs and fills in all the study specific information that is currently still missing.
% Several parameters are automatically retrieved from the DICOMs, some are retrieved from Data_Par.m just below. And the rest needs to be filled in.

% This still needs to be converted to files
for ii = 1:length(fList)
	switch (fList{ii})
		case 'Siemens_PCASL_volunteer'
			importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.LabelingPulseInterval = 0.4/1000;
			importStr{ii}.par.LabelingPulsesFlipAngle = 25;
			importStr{ii}.par.NumberSegments = 2;
			importStr{ii}.par.TotalAcquiredVolumes = [2 2];
			importStr{ii}.par.TotalReadoutTime = 0.0104;
			importStr{ii}.par.BackgroundSuppressionPulseTime = [0.85 0.1];
			importStr{ii}.par.BackgroundSuppressionNumberPulses = 2;

		case 'Siemens_PASL_multiTI'
			for cc = 1:10,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.NumberSegments = 2;
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 0;
			importStr{ii}.par.BackgroundSuppressionPulseTime = [0.85 0.1];
			importStr{ii}.par.BackgroundSuppressionNumberPulses = 2;

		case 'Philips_PCASL_3DGRASE_functional'
			importStr{ii}.par.Units = 'mL/100g/min';
			%importStr{ii}.par.ASLContext = bidsPar.strCbf;
			importStr{ii}.par.ASLContext = sprintf('%s\n',bidsPar.strCbf);
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.NumberSegments = 5;

		case 'Siemens_PASL_singleTI'
			importStr{ii}.par.ASLContext = sprintf('%s\n',bidsPar.strM0scan);
			for cc = 1:45,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 0.9;
			importStr{ii}.par.BolusCutOffTechnique = 'Q2TIPS';

		case 'Siemens_PCASL_2DEPI_AD'
			%importStr{ii}.par.ASLContext = '(Label+Control)*23';
			for cc = 1:46, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Siemens_PASL_3DGRASE_AD'
			importStr{ii}.par.LabelingType = 'PASL';
			%importStr{ii}.par.ASLContext = bidsPar.strDeltaM;
			importStr{ii}.par.ASLContext = sprintf('%s\n',bidsPar.strDeltaM);
			importStr{ii}.par.BolusCutOffFlag = false;
			
		case 'Siemens_PCASL_3DGRASE_AD'
			%importStr{ii}.par.ASLContext = 'M0+((Label+Control)*12)';
			importStr{ii}.par.ASLContext = sprintf('%s\n',bidsPar.strM0scan);
			for cc = 1:12, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Siemens_PCASL_2DEPI_pharma'
			%importStr{ii}.par.ASLContext = '(Label+Control)*44';
			for cc = 1:44, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Siemens_PASL_3DGRASE_AD2'
			%importStr{ii}.par.ASLContext = 'Label+Control';
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl);
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.LabelingLocationDescription = 'Labeling with FAIR';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 0;
			importStr{ii}.par.BolusCutOffTechnique = 'QUIPSSII';

		case 'Philips_PCASL_3DGRASE_volunteer'
			%importStr{ii}.par.ASLContext = '(M0*2)+((Label+Control)*7)';
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',bidsPar.strM0scan,bidsPar.strM0scan);
			for cc = 1:7, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Siemens_PCASL_3DGRASE_volunteer'
			%importStr{ii}.par.ASLContext = 'M0+((Label+Control)*15)';
			importStr{ii}.par.ASLContext = sprintf('%s\n',bidsPar.strM0scan);
			for cc = 1:15, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Siemens_PCASL_3DGRASE_vascular'
			%importStr{ii}.par.ASLContext = 'DeltaM*15';
			for cc = 1:90, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n',bidsPar.strDeltaM)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Philips_PCASL_2DEPI_pharma'
			%importStr{ii}.par.ASLContext = '(Label+Control)*38';
			for cc = 1:38, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Philips_PCASL_2DEPI_pharma2'
			%importStr{ii}.par.ASLContext = '(Label+Control)*32';
			for cc = 1:32, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case {'Philips_PCASL_2DEPI_Intera_volunteer','Philips_PCASL_2DEPI_Ingenia_volunteer'}
			%importStr{ii}.par.ASLContext = '(Label+Control)*75';
			for cc = 1:75, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case {'Philips_PCASL_2DEPI_dummyLL' 'Philips_PCASL_2DEPI_dummyMultiPLD' 'Philips_PCASL_2DEPI_dummyQUASAR'}
			importStr{ii}.par.LabelingType = 'PCASL';

		case {'Philips_PCASL_2DEPI_Bsup_AD1','Philips_PCASL_2DEPI_noBsup_AD'}
			%importStr{ii}.par.ASLContext = '(Label+Control)*40';
			for cc = 1:40, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case {'Philips_PCASL_2DEPI_Bsup_AD2','Philips_PCASL_2DEPI_Bsup_AD3'}
			%importStr{ii}.par.ASLContext = '(Label+Control)*30';
			for cc = 1:30, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Philips_PCASL_2DEPI_Bsup_AD4'
			importStr{ii}.par.LabelingType = 'PCASL';
			%importStr{ii}.par.ASLContext = bidsPar.strDeltaM;
			importStr{ii}.par.ASLContext = sprintf('%s\n',bidsPar.strDeltaM);
			importStr{ii}.par.LabelingEfficiency = 0.83;

		case {'Siemens_PASL_2DEPI_noBsup_AD2','Siemens_PASL_2DEPI_noBsup_AD'}
			%importStr{ii}.par.ASLContext = '(Label+Control)*31';(31*L+C)+L
			for cc = 1:31, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n',bidsPar.strLabel)];
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 0.6;
			importStr{ii}.par.BolusCutOffTechnique = 'Q2TIPS';
			importStr{ii}.par.LabelingSlabThickness = 80;

		case {'Siemens_PASL_3DGRASE_Prisma_Bsup_AD1','Siemens_PASL_3DGRASE_Prisma_Bsup_AD2'}
			%importStr{ii}.par.ASLContext = '(Label+Control)*10';
			for cc = 1:10, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 0;
			importStr{ii}.par.BolusCutOffTechnique = 'QUIPSSII';

		case 'Siemens_PASL_3DGRASE_Prisma_Bsup_AD3'
			%importStr{ii}.par.ASLContext = '(Label+Control)*2';
			for cc = 1:2, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = [0 200 400]/1000;
			importStr{ii}.par.BolusCutOffTechnique = 'QUIPSS';
			importStr{ii}.par.LabelingSlabThickness = 60;

		case {'Philips_PCASL_3DGRASE_R5.4_TopUp'}
			importStr{ii}.par.Units = 'mL/100g/min';
			%importStr{ii}.par.ASLContext = bidsPar.strCbf;
			importStr{ii}.par.ASLContext = sprintf('%s\n',bidsPar.strCbf);
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Philips_PCASL_2DEPI_volunteer3'
			for cc = 1:35, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strControl,bidsPar.strLabel)];end
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.AcquisitionVoxelSize = [3.75 3.75 5];

		case 'Siemens_PCASL_3DGRASE_volunteer2'
			for cc = 1:8, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strControl,bidsPar.strLabel)];end
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.AcquisitionVoxelSize = [3.4 3.4 4];

		case {'Philips_PCASL_2DEPI_glioma','Philips_PCASL_2DEPI_GBM'}
			%importStr{ii}.par.ASLContext = '(Control+Label)*30';
			for cc = 1:30, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strControl,bidsPar.strLabel)];end
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.LabelingPulseAverageGradient = 0.7;
			importStr{ii}.par.LabelingPulseMaximumGradient = 7;
			importStr{ii}.par.LabelingPulseFlipAngle = 23;
			importStr{ii}.par.LabelingPulseDuration = 0.0005;
			importStr{ii}.par.PCASLType = 'balanced';

		case {'Philips_PCASL_2DEPI_volunteer_1','Philips_PCASL_2DEPI_volunteer_2','Siemens_PCASL_2DEPI_volunteer'}
			%importStr{ii}.par.ASLContext = '(Control+Label)*70';
			for cc = 1:70, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strControl,bidsPar.strLabel)];end
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.LabelingLocationDescription = 'Fixed, 9 cm below ACPC';
			importStr{ii}.par.LabelingPulseAverageGradient = 0.6;
			importStr{ii}.par.LabelingPulseMaximumGradient = 6;
			importStr{ii}.par.LabelingPulseFlipAngle = 25;
			importStr{ii}.par.LabelingPulseInterval = 0.00124;
			importStr{ii}.par.PCASLType = 'balanced';

	end
	switch (fList{ii})
		case {'Philips_PCASL_2DEPI_volunteer_2','Siemens_PCASL_2DEPI_volunteer'}
			importStr{ii}.par.LabelingPulseInterval = 0.00115;
		case 'Philips_PCASL_3DGRASE_R5.4_TopUp'
			importStr{ii}.par.TotalReadoutTime = 0.0136;
	end

	%Check for Manufacturer and MRAcquisitionType in preference 
	%"Manufacturer" and "MRAcquisitionType"
	%in preference in the studyPar, but also in the DICOMs-JSONouput
	


	if ~isfield(importStr{ii}.par,'ReadoutSegments') && isfield(importStr{ii}.x,'NumberSegments')
		importStr{ii}.par.NumberSegments = importStr{ii}.x.NumberSegments;
	end

	% This can all go, will be in files
% 	% Labeling delays and durations
% 	if strcmpi(importStr{ii}.par.LabelingType,'PASL')
% 		%importStr{ii}.par.LabelingDuration = 0;% importStr{ii}.x.LabelingDuration           = 1.800;  % for PASL this is TI1
% 		importStr{ii}.par.PostLabelingDelay = importStr{ii}.x.InitialPostLabelDelay;
% 		if importStr{ii}.par.BolusCutOffFlag
% 			importStr{ii}.par.BolusCutOffDelayTime = importStr{ii}.par.BolusCutOffDelayTime + importStr{ii}.x.LabelingDuration;
% 		end
% 	else
% 		importStr{ii}.par.LabelingDuration = importStr{ii}.x.LabelingDuration;
% 		importStr{ii}.par.PostLabelingDelay = importStr{ii}.x.InitialPostLabelDelay;
% 	end

	
	
	end

	% Last round of edits - most parameters are filled in above, or (Background suppresion timing) prefilled with default values
	% But after this automatic prefilling, you might want to change a few parameters - this is done here
	% FILLIN
	% Either to change the automatically filled things above, or to supply further info about multi-PLD, vascular crushing, QUASAR etc.
	switch (fList{ii})
		case 'Siemens_PASL_multiTI'
			importStr{ii}.par.PostLabelingDelay = [300 300 600 600 900 900 1200 1200 1500 1500 1800 1800 2100 2100 2400 2400 2700 2700 3000 3000]/1000;

		case 'Siemens_PASL_singleTI'
			importStr{ii}.par.VascularCrushing = true;
			importStr{ii}.par.VascularCrushingVenc = 100;

		case 'Siemens_PCASL_3DGRASE_AD'
			importStr{ii}.par.LabelingDuration = [0 repmat(1800,[1,24])]/1000;
			importStr{ii}.par.VascularCrushing = true;
			importStr{ii}.par.VascularCrushingVenc = 10;

		case 'Philips_PCASL_2DEPI_dummyLL'
			%importStr{ii}.par.ASLContext = '(Label*15+Control*15)*5';
			for cc = 1:5
				for dd=1:15,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n',bidsPar.strLabel)];end
				for dd=1:15,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n',bidsPar.strControl)];end
			end
			importStr{ii}.par.FlipAngle = 25;
			importStr{ii}.par.LookLocker = true;
			importStr{ii}.par.PostLabelingDelay = repmat(250:250:3750,[1 10])/1000;

		case 'Philips_PCASL_2DEPI_dummyMultiPLD'
			%importStr{ii}.par.ASLContext = '(Label+Control)*75';
			for cc = 1:75, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strLabel,bidsPar.strControl)];end
			importStr{ii}.par.PostLabelingDelay = repmat([500 500 1000 1000 1500 1500 1800 1800 2200 2200],[1 15])/1000;

		case 'Philips_PCASL_2DEPI_dummyQUASAR'
			%importStr{ii}.par.ASLContext = '(Label*15+Control*15)*5';
			for cc = 1:5
				for dd=1:15,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n',bidsPar.strLabel)];end
				for dd=1:15,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n',bidsPar.strControl)];end
			end
			importStr{ii}.par.FlipAngle = [25*ones(1,120),12*ones(1,30)];
			importStr{ii}.par.LookLocker = true;
			importStr{ii}.par.VascularCrushing = true;
			importStr{ii}.par.VascularCrushingVenc = [zeros(1,60),10*ones(1,60),zeros(1,3)];
			importStr{ii}.par.PostLabelingDelay = repmat(250:250:3750,[1 10])/1000;

	end
%end




%% Export the fully anonymized datasets for public sharing
pthVec = {'GE_PCASL_3Dspiral_volunteer' 'Siemens_PCASL_3DGRASE_volunteer2' 'Philips_PCASL_2DEPI_volunteer3'};
for ii = 1:3
	xASL_Copy(fullfile(finalPath,pthVec{ii}),fullfile(anonymPath,pthVec{ii}));
	xASL_spm_deface(fullfile(anonymPath,pthVec{ii},'sub-Sub103','anat','sub-Sub103_T1w.nii'),true);
	gzip(fullfile(anonymPath,pthVec{ii},'sub-Sub103','anat','sub-Sub103_T1w.nii'));
	delete(fullfile(anonymPath,pthVec{ii},'sub-Sub103','anat','sub-Sub103_T1w.nii'));
end

pthVec = {'Siemens_PASL_multiTI','Siemens_PASL_singleTI','Siemens_PCASL_volunteer'};
for ii = 1:3
	xASL_Copy(fullfile(finalPath,pthVec{ii}),fullfile(anonymPath,pthVec{ii}));
	xASL_spm_deface(fullfile(anonymPath,pthVec{ii},'sub-Sub1','anat','sub-Sub1_T1w.nii'),true);
	gzip(fullfile(anonymPath,pthVec{ii},'sub-Sub1','anat','sub-Sub1_T1w.nii'));
	delete(fullfile(anonymPath,pthVec{ii},'sub-Sub1','anat','sub-Sub1_T1w.nii'));
end

%%
testSet = {'GE_PCASL_2DEPI_volunteer','GE_PCASL_3Dspiral_Product_pharma','GE_PCASL_3Dspiral_Product_volunteer','GE_PCASL_3Dspiral_volunteer','GE_PCASL_3Dspiral_WIP_pharma','GE_PCASL_3Dspiral_WIP_volunteer'};

for iTest = 1:length(testSet)
	[i,r] = xASL_bids_CompareStructures(['/pet/projekte/asl/data/BIDS/BIDS/' testSet{iTest} '/bids'],['/home/janpetr/tmp/comp/old/' testSet{iTest}]);
end
