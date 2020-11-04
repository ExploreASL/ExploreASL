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
for ii = 1:length(fList)
	

	% Defaults to be overwritten
	importStr{ii}.par = [];
	importStr{ii}.par.VascularCrushing = false;
	importStr{ii}.par.LabelingLocationDescription = 'Random description';
	importStr{ii}.par.LabelingDistance = 40;

	importStr{ii}.par.ASLContext = '';
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

		case 'GE_PCASL_3Dspiral_Product_pharma'
			importStr{ii}.par.Units = 'mL/100g/min';
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',bidsPar.strCbf,bidsPar.strM0scan);
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'GE_PCASL_3Dspiral_Product_volunteer'
			importStr{ii}.par.Units = 'mL/100g/min';
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',bidsPar.strM0scan,bidsPar.strCbf);
			importStr{ii}.par.LabelingType = 'PCASL';

		case {'GE_PCASL_3Dspiral_WIP_pharma',...
			  'Philips_PCASL_3DGRASE_R5.4_TopUp'}
			importStr{ii}.par.Units = 'mL/100g/min';
			%importStr{ii}.par.ASLContext = bidsPar.strCbf;
			importStr{ii}.par.ASLContext = sprintf('%s\n',bidsPar.strCbf);
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'GE_PCASL_3Dspiral_WIP_volunteer'
			importStr{ii}.par.Units = 'mL/100g/min';
			%importStr{ii}.par.ASLContext = 'CBF+M0';
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',bidsPar.strCbf,bidsPar.strM0scan);
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.AcquisitionVoxelSize = [4 4 8];
			importStr{ii}.par.LabelingPulseAverageGradient = 0.6;
			importStr{ii}.par.LabelingPulseMaximumGradient = 6;
			importStr{ii}.par.LabelingPulseFlipAngle = 18;
			importStr{ii}.par.LabelingPulseDuration = 0.0005;
			importStr{ii}.par.PCASLType = 'balanced';

		case 'GE_PCASL_3Dspiral_volunteer'
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',bidsPar.strM0scan,bidsPar.strDeltaM);
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.AcquisitionVoxelSize = [4 4 8];

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

		case {'GE_PCASL_2DEPI_volunteer','Philips_PCASL_2DEPI_volunteer_1','Philips_PCASL_2DEPI_volunteer_2','Siemens_PCASL_2DEPI_volunteer'}
			%importStr{ii}.par.ASLContext = '(Control+Label)*70';
			for cc = 1:70, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',bidsPar.strControl,bidsPar.strLabel)];end
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.LabelingLocationDescription = 'Fixed, 9 cm below ACPC';
			importStr{ii}.par.LabelingPulseAverageGradient = 0.6;
			importStr{ii}.par.LabelingPulseMaximumGradient = 6;
			importStr{ii}.par.LabelingPulseFlipAngle = 25;
			%importStr{ii}.par.AcquisitionDuration = 672;
			importStr{ii}.par.LabelingPulseInterval = 0.00124;
			importStr{ii}.par.PCASLType = 'balanced';

	end
	switch (fList{ii})
		case {'Philips_PCASL_2DEPI_volunteer_2','Siemens_PCASL_2DEPI_volunteer'}
			%importStr{ii}.par.AcquisitionDuration = 658;
			importStr{ii}.par.LabelingPulseInterval = 0.00115;
		case 'Philips_PCASL_3DGRASE_R5.4_TopUp'
			importStr{ii}.par.TotalReadoutTime = 0.0136;
	end

	% Process all the data and automatically fill in the missing parameters
	if strcmpi(importStr{ii}.x.MRAcquisitionType,'2D')
		importStr{ii}.par.PulseSequenceType = '2D_EPI';
	else
		if strcmpi(importStr{ii}.x.Manufacturer,'GE') || strcmpi(importStr{ii}.x.Manufacturer,'GE_WIP') || strcmpi(importStr{ii}.x.Manufacturer,'GE_product')
			importStr{ii}.par.PulseSequenceType = '3D_spiral';
		else
			importStr{ii}.par.PulseSequenceType = '3D_GRASE';
		end
	end

	%if ~isfield(importStr{ii}.par,'TotalAcquiredVolumes') && isfield(importStr{ii}.x,'NumberOfAverages') && (importStr{ii}.x.NumberOfAverages > 1)
	%	importStr{ii}.par.TotalAcquiredVolumes = importStr{ii}.x.NumberOfAverages;
	%end

	if ~isfield(importStr{ii}.par,'ReadoutSegments') && isfield(importStr{ii}.x,'NumberSegments')
		importStr{ii}.par.NumberSegments = importStr{ii}.x.NumberSegments;
	end

	% Labeling delays and durations
	if strcmpi(importStr{ii}.par.LabelingType,'PASL')
		%importStr{ii}.par.LabelingDuration = 0;% importStr{ii}.x.LabelingDuration           = 1.800;  % for PASL this is TI1
		importStr{ii}.par.PostLabelingDelay = importStr{ii}.x.InitialPostLabelDelay;
		if importStr{ii}.par.BolusCutOffFlag
			importStr{ii}.par.BolusCutOffDelayTime = importStr{ii}.par.BolusCutOffDelayTime + importStr{ii}.x.LabelingDuration;
		end
	else
		importStr{ii}.par.LabelingDuration = importStr{ii}.x.LabelingDuration;
		importStr{ii}.par.PostLabelingDelay = importStr{ii}.x.InitialPostLabelDelay;
	end

	if importStr{ii}.x.BackGrSupprPulses == 0
		importStr{ii}.par.BackgroundSuppression = false;
	else
		importStr{ii}.par.BackgroundSuppression = true;
		if ~isfield(importStr{ii}.par,'BackgroundSuppressionPulseTime') || isempty(importStr{ii}.par.BackgroundSuppressionPulseTime)
			switch (importStr{ii}.x.BackGrSupprPulses)
				case 2
					if importStr{ii}.par.PostLabelingDelay > 1.750
						importStr{ii}.par.BackgroundSuppressionPulseTime = [1.75 0.524];
						importStr{ii}.par.BackgroundSuppressionNumberPulses = 2;
					elseif importStr{ii}.par.PostLabelingDelay > 1.495
						importStr{ii}.par.BackgroundSuppressionPulseTime = [1.495 0.345];
						importStr{ii}.par.BackgroundSuppressionNumberPulses = 2;
					elseif importStr{ii}.par.PostLabelingDelay > 1.195
						warning('Backgrou suppresion not properly calculated');
						importStr{ii}.par.BackgroundSuppressionPulseTime = [1.195 0.245];
						importStr{ii}.par.BackgroundSuppressionNumberPulses = 2;
					else
						error('Pulses not fitting');
					end
				case 4
					if importStr{ii}.par.PostLabelingDelay > 1.510
						importStr{ii}.par.BackgroundSuppressionPulseTime = [1.510 0.875 0.375 0.095];
						importStr{ii}.par.BackgroundSuppressionNumberPulses = 4;
					else
						error('Pulses not fitting');
					end
				case 5
					if importStr{ii}.par.PostLabelingDelay > 1.510
						importStr{ii}.par.BackgroundSuppressionPulseTime = [(importStr{ii}.par.PostLabelingDelay+importStr{ii}.par.LabelingDuration+1) 1.510 0.875 0.375 0.095];
						importStr{ii}.par.BackgroundSuppressionNumberPulses = 5;
					else
						error('Pulses not fitting');
					end
				otherwise
					error('Unknown number of pulses');
			end
		else
			importStr{ii}.par.BackgroundSuppressionNumberPulses = length(importStr{ii}.par.BackgroundSuppressionPulseTime);
		end

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
end

% Create the structure with the correct field order
fieldOrderStruct = [];
for ii=1:length(bidsPar.listFieldOrder)
	fieldOrderStruct.(bidsPar.listFieldOrder{ii}) = '';
end

%% Go through all studies and check all the M0 and ASLs and modify the JSONs
% This step should be completely automatic, just taking the info filled above and using it to convert to full BIDS.
% No edits should be necessary unless there's something wrong that needs to be fixed.
for ii = 1:length(fList)
	% Make a copy of the par to the Flavors
	importStr{ii}.flavors = importStr{ii}.par;



	% Go through all subjects
	fSubs = xASL_adm_GetFileList(fullfile(outputPath,importStr{ii}.dirName,'analysis'),[],false,[],true);
	for jj = 1:length(fSubs)

		subLabel = xASL_adm_CorrectName(fSubs{jj},2);

		% Make a subject directory
		if ~exist(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel]),'dir')
			mkdir(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel]));
		end

		% Go throught the list of anat files
		for iiAnat = bidsPar.listAnatTypes
			% Check if it exists
			anatPath = '';
			if xASL_exist(fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj},[iiAnat{1},'.nii']),'file')
				anatPath = fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj},iiAnat{1});
			end

			if xASL_exist(fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj},[iiAnat{1} '_1'],[iiAnat{1},'.nii']),'file')
				anatPath = fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj},[iiAnat{1} '_1'],iiAnat{1});
			end

			if ~isempty(anatPath)

				if ~exist(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel],'anat'),'dir')
					mkdir(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel],'anat'));
				end

				% If yes, then copy the file
				xASL_Copy([anatPath '.nii'],fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel],'anat',...
					['sub-' subLabel '_' iiAnat{1} '.nii.gz']));

				% Load the JSON
				jsonAnat = spm_jsonread([anatPath,'.json']);
				jsonLocal = [];
				% Discard the fields
				% Copy all dicom ones
				for fn = fieldnames(jsonAnat)'
					% Fields to skip
					bCP = 1;
					for ll=1:length(bidsPar.listFieldsRemoveGeneral)
						if strcmp(bidsPar.listFieldsRemoveGeneral{ll},fn{1})
							bCP = 0;
						end
					end
					if bCP
						jsonLocal.(fn{1}) = jsonAnat.(fn{1});
					end
				end

				% Save the JSON
				jsonLocal = finalJsonCheck(jsonLocal,fieldOrderStruct,bidsPar.listRemoveIfEmpty);
				spm_jsonwrite(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel],'anat',['sub-' subLabel '_' iiAnat{1} '.json']),jsonLocal);
			end
		end

		fSes = xASL_adm_GetFileList(fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj}),'^ASL.+$',false,[],true);

		% Go through all sessions
		for kk = 1:length(fSes)


			% Make a subject directory
			if length(fSes)>1
				sesLabel = ['ses-' fSes{kk}(5:end)];
				sesLabelUnd = ['_' sesLabel];
				if ~exist(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel],sesLabel),'dir')
					mkdir(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel],sesLabel));
					mkdir(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel],sesLabel,'asl'));
				end
				inSesPath = fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj},fSes{kk});
				outSesPath = fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel],sesLabel);
			else
				sesLabel = '';
				sesLabelUnd = '';

				% Only one session - no session labeling
				if ~exist(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel]),'dir')
					mkdir(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel]));
					mkdir(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel],'perf'));
				end
				inSesPath = fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj},fSes{kk});
				outSesPath = fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel]);
			end

			% Check if there are multiple runs per session
			fRuns = xASL_adm_GetFileList(inSesPath,'^ASL4D_\d.nii+$',false,[],false);
			nSes = length(fRuns);

			for mm = 1:(max(nSes,1))
				if nSes
					aslLabel = ['ASL4D_' num2str(mm)];
					aslOutLabel = fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd '_run-' num2str(mm)]);
					aslOutLabelRelative = fullfile('perf',['sub-' subLabel sesLabelUnd '_run-' num2str(mm)]);
				else
					aslLabel = 'ASL4D';
					aslOutLabel = fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd]);
					aslOutLabelRelative = fullfile('perf',['sub-' subLabel sesLabelUnd]);
				end

				% Load the JSON
				%jsonLocal = xASL_import_json(fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj},fSes{kk},'ASL4D.json'));
				jsonDicom = spm_jsonread(fullfile(inSesPath,[aslLabel '.json']));
				if exist(fullfile(inSesPath,[aslLabel '_parms.mat']),'file')
					imParms = load(fullfile(inSesPath,[aslLabel '_parms.mat']));
				else
					imParms = [];
				end
				imNii = xASL_io_Nifti2Im(fullfile(inSesPath,[aslLabel '.nii']));

				rescaleParms = [];
				ParmsFields = {'RescaleSlope' 'RWVSlope'    'MRScaleSlope' 'RescaleIntercept'...
					'RescaleSlopeOriginal' 'RescaleSlope' 'MRScaleSlope' 'UsePhilipsFloatNotDisplayScaling' 'RWVSlope'};
				JSONFields  = {'PhilipsRescaleSlope'  'PhilipsRWVSlope' 'PhilipsScaleSlope' 'PhilipsRescaleIntercept'...
					'RescaleSlopeOriginal' 'RescaleSlope' 'MRScaleSlope' 'UsePhilipsFloatNotDisplayScaling' 'RWVSlope'};

				for pp = 1:length(ParmsFields)
					if isfield(jsonDicom,JSONFields{pp})
						rescaleParms.(ParmsFields{pp}) = jsonDicom.(JSONFields{pp});
					end
					if ~isempty(imParms)
						if isfield(imParms.parms,ParmsFields{pp})
							rescaleParms.(ParmsFields{pp}) = imParms.parms.(ParmsFields{pp});
						end
					end
				end
				if ~isempty(strfind(jsonDicom.Manufacturer,'Philips')) || ~isempty(strfind(jsonDicom.Manufacturer,'philips'))
					scaleFactor = xASL_adm_GetPhilipsScaling(rescaleParms,xASL_io_ReadNifti(fullfile(inSesPath,[aslLabel '.nii'])));
				else
					scaleFactor = 0;
				end

				if scaleFactor
					imNii = imNii .* scaleFactor;
					xASL_io_SaveNifti(fullfile(inSesPath,[aslLabel '.nii']),[aslOutLabel '_asl.nii.gz'],imNii,[],1,[]);
				elseif size(imNii,4) == 1
					% The fourth dimension is 1, so we have to write the file again, to make sure the
					xASL_io_SaveNifti(fullfile(inSesPath,[aslLabel '.nii']),[aslOutLabel '_asl.nii.gz'],imNii,[],1,[]);
				else
					% Copy the ASL
					xASL_Copy(fullfile(inSesPath,[aslLabel '.nii']),[aslOutLabel '_asl.nii.gz']);
				end


				% Copy the basic ones
				jsonLocal = importStr{ii}.par;

				% Copy all dicom ones
				for fn = fieldnames(jsonDicom)'
					% Fields to skip
					bCP = 1;
					for ll=1:length(bidsPar.listFieldsRemoveGeneral)
						if strcmp(bidsPar.listFieldsRemoveGeneral{ll},fn{1})
							bCP = 0;
						end
					end
					for ll=1:length(bidsPar.listFieldsRemoveASL)
						if strcmp(bidsPar.listFieldsRemoveASL{ll},fn{1})
							bCP = 0;
						end
					end
					if bCP
						jsonLocal.(fn{1}) = jsonDicom.(fn{1});
					end
				end

				% Check if BolusDuration field is present and not in conflict with the BolusCutoffDelayTime
				if isfield(jsonDicom,'BolusDuration')
					if ~isfield(importStr{ii}.par,'BolusCutOffDelayTime')
						warning('Bolus duration obtained from DICOM, but not correctly redefined.');
					elseif ~isequal(jsonDicom.BolusDuration,importStr{ii}.par.BolusCutOffDelayTime(1))
						warning('Bolus duration obtained from DICOM and the manuall defined one differ.');
					end
				end

				%jsonLocal.EchoTime = jsonDicom.EchoTime;
				%jsonLocal.MagneticFieldStrength = jsonDicom.MagneticFieldStrength;
				%jsonLocal.RepetitionTime = jsonDicom.RepetitionTime;
				%jsonLocal.Manufacturer = jsonDicom.Manufacturer;
				%jsonLocal.FlipAngle = jsonDicom.FlipAngle;

				% Free info about the sequence, now just the scanner type+software
				if isfield(jsonDicom,'ManufacturersModelName')
					jsonLocal.PulseSequenceDetails = jsonDicom.ManufacturersModelName;
				%else
				%	jsonLocal.PulseSequenceDetails = '';
				end
				if isfield(jsonDicom,'SoftwareVersions')
					if ~isempty(jsonLocal.PulseSequenceDetails)
						jsonLocal.PulseSequenceDetails = [jsonLocal.PulseSequenceDetails '-'];
					end
					jsonLocal.PulseSequenceDetails = [jsonLocal.PulseSequenceDetails jsonDicom.SoftwareVersions];
					%jsonLocal.MRSoftwareVersion = jsonDicom.SoftwareVersions;
				end

				% Fill in extra parameters based on the JSON from the data
				if importStr{ii}.par.PulseSequenceType(1) == '2'
					jsonLocal.SliceTiming = ((0:(size(imNii,3)-1))')*importStr{ii}.x.SliceReadoutTime;
				else
					if isfield(jsonLocal,'SliceTiming')
						jsonLocal = rmfield(jsonLocal,'SliceTiming');
					end
				end

				if isfield(jsonLocal,'EffectiveEchoSpacing')
					if jsonLocal.EffectiveEchoSpacing == 0
						jsonLocal = rmfield(jsonLocal,'EffectiveEchoSpacing');
					else
						jsonLocal.EffectiveEchoSpacing = abs(jsonLocal.EffectiveEchoSpacing);
					end
				end

				if isfield(jsonLocal,'TotalReadoutTime')
					if jsonLocal.TotalReadoutTime == 0
						jsonLocal = rmfield(jsonLocal,'TotalReadoutTime');
					else
						jsonLocal.TotalReadoutTime = abs(jsonLocal.TotalReadoutTime);
					end
				end

				if isfield(importStr{ii}.x,'RepetitionTime')
					jsonLocal.RepetitionTime = importStr{ii}.x.RepetitionTime;
				end

				% Check if TR is a vector - replace by the maximum then
				if length(jsonLocal.RepetitionTime) > 1
					jsonLocal.RepetitionTime = max(jsonLocal.RepetitionTime);
					warning('TR was a vector. Taking the maximum only.');
				end

				% Fill in the number of averages
				%ppStr = importStr{ii}.dirName;
				%if isfield(importStr{ii}.par,'TotalAcquiredVolumes')
				%	ppStr = [ppStr ' -' num2str(max(importStr{ii}.par.TotalAcquiredVolumes(:))) '-'];
				%else
				%	ppStr = [ppStr ' -.-'];
				%end

				%if isfield(imParms,'parms') && isfield(imParms.parms, 'NumberOfAverages')  && (max(imParms.parms.NumberOfAverages) > 1)
				%	ppStr = [ppStr ' -' num2str(max(imParms.parms.NumberOfAverages)) '-'];
				%else
				%	ppStr = [ppStr ' -.-'];
				%end

				% Import the number of averages
				if isfield(imParms,'parms') && isfield(imParms.parms,'NumberOfAverages') && (max(imParms.parms.NumberOfAverages) > 1)
					if isfield(importStr{ii}.par,'TotalAcquiredVolumes')
						if max(imParms.parms.NumberOfAverages) ~= importStr{ii}.par.TotalAcquiredVolumes
							warning('Discrepancy in the number of averages');
						end
					else
						%importStr{ii}.par.TotalAcquiredVolumes = max(imParms.parms.NumberOfAverages);
					end
				end

				%if isfield(importStr{ii}.par,'TotalAcquiredVolumes')
				%	ppStr = [ppStr ' -' num2str(max(importStr{ii}.par.TotalAcquiredVolumes(:))) '-'];
				%else
				%	ppStr = [ppStr ' -.-'];
				%end
				%fprintf('%s\n',ppStr);

				% Type of an M0 image
				bJsonLocalM0isFile = 0;
				if strcmpi(importStr{ii}.x.M0,'separate_scan')
					if isfield(importStr{ii}.x,'M0PositionInASL4D') && (max(importStr{ii}.x.M0PositionInASL4D(:))>0)
						jsonLocal.M0 = true;
					elseif xASL_exist(fullfile(inSesPath,'M0.nii'))
						if length(fSes)>1
							%jsonLocal.M0 = fullfile(importStr{ii}.dirName,['sub-' subLabel],['ses-' sesLabel],'asl',['sub-' subLabel sesLabelUnd '_' bidsPar.strM0scan '.nii.gz']);
							jsonLocal.M0 = fullfile('perf',['sub-' subLabel sesLabelUnd]);
							bJsonLocalM0isFile = 1;
						else
							%jsonLocal.M0 = fullfile(importStr{ii}.dirName,['sub-' subLabel],'asl',['sub-' subLabel sesLabelUnd '_M0Scan.nii.gz']);
							jsonLocal.M0 = fullfile('perf',['sub-' subLabel sesLabelUnd]);
							bJsonLocalM0isFile = 1;
						end
					else
						if ~isempty(strfind(importStr{ii}.par.ASLContext,bidsPar.strM0scan))
							jsonLocal.M0 = true;
						else
							jsonLocal.M0 = false;
						end
					end
				else
					if strcmpi(importStr{ii}.x.M0,'UseControlAsM0')
						jsonLocal.M0 = false;
					else
						if strcmpi(importStr{ii}.x.M0,'no_background_suppression')
							jsonLocal.M0 = false;
						else
							jsonLocal.M0 = importStr{ii}.x.M0;
						end
					end
				end

				% Copy some things from the local JSON to the flavors
				jsonToFlavors = {'PulseSequenceType' 'PulseSequenceDetails' 'SliceTiming' 'ASLContext' 'LabelingType' 'LabelingDuration' 'PostLabelingDelay' 'BackgroundSuppression'...
					'BackgroundSuppressionNumberPulses' 'BackgroundSuppressionPulseTime' 'M0' 'AcquisitionVoxelSize'...
					'LookLocker' 'LabelingEfficiency' 'BolusCutOffFlag' 'BolusCutOffDelayTime' 'BolusCutOffTechnique'};
				for ll=1:length(jsonToFlavors)
					if isfield(jsonLocal,jsonToFlavors{ll})
						importStr{ii}.flavors.(jsonToFlavors{ll}) = jsonLocal.(jsonToFlavors{ll});
					end
				end

				% Remove the AslContext field and save it as a separate file
				fContext = fopen([aslOutLabel '_' bidsPar.strAslContext '.tsv'],'w+');
				fwrite(fContext,sprintf('volume_type\n'));
				fwrite(fContext,jsonLocal.ASLContext);
				fclose(fContext);

				jsonLocal = rmfield(jsonLocal,'ASLContext');

				if mm == 1
					for nn = 1:2
						if nn == 1
							nnStrIn = '';
							if xASL_exist(fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj},fSes{kk},'M0PERev.nii'))
								nnStrOut = '_dir-ap';

								tagPhaseEncodingDirection = 'j-';
								jsonLocal.PhaseEncodingDirection = 'j-';
								tagIntendedFor = [];
								tagTotalReadoutTime = importStr{ii}.par.TotalReadoutTime;

								if bJsonLocalM0isFile
									%jsonLocal.M0 = [jsonLocal.M0 nnStrOut '.nii.gz,' jsonLocal.M0 '_dir-pa' '.nii.gz'];
									jsonLocal.M0 = [jsonLocal.M0 nnStrOut '_' bidsPar.strM0scan '.nii.gz'];
								end
							else
								if bJsonLocalM0isFile
									jsonLocal.M0 = [jsonLocal.M0 '_' bidsPar.strM0scan '.nii.gz'];
								end
								nnStrOut = '';
								tagPhaseEncodingDirection = [];
								tagIntendedFor = [];
								tagTotalReadoutTime = [];
							end
						else
							nnStrIn = 'PERev';
							nnStrOut = '_dir-pa';
							tagPhaseEncodingDirection = 'j';
							tagIntendedFor = fullfile('perf',['sub-' subLabel sesLabelUnd '_dir-ap' '_' bidsPar.strM0scan '.nii.gz']);

							if isfield(importStr{ii}.par,'TotalReadoutTime')
								tagTotalReadoutTime = importStr{ii}.par.TotalReadoutTime;
							else
								tagTotalReadoutTime = [];
							end
						end
						% If M0, then copy M0 and add ASL path to the IntendedFor
						if xASL_exist(fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj},fSes{kk},['M0' nnStrIn '.nii']))
							jsonM0 = spm_jsonread(fullfile(inSesPath,['M0' nnStrIn '.json']));
							imM0   = xASL_io_Nifti2Im(fullfile(inSesPath,['M0' nnStrIn '.json']));
							if exist(fullfile(inSesPath,['M0' nnStrIn '_parms.mat']),'file')
								imParmsM0 = load(fullfile(inSesPath,['M0' nnStrIn '_parms.mat']));
							else
								imParmsM0 = [];
							end

							rescaleParms = [];
							ParmsFields = {'RescaleSlope' 'RWVSlope'    'MRScaleSlope' 'RescaleIntercept'...
								'RescaleSlopeOriginal' 'RescaleSlope' 'MRScaleSlope' 'UsePhilipsFloatNotDisplayScaling' 'RWVSlope'};
							JSONFields  = {'PhilipsRescaleSlope'  'PhilipsRWVSlope' 'PhilipsScaleSlope' 'PhilipsRescaleIntercept'...
								'RescaleSlopeOriginal' 'RescaleSlope' 'MRScaleSlope' 'UsePhilipsFloatNotDisplayScaling' 'RWVSlope'};
							for pp = 1:length(ParmsFields)
								if isfield(jsonM0,JSONFields{pp})
									rescaleParms.(ParmsFields{pp}) = jsonM0.(JSONFields{pp});
								end
								if ~isempty(imParmsM0)
									if isfield(imParmsM0.parms,ParmsFields{pp})
										rescaleParms.(ParmsFields{pp}) = imParmsM0.parms.(ParmsFields{pp});
									end
								end
							end


							if ~isempty(strfind(jsonDicom.Manufacturer,'Philips')) || ~isempty(strfind(jsonDicom.Manufacturer,'philips'))
								scaleFactor = xASL_adm_GetPhilipsScaling(rescaleParms,xASL_io_ReadNifti(fullfile(inSesPath,['M0' nnStrIn '.nii'])));
							else
								scaleFactor = 0;
							end

							if scaleFactor
								imM0 = imM0 .* scaleFactor;
							end

							jsonM0Write = [];
							% Copy all dicom ones
							for fn = fieldnames(jsonM0)'
								% Fields to skip
								bCP = 1;
								for ll=1:length(bidsPar.listFieldsRemoveGeneral)
									if strcmp(bidsPar.listFieldsRemoveGeneral{ll},fn{1})
										bCP = 0;
									end
								end
								for ll=1:length(bidsPar.listFieldsRemoveASL)
									if strcmp(bidsPar.listFieldsRemoveASL{ll},fn{1})
										bCP = 0;
									end
								end
								if bCP
									jsonM0Write.(fn{1}) = jsonM0.(fn{1});
								end
							end

							if isfield(jsonLocal,'SliceTiming')
								% Issue a warning if the SliceTiming was already existing for M0, but still overwrite with ASL one
								if isfield(jsonM0Write,'SliceTiming')
									warning('SliceTiming already existed for M0, overwriting with ASL');
								end

								if size(imNii,3) == size(imM0,3)
									% Either copy if the save number of slices in M0 as in ASL
									jsonM0Write.SliceTiming = jsonLocal.SliceTiming;
								else
									% Or recalculate for M0 if the number of slices differ
									jsonM0Write.SliceTiming = ((0:(size(imM0,3)-1))')*importStr{ii}.x.SliceReadoutTime;
								end
							else
								if isfield(jsonM0Write,'SliceTiming')
									jsonM0Write = rmfield(jsonM0Write,'SliceTiming');
									warning('Removing pre-existing SliceTiming from M0, as there was no SliceTiming for ASL');
								end
							end

							if isfield(importStr{ii}.x,'RepetitionTime')
								jsonM0Write.RepetitionTime = importStr{ii}.x.RepetitionTime;
							else
								jsonM0Write.RepetitionTime = jsonM0.RepetitionTime;
							end

							jsonM0Write.IntendedFor = [aslOutLabelRelative '_asl.nii.gz'];

							if ~isempty(tagPhaseEncodingDirection)
								jsonM0Write.PhaseEncodingDirection = tagPhaseEncodingDirection;
							end

							if isfield(jsonM0Write,'EffectiveEchoSpacing')
								if jsonM0Write.EffectiveEchoSpacing == 0
									jsonM0Write = rmfield(jsonM0Write,'EffectiveEchoSpacing');
								else
									jsonM0Write.EffectiveEchoSpacing = abs(jsonM0Write.EffectiveEchoSpacing);
								end
							end

							if isfield(jsonM0Write,'TotalReadoutTime')
								if jsonM0Write.TotalReadoutTime == 0
									jsonM0Write = rmfield(jsonM0Write,'TotalReadoutTime');
								else
									jsonM0Write.TotalReadoutTime = abs(jsonM0Write.TotalReadoutTime);
								end
							end

							if ~isempty(tagIntendedFor)
								jsonM0Write.IntendedFor = tagIntendedFor;
							end

							if ~isempty(tagTotalReadoutTime)
								jsonM0Write.TotalReadoutTime = tagTotalReadoutTime;
							end

							if nn == 2 && ~exist(fullfile(outSesPath,'fmap'),'dir')
								mkdir(fullfile(outSesPath,'fmap'));
							end

							% if scaling modified then save instead of copy
							if scaleFactor || size(imM0,4) == 1
								if nn == 1
									xASL_io_SaveNifti(fullfile(inSesPath,['M0' nnStrIn '.nii']),fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.nii.gz']),imM0,[],1,[]);
								else
									xASL_io_SaveNifti(fullfile(inSesPath,['M0' nnStrIn '.nii']),fullfile(outSesPath,'fmap',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.nii.gz']),imM0,[],1,[]);
								end
							else
								% Copy the M0
								if nn == 1
									xASL_Copy(fullfile(inSesPath,['M0' nnStrIn '.nii']),...
										fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.nii.gz']));
								else
									xASL_Copy(fullfile(inSesPath,['M0' nnStrIn '.nii']),...
										fullfile(outSesPath,'fmap',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.nii.gz']));
								end
							end
							% Save JSON to new dir
							jsonM0Write = finalJsonCheck(jsonM0Write,fieldOrderStruct,bidsPar.listRemoveIfEmpty);
							if nn == 1
								spm_jsonwrite(fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.json']),jsonM0Write);
							else
								spm_jsonwrite(fullfile(outSesPath,'fmap',['sub-' subLabel sesLabelUnd nnStrOut '_' bidsPar.strM0scan '.json']),jsonM0Write);
							end
						end
					end
				else
					if bJsonLocalM0isFile
						jsonLocal.M0 = [jsonLocal.M0 '.nii.gz'];
					end
				end
				% Save JSON to new dir
				jsonLocal = finalJsonCheck(jsonLocal,fieldOrderStruct,bidsPar.listRemoveIfEmpty);
				spm_jsonwrite([aslOutLabel '_asl.json'],jsonLocal);

			end
		end
	end
end

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
function jsonOut = finalJsonCheck(jsonIn,fieldOrderStruct,listRemoveIfEmpty)
jsonOut = jsonIn;
if isfield(jsonOut,'CoilString')
	switch (jsonOut.Manufacturer)
		case 'Philips'
			jsonOut.ReceiveCoilName = jsonOut.CoilString;
		case 'GE'
			jsonOut.ReceiveCoilName = jsonOut.CoilString;
		case 'Siemens'
			jsonOut.ReceiveCoilActiveElements = jsonOut.CoilString;
		otherwise
			error('Unknown manufacturer')
	end
	jsonOut = rmfield(jsonOut,'CoilString');
end

jsonOut = rmfield(jsonOut,'NumberOfAverages');

if isfield(jsonOut,'NumberSegments')
	jsonOut.NumberShots = jsonOut.NumberSegments;
	jsonOut = rmfield(jsonOut,'NumberSegments');
end

if isfield(jsonOut,'PhaseEncodingAxis')% && strcmpi(jsonOut.Manufacturer,'Philips')
	if ~isfield(jsonOut,'PhaseEncodingDirection')
		jsonOut.PhaseEncodingDirection = jsonOut.PhaseEncodingAxis;
	end
	jsonOut = rmfield(jsonOut,'PhaseEncodingAxis');
end

for ii = 1:length(listRemoveIfEmpty)
	if isfield(jsonOut,listRemoveIfEmpty{ii})
		if isempty(jsonOut.(listRemoveIfEmpty{ii}))
			jsonOut = rmfield(jsonOut,listRemoveIfEmpty{ii});
		end
	end
end

% And sort the fields
jsonOut = xASL_adm_OrderFields(jsonOut,fieldOrderStruct);
end
