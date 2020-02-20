%% Load the zeroed-out DICOMs the BIDS format and fill in all the missing info to the JSON
ExploreASL_Initialize([],false);
basePath = '/home/janpetr/tmp/BIDS';
parPath  = '/home/janpetr/tmp/BIDSpar';
outputPath = '/home/janpetr/tmp/BIDSout';
finalPath = '/home/janpetr/tmp/BIDSfinal';
lRMFields = {'InstitutionName' 'InstitutionalDepartmentName' 'InstitutionAddress' 'DeviceSerialNumber' 'StationName' 'ProcedureStepDescription' 'SeriesDescription' 'ProtocolName'}; % Fields to exclude
lAnat = {'T1' 'T2' 'FLAIR'}; % A list of anatomical scans to include
%% Load the list of the directories
fList = xASL_adm_GetFileList(basePath,[],false,[],true);

%% Copy all to the 'raw' subfolder

for ii = 1:length(fList)
	fFiles = xASL_adm_GetFileList(fullfile(basePath,fList{ii}),[],false,[],false);
	fDirs = xASL_adm_GetFileList(fullfile(basePath,fList{ii}),[],false,[],true);
	for jj = 1:length(fFiles)
		xASL_Move(fullfile(basePath,fList{ii},fFiles{jj}),fullfile(basePath,fList{ii},'raw',fFiles{jj}));
	end
	for jj = 1:length(fDirs)
		if ~strcmp(fDirs{jj},'raw')
			xASL_Move(fullfile(basePath,fList{ii},fDirs{jj}),fullfile(basePath,fList{ii},'raw',fDirs{jj}));
		end
	end
end	

%% Specific handling
mkdir([basePath '/Philips_PCASL_3DGRASE_Divers/raw/Patient1/ASL_Session1/dat/dat/']);
system(['mv ' basePath '/Philips_PCASL_3DGRASE_Divers/raw/* ' basePath '/Philips_PCASL_3DGRASE_Divers/raw/Patient1/ASL_Session1/dat/dat/'])

%% Specify study-parameters for import
importStr = [];
for ii = 1:length(fList)
	importStr{ii}.dirName = fList{ii};
	switch (fList{ii})
		case 'Philips_PCASL_3DGRASE_Divers' 
			importStr{ii}.configName = 'Divers_Bonn';
		case 'Philips_PCASL_2DEPI_Frontier'
			importStr{ii}.configName = 'FRONTIER';
		case 'Philips_PCASL_2DEPI_Chili'
			importStr{ii}.bLoadConfig = false;
			importStr{ii}.imPar.folderHierarchy = { '^(Sub-\d{4})$' '^(ASL|T1)$' '^DICOM'};
			importStr{ii}.imPar.tokenOrdering = [ 1 0 2];
			importStr{ii}.imPar.tokenSessionAliases = { '' };
			importStr{ii}.imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^T1$', 'T1'};
			importStr{ii}.imPar.bMatchDirectories = true;
		case {'GE_PCASL_3Dspiral_Product_22q11', 'GE_PCASL_3Dspiral_WIP_Oslo_AntiPsychotics_Old','GE_PCASL_2DEPI_stripped_3CV','Philips_PCASL_2DEPI_stripped_3CV1','Philips_PCASL_2DEPI_stripped_3CV2',...
			  'Siemens_PCASL_2DEPI_stripped_3CV','GE_PCASL_3Dspiral_Product_GE', 'GE_PCASL_3Dspiral_WIP_KCL_INOX','Philips_PCASL_2DEPI_Bsup_EPAD1','Philips_PCASL_2DEPI_Bsup_EPAD2',...
			  'Philips_PCASL_2DEPI_Bsup_EPAD3','Siemens_PASL_2DEPI_noBsup2_EPAD','Siemens_PASL_2DEPI_noBsup_EPAD','Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD','Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD2',...
			  'Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD3','Philips_PCASL_2DEPI_Achieva_Bsup_GENFI','Philips_PCASL_2DEPI_Achieva_noBsup_GENFI','Philips_PCASL_2DEPI_BioCog_Old',...
			  'Philips_PCASL_2DEPI_CP_Tavi_HC','Philips_PCASL_2DEPI_intera_FIND','Siemens_PCASL_3DGRASE_Sleep_Oslo_trial','Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38',...
			  'Philips_PCASL_2DEPI_Ingenia_FIND','Philips_PCASL_3DGRASE_Dent_example','Siemens_PASL_3DGRASE_APGEM_1','Siemens_PCASL_2DEPI_BioCog','Siemens_PCASL_2DEPI_Harmy_recombine_ASLscans',...
			  'Siemens_PCASL_3DGRASE_failed_APGEM2','Siemens_PASL_3DGRASE_GENFI','Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar',...
			  'Philips_PCASL_2DEPI_intera_FIND_LL','Philips_PCASL_2DEPI_intera_FIND_multiPLD','Philips_PCASL_2DEPI_intera_FIND_QUASAR'}
			importStr{ii}.imPar.folderHierarchy = { '^(.)+$' '^(ASL|T1w|M0|T2|FLAIR)$' };
			importStr{ii}.imPar.tokenOrdering = [ 1 0 2];
			importStr{ii}.imPar.tokenSessionAliases = { '', ''};
			importStr{ii}.imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^T1w$', 'T1';'^M0$', 'M0';'^T2$', 'T2';'^FLAIR$' , 'FLAIR'};
			importStr{ii}.imPar.bMatchDirectories = true;
			importStr{ii}.bLoadConfig = false;
			
			%case
			%importStr{ii}.configName = '';
			%importStr{ii}.bLoadConfig = false;
	end
end	

%% Go through all studies and import them
for ii = 1:length(fList)
	% Load the study parameters
	if ~isfield(importStr{ii},'bLoadConfig') || importStr{ii}.bLoadConfig
		imPar = ExploreASL_ImportConfig(fullfile(basePath,importStr{ii}.configName));
	else
		importStr{ii}.imPar.AnalysisRoot = basePath;
		importStr{ii}.imPar.RawRoot = basePath;
		
		imPar = importStr{ii}.imPar;
	end
	
	% Change the studyID to the new one
	imPar.studyID = importStr{ii}.dirName;
	
	% Define the output directory
	imPar.AnalysisRoot = outputPath;
	
	% Import the whole session to JSON and NIFTI
	ExploreASL_Import(imPar,false, true, false, true, false);
end

system(['rm ' outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0_a_2.json']);
system(['rm ' outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0_a_2.nii']);
system(['mv ' outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0_1.nii ' outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0.nii']);
system(['mv ' outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0_1.json ' outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0.json']);

system(['mv ' outputPath '/GE_PCASL_3Dspiral_Product_GE/analysis/Sub1/ASL_1/ASL4D_1.nii ' outputPath '/GE_PCASL_3Dspiral_Product_GE/analysis/Sub1/ASL_1/ASL4D.nii']);
system(['mv ' outputPath '/GE_PCASL_3Dspiral_Product_GE/analysis/Sub1/ASL_1/ASL4D_1.json ' outputPath '/GE_PCASL_3Dspiral_Product_GE/analysis/Sub1/ASL_1/ASL4D.json']);
system(['mv ' outputPath '/GE_PCASL_3Dspiral_Product_GE/analysis/Sub1/ASL_1/ASL4D_a_2.nii ' outputPath '/GE_PCASL_3Dspiral_Product_GE/analysis/Sub1/ASL_1/M0.nii']);
system(['mv ' outputPath '/GE_PCASL_3Dspiral_Product_GE/analysis/Sub1/ASL_1/ASL4D_a_2.json ' outputPath '/GE_PCASL_3Dspiral_Product_GE/analysis/Sub1/ASL_1/M0.json']);

system(['mv ' outputPath '/Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38/analysis/Sub2/ASL_1/ASL4D_10.nii ' outputPath '/Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38/analysis/Sub2/ASL_1/ASL4D.nii']);
system(['mv ' outputPath '/Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38/analysis/Sub2/ASL_1/ASL4D_10.json ' outputPath '/Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38/analysis/Sub2/ASL_1/ASL4D.json']);
system(['rm ' outputPath '/Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38/analysis/Sub2/ASL_1/ASL4D_*']);

system(['mv ' outputPath '/Siemens_PCASL_2DEPI_Harmy_recombine_ASLscans/analysis/Sub1/ASL_1/ASL4D_3_1.json ' outputPath '/Siemens_PCASL_2DEPI_Harmy_recombine_ASLscans/analysis/Sub1/ASL_1/ASL4D_1.json']);
system(['mv ' outputPath '/Siemens_PCASL_2DEPI_Harmy_recombine_ASLscans/analysis/Sub1/ASL_1/ASL4D_3_1.nii ' outputPath '/Siemens_PCASL_2DEPI_Harmy_recombine_ASLscans/analysis/Sub1/ASL_1/ASL4D_1.nii']);

system(['mv ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0_1.nii ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0.nii']);
system(['mv ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0_1.json ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0.json']);
system(['mv ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0_1_2.nii ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0PERev.nii']);
system(['mv ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0_1_2.json ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0PERev.json']);

%% Specify the missing study parameters
for ii = 1:length(fList)
	% Load the data-par
	
	clear x;
	dataParPath = fullfile(parPath,importStr{ii}.dirName);
	
	if exist([dataParPath,'.mat'],'file')
        TempVar = load([dataParPath,'.mat']);
        FieldN = fields(TempVar);
        x = TempVar.(FieldN{1});
    elseif exist([dataParPath,'.json'],'file')
        % JSON import
        x = xASL_import_json([dataParPath,'.json']);
    elseif exist([dataParPath,'.m'],'file')
		PathJSON = xASL_init_ConvertM2JSON([dataParPath,'.m']); % convert .m to .json
		x = xASL_import_json(PathJSON);
		xASL_delete([dataParPath,'.m']);
	end
	
	% Save the x-structure for the study
	importStr{ii}.x = x;

	% Defaults to be overwritten
	importStr{ii}.par = [];
	importStr{ii}.par.VascularCrushing = 'false';
	importStr{ii}.par.LabelingSlabLocation = 'Random description';
	importStr{ii}.par.LabelingOrientation = 'Random description';
	importStr{ii}.par.LabelingDistance = 40;
	
	% Units  mL/100g/min
	% ASLContext
	switch (fList{ii})
		case 'Siemens_PCASL_2DEPI_Harmy_recombine_ASLscans'
			importStr{ii}.par.ASLContext = '(Label+Control)*23';
			importStr{ii}.par.LabelingType = 'PCASL';
			
		case 'Siemens_PCASL_3DGRASE_failed_APGEM2'
			importStr{ii}.par.ASLContext = 'M0+((Label+Control)*12)';
			importStr{ii}.par.LabelingType = 'PCASL';
		  
		case 'Siemens_PASL_3DGRASE_GENFI'
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.ASLContext = 'DeltaM';
			importStr{ii}.par.BolusCutOffFlag = false;
		  
		case 'Siemens_PCASL_2DEPI_BioCog'
			importStr{ii}.par.ASLContext = '(Label+Control)*44';
			importStr{ii}.par.LabelingType = 'PCASL';
		  
		case 'Siemens_PASL_3DGRASE_APGEM_1'
			importStr{ii}.par.ASLContext = 'Label+Control';
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.LabelingSlabLocation = 'FAIR';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 700;
			importStr{ii}.par.BolusCutOffTechnique = 'QUIPSSII';
		  
		case 'Philips_PCASL_3DGRASE_Dent_example'
			importStr{ii}.par.ASLContext = '(M0*2)+((Label+Control)*7)';
			importStr{ii}.par.LabelingType = 'PCASL';
			
		case 'Siemens_PCASL_3DGRASE_Sleep_Oslo_trial'
			importStr{ii}.par.ASLContext = 'M0+((Label+Control)*15)';
			importStr{ii}.par.LabelingType = 'PCASL';
			
		case 'Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38'
			importStr{ii}.par.ASLContext = 'DeltaM*15';
			importStr{ii}.par.LabelingType = 'PCASL';
			
		case 'Philips_PCASL_2DEPI_BioCog_Old'
			importStr{ii}.par.ASLContext = '(Label+Control)*38';
			importStr{ii}.par.LabelingType = 'PCASL';
			
		case 'Philips_PCASL_2DEPI_CP_Tavi_HC'		
			importStr{ii}.par.ASLContext = '(Label+Control)*32';
			importStr{ii}.par.LabelingType = 'PCASL';
			
		case {'Philips_PCASL_2DEPI_intera_FIND','Philips_PCASL_2DEPI_Ingenia_FIND','Philips_PCASL_2DEPI_intera_FIND_LL','Philips_PCASL_2DEPI_intera_FIND_multiPLD','Philips_PCASL_2DEPI_intera_FIND_QUASAR'}
			importStr{ii}.par.ASLContext = '(Label+Control)*75';
			importStr{ii}.par.LabelingType = 'PCASL';
						
		case {'Philips_PCASL_2DEPI_Achieva_Bsup_GENFI','Philips_PCASL_2DEPI_Achieva_noBsup_GENFI'}
			importStr{ii}.par.ASLContext = '(Label+Control)*40';
			importStr{ii}.par.LabelingType = 'PCASL';
			
		case {'Philips_PCASL_2DEPI_Bsup_EPAD1','Philips_PCASL_2DEPI_Bsup_EPAD2'}
			importStr{ii}.par.ASLContext = '(Label+Control)*30';
			importStr{ii}.par.LabelingType = 'PCASL';
			
		case 'Philips_PCASL_2DEPI_Bsup_EPAD3'
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.ASLContext = 'DeltaM';
			importStr{ii}.par.LabelingEfficiency = 0.83;
			
		case {'Siemens_PASL_2DEPI_noBsup2_EPAD','Siemens_PASL_2DEPI_noBsup_EPAD'}
			importStr{ii}.par.ASLContext = '(Label+Control)*31';
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 600;
			importStr{ii}.par.BolusCutOffTechnique = 'Q2TIPS';
			importStr{ii}.par.LabelingSlabThickness = 80;
			
		case {'Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD','Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD2'}
			importStr{ii}.par.ASLContext = '(Label+Control)*10';
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 0;
			importStr{ii}.par.BolusCutOffTechnique = 'QUIPSSII';
			
		case 'Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD3'
			importStr{ii}.par.ASLContext = '(Label+Control)*2';
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = [0 200 400];
			importStr{ii}.par.BolusCutOffTechnique = 'QUIPSS';
			importStr{ii}.par.LabelingSlabThickness = 60;
			
		case {'Philips_PCASL_3DGRASE_Divers','GE_PCASL_3Dspiral_Product_22q11', 'GE_PCASL_3Dspiral_WIP_Oslo_AntiPsychotics_Old','GE_PCASL_3Dspiral_Product_GE',...
			  'Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar'}
			importStr{ii}.par.Units = 'mL/100g/min';
			importStr{ii}.par.ASLContext = 'CBF';
			importStr{ii}.par.LabelingType = 'PCASL';
			
		case 'GE_PCASL_3Dspiral_WIP_KCL_INOX'
			importStr{ii}.par.Units = 'mL/100g/min';
			importStr{ii}.par.ASLContext = 'CBF+M0';
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.AcquisitionResolution = [4 4 8];
			importStr{ii}.par.AverageLabelingGradient = 0.6;
			importStr{ii}.par.SliceSelectiveLabelingGradient = 6;
			importStr{ii}.par.FlipAngleOfB1LabelingPulses = 18;
			importStr{ii}.par.PulseDuration = 0.0005;
			importStr{ii}.par.PCASLType = 'balanced';
			
		case {'Philips_PCASL_2DEPI_Frontier','Philips_PCASL_2DEPI_Chili'}
			importStr{ii}.par.ASLContext = '(Control+Label)*30';
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.AverageLabelingGradient = 0.7;
			importStr{ii}.par.SliceSelectiveLabelingGradient = 7;
			importStr{ii}.par.FlipAngleOfB1LabelingPulses = 23;
			importStr{ii}.par.PulseDuration = 0.0005;
			importStr{ii}.par.PCASLType = 'balanced';
							
		case {'GE_PCASL_2DEPI_stripped_3CV','Philips_PCASL_2DEPI_stripped_3CV1','Philips_PCASL_2DEPI_stripped_3CV2','Siemens_PCASL_2DEPI_stripped_3CV'}
			importStr{ii}.par.ASLContext = '(Control+Label)*70';
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.LabelingSlabLocation = 'Fixed, 9 cm below ACPC';
			importStr{ii}.par.AverageLabelingGradient = 0.6;
			importStr{ii}.par.SliceSelectiveLabelingGradient = 6;
			importStr{ii}.par.FlipAngleOfB1LabelingPulses = 25;
			importStr{ii}.par.AcquisitionDuration = 672;
			importStr{ii}.par.InterPulseSpacing = 0.00124;
			importStr{ii}.par.PCASLType = 'balanced';
			
			%case
			%importStr{ii}.ASLContext
			%importStr{ii}.bLoadConfig = false;
	end
	switch (fList{ii})
		case {'Philips_PCASL_2DEPI_stripped_3CV2','Siemens_PCASL_2DEPI_stripped_3CV'}
			importStr{ii}.par.AcquisitionDuration = 658;
			importStr{ii}.par.InterPulseSpacing = 0.00115;
	end
	
		
	
	% Process all the data and automatically fill in the missing parameters
	if strcmp(importStr{ii}.x.readout_dim,'2D')
		importStr{ii}.par.PulseSequenceType = '2D_EPI';
	else
		if strcmp(importStr{ii}.x.Vendor,'GE') || strcmp(importStr{ii}.x.Vendor,'GE_WIP') || strcmp(importStr{ii}.x.Vendor,'GE_product')
			importStr{ii}.par.PulseSequenceType = '3D_spiral';
		else
			importStr{ii}.par.PulseSequenceType = '3D_GRASE';
		end
	end
	
	% Labeling delays and durations
	if strcmp(importStr{ii}.par.LabelingType,'PASL')
		importStr{ii}.par.LabelingDuration = 0;% x.Q.LabelingDuration           = 1800;  % for PASL this is TI1
		importStr{ii}.par.InitialPostLabelDelay = importStr{ii}.x.Q.Initial_PLD;
		if importStr{ii}.par.BolusCutOffFlag 
			importStr{ii}.par.BolusCutOffTimingSequence = x.Q.LabelingDuration;
		end		
	else
		importStr{ii}.par.LabelingDuration = importStr{ii}.x.Q.LabelingDuration;
		importStr{ii}.par.InitialPostLabelDelay = importStr{ii}.x.Q.Initial_PLD;
	end
		
	if importStr{ii}.x.Q.BackGrSupprPulses == 0
		importStr{ii}.par.BackgroundSuppression = 'false';
	else
		importStr{ii}.par.BackgroundSuppression = 'true';
		switch (importStr{ii}.x.Q.BackGrSupprPulses)
			case 2
				if importStr{ii}.par.InitialPostLabelDelay > 1.750
					importStr{ii}.par.BackgroundSuppressionPulseTime = [1.75 0.524];
				elseif importStr{ii}.par.InitialPostLabelDelay > 1.495
					importStr{ii}.par.BackgroundSuppressionPulseTime = [1.495 0.345];
				else
					error('Pulses not fitting');
				end
			case 4
				if importStr{ii}.par.InitialPostLabelDelay > 1.510
					importStr{ii}.par.BackgroundSuppressionPulseTime = [1.510 0.875 0.375 0.095];
				else
					error('Pulses not fitting');
				end
			case 5
				if importStr{ii}.par.InitialPostLabelDelay > 1.510
					importStr{ii}.par.BackgroundSuppressionPulseTime = [importStr{ii}.par.InitialPostLabelDelay+importStr{ii}.par.LabelingDuration+1 1.510 0.875 0.375 0.095];
				else
					error('Pulses not fitting');
				end
			otherwise
				error('Unknown number of pulses');
		end
		
	end
	
	% Last round of edits
	switch (fList{ii})
		case 'Siemens_PCASL_3DGRASE_failed_APGEM2'
			%importStr{ii}.par.ASLContext = 'M0+((Label+Control)*12)';
			importStr{ii}.par.LabelingDuration = [0 repmat(1800,[1,24])];
			importStr{ii}.par.VascularCrushing = true;
			importStr{ii}.par.VascularCrushingVenc = 10;
		
		case 'Philips_PCASL_2DEPI_intera_FIND_LL'
			importStr{ii}.par.ASLContext = '(Label*15+Control*15)*5';
			importStr{ii}.par.FlipAngle = 25;
			importStr{ii}.par.LookLocker = true;
			importStr{ii}.par.InitialPostLabelDelay = repmat(250:250:3750,[1 10]);
		
		case 'Philips_PCASL_2DEPI_intera_FIND_multiPLD'
			importStr{ii}.par.ASLContext = '(Label+Control)*75';
			InitialPostLabelDelay = repmat([500 500 1000 1000 1500 1500 1800 1800 2200 2200],[1 15]);
			
		case 'Philips_PCASL_2DEPI_intera_FIND_QUASAR'
			importStr{ii}.par.ASLContext = '(Label*15+Control*15)*5';
			importStr{ii}.par.FlipAngle = [25*ones(1,120),12*ones(1,30)];
			importStr{ii}.par.LookLocker = true;
			importStr{ii}.par.VascularCrushing = true;
			importStr{ii}.par.VascularCrushingVenc = [zeros(1,60),10*ones(1,60),zeros(1,3)];
			importStr{ii}.par.InitialPostLabelDelay = repmat(250:250:3750,[1 10]);
			
	end
end

%% Go through all studies and check all the M0 and ASLs and modify the JSONs
for ii = 1:length(fList)
	% Make a copy of the par to the Flavors
	importStr{ii}.flavors = importStr{ii}.par;
	
	% Make an output directory
	if ~exist(fullfile(finalPath,importStr{ii}.dirName),'dir')
		mkdir(fullfile(finalPath,importStr{ii}.dirName));
	end
	
	% Go through all subjects
	fSubs = xASL_adm_GetFileList(fullfile(outputPath,importStr{ii}.dirName,'analysis'),[],false,[],true);
	for jj = 1:length(fSubs)
		
		subLabel = xASL_adm_CorrectName(fSubs{jj},2);
		
		% Make a subject directory
		if ~exist(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel]),'dir')
			mkdir(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel]));
		end
			
		% Go throught the list of anat files
		for iiAnat = lAnat
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
					for ll=1:length(lRMFields)
						if strcmp(lRMFields{ll},fn{1})
							bCP = 0;
						end
					end
					if bCP
						jsonLocal.(fn{1}) = jsonAnat.(fn{1});
					end
				end
				
				% Save the JSON
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
					mkdir(fullfile(finalPath,importStr{ii}.dirName,['sub-' subLabel],'asl'));
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
					aslOutLabel = fullfile(outSesPath,'asl',...
					['sub-' subLabel sesLabelUnd '_run-' num2str(mm)]);
				else
					aslLabel = 'ASL4D';
					aslOutLabel = fullfile(outSesPath,'asl',...
					['sub-' subLabel sesLabelUnd]);
				end
				
				% Copy the ASL
				xASL_Copy(fullfile(inSesPath,[aslLabel '.nii']),[aslOutLabel '_asl.nii.gz']);
			
				% Load the JSON
				%jsonLocal = xASL_import_json(fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj},fSes{kk},'ASL4D.json'));
				jsonDicom = spm_jsonread(fullfile(inSesPath,[aslLabel '.json']));
				imNii = xASL_io_Nifti2Im(fullfile(inSesPath,[aslLabel '.nii']));
				
				% Copy the basic ones
				jsonLocal = importStr{ii}.par;
				
				% Copy all dicom ones
				for fn = fieldnames(jsonDicom)'
					% Fields to skip
					bCP = 1;
					for ll=1:length(lRMFields)
						if strcmp(lRMFields{ll},fn{1})
							bCP = 0;
						end
					end
					if bCP
						jsonLocal.(fn{1}) = jsonDicom.(fn{1});
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
				else
					jsonLocal.PulseSequenceDetails = '';
				end
				if isfield(jsonDicom,'SoftwareVersions')
					if ~isempty(jsonLocal.PulseSequenceDetails)
						jsonLocal.PulseSequenceDetails = [jsonLocal.PulseSequenceDetails '-'];
					end
					jsonLocal.PulseSequenceDetails = [jsonLocal.PulseSequenceDetails jsonDicom.SoftwareVersions];
					jsonLocal.MRSoftwareVersion = jsonDicom.SoftwareVersions;
				end
				
				% Fill in extra parameters based on the JSON from the data
				if importStr{ii}.par.PulseSequenceType(1) == '2'
					jsonLocal.SliceTiming = ((0:(size(imNii,3)-1))')*importStr{ii}.x.Q.SliceReadoutTime;
				end
				
				% Type of an M0 image
				if strcmp(importStr{ii}.x.M0,'separate_scan')
					if isfield(importStr{ii}.x,'M0PositionInASL4D') && (max(importStr{ii}.x.M0PositionInASL4D(:))>0)
						jsonLocal.M0 = 'true';
					elseif exist(fullfile(inSesPath,'M0.nii'),'file')
						if length(fSes)>1
							jsonLocal.M0 = fullfile(importStr{ii}.dirName,['sub-' subLabel],['ses-' sesLabel],'asl',['sub-' subLabel sesLabelUnd '_M0Scan.nii.gz']);
						else
							jsonLocal.M0 = fullfile(importStr{ii}.dirName,['sub-' subLabel],'asl',['sub-' subLabel sesLabelUnd '_M0Scan.nii.gz']);
						end
					else
						if ~isempty(strfind(importStr{ii}.par.ASLContext,'M0'))
							jsonLocal.M0 = 'true';
						else
							jsonLocal.M0 = 'false';
						end
					end
				else
					if strcmp(importStr{ii}.x.M0,'UseControlAsM0')
						jsonLocal.M0 = 'false';
					else
						if strcmp(importStr{ii}.x.M0,'no_background_suppression')
						else
							jsonLocal.M0 = importStr{ii}.x.M0;
						end
					end
				end
				
				% Copy some things from the local JSON to the flavors
				jsonToFlavors = {'PulseSequenceType' 'PulseSequenceDetails' 'SliceTiming' 'ASLContext' 'LabelingType' 'LabelingDuration' 'InitialPostLabelDelay' 'BackgroundSuppression' 'M0' 'AcquisitionResolution'...
					'LookLocker' 'LabelingEfficiency' 'BolusCutOffFlag' 'BolusCutOffTimingSequence' 'BolusCutOffDelayTime' 'BolusCutOffTechnique'};
				for ll=1:length(jsonToFlavors)
					if isfield(jsonLocal,jsonToFlavors{ll})
						importStr{ii}.flavors.(jsonToFlavors{ll}) = jsonLocal.(jsonToFlavors{ll});
					end
				end
				
				% Save JSON to new dir
				spm_jsonwrite([aslOutLabel '_asl.json'],jsonLocal);
						
				if mm == 1
					for nn = 1:2
						if nn == 1
							nnStr = '';
						else
							nnStr = 'PERev';
						end
						% If M0, then copy M0 and add ASL path to the IntendedFor
						if exist(fullfile(outputPath,importStr{ii}.dirName,'analysis',fSubs{jj},fSes{kk},['M0' nnStr '.nii']),'file')
							jsonM0 = spm_jsonread(fullfile(inSesPath,['M0' nnStr '.json']));
							
							% Copy all dicom ones
							for fn = fieldnames(jsonM0)'
								% Fields to skip
								bCP = 1;
								for ll=1:length(lRMFields)
									if strcmp(lRMFields{ll},fn{1})
										bCP = 0;
									end
								end
								if bCP
									jsonM0Write.(fn{1}) = jsonM0.(fn{1});
								end
							end
							jsonM0Write.RepetitionTime = jsonM0.RepetitionTime;
							jsonM0Write.IntendedFor = [aslOutLabel '_asl.nii.gz'];
							% Copy the M0
							xASL_Copy(fullfile(inSesPath,['M0' nnStr '.nii']),...
								fullfile(outSesPath,'asl',['sub-' subLabel sesLabelUnd '_M0Scan' nnStr '.nii.gz']));
							% Save JSON to new dir
							spm_jsonwrite(fullfile(outSesPath,'asl',['sub-' subLabel sesLabelUnd '_M0Scan' nnStr '.json']),jsonM0Write);
						end
					end
				end
			end
		end
	end
end

	

%% Export for ASL JSON flavours
fCSVOut = fopen(fullfile(finalPath,'flavours.csv'),'w');

fieldsCSV = {'PulseSequenceType' 'PulseSequenceDetails' 'SliceTiming' 'ASLContext' 'LabelingType' 'LabelingDuration' 'InitialPostLabelDelay' 'BackgroundSuppression' 'M0' 'AcquisitionResolution'...
'LookLocker' 'LabelingEfficiency' 'BolusCutOffFlag' 'BolusCutOffTimingSequence' 'BolusCutOffDelayTime' 'BolusCutOffTechnique'};

for ii = 1:length(fList)
	% New line for a new study
	if ii>1
		fwrite(fCSVOut,sprintf('\n'));
	end
	% Write the name of the study
	fwrite(fCSVOut,importStr{ii}.dirName);
	
	% For all possible fields
	for jj = 1:length(fieldsCSV)
		fwrite(fCSVOut,'; ');
		if isfield(importStr{ii}.flavors,fieldsCSV{jj})
			if isnumeric(importStr{ii}.flavors.(fieldsCSV{jj}))
				otStr = xASL_num2str(importStr{ii}.flavors.(fieldsCSV{jj}));
				if size(otStr,1)>1
					otStr(:,end+1) = ',';
				end
				otStr = otStr';
				fwrite(fCSVOut,otStr(:)');
			else
				fwrite(fCSVOut,importStr{ii}.flavors.(fieldsCSV{jj}));
			end
		else
			fwrite(fCSVOut,' ');
		end
	end
end
% Close file
fclose(fCSVOut);

