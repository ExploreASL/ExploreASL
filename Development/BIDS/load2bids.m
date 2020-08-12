%% Load the zeroed-out DICOMs the BIDS format and fill in all the missing info to the JSON
ExploreASL_Initialize([],false);

% FILLIN
baseDir = '/pet/projekte/asl/data/BIDS/';
basePath = [baseDir 'BIDS']; % The raw dicoms should be here in /StudyName/raw/...
parPath  = [baseDir 'BIDSpar']; % Put it your data_par.m file to load your parameters with name StudyName.m or StudyName.json (it converts all to JSON on the first run)
outputPath = [baseDir 'BIDSout']; % The directory to save the DICOM2NII import
finalPath = [baseDir 'BIDSfinal']; % Takes files in NIFTI+JSON from outputPath and saves the complete BIDS format to finalPath
anonymPath = [baseDir 'BIDSanonymized']; % Takes files in NIFTI+JSON from outputPath and saves the complete BIDS format to finalPath

lRMFields = {'InstitutionName' 'InstitutionalDepartmentName' 'InstitutionAddress' 'DeviceSerialNumber' 'StationName' 'ProcedureStepDescription' 'SeriesDescription' 'ProtocolName'...
	         'PhilipsRescaleSlope'  'PhilipsRWVSlope' 'PhilipsScaleSlope' 'PhilipsRescaleIntercept' 'UsePhilipsFloatNotDisplayScaling'...
			 'RescaleSlopeOriginal' 'RescaleSlope'    'MRScaleSlope'      'RescaleIntercept',... % Fields to exclude
			 'Modality', 'ImagingFrequency', 'PatientPosition','MRAcquisitionType','ImageType','PhaseEncodingPolarityGE',... % Additional fields to remove by Patricia
			 'SeriesNumber','AcquisitionTime','AcquisitionNumber','SliceThickness','SpacingBetweenSlices','SAR',...
			 'PercentPhaseFOV','AcquisitionMatrixPE','ReconMatrixPE','PixelBandwidth','ImageOrientationPatientDICOM',...
			 'InPlanePhaseEncodingDirectionDICOM','ConversionSoftware','ConversionSoftwareVersion','AcquisitionMatrix',...
			 'EchoTrainLength','PhaseEncodingSteps','BodyPartExamined','ShimSetting','TxRefAmp','PhaseResolution',...
			 'RefLinesPE','BandwidthPerPixelPhaseEncode','ImageComments','ConsistencyInfo','WipMemBlock','Interpolation2D',...
			 'SaturationStopTime','BaseResolution','DerivedVendorReportedEchoSpacing','RawImage'}; % Fields to exclude
lRMFieldsASL = {'InversionTime'}; % Fields to exclude from ASL only
% The correct order of fields in the JSON
fieldOrder = {'Manufacturer','ManufacturersModelName','DeviceSerialNumber','StationName','SoftwareVersions','HardcopyDeviceSoftwareVersion',...%BIDS fields
	          'MagneticFieldStrength','ReceiveCoilName','ReceiveCoilActiveElements','GradientSetType','MRTransmitCoilSequence','MatrixCoilMode',...
			  'CoilCombinationMethod','PulseSequenceType','ScanningSequence','SequenceVariant','ScanOptions','SequenceName','PulseSequenceDetails',...
			  'NonlinearGradientCorrection','NumberShots','ParallelReductionFactorInPlane','ParallellAcquisitionTechnique','PartialFourier',...
			  'PartialFourierDirection','PhaseEncodingDirection','EffectiveEchoSpacing','TotalReadoutTime','EchoTime','InversionTime','SliceTiming',...
			  'SliceEncodingDirection','DwellTime','FlipAngle','MultibandAccelerationFactor','NegativeContrast','AnatomicalLandmarkCoordinates',...
			  'RepetitionTime','VolumeTiming','TaskName','Units',... % And ASL fields
			  'LabelingType','PostLabelingDelay','BackgroundSuppression','M0','VascularCrushing','AcquisitionVoxelSize','TotalAcquiredVolumes',...
			  'BackgroundSuppressionNumberPulses','BackgroundSuppressionPulseTime','VascularCrushingVenc','LabelingLocationDescription',...
			  'LabelingOrientation','LabelingDistance','LookLocker','LabelingEfficiency','PCASLType','CASLType','LabelingDuration',...
			  'LabelingPulseAverageGradient','LabelingPulseMaximumGradient','LabelingPulseAverageB1','LabelingPulseDuration','LabelingPulseFlipAngle',...
			  'LabelingPulseInterval','PASLType','LabelingSlabThickness','BolusCutOffFlag','BolusCutOffTimingSequence','BolusCutOffDelayTime',...
			  'BolusCutOffTechnique'};
removeEmptyFields = {'EffectiveEchoSpacing','TotalReadoutTime'};
lAnat = {'T1w' 'T2w' 'FLAIR'}; % A list of anatomical scans to include
aslcontextStr = 'aslcontext';
labelStr = 'label';
controlStr = 'control';
m0scanStr = 'm0scan';
cbfStr = 'cbf';
deltamStr = 'deltam';
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

%% Specific handling - you can manually copy/move some files if needed to get them to the correct structure
mkdir([basePath '/Philips_PCASL_3DGRASE_Divers/raw/Patient1/ASL_Session1/dat/dat/']);
system(['mv ' basePath '/Philips_PCASL_3DGRASE_Divers/raw/* ' basePath '/Philips_PCASL_3DGRASE_Divers/raw/Patient1/ASL_Session1/dat/dat/'])

%% Specify study-parameters for import
% It either loads this information from the ExploreASL_ImportConfig.m based on the name of your study, but it is much easier to fill in your information directly here
importStr = [];
for ii = 1:length(fList)
	importStr{ii}.dirName = fList{ii};
	switch (fList{ii})
		case {'Siemens_PASL_singleTI_GIFMI', 'Siemens_PCASL_GIFMI', 'Siemens_PASL_multiTI_GIFMI'}
			importStr{ii}.imPar.folderHierarchy = { '^(.)+$' ['^(ASL|T1w|M0|M0-PA|ASL_NS|ASL_SS|NS_TI0300|NS_TI0600|NS_TI0900|NS_TI1200|'...
				                                              'NS_TI1500|NS_TI1800|NS_TI2100|NS_TI2400|NS_TI2700|NS_TI3000|SS_TI0300|SS_TI0600|'...
															  'SS_TI0900|SS_TI1200|SS_TI1500|SS_TI1800|SS_TI2100|SS_TI2400|SS_TI2700|SS_TI3000)$']};
			importStr{ii}.imPar.tokenOrdering = [ 1 0 2];
			importStr{ii}.imPar.tokenSessionAliases = { '', ''};
			importStr{ii}.imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^T1w$', 'T1w';'^M0$', 'M0';'^M0-PA$', 'M0_2';'^ASL_NS$', 'ASL4D_NS';'^ASL_SS$', 'ASL4D_SS';...
				                                     '^NS_TI0300$', 'ASL4D_NS_300'; '^NS_TI0600$', 'ASL4D_NS_600'; '^NS_TI0900$', 'ASL4D_NS_900'; '^NS_TI1200$', 'ASL4D_NS_1200';'^NS_TI1500$', 'ASL4D_NS_1500';...
													 '^NS_TI1800$', 'ASL4D_NS_1800';'^NS_TI2100$', 'ASL4D_NS_2100';'^NS_TI2400$', 'ASL4D_NS_2400';'^NS_TI2700$', 'ASL4D_NS_2700';'^NS_TI3000$', 'ASL4D_NS_3000';...
													 '^SS_TI0300$', 'ASL4D_SS_300'; '^SS_TI0600$', 'ASL4D_SS_600'; '^SS_TI0900$', 'ASL4D_SS_900'; '^SS_TI1200$', 'ASL4D_SS_1200';'^SS_TI1500$', 'ASL4D_SS_1500';...
													 '^SS_TI1800$', 'ASL4D_SS_1800';'^SS_TI2100$', 'ASL4D_SS_2100';'^SS_TI2400$', 'ASL4D_SS_2400';'^SS_TI2700$', 'ASL4D_SS_2700';'^SS_TI3000$', 'ASL4D_SS_3000';};
			importStr{ii}.imPar.bMatchDirectories = true;
			importStr{ii}.bLoadConfig = false;
		case 'Philips_PCASL_3DGRASE_Divers'
			%importStr{ii}.configName = 'Divers_Bonn';
			importStr{ii}.bLoadConfig = false;
			importStr{ii}.imPar.folderHierarchy = { '^(Patient\d{1})$' '^(ASL)_(Session1|Session2|Session3|Session4|Session5|Session6|Session7)$' '^.*$' '^.*$'};
			importStr{ii}.imPar.tokenOrdering = [ 1 3 2];
			importStr{ii}.imPar.tokenSessionAliases = {'Session1', 'ASL_1' ; 'Session2', 'ASL_2' ; 'Session3', 'ASL_3' ; 'Session4', 'ASL_4' ; 'Session5', 'ASL_5' ; 'Session6', 'ASL_6' ; 'Session7', 'ASL_7'};
			importStr{ii}.imPar.tokenScanAliases = {'^ASL$', 'ASL4D';'^t1$', 'T1w'};
			importStr{ii}.imPar.bMatchDirectories = true;
			importStr{ii}.imPar.dcmwildcard = '*.';
		case 'Philips_PCASL_2DEPI_Frontier'
			%importStr{ii}.configName = 'FRONTIER';
			importStr{ii}.bLoadConfig = false;
			importStr{ii}.imPar.folderHierarchy = {'^(P\d{2})$' '^(ASL|DSC|M0)$'};
			importStr{ii}.imPar.tokenOrdering = [ 1 0 2];
			importStr{ii}.imPar.tokenSessionAliases = {};
			importStr{ii}.imPar.tokenScanAliases = {'^ASL$', 'ASL4D';'^DSC$','DSC4D';'^M0$','M0'};
			importStr{ii}.imPar.bMatchDirectories = true;
		case 'Philips_PCASL_2DEPI_Chili'
			importStr{ii}.bLoadConfig = false;
			importStr{ii}.imPar.folderHierarchy = { '^(Sub-\d{4})$' '^(ASL|T1|M0)$' '^DICOM'};
			importStr{ii}.imPar.tokenOrdering = [ 1 0 2];
			importStr{ii}.imPar.tokenSessionAliases = { '' };
			importStr{ii}.imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^T1$', 'T1w';'^M0$', 'M0'};
			importStr{ii}.imPar.bMatchDirectories = true;
		case {'GE_PCASL_3Dspiral_Product_22q11', 'GE_PCASL_3Dspiral_WIP_Oslo_AntiPsychotics_Old','GE_PCASL_2DEPI_stripped_3CV','Philips_PCASL_2DEPI_stripped_3CV1','Philips_PCASL_2DEPI_stripped_3CV2',...
			  'Siemens_PCASL_2DEPI_stripped_3CV','GE_PCASL_3Dspiral_Product_GE', 'GE_PCASL_3Dspiral_WIP_KCL_INOX','Philips_PCASL_2DEPI_Bsup_EPAD1','Philips_PCASL_2DEPI_Bsup_EPAD2',...
			  'Philips_PCASL_2DEPI_Bsup_EPAD3','Siemens_PASL_2DEPI_noBsup2_EPAD','Siemens_PASL_2DEPI_noBsup_EPAD','Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD','Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD2',...
			  'Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD3','Philips_PCASL_2DEPI_Achieva_Bsup_GENFI','Philips_PCASL_2DEPI_Achieva_noBsup_GENFI','Philips_PCASL_2DEPI_BioCog_Old',...
			  'Philips_PCASL_2DEPI_CP_Tavi_HC','Philips_PCASL_2DEPI_intera_FIND','Siemens_PCASL_3DGRASE_Sleep_Oslo_trial','Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38',...
			  'Philips_PCASL_2DEPI_Ingenia_FIND','Philips_PCASL_3DGRASE_Dent_example','Siemens_PASL_3DGRASE_APGEM_1','Siemens_PCASL_2DEPI_BioCog','Siemens_PCASL_2DEPI_Harmy_recombine_ASLscans',...
			  'Siemens_PCASL_3DGRASE_failed_APGEM2','Siemens_PASL_3DGRASE_GENFI','Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar',...
			  'Philips_PCASL_2DEPI_intera_FIND_LL','Philips_PCASL_2DEPI_intera_FIND_multiPLD','Philips_PCASL_2DEPI_intera_FIND_QUASAR','GE_PCASL_3Dspiral_UCL','Philips_PCASL_2DEPI_UCL','Siemens_PCASL_3DGRASE_UCL'}
			importStr{ii}.imPar.folderHierarchy = { '^(.)+$' '^(ASL|T1w|M0|T2|FLAIR)$' };
			importStr{ii}.imPar.tokenOrdering = [ 1 0 2];
			importStr{ii}.imPar.tokenSessionAliases = { '', ''};
			importStr{ii}.imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^T1w$', 'T1w';'^M0$', 'M0';'^T2$', 'T2w';'^FLAIR$' , 'FLAIR'};
			importStr{ii}.imPar.bMatchDirectories = true;
			importStr{ii}.bLoadConfig = false;
		case {'StudyName'}
			% FILLIN
			% The standard information from the ImportConfig
			importStr{ii}.imPar.folderHierarchy = { '^(.)+$' '^(ASL|T1w|M0|T2|FLAIR)$' };
			importStr{ii}.imPar.tokenOrdering = [ 1 0 2];
			importStr{ii}.imPar.tokenSessionAliases = { '', ''};
			importStr{ii}.imPar.tokenScanAliases = { '^ASL$', 'ASL4D';'^T1w$', 'T1w';'^M0$', 'M0';'^T2$', 'T2w';'^FLAIR$' , 'FLAIR'};
			importStr{ii}.imPar.bMatchDirectories = true;
			% Specify to false that this info should be used instead of loading it from the ExploreASL_ImportConfig.m
			importStr{ii}.bLoadConfig = false;
	end
end

%% Go through all studies and import them
% This simply runs the ExploreASL_Import
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

system(['rm ' outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0_1.json']);
system(['rm ' outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0_1.nii']);
system(['mv ' outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0_a_2.nii ' outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0.nii']);
%xASL_Move([outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0_a_2.json'],[outputPath '/GE_PCASL_3Dspiral_Product_22q11/analysis/11/ASL_1/M0.json']);

system(['rm ' outputPath '/Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38/analysis/Sub2/ASL_1/ASL4D_0170*']);
system(['rm ' outputPath '/Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38/analysis/Sub2/ASL_1/ASL4D_2460*']);
nii_files = xASL_adm_GetFileList([outputPath '/Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38/analysis/Sub2/ASL_1/'],'^*.nii$','FPList',[],false);
nii_files = xASL_bids_MergeNifti(nii_files, 'ASL');

system(['mv ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0_1.nii ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0.nii']);
system(['mv ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0_1.json ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0.json']);
system(['mv ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0_1_2.nii ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0PERev.nii']);
system(['mv ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0_1_2.json ' outputPath '/Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar/analysis/Sub1/ASL_1/M0PERev.json']);

if xASL_exist([outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS.nii'])
	system(['rm ' outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS.json']);
	system(['mv ' outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/ASL4D_SS.json ' outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/ASL4D.json']);
	imNS = xASL_io_Nifti2Im([outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS.nii']);
	imSS = xASL_io_Nifti2Im([outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/ASL4D_SS.nii']);
	imNS(:,:,:,2) = imSS;
	xASL_io_SaveNifti([outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS.nii'],[outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/ASL4D.nii'],imNS/10,[],1,[]);
	system(['rm ' outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS.nii']);
	system(['rm ' outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/ASL4D_SS.nii']);
	system(['mv ' outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/M0_2.json ' outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/M0PERev.json']);
	system(['mv ' outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/M0_2.nii ' outputPath '/Siemens_PCASL_GIFMI/analysis/Sub1/ASL_1/M0PERev.nii']);
end

if xASL_exist([outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS_300.nii'])
	mTIvec = [300,600,900,1200,1500,1800,2100,2400,2700,3000];
	for ii = 1:length(mTIvec)
		if ii>1
			system(['rm ' outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(ii)) '.json']);
			system(['rm ' outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_SS_' num2str(mTIvec(ii)) '.json']);
			imNSSS(:,:,:,2*(ii-1)+1) = xASL_io_Nifti2Im([outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(ii)) '.nii']);
			imNSSS(:,:,:,2*(ii-1)+2) = xASL_io_Nifti2Im([outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_SS_' num2str(mTIvec(ii)) '.nii']);
		else
			system(['mv ' outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(ii)) '.json ' outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D.json']);
			system(['rm ' outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_SS_' num2str(mTIvec(ii)) '.json']);
			imNSSS = xASL_io_Nifti2Im([outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(ii)) '.nii']);
			imNSSS(:,:,:,2) = xASL_io_Nifti2Im([outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_SS_' num2str(mTIvec(ii)) '.nii']);
		end
	end
	xASL_io_SaveNifti([outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(1)) '.nii'],[outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D.nii'],imNSSS/10,[],1,[]);
	for ii = 1:length(mTIvec)
		system(['rm ' outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_NS_' num2str(mTIvec(ii)) '.nii']);
		system(['rm ' outputPath '/Siemens_PASL_multiTI_GIFMI/analysis/Sub1/ASL_1/ASL4D_SS_' num2str(mTIvec(ii)) '.nii']);
	end
end

%% Specify the missing study parameters
% This is the most important part - it takes the imported JSONs and fills in all the study specific information that is currently still missing.
% Several parameters are automatically retrieved from the DICOMs, some are retrieved from Data_Par.m just below. And the rest needs to be filled in.
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
	
	% Fill in the generic text for the JSON 
	descriptionJSON{ii}.Name = fList{ii};
	descriptionJSON{ii}.BIDSVersion = '1.5.0';
	descriptionJSON{ii}.DatasetType = 'raw';
	descriptionJSON{ii}.License = 'RandomText';
	descriptionJSON{ii}.Authors = {'RandomText'};
	descriptionJSON{ii}.Acknowledgements = 'RandomText';
	descriptionJSON{ii}.HowToAcknowledge = 'Please cite this paper: https://www.ncbi.nlm.nih.gov/pubmed/001012092119281';
	descriptionJSON{ii}.Funding = {'RandomText'};
	descriptionJSON{ii}.EthicsApprovals = {'RandomText'};
	descriptionJSON{ii}.ReferencesAndLinks = {'RandomText'};
	descriptionJSON{ii}.DatasetDOI = 'RandomText';
	
	% Save the x-structure for the study
	%importStr{ii}.x = x;
	importStr{ii}.x = xASL_bids_parms2BIDS(x, [], 1, []);

	% Defaults to be overwritten
	importStr{ii}.par = [];
	importStr{ii}.par.VascularCrushing = false;
	importStr{ii}.par.LabelingLocationDescription = 'Random description';
	%importStr{ii}.par.LabelingOrientation = 'Random description';
	importStr{ii}.par.LabelingDistance = 40;

	importStr{ii}.par.ASLContext = '';
	switch (fList{ii})
		% FILLIN
		% - the basic information (ASLContext, LabelingType) needs to be filled for all
		%   - See 'Siemens_PCASL_2DEPI_Harmy_recombine_ASLscans' and 'Siemens_PCASL_3DGRASE_failed_APGEM2' for different types of Control/label order.
		%   - For PASL examples, see 'Siemens_PASL_3DGRASE_APGEM_1' or 'Siemens_PASL_2DEPI_noBsup2_EPAD'
		%   - Complete info for pCASL is in 'GE_PCASL_2DEPI_stripped_3CV'
		case 'StudyName'
			%importStr{ii}.par.ASLContext = '(Label+Control)*23';
			for cc = 1:23, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Siemens_PCASL_GIFMI'
			importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.LabelingPulseInterval = 0.4/1000;
			importStr{ii}.par.LabelingPulsesFlipAngle = 25;
			importStr{ii}.par.NumberSegments = 2;
			importStr{ii}.par.TotalAcquiredVolumes = [2 2];
			importStr{ii}.par.TotalReadoutTime = 0.0104;
			importStr{ii}.par.BackgroundSuppressionPulseTime = [0.85 0.1];
			importStr{ii}.par.BackgroundSuppressionNumberPulses = 2;
			
		case 'Siemens_PASL_multiTI_GIFMI'
			for cc = 1:10,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.NumberSegments = 2;
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 0;
			importStr{ii}.par.BackgroundSuppressionPulseTime = [0.85 0.1];
			importStr{ii}.par.BackgroundSuppressionNumberPulses = 2;
				
		case 'Siemens_PASL_singleTI_GIFMI'
			importStr{ii}.par.ASLContext = sprintf('%s\n',m0scanStr);
			for cc = 1:45,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 0.9;
			importStr{ii}.par.BolusCutOffTechnique = 'Q2TIPS';
			
		case 'Siemens_PCASL_2DEPI_Harmy_recombine_ASLscans'
			%importStr{ii}.par.ASLContext = '(Label+Control)*23';
			for cc = 1:46, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Siemens_PCASL_3DGRASE_failed_APGEM2'
			%importStr{ii}.par.ASLContext = 'M0+((Label+Control)*12)';
			importStr{ii}.par.ASLContext = sprintf('%s\n',m0scanStr);
			for cc = 1:12, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Siemens_PASL_3DGRASE_GENFI'
			importStr{ii}.par.LabelingType = 'PASL';
			%importStr{ii}.par.ASLContext = deltamStr;
			importStr{ii}.par.ASLContext = sprintf('%s\n',deltamStr);
			importStr{ii}.par.BolusCutOffFlag = false;

		case 'Siemens_PCASL_2DEPI_BioCog'
			%importStr{ii}.par.ASLContext = '(Label+Control)*44';
			for cc = 1:44, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Siemens_PASL_3DGRASE_APGEM_1'
			%importStr{ii}.par.ASLContext = 'Label+Control';
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',labelStr,controlStr);
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.LabelingLocationDescription = 'Labeling with FAIR';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 0;
			importStr{ii}.par.BolusCutOffTechnique = 'QUIPSSII';

		case 'Philips_PCASL_3DGRASE_Dent_example'
			%importStr{ii}.par.ASLContext = '(M0*2)+((Label+Control)*7)';
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',m0scanStr,m0scanStr);
			for cc = 1:7, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Siemens_PCASL_3DGRASE_Sleep_Oslo_trial'
			%importStr{ii}.par.ASLContext = 'M0+((Label+Control)*15)';
			importStr{ii}.par.ASLContext = sprintf('%s\n',m0scanStr);
			for cc = 1:15, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Siemens_PCASL_3DGRASE_RUNDMCSI_1774_asl_W38'
			%importStr{ii}.par.ASLContext = 'DeltaM*15';
			for cc = 1:92, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n',deltamStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Philips_PCASL_2DEPI_BioCog_Old'
			%importStr{ii}.par.ASLContext = '(Label+Control)*38';
			for cc = 1:38, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Philips_PCASL_2DEPI_CP_Tavi_HC'
			%importStr{ii}.par.ASLContext = '(Label+Control)*32';
			for cc = 1:32, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case {'Philips_PCASL_2DEPI_intera_FIND','Philips_PCASL_2DEPI_Ingenia_FIND','Philips_PCASL_2DEPI_intera_FIND_LL','Philips_PCASL_2DEPI_intera_FIND_multiPLD','Philips_PCASL_2DEPI_intera_FIND_QUASAR'}
			%importStr{ii}.par.ASLContext = '(Label+Control)*75';
			for cc = 1:75, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case {'Philips_PCASL_2DEPI_Achieva_Bsup_GENFI','Philips_PCASL_2DEPI_Achieva_noBsup_GENFI'}
			%importStr{ii}.par.ASLContext = '(Label+Control)*40';
			for cc = 1:40, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case {'Philips_PCASL_2DEPI_Bsup_EPAD1','Philips_PCASL_2DEPI_Bsup_EPAD2'}
			%importStr{ii}.par.ASLContext = '(Label+Control)*30';
			for cc = 1:30, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'Philips_PCASL_2DEPI_Bsup_EPAD3'
			importStr{ii}.par.LabelingType = 'PCASL';
			%importStr{ii}.par.ASLContext = deltamStr;
			importStr{ii}.par.ASLContext = sprintf('%s\n',deltamStr);
			importStr{ii}.par.LabelingEfficiency = 0.83;

		case {'Siemens_PASL_2DEPI_noBsup2_EPAD','Siemens_PASL_2DEPI_noBsup_EPAD'}
			%importStr{ii}.par.ASLContext = '(Label+Control)*31';
			for cc = 1:31, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 0.6;
			importStr{ii}.par.BolusCutOffTechnique = 'Q2TIPS';
			importStr{ii}.par.LabelingSlabThickness = 80;

		case {'Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD','Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD2'}
			%importStr{ii}.par.ASLContext = '(Label+Control)*10';
			for cc = 1:10, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = 0;
			importStr{ii}.par.BolusCutOffTechnique = 'QUIPSSII';

		case 'Siemens_PASL_3DGRASE_Prisma_Bsup_EPAD3'
			%importStr{ii}.par.ASLContext = '(Label+Control)*2';
			for cc = 1:2, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.LabelingType = 'PASL';
			importStr{ii}.par.BolusCutOffFlag = true;
			importStr{ii}.par.BolusCutOffDelayTime = [0 200 400]/1000;
			importStr{ii}.par.BolusCutOffTechnique = 'QUIPSS';
			importStr{ii}.par.LabelingSlabThickness = 60;

		case 'GE_PCASL_3Dspiral_Product_22q11'
			importStr{ii}.par.Units = 'mL/100g/min';
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',cbfStr,m0scanStr);
			importStr{ii}.par.LabelingType = 'PCASL';
			
		case 'GE_PCASL_3Dspiral_Product_GE'
			importStr{ii}.par.Units = 'mL/100g/min';
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',m0scanStr,cbfStr);
			importStr{ii}.par.LabelingType = 'PCASL';
		
		case {'Philips_PCASL_3DGRASE_Divers', 'GE_PCASL_3Dspiral_WIP_Oslo_AntiPsychotics_Old',...
			  'Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar'}
			importStr{ii}.par.Units = 'mL/100g/min';
			%importStr{ii}.par.ASLContext = cbfStr;
			importStr{ii}.par.ASLContext = sprintf('%s\n',cbfStr);
			importStr{ii}.par.LabelingType = 'PCASL';

		case 'GE_PCASL_3Dspiral_WIP_KCL_INOX'
			importStr{ii}.par.Units = 'mL/100g/min';
			%importStr{ii}.par.ASLContext = 'CBF+M0';
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',cbfStr,m0scanStr);
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.AcquisitionVoxelSize = [4 4 8];
			importStr{ii}.par.LabelingPulseAverageGradient = 0.6;
			importStr{ii}.par.LabelingPulseMaximumGradient = 6;
			importStr{ii}.par.LabelingPulseFlipAngle = 18;
			importStr{ii}.par.LabelingPulseDuration = 0.0005;
			importStr{ii}.par.PCASLType = 'balanced';

		case 'GE_PCASL_3Dspiral_UCL'
			importStr{ii}.par.ASLContext = sprintf('%s\n%s\n',m0scanStr,deltamStr);
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.AcquisitionVoxelSize = [4 4 8];
			
		case 'Philips_PCASL_2DEPI_UCL'
			for cc = 1:35, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',controlStr,labelStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';	
			importStr{ii}.par.AcquisitionVoxelSize = [3.75 3.75 5];
		
		case 'Siemens_PCASL_3DGRASE_UCL'
			for cc = 1:8, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',controlStr,labelStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';	
			importStr{ii}.par.AcquisitionVoxelSize = [3.4 3.4 4];
			
		case {'Philips_PCASL_2DEPI_Frontier','Philips_PCASL_2DEPI_Chili'}
			%importStr{ii}.par.ASLContext = '(Control+Label)*30';
			for cc = 1:30, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',controlStr,labelStr)];end
			importStr{ii}.par.LabelingType = 'PCASL';
			importStr{ii}.par.LabelingPulseAverageGradient = 0.7;
			importStr{ii}.par.LabelingPulseMaximumGradient = 7;
			importStr{ii}.par.LabelingPulseFlipAngle = 23;
			importStr{ii}.par.LabelingPulseDuration = 0.0005;
			importStr{ii}.par.PCASLType = 'balanced';

		case {'GE_PCASL_2DEPI_stripped_3CV','Philips_PCASL_2DEPI_stripped_3CV1','Philips_PCASL_2DEPI_stripped_3CV2','Siemens_PCASL_2DEPI_stripped_3CV'}
			%importStr{ii}.par.ASLContext = '(Control+Label)*70';
			for cc = 1:70, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',controlStr,labelStr)];end
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
		case {'Philips_PCASL_2DEPI_stripped_3CV2','Siemens_PCASL_2DEPI_stripped_3CV'}
			%importStr{ii}.par.AcquisitionDuration = 658;
			importStr{ii}.par.LabelingPulseInterval = 0.00115;
		case 'Philips_PCASL_3DGRASE_R5.4_PlusTopUp_TestKoen_FatSat_noDataPar'
			importStr{ii}.par.TotalReadoutTime = 0.0136;
	end

	% Process all the data and automatically fill in the missing parameters
	if strcmp(importStr{ii}.x.MRAcquisitionType,'2D')
		importStr{ii}.par.PulseSequenceType = '2D_EPI';
	else
		if strcmp(importStr{ii}.x.Manufacturer,'GE') || strcmp(importStr{ii}.x.Manufacturer,'GE_WIP') || strcmp(importStr{ii}.x.Manufacturer,'GE_product')
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
	if strcmp(importStr{ii}.par.LabelingType,'PASL')
		%importStr{ii}.par.LabelingDuration = 0;% importStr{ii}.x.LabelingDuration           = 1.800;  % for PASL this is TI1
		importStr{ii}.par.PostLabelingDelay = importStr{ii}.x.InitialPostLabelDelay;
		if importStr{ii}.par.BolusCutOffFlag
			importStr{ii}.par.BolusCutOffTimingSequence = importStr{ii}.x.LabelingDuration;
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
						warning('Not properly calculated');
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
		case 'Siemens_PASL_multiTI_GIFMI'
			importStr{ii}.par.PostLabelingDelay = [300 300 600 600 900 900 1200 1200 1500 1500 1800 1800 2100 2100 2400 2400 2700 2700 3000 3000]/1000;
			
		case 'Siemens_PASL_singleTI_GIFMI'
			importStr{ii}.par.VascularCrushing = true;
			importStr{ii}.par.VascularCrushingVenc = 100;
			
		case 'Siemens_PCASL_3DGRASE_failed_APGEM2'
			importStr{ii}.par.LabelingDuration = [0 repmat(1800,[1,24])]/1000;
			importStr{ii}.par.VascularCrushing = true;
			importStr{ii}.par.VascularCrushingVenc = 10;

		case 'Philips_PCASL_2DEPI_intera_FIND_LL'
			%importStr{ii}.par.ASLContext = '(Label*15+Control*15)*5';
			for cc = 1:5
				for dd=1:15,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n',labelStr)];end
				for dd=1:15,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n',controlStr)];end
			end
			importStr{ii}.par.FlipAngle = 25;
			importStr{ii}.par.LookLocker = true;
			importStr{ii}.par.PostLabelingDelay = repmat(250:250:3750,[1 10])/1000;

		case 'Philips_PCASL_2DEPI_intera_FIND_multiPLD'
			%importStr{ii}.par.ASLContext = '(Label+Control)*75';
			for cc = 1:75, importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n%s\n',labelStr,controlStr)];end
			importStr{ii}.par.PostLabelingDelay = repmat([500 500 1000 1000 1500 1500 1800 1800 2200 2200],[1 15])/1000;

		case 'Philips_PCASL_2DEPI_intera_FIND_QUASAR'
			%importStr{ii}.par.ASLContext = '(Label*15+Control*15)*5';
			for cc = 1:5
				for dd=1:15,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n',labelStr)];end
				for dd=1:15,importStr{ii}.par.ASLContext = [importStr{ii}.par.ASLContext sprintf('%s\n',controlStr)];end
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
for ii=1:length(fieldOrder)
	fieldOrderStruct.(fieldOrder{ii}) = '';
end

%% Go through all studies and check all the M0 and ASLs and modify the JSONs
% This step should be completely automatic, just taking the info filled above and using it to convert to full BIDS.
% No edits should be necessary unless there's something wrong that needs to be fixed.
for ii = 1:length(fList)
	% Make a copy of the par to the Flavors
	importStr{ii}.flavors = importStr{ii}.par;

	% Make an output directory
	if ~exist(fullfile(finalPath,importStr{ii}.dirName),'dir')
		mkdir(fullfile(finalPath,importStr{ii}.dirName));
	end
	
	spm_jsonwrite(fullfile(finalPath,importStr{ii}.dirName,'dataset_description.json'),descriptionJSON{ii});

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
				jsonLocal = finalJsonCheck(jsonLocal,fieldOrderStruct,removeEmptyFields);
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
					for ll=1:length(lRMFields)
						if strcmp(lRMFields{ll},fn{1})
							bCP = 0;
						end
					end
					for ll=1:length(lRMFieldsASL)
						if strcmp(lRMFieldsASL{ll},fn{1})
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
				if strcmp(importStr{ii}.x.M0,'separate_scan')
					if isfield(importStr{ii}.x,'M0PositionInASL4D') && (max(importStr{ii}.x.M0PositionInASL4D(:))>0)
						jsonLocal.M0 = true;
					elseif xASL_exist(fullfile(inSesPath,'M0.nii'))
						if length(fSes)>1
							%jsonLocal.M0 = fullfile(importStr{ii}.dirName,['sub-' subLabel],['ses-' sesLabel],'asl',['sub-' subLabel sesLabelUnd '_' m0scanStr '.nii.gz']);
							jsonLocal.M0 = fullfile('perf',['sub-' subLabel sesLabelUnd]);
							bJsonLocalM0isFile = 1;
						else
							%jsonLocal.M0 = fullfile(importStr{ii}.dirName,['sub-' subLabel],'asl',['sub-' subLabel sesLabelUnd '_M0Scan.nii.gz']);
							jsonLocal.M0 = fullfile('perf',['sub-' subLabel sesLabelUnd]);
							bJsonLocalM0isFile = 1;
						end
					else
						if ~isempty(strfind(importStr{ii}.par.ASLContext,m0scanStr))
							jsonLocal.M0 = true;
						else
							jsonLocal.M0 = false;
						end
					end
				else
					if strcmp(importStr{ii}.x.M0,'UseControlAsM0')
						jsonLocal.M0 = false;
					else
						if strcmp(importStr{ii}.x.M0,'no_background_suppression')
							jsonLocal.M0 = false;
						else
							jsonLocal.M0 = importStr{ii}.x.M0;
						end
					end
				end

				% Copy some things from the local JSON to the flavors
				jsonToFlavors = {'PulseSequenceType' 'PulseSequenceDetails' 'SliceTiming' 'ASLContext' 'LabelingType' 'LabelingDuration' 'PostLabelingDelay' 'BackgroundSuppression'...
					'BackgroundSuppressionNumberPulses' 'BackgroundSuppressionPulseTime' 'M0' 'AcquisitionVoxelSize'...
					'LookLocker' 'LabelingEfficiency' 'BolusCutOffFlag' 'BolusCutOffTimingSequence' 'BolusCutOffDelayTime' 'BolusCutOffTechnique'};
				for ll=1:length(jsonToFlavors)
					if isfield(jsonLocal,jsonToFlavors{ll})
						importStr{ii}.flavors.(jsonToFlavors{ll}) = jsonLocal.(jsonToFlavors{ll});
					end
				end
				
				% Remove the AslContext field and save it as a separate file
				fContext = fopen([aslOutLabel '_' aslcontextStr '.tsv'],'w+');
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
									jsonLocal.M0 = [jsonLocal.M0 nnStrOut '_' m0scanStr '.nii.gz'];
								end
							else
								if bJsonLocalM0isFile
									jsonLocal.M0 = [jsonLocal.M0 '_' m0scanStr '.nii.gz'];
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
							tagIntendedFor = fullfile('perf',['sub-' subLabel sesLabelUnd '_dir-ap' '_' m0scanStr '.nii.gz']);

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
								for ll=1:length(lRMFields)
									if strcmp(lRMFields{ll},fn{1})
										bCP = 0;
									end
								end
								for ll=1:length(lRMFieldsASL)
									if strcmp(lRMFieldsASL{ll},fn{1})
										bCP = 0;
									end
								end
								if bCP
									jsonM0Write.(fn{1}) = jsonM0.(fn{1});
								end
							end
							
							if isfield(jsonLocal,'SliceTiming')
								jsonM0Write.SliceTiming = jsonLocal.SliceTiming;
							else
								if isfield(jsonM0Write,'SliceTiming')
									jsonM0Write = rmfield(jsonM0Write,'SliceTiming');
								end
							end
							
							jsonM0Write.RepetitionTime = jsonM0.RepetitionTime;
							jsonM0Write.IntendedFor = [aslOutLabelRelative '_asl.nii.gz'];
							
							if ~isempty(tagPhaseEncodingDirection)
								jsonM0Write.PhaseEncodingDirection = tagPhaseEncodingDirection;
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
									xASL_io_SaveNifti(fullfile(inSesPath,['M0' nnStrIn '.nii']),fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd nnStrOut '_' m0scanStr '.nii.gz']),imM0,[],1,[]);
								else
									xASL_io_SaveNifti(fullfile(inSesPath,['M0' nnStrIn '.nii']),fullfile(outSesPath,'fmap',['sub-' subLabel sesLabelUnd nnStrOut '_' m0scanStr '.nii.gz']),imM0,[],1,[]);
								end
							else
								% Copy the M0
								if nn == 1
									xASL_Copy(fullfile(inSesPath,['M0' nnStrIn '.nii']),...
										fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd nnStrOut '_' m0scanStr '.nii.gz']));
								else
									xASL_Copy(fullfile(inSesPath,['M0' nnStrIn '.nii']),...
										fullfile(outSesPath,'fmap',['sub-' subLabel sesLabelUnd nnStrOut '_' m0scanStr '.nii.gz']));
								end
							end
							% Save JSON to new dir
							jsonM0Write = finalJsonCheck(jsonM0Write,fieldOrderStruct,removeEmptyFields);
							if nn == 1
								spm_jsonwrite(fullfile(outSesPath,'perf',['sub-' subLabel sesLabelUnd nnStrOut '_' m0scanStr '.json']),jsonM0Write);
							else
								spm_jsonwrite(fullfile(outSesPath,'fmap',['sub-' subLabel sesLabelUnd nnStrOut '_' m0scanStr '.json']),jsonM0Write);
							end
						end
					end
				else
					if bJsonLocalM0isFile
						jsonLocal.M0 = [jsonLocal.M0 '.nii.gz'];
					end
				end
				% Save JSON to new dir
				jsonLocal = finalJsonCheck(jsonLocal,fieldOrderStruct,removeEmptyFields);
				spm_jsonwrite([aslOutLabel '_asl.json'],jsonLocal);
				
			end
		end
	end
end

%% Export the fully anonymized datasets for public sharing
pthVec = {'GE_PCASL_3Dspiral_UCL' 'Siemens_PCASL_3DGRASE_UCL' 'Philips_PCASL_2DEPI_UCL'};
for ii = 1:3
	xASL_Copy(fullfile(finalPath,pthVec{ii}),fullfile(anonymPath,pthVec{ii}));
	xASL_spm_deface(fullfile(anonymPath,pthVec{ii},'sub-Sub103','anat','sub-Sub103_T1w.nii'),true);
	gzip(fullfile(anonymPath,pthVec{ii},'sub-Sub103','anat','sub-Sub103_T1w.nii'));
	delete(fullfile(anonymPath,pthVec{ii},'sub-Sub103','anat','sub-Sub103_T1w.nii'));
end

pthVec = {'Siemens_PASL_multiTI_GIFMI','Siemens_PASL_singleTI_GIFMI','Siemens_PCASL_GIFMI'};
for ii = 1:3
	xASL_Copy(fullfile(finalPath,pthVec{ii}),fullfile(anonymPath,pthVec{ii}));
	xASL_spm_deface(fullfile(anonymPath,pthVec{ii},'sub-Sub1','anat','sub-Sub1_T1w.nii'),true);
	gzip(fullfile(anonymPath,pthVec{ii},'sub-Sub1','anat','sub-Sub1_T1w.nii'));
	delete(fullfile(anonymPath,pthVec{ii},'sub-Sub1','anat','sub-Sub1_T1w.nii'));
end

%% Export for ASL JSON flavours
% This do not need to be run for BIDS conversion - it only save the info for all processed datasets.
fCSVOut = fopen(fullfile(finalPath,'flavours.csv'),'w');

fieldsCSV = {'PulseSequenceType' 'PulseSequenceDetails' 'SliceTiming' 'ASLContext' 'LabelingType' 'LabelingDuration' 'PostLabelingDelay' 'BackgroundSuppression' 'M0' 'AcquisitionVoxelSize'...
'BackgroundSuppressionNumberPulses' 'BackgroundSuppressionPulseTime' 'LookLocker' 'LabelingEfficiency' 'BolusCutOffFlag' 'BolusCutOffTimingSequence' 'BolusCutOffDelayTime' 'BolusCutOffTechnique'};

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

%%
function jsonOut = finalJsonCheck(jsonIn,fieldOrderStruct,removeEmptyFields)
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

if isfield(jsonOut,'PhaseEncodingAxis') && strcmp(jsonOut.Manufacturer,'Philips')
	jsonOut.PhaseEncodingDirection = jsonOut.PhaseEncodingAxis;
	jsonOut = rmfield(jsonOut,'PhaseEncodingAxis');
end

for ii = 1:length(removeEmptyFields)
	if isfield(jsonOut,removeEmptyFields{ii})
		if isempty(jsonOut.(removeEmptyFields{ii}))
			jsonOut = rmfield(jsonOut,removeEmptyFields{ii});
		end
	end
end

% And sort the fields
jsonOut = xASL_adm_OrderFields(jsonOut,fieldOrderStruct);
end