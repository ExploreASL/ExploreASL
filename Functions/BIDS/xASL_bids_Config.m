function bidsPar = xASL_bids_Config()
%xASL_bids_Config Creates several config variables specifying on how to convert and save Dicom to BIDS
%
% FORMAT: bidsPar = xASL_bids_Config()
%
% INPUT:
%   n/a
%
% OUTPUT: bidsPar - a structure containing all the relevant lists
%             listFieldsRemoveGeneral - a list of fields to be removed from all BIDS JSON files
%             listFieldsRemoveASL     - A list of fields to remove the ASL-BIDS JSON files only
%             listFieldOrder          - Gives the correct order of fields to be saved in JSON so that it corresponds to the BIDS definition
%             listRemoveIfEmpty       - List of fields to be removed if they are empty
%             listAnatTypes           - A list of anatomical scan-types to include
%             string*                 - Defined strings for certain ASL-BIDS keywords (stringAslContext, stringLabel, stringControl, stringM0scan, stringCbf, stringDeltaM)
%             datasetDescription      - A list of req, rec, and opt fields to be put to dataset_description.json - they all need to be predefined by the user
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Creates several structures necessary for configuring the DICOM to BIDS conversion and saving of BIDS JSON files and NII structure.
%
% EXAMPLE:     n/a
% __________________________________
% Copyright 2015-2023 ExploreASL

% BIDS version is hard-coded
bidsPar.BIDSVersion = '1.6.0';

% a list of fields to be removed from all BIDS JSON files
bidsPar.listFieldsRemoveGeneral = {'ProcedureStepDescription' 'SeriesDescription' 'ProtocolName'... % Anonymization
	'PhilipsRescaleSlope' 'PhilipsRWVSlope' 'PhilipsScaleSlope' 'PhilipsRescaleIntercept' 'UsePhilipsFloatNotDisplayScaling' 'RWVSlope',...
	'PhilipsRWVIntercept' 'RescaleSlopeOriginal' 'RescaleSlope' 'MRScaleSlope' 'RescaleIntercept', 'RWVIntercept',... % Scales slopes have to be applied and removed
	'Modality', 'ImagingFrequency', 'PatientPosition','ImageType','PhaseEncodingPolarityGE','EstimatedEffectiveEchoSpacing','EstimatedTotalReadoutTime',...
	'SeriesNumber','AcquisitionTime','AcquisitionNumber','SliceThickness','SpacingBetweenSlices','SAR','InternalPulseSequenceName','PhaseEncodingStepsNoPartialFourier', 'PartialFourierEnabled',...
	'Rows' 'Columns' 'InPlanePhaseEncodingDirection' 'NumberOfTemporalPositions','GELabelingType', 'PulseSequenceName','AcquisitionContrast',...
	'PhilipsNumberTemporalScans','TemporalPositionIdentifier','PhilipsLabelControl','NumberOfTemporalPositions','PhoenixProtocol','PhoenixAnalyzed','SiemensSliceTime',...
	'PercentPhaseFOV','AcquisitionMatrixPE','ReconMatrixPE','PixelBandwidth','ImageOrientationPatientDICOM','ComplexImageComponent',...
	'InPlanePhaseEncodingDirectionDICOM','ConversionSoftware','ConversionSoftwareVersion','AcquisitionMatrix','GELabelingDuration',...
	'EchoTrainLength','PhaseEncodingSteps','BodyPartExamined','ShimSetting','TxRefAmp','PhaseResolution','InstanceNumber', 'WaterFatShift',...
	'RefLinesPE','BandwidthPerPixelPhaseEncode','ImageComments','ConsistencyInfo','WipMemBlock','Interpolation2D','TotalAcquiredVolumes', 'NumberOfExcitations', 'PrescanReuseString', ...
	'SaturationStopTime','BaseResolution','DerivedVendorReportedEchoSpacing','RawImage', 'T1', 'PercentSampling', 'ImageOrientationText', 'PhaseOversampling','BolusDuration','scaleFactor','EchoNumber'}; % Fields to exclude as not defined in BIDS

% A list of fields to remove the ASL.json files as the field is not defined in ASL-BIDS 
bidsPar.listFieldsRemoveASL = {'SliceReadoutTime','RepetitionTime','InversionTime','LabelOffset','PostLabelDelay','NumRFBlocks',...
	'GELabelingDuration','RFGap','MeanGzx10','PhiAdjust','M0','LabelingType','ScanType','SequenceType'}; 

% A list of fields to remove the M0.json files as the field is not defined in ASL-BIDS
bidsPar.listFieldsRemoveM0 = {'TotalAcquiredPairs', 'ArterialSpinLabelingType','LabelingDuration','PostLabelingDelay','LabelingLocationDescription','LabelingOrientation','LabelingDistance','BackgroundSuppression','BackgroundSuppressionNumberPulses','BackgroundSuppressionPulseTime',...
	                          'LabelingEfficiency','LabelingDuration','PCASLType','CASLType','LabelingPulseAverageGradient','LabelingPulseMaximumGradient','LabelingPulseAverageB1','LabelingPulseDuration','LabelingPulseFlipAngle','LabelingPulseInterval','BolusCutOffFlag','PASLType',...
							  'LabelingSlabThickness','BolusCutOffDelayTime','BolusCutOffTechnique',};

% A list of fields to remove the all non-ASL, non-M0 JSON files as the field is not defined in BIDS
bidsPar.listFieldsRemoveNonASL = {'TotalAcquiredPairs'}; 

% Gives the correct order of fields to be saved in JSON so that it corresponds to the BIDS definition
% This is not mandatory, just makes things more accessible
bidsPar.listFieldOrder = {'Manufacturer','ManufacturersModelName','DeviceSerialNumber','StationName','SoftwareVersions','HardcopyDeviceSoftwareVersion',...%BIDS fields
	          'MagneticFieldStrength','ReceiveCoilName','ReceiveCoilActiveElements','GradientSetType','MRTransmitCoilSequence','MatrixCoilMode',...
			  'CoilCombinationMethod','PulseSequenceType','ScanningSequence','SequenceVariant','ScanOptions','SequenceName','PulseSequenceDetails',...
			  'NonlinearGradientCorrection','MRAcquisitionType','MTState','MTOffsetFrequency','MTPulseBandwidth','MTNumberOfPulses',...
			  'MTPulseShape','MTPulseDuration','SpoilingState','SpoilingType','SpoilingRFPhaseIncrement','SpoilingGradientMoment',...
			  'SpoilingGradientDuration','NumberShots','ParallelReductionFactorInPlane','ParallelReductionOutOfPlane', 'ParallelAcquisitionTechnique','PartialFourier',...
			  'PartialFourierDirection','PhaseEncodingDirection','EffectiveEchoSpacing','TotalReadoutTime','EchoTime','MixingTime','InversionTime',...
			  'SliceTiming','SliceEncodingDirection','DwellTime','FlipAngle','NegativeContrast','MultibandAccelerationFactor','AnatomicalLandmarkCoordinates',...
			  'InstitutionName','InstitutionalDepartmentName','InstitutionAddress',...
			  'ContrastBolusIngredient','RepetitionTimeExcitation','RepetitionTimePreparation',...
			  'RepetitionTime','VolumeTiming','TaskName','NumberOfVolumesDiscardedByScanner',...
			  'NumberOfVolumesDiscardedByUser','Units','DelayTime','AcquisitionDuration','DelayAfterTrigger',... 
			  'Instructions','TaskDescription','CogAtlasID','CogPOID','MultipartID',...% And ASL-BIDS fields
			  'ArterialSpinLabelingType','PostLabelingDelay','BackgroundSuppression','M0Type',... % ASL REQUIRED FIELDS
			  'VascularCrushing','AcquisitionVoxelSize','M0Estimate','TotalAcquiredPairs',...
			  'BackgroundSuppressionNumberPulses','BackgroundSuppressionPulseTime','VascularCrushingVENC','LabelingLocationDescription',...
			  'LabelingOrientation','LabelingDistance','LookLocker','LabelingEfficiency','LabelingDuration','PCASLType','CASLType',...
			  'LabelingPulseAverageGradient','LabelingPulseMaximumGradient','LabelingPulseAverageB1','LabelingPulseDuration','LabelingPulseFlipAngle',...
			  'LabelingPulseInterval','BolusCutOffFlag','PASLType','LabelingSlabThickness','BolusCutOffDelayTime',...
			  'BolusCutOffTechnique','EchoTime1','EchoTime2','EchoTime'};

% Definition of ASL fields		  
bidsPar.ASLfields.Required = {'TotalAcquiredPairs','RepetitionTimePreparation','EchoTime','MagneticFieldStrength','MRAcquisitionType','ArterialSpinLabelingType','PostLabelingDelay',...
	'BackgroundSuppression','M0Type'};		  
bidsPar.ASLfields.Recommended = {'PulseSequenceDetails','Manufacturer','SliceEncodingDirection','FlipAngle','VascularCrushing',...
	'AcquisitionVoxelSize','LabelingLocationDescription','LabelingOrientation','LabelingDistance',};		  
bidsPar.ASLfields.Optional = {'LookLocker','LabelingEfficiency','Units','PulseSequenceType'};		  
bidsPar.M0fields.Required = {'RepetitionTimePreparation','EchoTime','MagneticFieldStrength','MRAcquisitionType'};

% Defined strings for certain ASL-BIDS keywords
bidsPar.stringASL = 'asl';
bidsPar.stringAslContext = 'aslcontext';
bidsPar.stringControl = 'control';
bidsPar.stringLabel = 'label';
bidsPar.stringM0scan = 'm0scan';
bidsPar.stringCbf = 'cbf';
bidsPar.stringDeltaM = 'deltam';
bidsPar.stringM0Separate = 'Separate';
bidsPar.stringM0Included = 'Included';
bidsPar.stringM0Estimate = 'Estimate';
bidsPar.stringM0Absent   = 'Absent';
bidsPar.stringPerfusion = 'perf';
bidsPar.stringFmap = 'fmap';

bidsPar.sidecarName = {'.json' '_aslcontext.tsv' '_labeling.jpg'};
bidsPar.sidecarRequired =[1 0 0];
bidsPar.sidecarTypeSpecific = {'no' 'asl' 'asl'};
bidsPar.sidecarSuffixType = [1 0 0]; % specifies if the sidecar suffix keeps the scantype (e.g. yes for *_asl.json, not for *_aslcontext

% Conditional dependencies
% RepetitionTime and VolumeTiming are mutually exclusive
bidsPar.ASLCondition{1}.field = 'RepetitionTime';
bidsPar.ASLCondition{1}.value = '';
bidsPar.ASLCondition{1}.RequiredFilled = {};
bidsPar.ASLCondition{1}.RequiredEmpty = {'VolumeTiming'};
bidsPar.ASLCondition{1}.RecommendedFilled = {};

bidsPar.ASLCondition{2}.field = 'VolumeTiming';
bidsPar.ASLCondition{2}.value = '';
bidsPar.ASLCondition{2}.RequiredFilled = {'AcquisitionDuration','SliceTiming'};
bidsPar.ASLCondition{2}.RequiredEmpty = {'RepetitionTime'};
bidsPar.ASLCondition{2}.RecommendedFilled = {};

bidsPar.ASLCondition{3}.field = 'MRAcquisitionType';
bidsPar.ASLCondition{3}.value = '2D';
bidsPar.ASLCondition{3}.RequiredFilled = {'SliceTiming'};
bidsPar.ASLCondition{3}.RequiredEmpty = {};
bidsPar.ASLCondition{3}.RecommendedFilled = {};

bidsPar.ASLCondition{4}.field = 'LookLocker';
bidsPar.ASLCondition{4}.value = true;
bidsPar.ASLCondition{4}.RequiredFilled = {'FlipAngle'};
bidsPar.ASLCondition{4}.RequiredEmpty = {};
bidsPar.ASLCondition{4}.RecommendedFilled = {};

bidsPar.ASLCondition{5}.field = 'BackgroundSuppression';
bidsPar.ASLCondition{5}.value = true;
bidsPar.ASLCondition{5}.RequiredFilled = {};
bidsPar.ASLCondition{5}.RequiredEmpty = {};
bidsPar.ASLCondition{5}.RecommendedFilled = {'BackgroundSuppressionNumberPulses','BackgroundSuppressionPulseTime'};

bidsPar.ASLCondition{6}.field = 'VascularCrushing';
bidsPar.ASLCondition{6}.value = true;
bidsPar.ASLCondition{6}.RequiredFilled = {};
bidsPar.ASLCondition{6}.RequiredEmpty = {};
bidsPar.ASLCondition{6}.RecommendedFilled = {'VascularCrushingVENC'};

bidsPar.ASLCondition{7}.field = 'ArterialSpinLabelingType';
bidsPar.ASLCondition{7}.value = 'PCASL';
bidsPar.ASLCondition{7}.RequiredFilled = {'LabelingDuration'};
bidsPar.ASLCondition{7}.RequiredEmpty = {'CASLType','PASLType','LabelingSlabThickness','BolusCutOffFlag','BolusCutOffDelayTime',...
	'BolusCutOffTechnique'};
bidsPar.ASLCondition{7}.RecommendedFilled = {'PCASLType','LabelingPulseAverageGradient','LabelingPulseMaximumGradient',...
	'LabelingPulseAverageB1','LabelingPulseDuration','LabelingPulseFlipAngle','LabelingPulseInterval'};

bidsPar.ASLCondition{8}.field = 'ArterialSpinLabelingType';
bidsPar.ASLCondition{8}.value = '^CASL';
bidsPar.ASLCondition{8}.RequiredFilled = {'LabelingDuration'};
bidsPar.ASLCondition{8}.RequiredEmpty = {'PCASLType','PASLType','LabelingSlabThickness','BolusCutOffFlag','BolusCutOffDelayTime',...
	'BolusCutOffTechnique'};
bidsPar.ASLCondition{8}.RecommendedFilled = {'CASLType','LabelingPulseAverageGradient','LabelingPulseMaximumGradient',...
	'LabelingPulseAverageB1','LabelingPulseDuration','LabelingPulseFlipAngle','LabelingPulseInterval'};

bidsPar.ASLCondition{9}.field = 'ArterialSpinLabelingType';
bidsPar.ASLCondition{9}.value = 'PASL';
bidsPar.ASLCondition{9}.RequiredFilled = {'BolusCutOffFlag'};
bidsPar.ASLCondition{9}.RequiredEmpty = {'PCASLType','CASLType','LabelingPulseAverageGradient','LabelingPulseMaximumGradient',...
	'LabelingPulseAverageB1','LabelingPulseDuration','LabelingPulseFlipAngle','LabelingPulseInterval','LabelingDuration'};
bidsPar.ASLCondition{9}.RecommendedFilled = {'PASLType','LabelingSlabThickness'};

bidsPar.ASLCondition{10}.field = 'BolusCutOffFlag';
bidsPar.ASLCondition{10}.value = true;
bidsPar.ASLCondition{10}.RequiredFilled = {'BolusCutOffDelayTime','BolusCutOffTechnique'};
bidsPar.ASLCondition{10}.RequiredEmpty = {};
bidsPar.ASLCondition{10}.RecommendedFilled = {};

bidsPar.ASLCondition{11}.field = 'M0Type';
bidsPar.ASLCondition{11}.value = bidsPar.stringM0Estimate;
bidsPar.ASLCondition{11}.RequiredFilled = {'M0Estimate'};
bidsPar.ASLCondition{11}.RequiredEmpty = {};
bidsPar.ASLCondition{11}.RecommendedFilled = {};

bidsPar.ASLCondition{12}.field = 'BolusCutOffFlag';
bidsPar.ASLCondition{12}.value = false;
bidsPar.ASLCondition{12}.RequiredFilled = {};
bidsPar.ASLCondition{12}.RequiredEmpty = {'BolusCutOffDelayTime','BolusCutOffTechnique'};
bidsPar.ASLCondition{12}.RecommendedFilled = {};

% List of fields to be removed if they are empty
bidsPar.listRemoveIfEmpty = {'EffectiveEchoSpacing','TotalReadoutTime'};

% A list of anatomical scan-types to include
bidsPar.listAnatTypes = {'T1' 'T1w' 'T2w' 'T1c' 'FLAIR'}; 

% A list of directories for BIDS to Legacy conversion
bidsPar.BIDS2LegacyFolderConfiguration =...
	{'type'   'modality' 'foldernames' 'filenames'      'run_locations'   'run_1_index';...
	 'T1w'    'anat'     ''            'T1'             'file'            false;...
	 'T2w'    'anat'     ''            'T2'             'file'            false;...
	 'T1'     'anat'     ''            'T1'             'file'            false;...
	 'T1c'    'anat'     ''            'T1c'            'file'            false;...
	 'FLAIR'  'anat'     ''            'FLAIR'          'file'            false;...
	 'asl'    'perf'     'ASL'         'ASL4D'          'folder'          true;...
	 'm0scan' 'perf'     'ASL'         'M0'             'folder'          true;...
	 'm0scan' 'fmap'     'ASL'         'M0_RevPE'       'folder'          true;};

% Definition of the dataset_description.json fields
% https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
bidsPar.datasetDescription.filename = 'dataset_description';
bidsPar.datasetDescription.Required = {'Name' 'BIDSVersion'};
bidsPar.datasetDescription.Recommended = {'HEDVersion' 'DatasetType' 'License'};
bidsPar.datasetDescription.Optional = {'Authors' 'Acknowledgements' 'HowToAcknowledge' 'Funding'...
	                                   'EthicsApprovals' 'ReferencesAndLinks' 'DatasetDOI'};

end
