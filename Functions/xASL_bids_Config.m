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
%             str*                    - Defined strings for certain ASL-BIDS keywords (strAslContext, strLabel, strControl, strM0scan, strCbf, strDeltaM)
%             datasetDescription      - A list of req, rec, and opt fields to be put to dataset_description.json - they all need to be predefined by the user
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Creates several structures necessary for configuring the DICOM to BIDS conversion and saving of BIDS JSON files and NII structure.
% __________________________________
% Copyright 2015-2020 ExploreASL

% BIDS version is hard-coded
bidsPar.BIDSVersion = '1.5.0';

% a list of fields to be removed from all BIDS JSON files
bidsPar.listFieldsRemoveGeneral = {'InstitutionName' 'InstitutionalDepartmentName' 'InstitutionAddress' 'DeviceSerialNumber' 'StationName'...
	'ProcedureStepDescription' 'SeriesDescription' 'ProtocolName'... % Anonymization
	'PhilipsRescaleSlope' 'PhilipsRWVSlope' 'PhilipsScaleSlope' 'PhilipsRescaleIntercept' 'UsePhilipsFloatNotDisplayScaling' 'RWVSlope',...
	'PhilipsRWVIntercept' 'RescaleSlopeOriginal' 'RescaleSlope' 'MRScaleSlope' 'RescaleIntercept', 'RWVIntercept',... % Scales slopes have to be applied and removed
	'Modality', 'ImagingFrequency', 'PatientPosition','MRAcquisitionType','ImageType','PhaseEncodingPolarityGE',...
	'SeriesNumber','AcquisitionTime','AcquisitionNumber','SliceThickness','SpacingBetweenSlices','SAR',...
	'PercentPhaseFOV','AcquisitionMatrixPE','ReconMatrixPE','PixelBandwidth','ImageOrientationPatientDICOM',...
	'InPlanePhaseEncodingDirectionDICOM','ConversionSoftware','ConversionSoftwareVersion','AcquisitionMatrix',...
	'EchoTrainLength','PhaseEncodingSteps','BodyPartExamined','ShimSetting','TxRefAmp','PhaseResolution',...
	'RefLinesPE','BandwidthPerPixelPhaseEncode','ImageComments','ConsistencyInfo','WipMemBlock','Interpolation2D',...
	'SaturationStopTime','BaseResolution','DerivedVendorReportedEchoSpacing','RawImage','PhaseOversampling','BolusDuration'}; % Fields to exclude as not defined in BIDS

% A list of fields to remove the ASL-BIDS JSON files only
bidsPar.listFieldsRemoveASL = {'SliceReadoutTime','InversionTime','LabelOffset','PostLabelDelay','NumRFBlocks',...
	'GELabelingDuration','RFGap','MeanGzx10','PhiAdjust'}; % Fields to exclude from ASL only

% Gives the correct order of fields to be saved in JSON so that it corresponds to the BIDS definition
% This is not mandatory, just makes things more accessible
bidsPar.listFieldOrder = {'Manufacturer','ManufacturersModelName','DeviceSerialNumber','StationName','SoftwareVersions','HardcopyDeviceSoftwareVersion',...%BIDS fields
	          'MagneticFieldStrength','ReceiveCoilName','ReceiveCoilActiveElements','GradientSetType','MRTransmitCoilSequence','MatrixCoilMode',...
			  'CoilCombinationMethod','PulseSequenceType','ScanningSequence','SequenceVariant','ScanOptions','SequenceName','PulseSequenceDetails',...
			  'NonlinearGradientCorrection','NumberShots','ParallelReductionFactorInPlane','ParallellAcquisitionTechnique','PartialFourier',...
			  'PartialFourierDirection','PhaseEncodingDirection','EffectiveEchoSpacing','TotalReadoutTime','EchoTime','InversionTime','SliceTiming',...
			  'SliceEncodingDirection','DwellTime','FlipAngle','MultibandAccelerationFactor','NegativeContrast','AnatomicalLandmarkCoordinates',...
			  'RepetitionTime','VolumeTiming','TaskName','Units',... % And ASL-BIDS fields
			  'LabelingType','PostLabelingDelay','BackgroundSuppression','M0','VascularCrushing','AcquisitionVoxelSize','TotalAcquiredVolumes',...
			  'BackgroundSuppressionNumberPulses','BackgroundSuppressionPulseTime','VascularCrushingVenc','LabelingLocationDescription',...
			  'LabelingOrientation','LabelingDistance','LookLocker','LabelingEfficiency','PCASLType','CASLType','LabelingDuration',...
			  'LabelingPulseAverageGradient','LabelingPulseMaximumGradient','LabelingPulseAverageB1','LabelingPulseDuration','LabelingPulseFlipAngle',...
			  'LabelingPulseInterval','PASLType','LabelingSlabThickness','BolusCutOffFlag','BolusCutOffDelayTime',...
			  'BolusCutOffTechnique'};

% Definition of ASL fields		  
bidsPar.ASLfields.Required = {'EchoTime','MagneticFieldStrength','PulseSequenceType','LabelingType','PostLabelingDelay',...
	'BackgroundSuppression','M0'};		  
bidsPar.ASLfields.Recommended = {'PulseSequenceDetails','Manufacturer','SliceEncodingDirection','FlipAngle','VascularCrushing',...
	'AcquisitionVoxelSize','LabelingLocationDescription','LabelingOrientation','LabelingDistance',};		  
bidsPar.ASLfields.Optional = {'TotalAcquiredVolumes','LookLocker','LabelingEfficiency','Units'};		  

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

bidsPar.ASLCondition{3}.field = 'PulseSequenceType';
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
bidsPar.ASLCondition{6}.RecommendedFilled = {'VascularCrushingVenc'};

bidsPar.ASLCondition{7}.field = 'LabelingType';
bidsPar.ASLCondition{7}.value = 'PCASL';
bidsPar.ASLCondition{7}.RequiredFilled = {'LabelingDuration'};
bidsPar.ASLCondition{7}.RequiredEmpty = {'CASLType','PASLType','LabelingSlabThickness','BolusCutOffFlag','BolusCutOffDelayTime',...
	'BolusCutOffTechnique'};
bidsPar.ASLCondition{7}.RecommendedFilled = {'PCASLType','LabelingPulseAverageGradient','LabelingPulseMaximumGradient',...
	'LabelingPulseAverageB1','LabelingPulseDuration','LabelingPulseFlipAngle','LabelingPulseInterval'};

bidsPar.ASLCondition{8}.field = 'LabelingType';
bidsPar.ASLCondition{8}.value = '^CASL';
bidsPar.ASLCondition{8}.RequiredFilled = {'LabelingDuration'};
bidsPar.ASLCondition{8}.RequiredEmpty = {'PCASLType','PASLType','LabelingSlabThickness','BolusCutOffFlag','BolusCutOffDelayTime',...
	'BolusCutOffTechnique'};
bidsPar.ASLCondition{8}.RecommendedFilled = {'CASLType','LabelingPulseAverageGradient','LabelingPulseMaximumGradient',...
	'LabelingPulseAverageB1','LabelingPulseDuration','LabelingPulseFlipAngle','LabelingPulseInterval'};

bidsPar.ASLCondition{9}.field = 'LabelingType';
bidsPar.ASLCondition{9}.value = 'PASL';
bidsPar.ASLCondition{9}.RequiredFilled = {'BolusCutOffFlag'};
bidsPar.ASLCondition{9}.RequiredEmpty = {'PCASLType','CASLType','LabelingPulseAverageGradient','LabelingPulseMaximumGradient',...
	'LabelingPulseAverageB1','LabelingPulseDuration','LabelingPulseFlipAngle','LabelingPulseInterval'};
bidsPar.ASLCondition{9}.RecommendedFilled = {'PASLType','LabelingSlabThickness'};

bidsPar.ASLCondition{10}.field = 'BolusCutOffFlag';
bidsPar.ASLCondition{10}.value = true;
bidsPar.ASLCondition{10}.RequiredFilled = {};
bidsPar.ASLCondition{10}.RequiredEmpty = {};
bidsPar.ASLCondition{10}.RecommendedFilled = {'BolusCutOffDelayTime','BolusCutOffTechnique'};

% List of fields to be removed if they are empty
bidsPar.listRemoveIfEmpty = {'EffectiveEchoSpacing','TotalReadoutTime'};

% A list of anatomical scan-types to include
bidsPar.listAnatTypes = {'T1w' 'T2w' 'FLAIR'}; 

% Defined strings for certain ASL-BIDS keywords
bidsPar.strAslContext = 'aslcontext';
bidsPar.strLabel = 'label';
bidsPar.strControl = 'control';
bidsPar.strM0scan = 'm0scan';
bidsPar.strCbf = 'cbf';
bidsPar.strDeltaM = 'deltam';		  

% Definition of the dataset_description.json fields
% https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
bidsPar.datasetDescription.filename = 'dataset_description';
bidsPar.datasetDescription.Required = {'Name' 'BIDSVersion'};
bidsPar.datasetDescription.Recommended = {'HEDVersion' 'DatasetType' 'License'};
bidsPar.datasetDescription.Optional = {'Authors' 'Acknowledgements' 'HowToAcknowledge' 'Funding'...
	                                   'EthicsApprovals' 'ReferencesAndLinks' 'DatasetDOI'};

end