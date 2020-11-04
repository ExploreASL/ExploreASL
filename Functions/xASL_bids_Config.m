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
	'PhilipsRescaleSlope'  'PhilipsRWVSlope' 'PhilipsScaleSlope' 'PhilipsRescaleIntercept' 'UsePhilipsFloatNotDisplayScaling' 'RWVSlope' ...
	'PhilipsRWVIntercept' 'RescaleSlopeOriginal' 'RescaleSlope'    'MRScaleSlope'      'RescaleIntercept',... % Scales slopes have to be applied and removed
	'Modality', 'ImagingFrequency', 'PatientPosition','MRAcquisitionType','ImageType','PhaseEncodingPolarityGE',...
	'SeriesNumber','AcquisitionTime','AcquisitionNumber','SliceThickness','SpacingBetweenSlices','SAR',...
	'PercentPhaseFOV','AcquisitionMatrixPE','ReconMatrixPE','PixelBandwidth','ImageOrientationPatientDICOM',...
	'InPlanePhaseEncodingDirectionDICOM','ConversionSoftware','ConversionSoftwareVersion','AcquisitionMatrix',...
	'EchoTrainLength','PhaseEncodingSteps','BodyPartExamined','ShimSetting','TxRefAmp','PhaseResolution',...
	'RefLinesPE','BandwidthPerPixelPhaseEncode','ImageComments','ConsistencyInfo','WipMemBlock','Interpolation2D',...
	'SaturationStopTime','BaseResolution','DerivedVendorReportedEchoSpacing','RawImage','PhaseOversampling','BolusDuration'}; % Fields to exclude as not defined in BIDS

% A list of fields to remove the ASL-BIDS JSON files only
bidsPar.listFieldsRemoveASL = {'InversionTime','LabelOffset','PostLabelDelay','NumRFBlocks','RFGap','MeanGzx10','PhiAdjust'}; % Fields to exclude from ASL only

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
bidsPar.ASLfields.Required = {};		  
bidsPar.ASLfields.Recommended = {};		  
bidsPar.ASLfields.Optional = {};		  

% Conditional dependencies
bidsPar.ASLCondition{1}.field = '';
bidsPar.ASLCondition{1}.value = '';
bidsPar.ASLCondition{1}.RequiredFilled = {};
bidsPar.ASLCondition{1}.RequiredEmpty = {};
bidsPar.ASLCondition{1}.RecommendedFilled = {};

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