function header = xASL_io_DcmtkRead(filepath, bPixel, bTryDCMTK, bSkipNonDicoms)
% ------------------------------------------------------------------------------------------------------------------------------------------------------
% SHORT DESCRIPTION: Reads DICOM headers using DCMTK
%
% FORMAT: header = xASL_io_DcmtkRead(filepath[, bPixel, bTryDCMTK, bSkipNonDicoms])
%
% INPUT:
%         filepath          - full path to the DICOM file (REQUIRED, STRING)
%         bPixel            - read pixel data (OPTIONAL, BOOLEAN, DEFAULT = false)
%         bTryDCMTK         - try using DCMTK, if false then directly use SPM (OPTIONAL, BOOLEAN, DEFAULT = true)
%         bSkipNonDicoms    - catch any errors reading a file and instead of crashing the code, skip the file
%                             this is useful e.g. when parsing multiple folders that are partly filled with DICOMs
%                             e.g., with xASL_adm_SortDicomToFolders
%                             as DICOM files do not necessarily have a unique extension
%                             (OPTIONAL, DEFAULT=false)
% OUTPUT:
%         header (structure) - structure containing parsed DICOM header
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Calls the MEX function that uses DCMTK library to read the DICOM header.
%               To change which parameters are read and their names - the MEX file needs to be edited.
%               This function also corrects formating of certain parameters.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      header = xASL_io_DcmtkRead('/tmp/file.dcm', [], 1);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:
% __________________________________
% Copyright 2015-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

if nargin < 2 || isempty(bPixel)
	bPixel = false;
end
if nargin < 3 || isempty(bTryDCMTK)
	bTryDCMTK = true;
end
if nargin < 4 || isempty(bSkipNonDicoms)
    bSkipNonDicoms = false;
end
if bTryDCMTK
	try
		% Read using DCMTK
		header = xASL_mex_DcmtkRead(filepath, double(bPixel));
	catch ME
        fprintf('%s%s\n','Warning, reading a DICOM with DCMTK failed. Will continue with SPM. Warning message: ', ME.message);
		% Read using SPM routines
		header = spm_dicom_headers(filepath);
        % Verify that this is a DICOM file
        if xASL_io_DcmtkRead_bSkipNonDicoms(bSkipNonDicoms, header, filepath), return; end
		header = xASL_io_DcmtkRead_TrimSPM(header{1});
	end
else
	% Read using SPM routines
	header = spm_dicom_headers(filepath);
    % Verify that this is a DICOM file
    if xASL_io_DcmtkRead_bSkipNonDicoms(bSkipNonDicoms, header, filepath), return; end
	header = xASL_io_DcmtkRead_TrimSPM(header{1});
end
% Correct the acquisition time formatting from string to seconds
if ~isempty(header.AcquisitionTime) && ~isnumeric(header.AcquisitionTime)
	% The old import uses probably str2double, an looses some precission by that - but that won't be more than seconds...
	header.AcquisitionTime = str2num(header.AcquisitionTime);
	%acqTim = str2double(header.AcquisitionTime);
	%atH   = round(acqTim,-4);
	%acqTim = acqTim - atH;
	%atMin = round(acqTim,-2);
	%acqTim = acqTim - atMin;
	%atSec = round(acqTim,0);
	%atMs  = acqTim - atSec;
	%header.AcquisitionTime = atMs + atSec + atMin/100*60 + atH/10000*3600;
end
% Correct the seies time formatting from string to seconds
if ~isempty(header.SeriesTime) && ~isnumeric(header.SeriesTime)
	header.SeriesTime = str2num(header.SeriesTime);
end
% Convert the MRScaleSlope - i.e. the Private_2005_100e tag
listConvert = {'MRScaleSlope' 'PhilipsNumberTemporalScans' 'PhilipsLabelControl' 'PhoenixProtocol'};
typeConvert = {'float' 'decimal' 'char' 'char'};
endianConvert = [0 1 1 1];
% If the field exists and is still in a string format
for iField = 1:length(listConvert)
	if isfield(header,listConvert{iField}) && ~isempty(header.(listConvert{iField})) && ischar(header.(listConvert{iField}))
		% FL
		% First try normal conversion to a double precision
		num = str2double(header.(listConvert{iField}));
		if isnan(num)
			% Conversion failed - it is not a string, but it is probably given in hex
			num = xASL_adm_Hex2Num(header.(listConvert{iField}),typeConvert{iField},endianConvert(iField));
			if ~isnan(num)
				header.(listConvert{iField}) = num;
			end
		else
			% Conversion went fine, convert to single
			header.(listConvert{iField}) = single(num);
		end
	end
end
%% Philips fields
% Convert the RescaleSlopeOriginal - i.e. the Private_2005_140a tag
% If the field exists and is still in a string format
if isfield(header,'RescaleSlopeOriginal') && ~isempty(header.RescaleSlopeOriginal) && ischar(header.RescaleSlopeOriginal)
	% DS
	% First try normal conversion to a double precision
	num = str2double(header.RescaleSlopeOriginal);
	if isnan(num)
		% Conversion failed - it is not a string, but it is probably given in hex
		num = xASL_adm_Hex2Num(header.RescaleSlopeOriginal,'decimal',1);
		if ~isnan(num)
			header.RescaleSlopeOriginal = num;
		end
	else
		header.RescaleSlopeOriginal = num;
	end
end
if ~isempty(header.EffectiveEchoSpacing) 
	% Unknown format
	header.EffectiveEchoSpacing = str2double(header.EffectiveEchoSpacing);
end
if ~isempty(header.MRSeriesWaterFatShift) 
	% FL
	% First try normal conversion to a double precision
	num = str2double(header.MRSeriesWaterFatShift);
	if isnan(num)
		% Conversion failed - it is not a string, but it is probably given in hex
		num = xASL_adm_Hex2Num(header.MRSeriesWaterFatShift,'float',0);
		if ~isnan(num)
			header.MRSeriesWaterFatShift = num;
		end
	else
		% Conversion went fine, convert to single
		header.MRSeriesWaterFatShift = single(num);
	end
end
if ~isempty(header.MRSeriesEPIFactor) 
	% SL
	num = str2double(header.MRSeriesEPIFactor);
	if isnan(num)
		% Conversion failed - it is not a string, but it is probably given in hex
		num = xASL_adm_Hex2Num(header.MRSeriesEPIFactor,'sint',0);
		if ~isnan(num)
			header.MRSeriesEPIFactor = num;
		end
	else
		% Conversion went fine, convert to single
		header.MRSeriesEPIFactor = single(num);
	end
end
%% GE fields
if ~isempty(header.AssetRFactor) 
	% Unknown format
	header.AssetRFactor = str2double(header.AssetRFactor);
end
%% Siemens fields
if ~isempty(header.BandwidthPerPixelPhaseEncode) 
	% Unknown format
	header.BandwidthPerPixelPhaseEncode = str2double(header.BandwidthPerPixelPhaseEncode);
end
	
% if isfield(header,'Manufacturer') && ~isempty(strfind(lower(header.Manufacturer),'philips'))
% 	if isfield(header,'Private_2005_100e')
% 		% The old import script uses single, but looses some precision by that.
% 		% header.MRScaleSlope = single(str2double(header.Private_2005_100e));
% 		if ischar(header.Private_2005_100e) || isstring(header.Private_2005_100e) % note that isstring crashes before 2016b
% 			header.MRScaleSlope = str2double(header.Private_2005_100e);
% 		else
% 			header.MRScaleSlope = double(header.Private_2005_100e);
% 		end
% 		rmfield(header,'Private_2005_100e');
% 	end
% end
end
%-------------------------------------------------------------------------------------------------------------
% Rename and sort some of the fields
% And remove those fields that are not necessary
function header = xASL_io_DcmtkRead_TrimSPM(header)
% Those were already defined in SPM, so we have to redefine them here
if isfield(header,'PhaseNumber')
	header.PhilipsNumberTemporalScans = header.PhaseNumber;
end
if isfield(header,'EPIFactor')
	header.MRSeriesEPIFactor = header.EPIFactor;
end
if isfield(header,'WaterFatShift')
	header.MRSeriesWaterFatShift = header.WaterFatShift;
end
if isfield(header,'CSASeriesHeaderInfo')
	header.PhoenixProtocol = char(header.CSASeriesHeaderInfo);
end
% Read the nested fields
if isfield(header,'SharedFunctionalGroupsSequence')
	if iscell(header.SharedFunctionalGroupsSequence)
		timingItem = header.SharedFunctionalGroupsSequence{1}.MRTimingAndRelatedParametersSequence{1};
	else
		timingItem = header.SharedFunctionalGroupsSequence.Item_1.MRTimingAndRelatedParametersSequence.Item_1;
	end
	if isfield(timingItem,'RepetitionTime')
		header.RepetitionTime = timingItem.RepetitionTime;
	end
end
%% Philips fields
if isfield(header, 'PerFrameFunctionalGroupsSequence')
	privateItem = [];
	if iscell(header.PerFrameFunctionalGroupsSequence)
		echoItem   = header.PerFrameFunctionalGroupsSequence{1}.MREchoSequence{1};
		pixelItem  = header.PerFrameFunctionalGroupsSequence{1}.PixelValueTransformationSequence{1};
		if isfield(header.PerFrameFunctionalGroupsSequence{1}, 'Private_2005_140f')
			privateItem = header.PerFrameFunctionalGroupsSequence{1}.Private_2005_140f{1};
		end
	else
		echoItem   = header.PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1;
		pixelItem  = header.PerFrameFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1;
		if isfield(header.PerFrameFunctionalGroupsSequence.Item_1, 'Private_2005_140f')
			privateItem = header.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1;
		end
	end
	if isfield(echoItem, 'EffectiveEchoTime')
		header.EchoTime = echoItem.EffectiveEchoTime;
	end
		
	if isfield(pixelItem, 'RescaleSlope')
		header.RescaleSlope = pixelItem.RescaleSlope;
	end
	if isfield(pixelItem, 'RescaleIntercept')
		header.RescaleIntercept = pixelItem.RescaleIntercept;
	end
		
	if ~isempty(privateItem)
		if isfield(privateItem, 'NumberOfTemporalPositions')
			header.NumberOfTemporalPositions = privateItem.NumberOfTemporalPositions;
		end
		if isfield(privateItem, 'TemporalPositionIdentifier')
			header.TemporalPositionIdentifier = privateItem.TemporalPositionIdentifier;
		end
		if isfield(privateItem, 'PhilipsNumberTemporalScans')
			header.PhilipsNumberTemporalScans = privateItem.PhilipsNumberTemporalScans;
		end
		if isfield(privateItem, 'PhilipsLabelControl')
			header.PhilipsLabelControl = privateItem.PhilipsLabelControl;
		end
		if isfield(privateItem, 'RescaleSlopeOriginal')
			header.RescaleSlopeOriginal = privateItem.RescaleSlopeOriginal;
		end
		if isfield(privateItem, 'MRScaleSlope')
			header.MRScaleSlope = privateItem.MRScaleSlope;
		end
		if isfield(privateItem, 'NumberOfAverages')
			header.NumberOfAverages = privateItem.NumberOfAverages;
		end
	end
end
if isfield(header,'RealWorldValueMappingSequence')
	rwItem = header.RealWorldValueMappingSequence;
	if isfield(rwItem,'RealWorldValueIntercept')
		header.RWVIntercept = rwItem.RealWorldValueIntercept;
	end
	if isfield(rwItem,'RealWorldValueSlope')
		header.RWVSlope = rwItem.RealWorldValueSlope;
	end
end
if isfield(header, 'Private_0043_192c')
	header.EffectiveEchoSpacing = header.Private_0043_192c;
end
if isfield(header, 'Private_2001_1022')
	header.MRSeriesWaterFatShift = header.Private_2001_1022;
end
if isfield(header, 'Private_2001_1013')
	header.MRSeriesEPIFactor = header.Private_2001_1013;
end
if isfield(header, 'Private_0019_1028')
	header.BandwidthPerPixelPhaseEncode = header.Private_0019_1028;
end
%% GE fields
if isfield(header, 'Private_0043_1083')
	header.AssetRFactor = header.Private_0043_1083;
end
if isfield(header, 'Private_0019_109c')
	header.GELabelingType = header.Private_0019_109c;
end
if isfield(header, 'Private_0043_10a5')
	header.GELabelingDuration = header.Private_0043_10a5;
end
%% Siemens fields
if isfield(header, 'Private_0029_1020')
	header.PhoenixProtocol = header.Private_0029_1020;
end
if isfield(header, 'Private_0019_1029')
	header.SiemensSliceTime = header.Private_0019_1029;
end
%% General field cleanup
% Keep only those fields below...
listFieldsKeep = {'RepetitionTime', 'EchoTime', 'RescaleSlope', 'RescaleIntercept',...
        'NumberOfTemporalPositions', 'NumberOfAverages', 'AcquisitionTime',...
        'PixelData', 'MediaStorageSOPClassUID', 'Manufacturer', 'ManufacturersModelName', 'MRScaleSlope',...
		'StudyDate', 'StudyInstanceUID', 'SeriesInstanceUID', 'ImageType',...
	    'SeriesDescription', 'ProtocolName', 'SeriesTime', 'AcquisitionDate', 'SeriesDate',...
		'AssetRFactor', 'EffectiveEchoSpacing', 'AcquisitionMatrix', 'MRSeriesWaterFatShift',...
	    'MRSeriesEPIFactor', 'BandwidthPerPixelPhaseEncode', 'InPlanePhaseEncodingDirection',...
	    'Rows', 'Columns', 'RescaleSlopeOriginal', 'RWVIntercept', 'RWVSlope',...
	    'AcquisitionContrast', 'ComplexImageComponent', 'GELabelingType', 'PulseSequenceName',...
		'InversionTime', 'GELabelingDuration', 'PhilipsNumberTemporalScans', 'ProtocolName', ...
		'PhilipsLabelControl', 'TemporalPositionIdentifier', 'PhoenixProtocol', 'SoftwareVersions',...
		'SiemensSliceTime', 'StudyID', 'SeriesNumber', 'AcquisitionNumber', 'InstanceNumber' };
headerOld = header;
header = struct;
for iField=1:length(listFieldsKeep)
	if isfield(headerOld,listFieldsKeep{iField})
		header.(listFieldsKeep{iField}) = headerOld.(listFieldsKeep{iField});
	else
		header.(listFieldsKeep{iField}) = [];
	end
end
end
function [isNotADicom] = xASL_io_DcmtkRead_bSkipNonDicoms(bSkipNonDicoms, header, filepath)
%xASL_io_DcmtkRead_bSkipNonDicoms Manage bSkipNonDicoms
%
%         bSkipNonDicoms    - catch any errors reading a file and instead of crashing the code, skip the file
%                             this is useful e.g. when parsing multiple folders that are partly filled with DICOMs
%                             e.g., with xASL_adm_SortDicomToFolders
%                             as DICOM files do not necessarily have a unique extension
%                             (OPTIONAL, DEFAULT=false)
isNotADicom = false;
if isempty(header)
    if bSkipNonDicoms
        fprintf('%s\n', ['Non-dicom: ' filepath]);
        fprintf('Skipping\n');
        isNotADicom = true;
    else
        error(['Non-dicom: ' filepath]);
    end
end
end