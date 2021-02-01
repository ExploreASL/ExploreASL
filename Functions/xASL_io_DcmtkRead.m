function header = xASL_io_DcmtkRead(filepath, bPixel)
% ------------------------------------------------------------------------------------------------------------------------------------------------------
% SHORT DESCRIPTION: Reads DICOM headers using DCMTK
%
% FORMAT: header = xASL_io_DcmtkRead(filepath, bPixel)
%
% INPUT:
%         filepath (string) - full path to the DICOM file
%         bPixel (bool) - read pixel data, default 0
% OUTPUT:
%         header (structure) - structure containing parsed DICOM header
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Calls the MEX function that uses DCMTK library to read the DICOM header.
%               To change which parameters are read and their names - the MEX file needs to be edited.
%               This function also corrects formating of certain parameters.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:
% __________________________________
% Copyright 2015-2020 ExploreASL
%
% 2018-12-17, Jan Petr
%
% $Id: ExploreASL_Import 2018-12-17 xASL ver 1.0.0 $

if nargin == 1
	bPixel = 0;
end

if ismac()
	% Read using SPM routines
	header = spm_dicom_headers(filepath);
	header = xASL_io_DcmtkRead_TrimSPM(header{1});
else
	% Read using DCMTK
	header = xASL_mex_DcmtkRead(filepath, bPixel);
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

if ~isempty(header.AssetRFactor) 
	% Unknown format
	header.AssetRFactor = str2double(header.AssetRFactor);
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
	timingItem = header.SharedFunctionalGroupsSequence.Item_1.MRTimingAndRelatedParametersSequence.Item_1;
	if isfield(timingItem,'RepetitionTime')
		header.RepetitionTime = timingItem.RepetitionTime;
	end
end

if isfield(header,'PerFrameFunctionalGroupsSequence')
	echoItem = header.PerFrameFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1;
	if isfield(echoItem,'EffectiveEchoTime')
		header.EchoTime = echoItem.EffectiveEchoTime;
	end
	
	pixelItem = header.PerFrameFunctionalGroupsSequence.Item_1.PixelValueTransformationSequence.Item_1;
	if isfield(pixelItem,'RescaleSlope')
		header.RescaleSlope = pixelItem.RescaleSlope;
	end
	if isfield(pixelItem,'RescaleIntercept')
		header.RescaleIntercept = pixelItem.RescaleIntercept;
	end
	
	privatItem = header.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1;
	
	if isfield(privatItem,'NumberOfTemporalPositions')
		header.NumberOfTemporalPositions = privatItem.NumberOfTemporalPositions;
	end
	if isfield(privatItem,'TemporalPositionIdentifier')
		header.TemporalPositionIdentifier = privatItem.TemporalPositionIdentifier;
	end
	if isfield(privatItem,'PhilipsNumberTemporalScans')
		header.PhilipsNumberTemporalScans = privatItem.PhilipsNumberTemporalScans;
	end
	if isfield(privatItem,'PhilipsLabelControl')
		header.PhilipsLabelControl = privatItem.PhilipsLabelControl;
	end
	if isfield(privatItem,'RescaleSlopeOriginal')
		header.RescaleSlopeOriginal = privatItem.RescaleSlopeOriginal;
	end
	if isfield(privatItem,'MRScaleSlope')
		header.MRScaleSlope = privatItem.MRScaleSlope;
	end
	if isfield(privatItem,'NumberOfAverages')
		header.NumberOfAverages = privatItem.NumberOfAverages;
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

% Keep only those fields below...
listFieldsKeep = {'RepetitionTime', 'EchoTime', 'RescaleSlope', 'RescaleIntercept',...
        'NumberOfTemporalPositions', 'NumberOfAverages', 'AcquisitionTime',...
        'PixelData', 'MediaStorageSOPClassUID', 'Manufacturer',  'MRScaleSlope',...
		'StudyDate', 'StudyInstanceUID', 'SeriesInstanceUID', 'ImageType',...
	    'SeriesDescription', 'ProtocolName', 'SeriesTime', 'AcquisitionDate', 'SeriesDate',...
		'AssetRFactor', 'EffectiveEchoSpacing', 'AcquisitionMatrix', 'MRSeriesWaterFatShift',...
	    'MRSeriesEPIFactor', 'BandwidthPerPixelPhaseEncode', 'InPlanePhaseEncodingDirection',...
	    'Rows', 'Columns', 'RescaleSlopeOriginal', 'RWVIntercept', 'RWVSlope',...
	    'AcquisitionContrast', 'ComplexImageComponent', 'GELabelingType', 'PulseSequenceName',...
		'InversionTime', 'GELabelingDuration', 'PhilipsNumberTemporalScans',...
		'PhilipsLabelControl', 'TemporalPositionIdentifier', 'PhoenixProtocol', 'SoftwareVersions',...
		'SiemensSliceTime'};

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
