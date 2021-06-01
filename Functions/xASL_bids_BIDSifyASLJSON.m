function jsonOut = xASL_bids_BIDSifyASLJSON(jsonIn, studyPar, headerASL)
%xASL_bids_BIDSifyASLJSON Goes through the JSON structure of an ASL file and makes sure that all the necessary conversions and checks to BIDS format are applied
%
% FORMAT: jsonOut = xASL_bids_BIDSifyASLJSON(jsonIn)
%
% INPUT:
%   jsonIn    - JSON with the input fields - from DICOMs (REQUIRED)
%   studyPar  - Manually defined parameters (REQUIRED)
%   headerASL - Header of a NIFTI file from the ASL (REQUIRED)
%
% OUTPUT: 
%   jsonOut   - ordered and checked JSON structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
% It makes all the conversions to a proper BIDS structure, checks the existence of all BIDS fields, removes superfluous fields, checks all the conditions and orderes
% the structure on the output. It works according to the normal BIDS, or ASL-BIDS definition
%
% 1. Obtain the dimensions of the ASL data
% 2. Take all the manually predefined fields from studyPar
% 3. Extract the scaling factors from the JSON header
% 4. Convert certain DICOM fields
% 5. Prioritize DICOM fields over the manually provided studyPar fields
% 6. Field check and name conversion
% 7. Check for Hadamard sequence
% 8. Merge data from the Phoenix protocol
% 9. Background suppression check
% 10. SliceTiming check
% 11. Check if length of vector fields match the number of volumes
% 12. Reformat ASLcontext field
% 13. Verify TotalAcquiredPairs against ASLContext
% 14. Final field check
%
% EXAMPLE: n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL

%% 1. Obtain the dimensions of the ASL data
dimASL = headerASL.dat.dim;
if length(dimASL) < 4
	dimASL(4) = 1;
end
	
%% 2. Take all the manually predefined fields from studyPar
jsonOut = studyPar;
	
%% 3. Extract the scaling factors from the JSON header
if ~isempty(regexpi(jsonIn.Manufacturer,'Philips'))
	jsonOut.scaleFactor = xASL_adm_GetPhilipsScaling(jsonIn, headerASL);
else
	jsonOut.scaleFactor = 0;
end

%% 4. Convert certain DICOM fields
% TotalAcquiredPairs - use only the maximum values from DICOM
if isfield(jsonIn,'TotalAcquiredPairs') && ~isempty(jsonIn.TotalAcquiredPairs) && length(jsonIn.TotalAcquiredPairs)>1
	jsonIn.TotalAcquiredPairs = max(jsonIn.TotalAcquiredPairs);
end

% In case of LookLocker and manually defined multiple flip angle, use this as a priority
if isfield(studyPar,'LookLocker') && ~isempty(studyPar.LookLocker) && studyPar.LookLocker
	if isfield(studyPar,'FlipAngle') && ~isempty(studyPar.FlipAngle) && length(studyPar.FlipAngle)>1
		if isfield(jsonIn,'FlipAngle') && ~isequal(jsonIn.FlipAngle,jsonOut.FlipAngle)
			jsonIn.FlipAngle = jsonOut.FlipAngle;
		end
	end
end

%% 5. Prioritize DICOM fields over the manually provided studyPar fields
% Overwrite differing fields with those from Dicom, but report all differences
strDifferentFields = '';
for fn = fieldnames(jsonIn)'
	if isfield(jsonOut,fn{1})
		% If the field is there, then report different fields
		if ~isequal(jsonOut.(fn{1}),jsonIn.(fn{1}))
			% Just if this is not only a different vector orientation
			if ~isnumeric(jsonOut.(fn{1})) ||...
					(size(jsonOut.(fn{1}),1)>1 && size(jsonOut.(fn{1}),1) >1) ||...
					~isequal((jsonOut.(fn{1}))',jsonIn.(fn{1}))
				strDifferentFields = [strDifferentFields ' ' fn{1}];
			end
			if strcmp(fn{1},'TotalAcquiredPairs')
				jsonIn.(fn{1}) = jsonOut.(fn{1});
			end
		end
	end
	% Prioritize the DICOM values in general case
	jsonOut.(fn{1}) = jsonIn.(fn{1});
end
% Report if certain fields were different as a warning
if ~isempty(strDifferentFields)
	warning('The following user defined and DICOM fields differ:\n %s \n',strDifferentFields);
end

%% 6. Field check and name conversion
% Check repetition time
if isfield(studyPar,'RepetitionTimePreparation')
	% RTP from studyPar has the highest priority
	jsonOut.RepetitionTimePreparation = studyPar.RepetitionTimePreparation;
elseif isfield(studyPar,'RepetitionTime')
	% RT from studyPar comes next
	jsonOut.RepetitionTimePreparation = studyPar.RepetitionTime;
elseif isfield(jsonIn,'RepetitionTime')
	% RT from the DICOM has the lowest priority
	jsonOut.RepetitionTimePreparation = jsonIn.RepetitionTime;
	% Check if TR is a vector - replace by the maximum then
	if length(jsonOut.RepetitionTimePreparation) > 1
		jsonOut.RepetitionTimePreparation = max(jsonOut.RepetitionTimePreparation);
		warning('TR was a vector. Taking the maximum only.');
	end
end

% Check echo time, for vectors
if isfield(jsonOut,'EchoTime') && length(jsonOut.EchoTime)>1
	% Remove zero entries
	jsonOut.EchoTime = jsonOut.EchoTime(jsonOut.EchoTime ~= 0);
end

% Convert the field name of the old LabelingType
if isfield(jsonOut,'LabelingType')
	if isfield(jsonOut,'ArterialSpinLabelingType') && ~isequal(jsonOut.LabelingType,jsonOut.ArterialSpinLabelingType)
		warning('Both LabelingType and ArterialSpinLabelingType are defined with a different value');
	else
		jsonOut.ArterialSpinLabelingType = jsonOut.LabelingType;
	end
end

% The Labeling defined in a private GE field has a priority
if isfield(jsonOut,'GELabelingDuration') && ~isempty(jsonOut.GELabelingDuration)
	if isfield(jsonOut,'LabelingDuration') && ~isequal(jsonOut.GELabelingDuration,jsonOut.LabelingDuration)
		warning('Labeling duration mismatch with GE private field.');
	end
	jsonOut.LabelingDuration = jsonOut.GELabelingDuration;
	
	% GELabelingDuration comes together with the PostLabelingDelay defined in the standard DICOM field called InversionTime
	if isfield(jsonOut,'InversionTime') && ~isempty(jsonOut.InversionTime)
		% Verify if this doesn't differ from the predefined file, but the DICOM field has priority
		if isfield(jsonOut,'PostLabelingDelay') && ~isequal(jsonOut.PostLabelingDelay,jsonOut.InversionTime)
			warning('PostLabelingDelay mismatch with the GE-DICOM value in Inversion time.');
		end
		jsonOut.PostLabelingDelay = jsonOut.InversionTime;
	end
end

% Free info about the sequence, now just the scanner type+software
if isfield(jsonIn,'ManufacturersModelName')
	jsonOut.PulseSequenceDetails = jsonIn.ManufacturersModelName;
end
if isfield(jsonIn,'SoftwareVersions')
	if ~isfield(jsonOut,'PulseSequenceDetails')
		jsonOut.PulseSequenceDetails = '';
	end
	if ~isempty(jsonOut.PulseSequenceDetails)
		jsonOut.PulseSequenceDetails = [jsonOut.PulseSequenceDetails '-'];
	end
	jsonOut.PulseSequenceDetails = [jsonOut.PulseSequenceDetails jsonIn.SoftwareVersions];
end

% Process all the data and automatically fill in the missing parameters
if ~isfield(jsonOut,'MRAcquisitionType')
	error('MRAcquisitionType has to be defined either in the data or studyPar');
end

if strcmpi(jsonOut.MRAcquisitionType,'2D')
	jsonOut.PulseSequenceType = '2D_EPI';
else
	if strcmpi(jsonOut.Manufacturer,'GE') || strcmpi(jsonOut.Manufacturer,'GE_WIP') || strcmpi(jsonOut.Manufacturer,'GE_product')
		jsonOut.PulseSequenceType = '3D_spiral';
	else
		jsonOut.PulseSequenceType = '3D_GRASE';
	end
end
    
%% 7. Check for Hadamard sequence
% FME sequence check
if isfield(jsonIn,'SeriesDescription')
	isHadamardFME = ~isempty(regexp(char(jsonIn.SeriesDescription),'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once'));
else
	isHadamardFME = false;
end
if isHadamardFME
	if isfield(jsonOut,'EchoTime') && isfield(jsonOut,'PostLabelingDelay')
		% From the import, the length of EchoTime should correspond to the number of volumes
		if length(jsonOut.EchoTime)~=length(jsonOut.PostLabelingDelay)
			% So here, we first make sure that each PLD is repeated for the whole block of echo-times
			numberTEs = length(uniquetol(jsonOut.EchoTime,0.001)); % Obtain the number of echo times
			repeatedPLDs = repmat(jsonOut.PostLabelingDelay(:),1,numberTEs)';
			repeatedPLDs = repeatedPLDs(:);
			
			if length(repeatedPLDs) > length(jsonOut.EchoTime) || mod(length(jsonOut.EchoTime),length(repeatedPLDs)) ~= 0
				warning('Did not succeed in repeating PLDs for each TE for Hadamard sequence import');
			end
			% Make sure that number of volumes can be divided by the repeated PLDs
			jsonOut.PostLabelingDelay = repeatedPLDs;
		end
	end
end

% Check for FME Hadamard sequences
if isfield(jsonOut,'SeriesDescription')
	isHadamardFME = ~isempty(regexp(char(jsonOut.SeriesDescription),'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once'));
	if isHadamardFME
		startDetails = regexp(char(jsonOut.SeriesDescription),'\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once');
		jsonOut.HadamardType = xASL_str2num(jsonOut.SeriesDescription(startDetails:startDetails+1),'auto');
		jsonOut.HadamardNumberPLD = xASL_str2num(jsonOut.SeriesDescription(startDetails+3:startDetails+4),'auto');
		jsonOut.HadamardNumberTE = xASL_str2num(jsonOut.SeriesDescription(startDetails+10:startDetails+11),'auto');
		fprintf('FME sequence, Hadamard-%d encoded images, %d PLDs, %d TEs\n', jsonOut.HadamardType, jsonOut.HadamardNumberPLD, jsonOut.HadamardNumberTE);
	end
else
	isHadamardFME = false;
end
	
%% 8. Merge data from the Phoenix protocol
% Check if the Phoenix protocol is present and parse it
if isfield(jsonIn,'PhoenixProtocol') && ~isempty(regexpi(jsonIn.Manufacturer,'Siemens'))
	PhoenixParsed = xASL_bids_PhoenixProtocolReader(jsonIn.PhoenixProtocol);
	jsonIn.PhoenixAnalyzed = xASL_bids_PhoenixProtocolAnalyzer(PhoenixParsed);
end
    
% Overwrite differing fields with those from the Phoenix protocol, but report all differences
if isfield(jsonIn,'PhoenixAnalyzed') && ~isempty(jsonIn.PhoenixAnalyzed)
	strDifferentFields = '';
	for fn = fieldnames(jsonIn.PhoenixAnalyzed)'
		if isfield(jsonOut,fn{1})
			% If the field is there, then report different fields
			if ~isequal(jsonOut.(fn{1}),jsonIn.PhoenixAnalyzed.(fn{1}))
				% Just if this is not only a different vector orientation
				if ~isnumeric(jsonOut.(fn{1})) ||...
						(size(jsonOut.(fn{1}),1)>1 && size(jsonOut.(fn{1}),1) >1) ||...
						~isequal((jsonOut.(fn{1}))',jsonIn.PhoenixAnalyzed.(fn{1}))
					strDifferentFields = [strDifferentFields ' ' fn{1}];
				end
			end
		end
		% Prioritize the Phoenix field values in the general case
		if ~(strcmp(fn{1},'EchoTime') && isHadamardFME)
			jsonOut.(fn{1}) = jsonIn.PhoenixAnalyzed.(fn{1});
		end
	end
	% Report if certain fields were different as a warning
	if ~isempty(strDifferentFields)
		warning('The following user-defined/DICOM fields and DICOM-Phoenix fields differ:\n %s \n',strDifferentFields);
	end
end

%% 9. Background suppression check
% BSup sanity check
if ~isfield(jsonOut,'BackgroundSuppression')
	warning('BackgroundSuppression field should be define in BIDS, setting the default to false');
	jsonOut.BackgroundSuppression = false;
end

if jsonOut.BackgroundSuppression == false
	% remove pulsenumbers and timings if BSup is OFF
	if isfield(jsonOut,'BackgroundSuppressionNumberPulses')
		jsonOut = rmfield(jsonOut,'BackgroundSuppressionNumberPulses');
	end
	if isfield(jsonOut,'BackgroundSuppressionPulseTime')
		jsonOut = rmfield(jsonOut,'BackgroundSuppressionPulseTime');
	end
else
	% If times are given, but not the number of pulses, then assign the length
	if isfield(jsonOut,'BackgroundSuppressionPulseTime')
		if isfield(jsonOut,'BackgroundSuppressionNumberPulses')
			if jsonOut.BackgroundSuppressionNumberPulses ~= length(jsonOut.BackgroundSuppressionPulseTime)
				fprintf('Warning: Number of pulses and their timings do not match.');
			end
		else
			jsonOut.BackgroundSuppressionNumberPulses = length(jsonOut.BackgroundSuppressionPulseTime);
		end
	end
end

%% 10. SliceTiming check
% SliceReadoutTime from the manual entry is prioritized
if isfield (studyPar,'SliceReadoutTime')
	if isfield(studyPar,'SliceTiming') && ~isequal(studyPar.SliceTiming,studyPar.SliceReadoutTime)
		sprintf('Warning, difference in SliceTiming and SliceReadoutTime');
	end
	studyPar.SliceTiming = studyPar.SliceReadoutTime;
end

% Fill in extra parameters based on the JSON from the data
if jsonOut.MRAcquisitionType(1) == '2'
	% Take the studyPar as a prior source of SliceTiming since this is mostly wrong in DICOM otherwise
	if isfield(studyPar,'SliceTiming')
		jsonOut.SliceTiming = studyPar.SliceTiming;
	end
	
	% The Siemens field is rather reliable though
	if isfield(jsonOut,'SiemensSliceTime') && ~isempty(jsonOut.SiemensSliceTime)
		jsonOut.SliceTiming = jsonOut.SiemensSliceTime;
	end
	
	% If the length of SliceTiming fits to the number of slices, do nothing
	if length(jsonOut.SliceTiming) ~= dimASL(3)
		% if the length of studyPar.sliceTiming is higher than 1 and the difference non-zero then use this
		if length(jsonOut.SliceTiming) > 1 && abs(jsonOut.SliceTiming(2)-jsonOut.SliceTiming(1)) > 0
			jsonOut.SliceTiming = jsonOut.SliceTiming(2)-jsonOut.SliceTiming(1);
		end
		
		if abs(jsonOut.SliceTiming) > 0
			jsonOut.SliceTiming = ((0:(dimASL(3)-1))')*jsonOut.SliceTiming;
		end
	end
else
	% 3D sequences should not have a SliceTiming or have it defined as zero
	if isfield(jsonOut,'SliceTiming')
		jsonOut = rmfield(jsonOut,'SliceTiming');
	end
end

%% 11. Check if length of vector fields match the number of volumes
% If Post-labeling delay or labeling duration, and other fields, is longer than 1, but shorter then number of volumes then repeat it
listFieldsRepeat = {'PostLabelingDelay', 'LabelingDuration','VascularCrushingVENC','FlipAngle','RepetitionTimePreparation'};
for iRepeat = 1:length(listFieldsRepeat)
	if isfield(jsonOut,(listFieldsRepeat{iRepeat})) && (length(jsonOut.(listFieldsRepeat{iRepeat})) > 1) && (dimASL(4) ~= length(jsonOut.(listFieldsRepeat{iRepeat})))
		if mod(dimASL(4),length(jsonOut.(listFieldsRepeat{iRepeat})))
			error('Cannot find a match between the %s and the 4th dimension of the NIFTI.\n',listFieldsRepeat{iRepeat});
		else
			jsonOut.(listFieldsRepeat{iRepeat}) = repmat(jsonOut.(listFieldsRepeat{iRepeat})(:),[dimASL(4)/length(jsonOut.(listFieldsRepeat{iRepeat})) 1]);
		end
	end
end

%% 12. Reformat ASLcontext field
% Remove ',' and ';' at the end
if (jsonOut.ASLContext(end) == ';') || (jsonOut.ASLContext(end) == ',')
	jsonOut.ASLContext = jsonOut.ASLContext(1:(end-1));
end

% Replace all ',' and ';' by \n
jsonOut.ASLContext = strrep(jsonOut.ASLContext, ' ','');
jsonOut.ASLContext = strrep(jsonOut.ASLContext, ';',',');
lengthASLContext = sum(jsonOut.ASLContext == ',')+1;
jsonOut.ASLContext = strrep(jsonOut.ASLContext, ',',sprintf('\n'));

% Check if the length is the same
if dimASL(4) ~= lengthASLContext
	% Check if we can simply repeat it
	if mod(dimASL(4),lengthASLContext)
		error('Cannot find a match between the ASLContext and the 4th dimension of the NIFTI');
	else
		numRepeat = dimASL(4)/lengthASLContext;
		tmpStr = jsonOut.ASLContext;
		for iRepeat = 2:numRepeat
			jsonOut.ASLContext = sprintf('%s\n%s',jsonOut.ASLContext,tmpStr);
		end
	end
end
jsonOut.ASLContext = sprintf('%s\n',jsonOut.ASLContext);

% If Post-labeling delay or labeling duration is longer than 1 then we need to verify this against
% the ASLContext - these fields should be 0 for m0scans
listFieldsZero = {'PostLabelingDelay', 'LabelingDuration'};

% Find m0scans indices in ASLContext
ASLContextCell = strsplit(jsonOut.ASLContext,'\n'); % Split to cells by line-end
ASLContextM0Index = regexp(ASLContextCell,'^m0scan'); % Find m0scans
ASLContextM0Index = cellfun(@(x)~isempty(x),ASLContextM0Index); % Create a vector out of it
ASLContextM0Index = ASLContextM0Index(1:dimASL(4)); % Remove the last empty field

% Go through all variables, check those that have length bigger than 1
for iRepeat = 1:length(listFieldsZero)
	if isfield(jsonOut,listFieldsZero{iRepeat}) && length(jsonOut.(listFieldsZero{iRepeat})) > 1
		% Make sure the vector is a row vector
		jsonOut.(listFieldsZero{iRepeat}) = jsonOut.(listFieldsZero{iRepeat})(:)';
		% Check the all m0scans have zeros
		if ~isequal(ASLContextM0Index,jsonOut.(listFieldsZero{iRepeat})==0)
			% If not, then set to zeros and report a warning
			jsonOut.(listFieldsZero{iRepeat})(ASLContextM0Index) = 0;
			warning(['Had to set non-zero values for m0scan to zero in ' listFieldsZero{iRepeat}]);
		end
	end
end
	
%% 13. Verify TotalAcquiredPairs against ASLContext
% Import the number of averages
if isfield(jsonOut,'NumberOfAverages') && (max(jsonOut.NumberOfAverages) > 1)
	if isfield(studyPar,'TotalAcquiredPairs')
		if max(jsonOut.NumberOfAverages) ~= studyPar.TotalAcquiredPairs
			warning('Discrepancy in the number of averages');
		end
	end
end

% First count the number of controls, labels, and deltaMs
ASLContextControlIndex = cellfun(@(x)~isempty(x),regexpi(ASLContextCell,'^control')); % Create a vector out of it
ASLContextControlIndex = ASLContextControlIndex(1:dimASL(4)); % Remove the last empty field
ASLContextLabelIndex = cellfun(@(x)~isempty(x),regexpi(ASLContextCell,'^label')); % Create a vector out of it
ASLContextLabelIndex = ASLContextLabelIndex(1:dimASL(4)); % Remove the last empty field
ASLContextDeltaMIndex = cellfun(@(x)~isempty(x),regexpi(ASLContextCell,'^deltam')); % Create a vector out of it
ASLContextDeltaMIndex = ASLContextDeltaMIndex(1:dimASL(4)); % Remove the last empty field

% If TotalAcquiredPairs is 1, but more control/label pairs od deltaMs are present, then set this to the correct
% number
if ~isfield(jsonOut,'TotalAcquiredPairs') || jsonOut.TotalAcquiredPairs == 1
	if sum(ASLContextDeltaMIndex) > 0
		if sum(ASLContextControlIndex) == 0 && sum(ASLContextControlIndex) == 0
			jsonOut.TotalAcquiredPairs = sum(ASLContextDeltaMIndex);
		else
			warning('Cannot calculate TotalAcquiredPairs when both controls and deltaMs are present');
		end
	elseif sum(ASLContextControlIndex) > 0
		if sum(ASLContextControlIndex) == sum(ASLContextLabelIndex)
			jsonOut.TotalAcquiredPairs = sum(ASLContextControlIndex);
		else
			warning('Cannot calculte TotalAcquiredPairs when control and label numbers differ');
		end
	end
end

%% 14. Final field check
if isfield(jsonOut,'BolusCutOffFlag') && jsonOut.BolusCutOffFlag
	if isfield(jsonOut,'BolusCutOffTechnique') && strcmpi(jsonOut.BolusCutOffTechnique,'Q2TIPS')
		if ~isfield(jsonOut,'BolusCutOffDelayTime') || length(jsonOut.BolusCutOffDelayTime)~=2
			warning('Q2TIPS BolusCutOff has to have 2 values defined');
		end
	end
end

end
