function jsonOut = xASL_bids_BIDSifyASLJSON(jsonIn, studyPar, headerASL)
%xASL_bids_BIDSifyASLJSON Goes through the JSON structure of an ASL file and makes sure that all the necessary conversions and checks to BIDS format are applied
%
% FORMAT: jsonOut = xASL_bids_BIDSifyASLJSON(jsonIn, studyPar, headerASL)
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
% 0. Admin
% 1. Obtain the dimensions of the ASL data
% 2. Take all the manually predefined fields from studyPar
% 3. Extract the scaling factors from the JSON header
% 4. Convert certain DICOM fields
% 5. Prioritize DICOM fields over the manually provided studyPar fields
% 6. Field check and name conversion
% 7. Check for time encoded sequence
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
% Copyright (c) 2015-2022 ExploreASL

%% 0. Admin
if nargin < 1 || isempty(jsonIn)
	error('Missing input parameter jsonIn');
end

if nargin < 2 || isempty(studyPar)
	error('Missing input parameter studyPar');
end

if nargin < 3 || isempty(headerASL)
	error('Missing input parameter headerASL');
end

%% 1. Obtain the dimensions of the ASL data
dimASL = headerASL.dat.dim;
if length(dimASL) < 4
	dimASL(4) = 1;
end
	
%% 2. Take all the manually predefined fields from studyPar
jsonOut = studyPar;

% Check if required fields exist in studyPar but not in jsonIn
jsonInMerged = xASL_bids_MergeStudyPar(jsonIn, studyPar, 'asl');
	
%% 3. Extract the scaling factors from the JSON header
if ~isempty(regexpi(jsonInMerged.Manufacturer, 'Philips'))
	jsonOut.scaleFactor = xASL_adm_GetPhilipsScaling(jsonInMerged, headerASL);
else
	jsonOut.scaleFactor = 1;
end

%% 4. Convert certain DICOM fields
% TotalAcquiredPairs - use only the maximum values from DICOM
if isfield(jsonInMerged,'TotalAcquiredPairs') && ~isempty(jsonInMerged.TotalAcquiredPairs) && length(jsonInMerged.TotalAcquiredPairs)>1
	jsonInMerged.TotalAcquiredPairs = max(jsonInMerged.TotalAcquiredPairs);
end

% In case of LookLocker and manually defined multiple flip angle, use this as a priority
if isfield(studyPar,'LookLocker') && ~isempty(studyPar.LookLocker) && studyPar.LookLocker
	if isfield(studyPar,'FlipAngle') && ~isempty(studyPar.FlipAngle) && length(studyPar.FlipAngle)>1
		if isfield(jsonInMerged,'FlipAngle') && ~isequal(jsonInMerged.FlipAngle, jsonOut.FlipAngle)
			jsonInMerged.FlipAngle = jsonOut.FlipAngle;
		end
	end
end

%% 5. Prioritize DICOM fields over the manually provided studyPar fields
% Overwrite differing fields with those from Dicom, but report all differences
strDifferentFields = '';
contentDICOM = '';
contentJSON = '';
nextN = 1;

for fn = fieldnames(jsonInMerged)'
	if isfield(jsonOut,fn{1})
		% If the field is there, then report different fields
		if ~isequal(jsonOut.(fn{1}),jsonInMerged.(fn{1}))
			% Just if this is not only a different vector orientation
			if ~isnumeric(jsonOut.(fn{1})) ||...
					size(jsonOut.(fn{1}),1)>1 && size(jsonOut.(fn{1}),1) >1 ||...
					~isequal((jsonOut.(fn{1}))',jsonInMerged.(fn{1}))
				strDifferentFields{nextN} = fn{1};
                contentDICOM{nextN} = jsonInMerged.(fn{1});
                contentJSON{nextN} = jsonOut.(fn{1});
            end
            nextN = nextN+1;
			if strcmp(fn{1},'TotalAcquiredPairs')
				jsonInMerged.(fn{1}) = jsonOut.(fn{1});
			end
		end
	end
	% Prioritize the DICOM values in general case
	jsonOut.(fn{1}) = jsonInMerged.(fn{1});
end
% Report if certain fields were different as a warning
if ~isempty(strDifferentFields)
    for iField=1:numel(strDifferentFields)
        warning([strDifferentFields{iField} ' differed between DICOM (' xASL_num2str(contentDICOM{iField})...
            ') & studyPar (' xASL_num2str(contentJSON{iField}) ')']);
        
        if strcmp(strDifferentFields{iField}, 'TotalAcquiredPairs')
            fprintf('%s\n', 'Using the studyPar value');
        else
            fprintf('%s\n', 'Using the DICOM value');
        end
        % in the future mention here what we chose, which should be the
        % user-specified value
    end
end

%% 6. Field check and name conversion
% Check repetition time
if isfield(studyPar,'RepetitionTimePreparation')
	% RTP from studyPar has the highest priority
	jsonOut.RepetitionTimePreparation = studyPar.RepetitionTimePreparation;
elseif isfield(studyPar,'RepetitionTime')
	% RT from studyPar comes next
	jsonOut.RepetitionTimePreparation = studyPar.RepetitionTime;
elseif isfield(jsonInMerged,'RepetitionTimePreparation')
	jsonOut.RepetitionTimePreparation = jsonInMerged.RepetitionTimePreparation;
elseif isfield(jsonInMerged,'RepetitionTime')
	% RT from the DICOM has the lowest priority
	jsonOut.RepetitionTimePreparation = jsonInMerged.RepetitionTime;
	% Check if TR is a vector - replace by the maximum then
	if length(jsonOut.RepetitionTimePreparation) > 1
		jsonOut.RepetitionTimePreparation = max(jsonOut.RepetitionTimePreparation);
		warning('TR was a vector. Taking the maximum only.');
	end
end

% In Philips and dcm2niix20220720, often RepetitionTimePreparation is assigned the time right after the last readout and 
% RepetitionTimeExcitation the length of the whole cycle - that needs to be fixed.
if strcmpi(jsonOut.Manufacturer,'Philips')
	if isfield(jsonOut,'RepetitionTimeExcitation') && jsonOut.RepetitionTimeExcitation >= jsonOut.RepetitionTimePreparation
		jsonOut.RepetitionTimePreparation = jsonOut.RepetitionTimeExcitation;
		jsonOut = rmfield(jsonOut,'RepetitionTimeExcitation');
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
		warning(['Both LabelingType (' jsonOut.LabelingType ') and ArterialSpinLabelingType (' jsonOut.ArterialSpinLabelingType ') are defined. Using ' jsonOut.LabelingType]);
	else
		jsonOut.ArterialSpinLabelingType = jsonOut.LabelingType;
	end
end

% The Labeling defined in a private GE field has a priority
if isfield(jsonOut,'GELabelingDuration') && ~isempty(jsonOut.GELabelingDuration)
	% Verify if this doesn't differ from the predefined file, but the DICOM field has priority
	if isfield(jsonOut,'LabelingDuration') && ~isequal(jsonOut.GELabelingDuration,jsonOut.LabelingDuration)
		% if the DICOM information is reasonable - less LDs than volumes, then we report a warning
		if dimASL(4)>=numel(jsonOut.GELabelingDuration)
			warning(['StudyPar Labeling duration (' xASL_num2str(jsonOut.LabelingDuration) ') and GE DICOM private field (' xASL_num2str(jsonOut.GELabelingDuration) ') differ. Using ' xASL_num2str(jsonOut.GELabelingDuration)]);
			jsonOut.LabelingDuration = jsonOut.GELabelingDuration;
		elseif dimASL(4)>=numel(unique(jsonOut.GELabelingDuration))
			tempLabelingDuration = unique(jsonOut.GELabelingDuration);
			if tempLabelingDuration(1) == 0
				tempLabelingDuration(1:end-1) = tempLabelingDuration(2:end);
				tempLabelingDuration(end) = 0;
			end
			jsonOut.LabelingDuration = tempLabelingDuration;
			warning(['StudyPar Labeling duration (' xASL_num2str(jsonOut.LabelingDuration) ') and GE DICOM private field (' xASL_num2str(jsonOut.GELabelingDuration) ') differ. Using ' xASL_num2str(jsonOut.LabelingDuration)]);
		else
			% Otherwise, the information from DICOM appears to be wrong (as is often the case for eASL multi-PLD)
			% and we thus use the provided information. We thus keep the LabelingDuration field untouched.
			warning(['StudyPar Labeling duration (' xASL_num2str(jsonOut.LabelingDuration) ') and GE DICOM private field (' xASL_num2str(jsonOut.GELabelingDuration) ') differ. Using ' xASL_num2str(jsonOut.LabelingDuration)]);
		end
	else
		% All is good and we use the DICOM field
		jsonOut.LabelingDuration = jsonOut.GELabelingDuration;
	end
	
	% GELabelingDuration comes together with the PostLabelingDelay defined in the standard DICOM field called InversionTime
	if isfield(jsonOut,'InversionTime') && ~isempty(jsonOut.InversionTime)
		% Verify if this doesn't differ from the predefined file, but the DICOM field has priority
		if isfield(jsonOut,'PostLabelingDelay') && ~isequal(jsonOut.PostLabelingDelay,jsonOut.InversionTime)
			% if the DICOM information is reasonable - less PLDs than volumes, then we report a warning
			if dimASL(4)>=numel(jsonOut.InversionTime)
				%warning(['StudyPar PostLabelingDelay (' xASL_num2str(jsonOut.PostLabelingDelay) ') and GE DICOM Inversion time (' xASL_num2str(jsonOut.InversionTime) ') differ. Using ' xASL_num2str(jsonOut.InversionTime)]);
				%jsonOut.PostLabelingDelay = jsonOut.InversionTime;
				warning(['StudyPar PostLabelingDelay (' xASL_num2str(jsonOut.PostLabelingDelay) ') and GE DICOM Inversion time (' xASL_num2str(jsonOut.InversionTime) ') differ. Using ' xASL_num2str(jsonOut.PostLabelingDelay)]);
			else
				% Otherwise, the information from DICOM appears to be wrong (as is often the case for eASL multi-PLD)
				% and we thus use the provided information. We thus keep the PostLabelingDelay field untouched.
			end
		else
			jsonOut.PostLabelingDelay = jsonOut.InversionTime;
		end
	end
end

% For GE and multi-PLD or single-PLD not defined in the GELabelingDurationField, we prefer LabelingDuration from study par due to issues with eASL
if strcmp(jsonInMerged.Manufacturer, 'GE') && isfield(studyPar,'LabelingDuration') && ~isequal(studyPar.LabelingDuration, jsonOut.LabelingDuration) &&...
	( length(jsonOut.LabelingDuration)>1 || ~isfield(jsonInMerged,'GELabelingDuration'))
	warning(['StudyPar Labeling duration (' xASL_num2str(studyPar.LabelingDuration) ') and DICOM LabelingDuration field (' xASL_num2str(jsonOut.LabelingDuration) ') differ. Using ' xASL_num2str(studyPar.LabelingDuration)]);
	jsonOut.LabelingDuration = studyPar.LabelingDuration;
end
	
% Free info about the sequence, now just the scanner type+software
if isfield(jsonInMerged,'ManufacturersModelName')
	jsonOut.PulseSequenceDetails = jsonInMerged.ManufacturersModelName;
end
if isfield(jsonInMerged,'SoftwareVersions')
	if ~isfield(jsonOut,'PulseSequenceDetails')
		jsonOut.PulseSequenceDetails = '';
	end
	if ~isempty(jsonOut.PulseSequenceDetails)
		jsonOut.PulseSequenceDetails = [jsonOut.PulseSequenceDetails '-'];
	end
	jsonOut.PulseSequenceDetails = [jsonOut.PulseSequenceDetails jsonInMerged.SoftwareVersions];
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
    
%% 7. Check for time encoded sequence
[jsonOut, bTimeEncoded, bTimeEncodedFME] = xASL_bids_BIDSifyCheckTimeEncoded(jsonInMerged, jsonOut, dimASL(4));

	
%% 8. Merge data from the Phoenix protocol
% Check if the Phoenix protocol is present and parse it
if isfield(jsonInMerged,'PhoenixProtocol') && ~isempty(regexpi(jsonInMerged.Manufacturer,'Siemens'))
	PhoenixParsed = xASL_bids_PhoenixProtocolReader(jsonInMerged.PhoenixProtocol);
	jsonInMerged.PhoenixAnalyzed = xASL_bids_PhoenixProtocolAnalyzer(PhoenixParsed);
end
    
% Overwrite differing fields with those from the Phoenix protocol, but report all differences
if isfield(jsonInMerged,'PhoenixAnalyzed') && ~isempty(jsonInMerged.PhoenixAnalyzed)
	strDifferentFields = '';
	for fn = fieldnames(jsonInMerged.PhoenixAnalyzed)'
		if isfield(jsonOut,fn{1})
			% If the field is there, then report different fields
			if ~isequal(jsonOut.(fn{1}),jsonInMerged.PhoenixAnalyzed.(fn{1}))
				% Just if this is not only a different vector orientation
				if ~isnumeric(jsonOut.(fn{1})) ||...
						(size(jsonOut.(fn{1}),1)>1 && size(jsonOut.(fn{1}),1) >1) ||...
						~isequal((jsonOut.(fn{1}))',jsonInMerged.PhoenixAnalyzed.(fn{1}))
					strDifferentFields = [strDifferentFields ' ' fn{1}];
				end
			end
		end
		% Prioritize the Phoenix field values in the general case
		if ~(strcmp(fn{1},'EchoTime') && bTimeEncodedFME)
			jsonOut.(fn{1}) = jsonInMerged.PhoenixAnalyzed.(fn{1});
		end
	end
	% Report if certain fields were different as a warning
	if ~isempty(strDifferentFields)
		fprintf('Warning: The following user-defined/DICOM fields and DICOM-Phoenix fields differ: %s\n',strDifferentFields);
	end
end

%% 9. Background suppression check
% BSup sanity check
if ~isfield(jsonOut,'BackgroundSuppression')
	warning('BackgroundSuppression field should be defined in BIDS, trying to fix this');
	if isfield(jsonOut,'BackgroundSuppressionNumberPulses')
		if jsonOut.BackgroundSuppressionNumberPulses > 0
			fprintf('%s\n', 'BackgroundSuppressionNumberPulses > 0 field detected, setting BackgroundSuppression to true');
			jsonOut.BackgroundSuppression = true;
		else
			fprintf('%s\n', 'BackgroundSuppressionNumberPulses == 0 field detected, setting BackgroundSuppression to false');
			jsonOut.BackgroundSuppression = false;
		end
    elseif strcmp(jsonOut.MRAcquisitionType, '3D') || ~isempty(strfind(jsonOut.PulseSequenceType, '3D'))
		fprintf('%s\n', '3D acquisition detected, setting BackgroundSuppression to true');
		jsonOut.BackgroundSuppression = true;        
    else
		fprintf('%s\n', 'Setting BackgroundSuppression to false');
		jsonOut.BackgroundSuppression = false;
	end
elseif strcmp(jsonOut.MRAcquisitionType, '3D') || ~isempty(strfind(jsonOut.PulseSequenceType, '3D'))
    if ~jsonOut.BackgroundSuppression
        warning('BackgroundSuppression set to off, which is unlikely for a 3D acquisition');
    end
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

% Fill in extra parameters based on the JSON from the data
if jsonOut.MRAcquisitionType(1) == '2'
    % Parse slicetiming for 2D acquisitions
    
	% Take the studyPar as a prior source of  or SliceReadoutTime since this is mostly wrong in DICOM otherwise
	if isfield(studyPar,'SliceTiming')
		jsonOut.SliceTiming = studyPar.SliceTiming;
	end
	
	if isfield (studyPar, 'SliceReadoutTime')
		if isfield(studyPar, 'SliceTiming') && ~isequal(studyPar.SliceTiming,studyPar.SliceReadoutTime)
			sprintf('Warning, difference in SliceTiming and SliceReadoutTime in studyPar, taking SliceReadoutTime');
		end
		jsonOut.SliceTiming = studyPar.SliceReadoutTime;
	end
		
	% The Siemens field is rather reliable though
	if isfield(jsonOut,'SiemensSliceTime') && ~isempty(jsonOut.SiemensSliceTime)
		jsonOut.SliceTiming = jsonOut.SiemensSliceTime;
	end
	
	% If the length of SliceTiming fits to the number of slices, do nothing
	if isfield(jsonOut,'SliceTiming') && (length(jsonOut.SliceTiming) ~= dimASL(3))
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
if ~isfield(jsonOut, 'ASLContext')
    error('Missing JSON field: ASLContext');
    % Providing info here, if ASLContext doesn't exist then the script
    % would crash below anyway. PM: would it be good to verify the
    % existence of a few obligatory fields?
end

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
		fprintf('%s\n',['ASLContext has a length of ' xASL_num2str(lengthASLContext)]);
        fprintf('%s\n',['NIfTI image has ' xASL_num2str(dimASL(4)) ' volumes']);
        error('These ASLContext and the 4th dimension of the NIfTI do not match');
	else
		NumberEchoTimes = numel(unique(jsonOut.EchoTime));
		
        if NumberEchoTimes > 1 && mod(dimASL(4), lengthASLContext * NumberEchoTimes) == 0
			% For multi-TE acquisition, we check if we can additionally first repeat the echo times
			% Save the original ASLContext and clean the output
			ASLContextTemp = jsonOut.ASLContext;
			jsonOut.ASLContext = [];
			
			% Find indices of all the line-ends
			indexLineEnd = strfind(ASLContextTemp, sprintf('\n'));
			tmpStr = sprintf('\n');
			
			% Find indices of the starting and ending of each vector
			indexStart = [1, indexLineEnd+length(tmpStr)];
			indexEnd   = [indexLineEnd-1, length(ASLContextTemp)];
			
			% Repeat each string
			for iContext = 1:length(indexStart)
				% Repeat it for each TE
				for iTE = 1:NumberEchoTimes
					if iContext == 1 && iTE == 1
						jsonOut.ASLContext = sprintf('%s', ASLContextTemp(indexStart(iContext): indexEnd(iContext)));
					else
						jsonOut.ASLContext = sprintf('%s\n%s', jsonOut.ASLContext, ASLContextTemp(indexStart(iContext): indexEnd(iContext)));
					end
				end
			end
			
			% The repeat the whole block to fill the repetitions
			numRepeat = dimASL(4)/lengthASLContext/NumberEchoTimes;
            tmpStr = jsonOut.ASLContext;
            for iRepeat = 2:numRepeat
                jsonOut.ASLContext = sprintf('%s\n%s',jsonOut.ASLContext,tmpStr);
            end
        else
            numRepeat = dimASL(4)/lengthASLContext;
            tmpStr = jsonOut.ASLContext;
            for iRepeat = 2:numRepeat
                jsonOut.ASLContext = sprintf('%s\n%s',jsonOut.ASLContext,tmpStr);
            end
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
		if sum(ASLContextM0Index .* (jsonOut.(listFieldsZero{iRepeat})~=0))
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
			warning('Discrepancy in the number of averages from DICOM and TotalAcquiredPairs specified in studyPar.json...');
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
	if sum(ASLContextDeltaMIndex) > 0 % In case deltaMs are present
		
		% If only deltaM and no C/L pairs are present
		if sum(ASLContextControlIndex) == 0 && sum(ASLContextLabelIndex) == 0
			nPLD = length(unique(jsonOut.PostLabelingDelay(jsonOut.PostLabelingDelay>0)));
			
			if (nPLD == 3 || nPLD == 7) && mod(sum(ASLContextDeltaMIndex), nPLD + 1) == 0 && regexpi(jsonOut.Manufacturer,'GE')
				% Multi-PLD GE, with Hadamard encoding
				% We need to divide the TotalAcquiredPairs by the number of PLDs+1 (the number of Hadamard phases)
				jsonOut.TotalAcquiredPairs = sum(ASLContextDeltaMIndex) / (nPLD+1);
			elseif mod(sum(ASLContextDeltaMIndex), nPLD) == 0
				jsonOut.TotalAcquiredPairs = sum(ASLContextDeltaMIndex) / nPLD;
			else
				jsonOut.TotalAcquiredPairs = sum(ASLContextDeltaMIndex);
			end
		else
			warning('Cannot calculate TotalAcquiredPairs when both controls and deltaMs are present...');
		end
	elseif sum(ASLContextControlIndex) > 0 % In case C/L are present, than based this on the number of C/L pairs
		if sum(ASLContextControlIndex) == sum(ASLContextLabelIndex)
			jsonOut.TotalAcquiredPairs = sum(ASLContextControlIndex);
        elseif bTimeEncoded == 0
			warning('Cannot calculte TotalAcquiredPairs when control and label numbers differ...');
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


