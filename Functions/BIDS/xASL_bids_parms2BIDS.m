function outParms = xASL_bids_parms2BIDS(inXasl, inBids, bOutBids, bPriorityBids)
% Takes the input parameters from xASL legacy format (inXasl) and BIDS format, merges them and converts to either xASL legacy or BIDS format.
% FORMAT: outBids = xASL_bids_parms2BIDS(inXasl[, inBids, bOutBids, priorityBids])
% 
% INPUT:
%   inXasl       - a structure with input parameters in the legacy xASL format (REQUIRED)
%   inBids       - a structure with input parameters in the BIDS format (OPTIONAL, DEFAULT = [])
%   bOutBids     - the output structures is in BIDS format (==1, default) or xASL format (==0) (OPTIONAL, DEFAULT = 1)
%   bPriorityBids - in case of conflicts, the BIDS input is preferred (==1, default), otherwise (==0), xASL is prefered (OPTIONAL, DEFAULT = 1)
%
% OUTPUT:
% outParms       - the merged output structure in the selected format
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions takes two parameter structures and merges them. At the same time, renames all fields
%              according to the output type (note that only some fields have two standardised names different between the two formats.
%              In case of duplicities, takes the field value from the preferred format. 
%              Also takes into account that the units in BIDS are s, but in xASL ms.
%              This function performs the following steps:
%
%              1. Define field names that need to be convert/renamed/merged
%              2. Convert XASL fields to the output format (BIDS or XASL)
%              3. Convert BIDS fields to the output format (BIDS or XASL)
%              4. Merge the BIDS and XASL fields, convert field values
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: outParms = xASL_bids_parms2BIDS(inXasl, inBids);
%          outParms = xASL_bids_parms2BIDS(inXasl, [], 1, 0);
% __________________________________
% Copyright 2015-2023 ExploreASL


%% Admin
if nargin < 1
	error('Need at least 1 input parameters.');
end
if nargin < 2
	inBids = [];
end
if isempty(inXasl) && isempty(inBids)
	error('At least one of the structures (XASL or BIDS) need to be assigned.');
end
if nargin < 3 || isempty(bOutBids)
	bOutBids = 1;
end
if nargin < 4 || isempty(bPriorityBids)
	bPriorityBids = 1;
end

%% ----------------------------------------------------------------------
%% 1) Define field names that need to be convert/renamed/merged

% Fields with these names need to have the time converted between XASL legacy and BIDS, and define their recommended range in ms
convertTimeFieldsXASL =       {'EchoTime' 'RepetitionTime' 'Initial_PLD' 'LabelingDuration' 'GELabelingDuration' 'InversionTime' 'SliceReadoutTime' 'BloodT1' 'T2' 'TissueT1' 'SiemensSliceTime' 'BackgroundSuppressionPulseTime'};
convertTimeFieldsRange =       [0.5        5                10            10                 10                   10              5                  100       10   100        5                  5;...% Minimum in ms
                                500        20000            10000         5000               5000                 5000            400                5000      500  5000       400                10000];% Maximum in ms   
convertTimeFieldsAllowOutliers=[0          0                1             1                  1                    0               0                  0         0    0          0                  0]; % For multiple values, don't print a warning if mean value is within the range
					  
% Fields that are entered under the subfield 'Q' for xASL on the output
xASLqFields = {'LabelingType' 'Initial_PLD' 'BackGrSupprPulses' 'LookLocker' 'FlipAngle' 'LabelingDuration' 'SliceReadoutTime' 'NumberOfAverages' 'BloodT1'...
	           'BackgroundSuppressionPulseTime' 'BackgroundSuppressionNumberPulses' 'TimeEncodedMatrixSize' 'NumberEchoTimes' 'TimeEncodedMatrixType' 'SoftwareVersions' 'TimeEncodedMatrix'};

% Some JSON fields need to be updated to fit the BIDS definition
% This is a one way process of changing the names from old to new
updateNamesBIDSold = {'PhilipsRescaleSlope' 'PhilipsRWVSlope' 'PhilipsScaleSlope' 'PhilipsRescaleIntercept' 'PhilipsRWVIntercept' 'BackGrSupprPulses'                 'TotalAcquiredVolumes'  'InitialPostLabelDelay'};
updateNamesBIDSnew = {'RescaleSlope'        'RWVSlope'        'MRScaleSlope'      'RescaleIntercept'        'RWVIntercept'        'BackgroundSuppressionNumberPulses' 'TotalAcquiredPairs'    'PostLabelingDelay'};

% These fields have different names in xASL and in BIDS
% They are therefore renamed depending on the type of output
changeNamesXASL = {'Vendor'       'readoutDim'        'Initial_PLD'         'LabelingType'              'RepetitionTime'            'NumberOfAverages'     'SliceReadoutTime' 'Sequence'};
changeNamesBIDS = {'Manufacturer' 'MRAcquisitionType' 'PostLabelingDelay'   'ArterialSpinLabelingType'  'RepetitionTimePreparation' 'TotalAcquiredPairs'   'SliceTiming'      'PulseSequenceType'};

%% ----------------------------------------------------------------------
%% 2) Convert XASL fields to the output format (BIDS or XASL legacy)

% Goes through all XASL fields
if ~isempty(inXasl)
    % Flatten inXasl structure (include x.Q, x.settings, x.dataset, x.modules, x.modules.structural, x.modules.asl, x.modules.population parameters)
	inXasl = xASL_wrp_Quantify_FlattenX(inXasl);
	
	% If output is in BIDS, then XASL fields need to be converted
	if bOutBids
		FieldsA = fields(inXasl);
		for iA = 1:length(FieldsA)
			
			% Convert the units for all time fields from ms to s
			for iT = find(strcmp(FieldsA{iA}, convertTimeFieldsXASL))
				% Convert only numeric fields
				% In certain cases (e.g. SliceReadoutTime='shortestTR'), a string is given which is then skipped and not converted
				if isnumeric(inXasl.(FieldsA{iA}))
					% For non-zero fields, check if they are within the predefined range
					if inXasl.(FieldsA{iA}) ~= 0
						% For SliceTime - compare the difference of the first two fields
						if strcmpi(FieldsA{iA},'SliceTiming') || strcmpi(FieldsA{iA},'SliceReadoutTime') || strcmpi(FieldsA{iA},'SiemensSliceTime')
							valueCheck = inXasl.(FieldsA{iA});
							if length(valueCheck) > 1
								valueCheck = valueCheck(2) - valueCheck(1);
							end
						else
							valueCheck = inXasl.(FieldsA{iA});
						end
						% If outside of the recommended range, then still convert, but issue a warning
						if max(valueCheck < convertTimeFieldsRange(1,iT)) || max(valueCheck > convertTimeFieldsRange(2,iT))
							warning(['Field ' FieldsA{iA} ' in xASL structure has a value ' xASL_num2str(valueCheck)...
								', which is outside of the recommended range <'...
								xASL_num2str(convertTimeFieldsRange(1,iT)) ',' xASL_num2str(convertTimeFieldsRange(2,iT)) '> ms.']);
						end
					end
					inXasl.(FieldsA{iA}) = inXasl.(FieldsA{iA})/1000;
				end
			end
			
			% Rename XASL fields to BIDS field according to the information in changeNamesXASL and changeNamesBIDS
			for iL = find(strcmp(FieldsA{iA},changeNamesXASL))
				inXasl.(changeNamesBIDS{iL}) = inXasl.(changeNamesXASL{iL});
				
				% Remove only if the name was different
				if ~strcmp(changeNamesBIDS{iL},changeNamesXASL{iL})
					inXasl = rmfield(inXasl,changeNamesXASL{iL});
				end
			end
			
			if isfield(inXasl,'LabelingType') && strcmpi(inXasl.LabelingType,'PASL') && isfield(inXasl,'LabelingDuration')
				inXasl.BolusCutOffDelayTime(1) = inXasl.LabelingDuration;
				inXasl = rmfield(inXasl,'LabelingDuration');
				inXasl.BolusCutOffFlag = true;
			end
			
			if isfield(inXasl,'BackgroundSuppressionNumberPulses') && inXasl.BackgroundSuppressionNumberPulses == 0
				inXasl.BackgroundSuppression = false;
				inXasl = rmfield(inXasl,'BackgroundSuppressionNumberPulses');
			end
			
		end
	end
end

%% ----------------------------------------------------------------------
%% 3) Convert BIDS fields to the output format (BIDS or XASL)

% Goes through all BIDS fields
if ~isempty(inBids)
	% Move all fields from Q to the main structure
	if isfield(inBids,'Q')
		FieldsQ = fields(inBids.Q);
		for iQ=1:length(FieldsQ)
			inBids.(FieldsQ{iQ}) = inBids.Q.(FieldsQ{iQ});
		end
		inBids = rmfield(inBids,'Q');
	end
	
	% Update all the old JSON fields to the respective BIDS field name
	FieldsA = fields(inBids);
	for iA = 1:length(FieldsA)
		for iL = find(strcmp(FieldsA{iA}, updateNamesBIDSold))
			inBids.(updateNamesBIDSnew{iL}) = inBids.(updateNamesBIDSold{iL});
			inBids = rmfield(inBids,updateNamesBIDSold{iL});
		end
	end
	
	% When the output is in XASL we need to convert from BIDS to XASL
	if bOutBids ~= 1
		% Convert M0Type and M0Estimate to xASL Legacy
		if isfield(inBids,'M0Type')
			switch(inBids.M0Type)
				case 'Estimate'
					inBids.M0 = inBids.M0Estimate;
					inBids = rmfield(inBids,'M0Estimate');
				case 'Separate'
					inBids.M0 = 'separate_scan';
				case 'Included'
					inBids.M0 = 'separate_scan';
				case 'Absent'
					inBids.M0 = 'UseControlAsM0';
			end
			
			% Remove the old field
			inBids = rmfield(inBids,'M0Type');
		end
		
		% Preconvert certain names upfront - so that the values can be converted s-ms in the step below
		if isfield(inBids,'ArterialSpinLabelingType') && strcmpi(inBids.ArterialSpinLabelingType,'PASL') && isfield(inBids,'BolusCutOffDelayTime')
			inBids.LabelingDuration = inBids.BolusCutOffDelayTime(1);
			inBids = rmfield(inBids,'BolusCutOffDelayTime');
		end
			
		if isfield(inBids,'BackgroundSuppression') && (inBids.BackgroundSuppression == false)
			inBids.BackgroundSuppressionNumberPulses = 0;
			inBids = rmfield(inBids,'BackgroundSuppression');
		end
		
		FieldsA = fields(inBids);
		for iA = 1:length(FieldsA)
			% Rename all listed fields to XASL variant
			FieldNameChanged = FieldsA{iA};
			for iL = find(strcmp(FieldsA{iA}, changeNamesBIDS))
				inBids.(changeNamesXASL{iL}) = inBids.(changeNamesBIDS{iL});
				inBids = rmfield(inBids, changeNamesBIDS{iL});
				% Update the name of the field after the change
				FieldNameChanged = changeNamesXASL{iL};
			end
						
			% Convert the units for all time fields from s to ms
			for iT = find(strcmp(FieldNameChanged, convertTimeFieldsXASL))
				% Convert only if the field is numeric
				% In certain cases (e.g. SliceReadoutTime='shortestTR') a string can be given, which needs to be skipped then
				if isnumeric(inBids.(FieldNameChanged))
					inBids.(FieldNameChanged) = inBids.(FieldNameChanged)*1000;
					
					% Check if the value is within the recommended range after conversion and issue a warning if not
					if max(inBids.(FieldNameChanged) ~= 0) % if any of the numeric values is not zero
						
						% For SliceTiming - compare the difference of the first two fields
						if strcmpi(FieldNameChanged,'SliceTiming') || strcmpi(FieldNameChanged,'SliceReadoutTime') || strcmpi(FieldNameChanged,'SiemensSliceTime')
							valueCheck = inBids.(FieldNameChanged);
							if length(valueCheck) > 1
								valueCheck = valueCheck(2) - valueCheck(1);
							end
						else
							valueCheck = inBids.(FieldNameChanged);
						end
						
						if min(valueCheck) < convertTimeFieldsRange(1,iT) || max(valueCheck) > convertTimeFieldsRange(2,iT)
							if numel(valueCheck)==1
								warningMessage = ['a value of ' xASL_num2str(valueCheck)];
							else
								warningMessage = ['multiple values, from ' xASL_num2str(min(valueCheck)) ' to ' xASL_num2str(max(valueCheck))];
							end
							
                            warningMessage = ['Field ' FieldNameChanged ' in xASL structure has ' warningMessage...
                                ', which is outside of the recommended range <'...
                                xASL_num2str(convertTimeFieldsRange(1,iT)) ',' xASL_num2str(convertTimeFieldsRange(2,iT)) '> ms.'];
							
							if numel(valueCheck) == 1 
								% For a single value, always print the warning
								warning(warningMessage);
							elseif mean(valueCheck) < convertTimeFieldsRange(1, iT) || mean(valueCheck) > convertTimeFieldsRange(2, iT)
								% For multiple values, check if the mean is also outside of the range
								warning([warningMessage ' And the mean value is also ouside of the range.'])
							else
								fprintf([warningMessage ' But the mean is within the range.\n'])
							end
						end
					end
				else
					% But a warning is issued as this is not according to the BIDS specification
					warning(['Field ' FieldNameChanged ' should have been converted from ms to s, but the conversion was skipped as it contains a non-numeric value: ' inBids.(FieldNameChanged) '. '...
						     'This is not according to the BIDS specification.']);
				end
			end
		end
	end
end

%% ----------------------------------------------------------------------
%% 4) Merge the BIDS and XASL fields, convert field values

% Performs merging of XASL and BIDS, prioritizing either BIDS or XASL
% Don't overwrite existing field with a zero, empty, nan or inf value
if bPriorityBids
	outParms = inXasl;
	if ~isempty(inBids)
		FieldsA = fields(inBids);
		for iA = 1:length(FieldsA)
			if ~isfield(outParms, (FieldsA{iA}))
				% If the field is not yet there, then copy it
				outParms.(FieldsA{iA}) = inBids.(FieldsA{iA});
			else
				% If the field is already there, then check if the field to overwrite is not contain NaN or full of zeros or empty strings
                % Do this separately for cell and non-cell arrays
                if iscell(inBids.(FieldsA{iA}))
                    if sum(~cellfun(@isempty, inBids.(FieldsA{iA})))>0
                        outParms.(FieldsA{iA}) = inBids.(FieldsA{iA});
                    end
                elseif isstruct(inBids.(FieldsA{iA}))
                    % Ignore structs within BIDS metadata for now
                else
                    if (sum(isnan(inBids.(FieldsA{iA})))==0) && (sum(inBids.(FieldsA{iA}) ~= 0) > 0) && (sum(~isinf(inBids.(FieldsA{iA})))>0) && (~isempty(inBids.(FieldsA{iA})))
                        outParms.(FieldsA{iA}) = inBids.(FieldsA{iA});
                    end
                end
			end
		end
	end
else
	outParms = inBids;
	FieldsA = fields(inXasl);
	for iA = 1:length(FieldsA)
		if ~isfield(outParms, (FieldsA{iA}))
			% If the field is not yet there, then copy it
			outParms.(FieldsA{iA}) = inXasl.(FieldsA{iA});
		else
			% If the field is already there, then check if the field to overwrite is not contain NaN or full of zeros or empty strings
			% Do this separately for cell and non-cell arrays
			if iscell(inXasl.(FieldsA{iA}))
				if sum(~cellfun(@isempty, inXasl.(FieldsA{iA})))>0
					outParms.(FieldsA{iA}) = inXasl.(FieldsA{iA});
				end
			else
				if (sum(isnan(inXasl.(FieldsA{iA})))==0) && (sum(inXasl.(FieldsA{iA}) ~= 0)>0) && (sum(~isinf(inXasl.(FieldsA{iA})))>0) && (~isempty(inXasl.(FieldsA{iA})))
					outParms.(FieldsA{iA}) = inXasl.(FieldsA{iA});
				end
			end
		end
	end
end

% Last conversions
FieldsA = fields(outParms);
for iA = 1:length(FieldsA)
	if strcmp(FieldsA{iA},'AcquisitionTime') && ischar(outParms.AcquisitionTime)
		tP = xASL_adm_CorrectName(outParms.AcquisitionTime, 2);
		tP = xASL_adm_ConvertNr2Time(xASL_adm_ConvertTime2Nr(tP));
		if length(tP) > 7
			outParms.AcquisitionTime = [tP(1:6) '.' tP(7:8)];
		else
			outParms.AcquisitionTime = tP;
		end
	end
end

% Remove zeros from PhoenixProtocol
if isfield(outParms,'PhoenixProtocol')
	nullChar = char(0);
	outParms.PhoenixProtocol = strrep(outParms.PhoenixProtocol,nullChar,'');
end

% If output is in XASL, then copy into Q subfield
if bOutBids ~= 1
	FieldsA = fields(outParms);
	for iA = 1:length(FieldsA)
		for iL = find(strcmp(FieldsA{iA},xASLqFields))
			outParms.Q.(xASLqFields{iL}) = outParms.(xASLqFields{iL});
			outParms = rmfield(outParms, xASLqFields{iL});
		end
	end
end

% If we convert to xASL format, we have to move fields to correct position in the x structure
if ~bOutBids
    outParms = xASL_io_CheckDeprecatedFieldsX(outParms);
end


end


% Flatten the x structure to a single level
function inXasl = xASL_wrp_Quantify_FlattenX(inXasl)

    % Flatten inXasl structure: move x.Q, x.settings, & x.dataset main level of x
    firstLevelStructs = {'Q', 'settings', 'dataset'};
    
    % Iterate over first level structs
    for iStruct = 1:numel(firstLevelStructs)
        if isfield(inXasl,firstLevelStructs{iStruct})
            % Get the current fields
            subFields = fields(inXasl.(firstLevelStructs{iStruct}));
            % Move fields to main level of x
            for iField=1:length(subFields)
                inXasl.(subFields{iField}) = inXasl.(firstLevelStructs{iStruct}).(subFields{iField});
            end
            % Remove current substruct
            inXasl = rmfield(inXasl,firstLevelStructs{iStruct});
        end
    end


end




