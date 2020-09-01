function outParms = xASL_bids_parms2BIDS(inXasl, inBids, bOutBids, bPriorityBids)
% Takes the input parameters from xASL format (parms) and BIDS format, merges them and converts to either xASL or BIDS.
% FORMAT: outBids = xASL_bids_parms2BIDS(inParms[, inBids, bOutBids, priorityBids])
% 
% INPUT:
%   inXasl       - a structure with input parameters in the xASL format (REQUIRED)
%   inBids       - a structure with input parameters in the BIDS format (OPTIONAL, DEFAULT = [])
%   bOutBids     - the output structures is in BIDS format (==1, default) or xASL format (==0) (OPTIONAL, DEFAULT = 1)
%   bPriorityBids - in case of conflicts, the BIDS input is preferred (==1, default), otherwise (==0), xASL is prefered (OPTIONAL, DEFAULT = 1)
%
% OUTPUT:
% outParms       - the merged output structure in the correct format
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
% Copyright 2015-2020 ExploreASL


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

% Fields with these names need to have the time converted between XASL and BIDS, and define their recommended range in ms
convertTimeFieldsXASL = {'EchoTime' 'RepetitionTime' 'Initial_PLD' 'LabelingDuration' 'SliceReadoutTime' 'BloodT1' 'T2' 'TissueT1'};
convertTimeFieldsRange = [0.5        5                10            5                  5                  100       10   100;...% Minimum in ms
                          500        20000            10000         5000               400                5000      500  5000];% Maximum in ms   
					  
% Fields that are entered under the subfield 'Q' for xASL on the output
xASLqFields = {'LabelingType' 'Initial_PLD' 'BackGrSupprPulses' 'LabelingDuration' 'SliceReadoutTime' 'NumberOfAverages' 'BloodT1'};

% Some JSON fields need to be updated to fit the BIDS definition
% This is a one way process of changing the names from old to new
updateNamesBIDSold = {'PhilipsRescaleSlope' 'PhilipsRWVSlope' 'PhilipsScaleSlope' 'PhilipsRescaleIntercept' 'PhilipsRWVIntercept'};
updateNamesBIDSnew = {'RescaleSlope'        'RWVSlope'        'MRScaleSlope'      'RescaleIntercept'        'RWVIntercept'};

% These fields have different names in xASL and in BIDS
% They are therefore renamed depending on the type of output
changeNamesXASL = {'Vendor'       'readout_dim'       'Initial_PLD'};
changeNamesBIDS = {'Manufacturer' 'MRAcquisitionType' 'InitialPostLabelDelay'};


%% ----------------------------------------------------------------------
%% 2) Convert XASL fields to the output format (BIDS or XASL)

% Goes through all XASL fields
if ~isempty(inXasl)
	% Move all fields from Q to the main structure
	if isfield(inXasl,'Q')
		FieldsQ = fields(inXasl.Q);
		for iQ=1:length(FieldsQ)
			inXasl.(FieldsQ{iQ}) = inXasl.Q.(FieldsQ{iQ});
		end
		inXasl = rmfield(inXasl,'Q');
	end
	
	% If output is in BIDS, then XASL fields need to be converted
	if bOutBids
		FieldsA = fields(inXasl);
		for iA = 1:length(FieldsA)
			
			% Convert the units for all time fields from ms to s
			for iT = find(strcmp(FieldsA{iA},convertTimeFieldsXASL))
				% For non-zero fields, check if they are within the predefined range
				if inXasl.(FieldsA{iA}) ~= 0
					% If outside of the recommended range, then still convert, but issue a warning
					if inXasl.(FieldsA{iA}) < convertTimeFieldsRange(1,iT) || inXasl.(FieldsA{iA}) > convertTimeFieldsRange(2,iT)
						warning(['Field ' FieldsA{iA} ' in xASL structure has a value ' num2str(inXasl.(FieldsA{iA}))...
							', which is outside of the recommended range <'...
							num2str(convertTimeFieldsRange(1,iT)) ',' num2str(convertTimeFieldsRange(2,iT)) '> ms.']);
					end
				end
				inXasl.(FieldsA{iA}) = inXasl.(FieldsA{iA})/1000;
			end
			
			% Rename XASL fields to BIDS field according to the information in changeNamesXASL and changeNamesBIDS
			for iL = find(strcmp(FieldsA{iA},changeNamesXASL))
				inXasl.(changeNamesBIDS{iL}) = inXasl.(changeNamesXASL{iL});
				inXasl = rmfield(inXasl,changeNamesXASL{iL});
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
		for iL = find(strcmp(FieldsA{iA},updateNamesBIDSold))
			inBids.(updateNamesBIDSnew{iL}) = inBids.(updateNamesBIDSold{iL});
			inBids = rmfield(inBids,updateNamesBIDSold{iL});
		end
	end
	
	% When the output is in XASL we need to convert from BIDS to XASL
	if bOutBids ~= 1
		FieldsA = fields(inBids);
		for iA = 1:length(FieldsA)
			% Rename all listed fields to XASL variant
			FieldNameChanged = FieldsA{iA};
			for iL = find(strcmp(FieldsA{iA},changeNamesBIDS))
				inBids.(changeNamesXASL{iL}) = inBids.(changeNamesBIDS{iL});
				inBids = rmfield(inBids,changeNamesBIDS{iL});
				% Update the name of the field after the change
				FieldNameChanged = changeNamesXASL{iL};
			end
			
			% Convert the units for all time fields from s to ms
			for iT = find(strcmp(FieldNameChanged,convertTimeFieldsXASL))
				inBids.(FieldNameChanged) = inBids.(FieldNameChanged)*1000;
				
				% Check if the value is within the recommended range after conversion and issue a warning if not
				if inBids.(FieldNameChanged) ~= 0
					if inBids.(FieldNameChanged) < convertTimeFieldsRange(1,iT) || inBids.(FieldNameChanged) > convertTimeFieldsRange(2,iT)
						warning(['Field ' FieldNameChanged ' in xASL structure has a value ' num2str(inBids.(FieldNameChanged))...
							', which is outside of the recommended range <'...
							num2str(convertTimeFieldsRange(1,iT)) ',' num2str(convertTimeFieldsRange(2,iT)) '> ms.']);
					end
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
			if (~isfield(outParms,(FieldsA{iA}))) ||...
					((sum(isnan(inBids.(FieldsA{iA})))==0) && (sum(inBids.(FieldsA{iA}) == 0)~=length(inBids.(FieldsA{iA}))) && (sum(isinf(inBids.(FieldsA{iA})))==0) && (~isempty(inBids.(FieldsA{iA}))))
				outParms.(FieldsA{iA}) = inBids.(FieldsA{iA});
			end
		end
	end
else
	outParms = inBids;
	FieldsA = fields(inXasl);
	for iA = 1:length(FieldsA)
		if (~isfield(outParms,(FieldsA{iA}))) ||...
				((sum(isnan(inXasl.(FieldsA{iA})))==0) && (sum(inXasl.(FieldsA{iA}) == 0)~=length(inXasl.(FieldsA{iA}))) && (sum(isinf(inXasl.(FieldsA{iA})))==0) && (~isempty(inXasl.(FieldsA{iA}))))
			outParms.(FieldsA{iA}) = inXasl.(FieldsA{iA});
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

% If output is in XASL, then copy into Q subfield
if bOutBids ~= 1
	FieldsA = fields(outParms);
	for iA = 1:length(FieldsA)
		for iL = find(strcmp(FieldsA{iA},xASLqFields))
			outParms.Q.(xASLqFields{iL}) = outParms.(xASLqFields{iL});
			outParms = rmfield(outParms,xASLqFields{iL});
		end
	end
end

end
	