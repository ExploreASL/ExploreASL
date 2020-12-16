function jsonOut = ExploreASL_bids_VendorFieldCheck(jsonIn)
%ExploreASL_bids_VendorFieldCheck Goes through the JSON structure before saving it and ensures that all 
% vendor specific fields are properly checked and renamed if necessary
%
% FORMAT: jsonOut = ExploreASL_bids_VendorFieldCheck(jsonIn,bIsASL)
%
% INPUT:
%   jsonIn  - JSON with the input fields (REQUIRED)
%
% OUTPUT: 
%   jsonOut - ordered and checked JSON structure
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:
% It checks all the JSON fields, make sure that they are renamed from vendor specific names to common BIDS names
%
% __________________________________
% Copyright 2015-2020 ExploreASL
jsonOut = struct;
jsonRemove = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rename certain fields from the vendor-name to BIDS name

% Rename the coil names for different manufacturers
if isfield(jsonIn,'CoilString')
	switch (jsonIn.Manufacturer)
		case 'Philips'
			jsonOut.ReceiveCoilName = jsonIn.CoilString;
		case 'GE'
			jsonOut.ReceiveCoilName = jsonIn.CoilString;
		case 'Siemens'
			jsonOut.ReceiveCoilActiveElements = jsonIn.CoilString;
		otherwise
			error('Unknown manufacturer')
	end
	jsonRemove.CoilString = '';
end

if isfield(jsonIn,'EffectiveEchoSpacing') 
	if isempty(jsonIn.EffectiveEchoSpacing) || sum(jsonIn.EffectiveEchoSpacing == 0) || sum(isnan(jsonIn.EffectiveEchoSpacing))
		jsonRemove.EffectiveEchoSpacing = '';
	else
		jsonIn.EffectiveEchoSpacing = abs(jsonIn.EffectiveEchoSpacing);
	end
end

if isfield(jsonIn,'TotalReadoutTime') 
	if isempty(jsonIn.TotalReadoutTime) || sum(jsonIn.TotalReadoutTime == 0) || sum(isnan(jsonIn.TotalReadoutTime))
		jsonRemove.TotalReadoutTime = '';
	else
		jsonIn.TotalReadoutTime = abs(jsonIn.TotalReadoutTime);
	end
end

% Rename the fields with number of segments
jsonRemove.NumberOfAverages = '';
if isfield(jsonIn,'NumberSegments')
	jsonOut.NumberShots = jsonIn.NumberSegments;
	jsonRemove.NumberSegments = '';
end

% Rename the phase encoding directions fields
if isfield(jsonIn,'PhaseEncodingAxis')
	if ~isfield(jsonIn,'PhaseEncodingDirection')
		jsonOut.PhaseEncodingDirection = jsonIn.PhaseEncodingAxis;
	end
	jsonRemove.PhaseEncodingAxis = '';
end

% Go through all input fields and copy to output, but skip those in jsonRemove
for nameField = fieldnames(jsonIn)'
	if ~isfield(jsonRemove, nameField{1})
		jsonOut.(nameField{1}) = jsonIn.(nameField{1});
	end
end


end