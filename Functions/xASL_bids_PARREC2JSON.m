function parms = xASL_bids_PARREC2JSON(pathPar, PathJSON)
%xASL_bids_PARREC2JSON Read relevant dicom fields from Philips PAR file and save them as JSON-file
%
% FORMAT: parms = xASL_adm_Par2Parms(pathPar, PathJSON)
%
% INPUT:
%   pathPar   - path to the philips PAR-file
%   PathJSON - path to the JSON file to add fields to (or create if non-existing)
%
% OUTPUT:
%   parms     - OPTIONAL - outputs the loaded/extracted parameters
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Opens the Philips type PAR file. Reads the relevant DICOM header fields and adds them to the .json sidecar file.
%
% EXAMPLE: parms = xASL_adm_Par2Parms('patient.par','ASL4D.json',true)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2019 ExploreASL

% Admin

if nargin<2 || isempty(pathPar) || isempty(PathJSON)
	error('Need the paths of .PAR and .json files');
end

% Parms already exist
if exist(PathJSON, 'file')
	if nargout > 0
		fprintf('Loading parameters from %s\n', PathJSON);
		X = load(PathJSON);
		parms = X.parms;
	end
	return;
end

% Recreate again
fprintf('Recreating parameter file: %s\n', PathJSON);
parms = struct();

% Parse header
hdr = xASL_adm_ParReadHeader(pathPar);

% Extract DICOM tags
parms.RepetitionTime                = hdr.RepetitionTime;
parms.EchoTime                      = unique([hdr.SliceInformation(:).EchoTime]);
parms.NumberOfTemporalPositions     = hdr.MaxNumberOfDynamics;
parms.MRScaleSlope                  = unique([hdr.SliceInformation(:).ScaleSlope]);
parms.RescaleSlopeOriginal          = unique([hdr.SliceInformation(:).RescaleSlope]);
parms.RescaleIntercept              = unique([hdr.SliceInformation(:).RescaleIntercept]);

% Save to JSON file
xASL_bids_InsertJSONFields(parms, PathJSON);


end
