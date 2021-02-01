function parms = xASL_bids_Par2JSON(pathPar, pathJSON, bRecreate)
%xASL_bids_Par2JSON Reads the relevant DICOM fields from Philips PAR-file and saves them in a JSON sidecar
%
% FORMAT: parms = xASL_bids_Par2JSON(pathPar, pathJSON[, bRecreate])
%
% INPUT:
%   pathPar   - path to the Philips PAR-file
%   pathJSON  - path to the output JSON sidecar
%   bRecreate - Recreate the parms file if already exists (OPTIONAL, DEFAULT false)
%
% OUTPUT:
%   parms     - OPTIONAL - outputs the loaded/extracted parameters
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Opens the Philips PAR file. Reads the relevant DICOM headers and saves them to JSON sidecar in a BIDS format.
%              Only recreates an existing file if bRecreate option is set to TRUE.
%
% EXAMPLE: parms = xASL_bids_Par2JSON('patient.par','ASL4D.json',true)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2020 ExploreASL

% Admin

if nargin<2 || isempty(pathPar) || isempty(pathJSON)
	error('Need the paths of .PAR and .JSON files.');
end

if nargin<3 || isempty(bRecreate)
    bRecreate = false;
end

% Parms already exist
if ~bRecreate && exist(pathJSON, 'file')
	if nargout > 0
		fprintf('Loading parameters from Philips PAR file: %s\n', pathJSON);
		parms = spm_jsonread(pathJSON);
	end
	return;
end

% Recreate again
fprintf('Recreating parameter file in a JSON sidecar: %s\n', pathJSON);
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

% Converts parameters from the legacy to the BIDS format
parms = xASL_bids_parms2BIDS(parms, [], 1, 0);
		
% Saves the JSON sidecar
spm_jsonwrite(pathJSON, parms);

return
