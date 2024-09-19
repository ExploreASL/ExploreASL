function parms = xASL_bids_Par2JSON(pathPar, pathJSON)
%xASL_bids_Par2JSON Reads the relevant DICOM fields from Philips PAR-file and adds them in a JSON sidecar
%
% FORMAT: parms = xASL_bids_Par2JSON(pathPar, pathJSON)
%
% INPUT:
%   pathPar   - path to the Philips PAR-file
%   pathJSON  - path to the output JSON sidecar or a list of files
%
% OUTPUT:
%   parms     - OPTIONAL - outputs the loaded/extracted parameters
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Opens the Philips PAR file. Reads the relevant DICOM headers and saves them to JSON sidecar in a BIDS format.
%              The JSON file is created automatically by the dcm2nii readout, so it always looks for this JSON file and 
%              add the same time reads the PAR file and adds further parameters to the JSON that were not identified by
%              the dcm2nii tool.
%
% EXAMPLE: parms = xASL_bids_Par2JSON('patient.par','ASL4D.json')
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2023 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


% Admin

if nargin<2 || isempty(pathPar) || isempty(pathJSON)
	error('Need the paths of .PAR and .JSON files.');
end

if iscell(pathJSON)
	nJSON = length(pathJSON);
else
	nJSON = 1;
end

for iJSON = 1:nJSON
	% Either load JSON or one from the list of JSONs
	if iscell(pathJSON)
		currentPathJSON = pathJSON{iJSON};
	else
		currentPathJSON = pathJSON;
	end

	% Parms already exist, so its content is loaded
	if exist(currentPathJSON, 'file')
		parms = xASL_io_ReadJson(currentPathJSON);
	end

	% Parse the PAR header
	hdr = xASL_adm_ParReadHeader(pathPar);

	% Extract DICOM tags
	if isfield(hdr,'InversionTime')
		parms.InversionTime = hdr.InversionTime;
	elseif isfield(parms,'InversionTime')
		parms.InversionTime = parms.InversionTime*1000;
	end
	parms.RepetitionTime                = hdr.RepetitionTime;
	parms.EchoTime                      = unique([hdr.SliceInformation(:).EchoTime]);
	parms.NumberOfTemporalPositions     = hdr.MaxNumberOfDynamics;
	parms.MRScaleSlope                  = unique([hdr.SliceInformation(:).ScaleSlope]);
	parms.RescaleSlopeOriginal          = unique([hdr.SliceInformation(:).RescaleSlope]);
	parms.RescaleIntercept              = unique([hdr.SliceInformation(:).RescaleIntercept]);
	switch (hdr.ScanMode)
		case '3D'
			parms.MRAcquisitionType             = '3D';
		case 'MS'
			parms.MRAcquisitionType             = '2D';
	end
	parms.VascularCrushing = logical(hdr.FlowCompensation);
	parms.AcquisitionDuration = hdr.ScanDuration;

	% Converts parameters from the legacy to the BIDS format
	parms = xASL_bids_parms2BIDS(parms, [], 1, 0);

	% To make sure that this is not removed for non-ASL sequences
	parms.RepetitionTime                = hdr.RepetitionTime/1000;

	% Saves the JSON sidecar
	xASL_io_WriteJson(currentPathJSON, parms);
end

return
