function parms = xASL_adm_Par2Parms(pathPar, pathParms, bRecreate)
% Read relevant dicom fields from PAR file and saves them as mat-file
%
% FORMAT: parms = xASL_adm_Par2Parms(pathPar, pathParms[, bRecreate])
%
% INPUT:
%   pathPar   - path to the philips PAR-file
%   pathParms - path to the output MAT-file
%   bRecreate - Recreate the parms file if already exists (OPTIONAL, DEFAULT false)
%
% OUTPUT:
%   parms     - OPTIONAL - outputs the loaded/extracted parameters
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Opens the Philips type PAR file. Reads the relevant DICOM headers and saves them to .MAT file.
%              Only recreates an existing file if bRecreate option is set to TRUE.
%
% EXAMPLE: parms = xASL_adm_Par2Parms('patient.par','ASL4D.mat',true)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2019 ExploreASL

% Admin

if nargin<2 || isempty(pathPar) || isempty(pathParms)
	error('Need the paths of .PAR and .MAT files.');
end

if nargin<3 || isempty(bRecreate)
    bRecreate = false;
end

% Parms already exist
if ~bRecreate && exist(pathParms, 'file')
	if nargout > 0
		fprintf('Loading parameters from %s\n',pathParms);
		X = load(pathParms);
		parms = X.parms;
	end
	return;
end

% Recreate again
fprintf('Recreating parameter file: %s\n',pathParms);
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

% Save to .MAT file
save(pathParms, 'parms');
return
