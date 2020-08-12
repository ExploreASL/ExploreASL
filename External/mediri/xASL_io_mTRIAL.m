function xASL_io_mTRIAL(root)
%xASL_dcm_Import Function for the docker version of ExploreASL.
%
% FORMAT:       xASL_io_mTRIAL(root);
%
% INPUT:        root - Path to DICOM dataset with BIDS(-like) structure.
%
% OUTPUT:       n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function is used for the DICOM import in the docker workflow:
%
%               1. A predefined data structure is converted into NIFTI/JSON format.
%               The function expects an 'incoming' folder, containing a
%               'raw' folder. The next levels are a 'subject' and a 'visist'
%               folder. Within the 'visit' folder are optimally four
%               subfolders: ASL, T1, FLAIR and M0.
%               2. The output of the ExploreASL_Import function is
%               restructured to enable the following step.
%
% Detailed description of the incoming data structure:
% 
%               - /incoming/dockerInterfaceParameters.json
%               - /incoming/raw/sub-###/visit-###/ASL
%               - /incoming/raw/sub-###/visit-###/T1
%               - /incoming/raw/sub-###/visit-###/FLAIR
%               - /incoming/raw/sub-###/visit-###/M0
%
% dockerInterfaceParameters.json: (Siemens e.g.)
%
% {
% "Sequence":          "3D GRASE",
% "LabelingType":      "PCASL",
% "readout_dim":       "2D",
% "BackGrSupprPulses": "5",
% "Initial_PLD":       "2000",
% "LabelingDuration":  "1800",
% "SliceReadoutTime":  "40"
% }
% 
% BIDS TEST DATA:
% 
% 'Siemens'
% data.x.M0 = 'UseControlAsM0'; % Should actually be: data.x.M0 = 'separate_scan'; % but separate M0 does not work on mTRIAL right now
% data.x.Sequence = '3D_GRASE';
% data.x.Q.LabelingType = "PCASL";
% data.x.Q.BackGrSupprPulses = 5;
% data.x.Q.LabelingDuration = 1800;
% data.x.Q.Initial_PLD = 2000;
%
% 'Philips'
% data.x.M0 = 'UseControlAsM0'; % Should actually be: data.x.M0 = 'separate_scan'; % but separate M0 does not work on mTRIAL right now
% data.x.Sequence = '2D_EPI';
% data.x.Q.LabelingType = "PCASL";
% data.x.Q.BackGrSupprPulses = 2;
% data.x.Q.LabelingDuration = 1800;
% data.x.Q.Initial_PLD = 2000;
%
% 'GE'
% data.x.M0 = 'UseControlAsM0';
% data.x.Sequence = '3D_spiral';
% data.x.Q.LabelingType = "PCASL";
% data.x.Q.BackGrSupprPulses = 5;
% data.x.Q.LabelingDuration = 2025;
% data.x.Q.Initial_PLD = 1450;
% 
%
% EXAMPLE:      xASL_io_mTRIAL('/opt/incoming/');
%
%               xASL_io_mTRIAL('C:\...\incoming');
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Workflow

% Define path to dockerInterfaceParameters
dockerInterfaceFile = fullfile(root,'dockerInterfaceParameters.json');

% Define path to DataParFile
DataParFile = fullfile(root,'data','DataParFile.json');

% Read in DICOM data and convert the data to the NIFTI/JSON structure
ExploreASL_Import(ExploreASL_ImportConfig(root),false,true);

% Make the structure readable for ExploreASL_Master

% Rename folder
movefile(fullfile(root,'analysis'),fullfile(root,'data'));

% Parameters
data.x.name = "incoming";
data.x.subject_regexp = "^sub$";
data.x.Quality = 1;
data.x.bNativeSpaceAnalysis = 1;
data.x.DELETETEMP = 1;

% Add Q field
data.x.Q = struct;
data.x.M0 = 'UseControlAsM0';

% Read docker interface file
x_temporary = xASL_import_json(dockerInterfaceFile);

% Reassign fields
if isfield(x_temporary,'Sequence'),             data.x.Sequence = x_temporary.Sequence; end
if isfield(x_temporary,'LabelingType'),         data.x.Q.LabelingType = x_temporary.LabelingType; end
if isfield(x_temporary,'readout_dim'),          data.x.readout_dim = x_temporary.readout_dim; end
if isfield(x_temporary,'BackGrSupprPulses'),    data.x.Q.BackGrSupprPulses = x_temporary.BackGrSupprPulses; end
if isfield(x_temporary,'Initial_PLD'),          data.x.Q.Initial_PLD = x_temporary.Initial_PLD; end
if isfield(x_temporary,'LabelingDuration'),     data.x.Q.LabelingDuration = x_temporary.LabelingDuration; end

% Fix Sequence name
switch data.x.Sequence
    case "3D Spiral"
        data.x.Sequence = '3D_spiral';
    case "3D GRASE"
        data.x.Sequence = '3D_GRASE';
    case "2D EPI"
        data.x.Sequence = '2D_EPI';
end

% Get vendor & Acquisition
pathASL4D = fullfile(root,'data','sub','ASL_1','ASL4D.json');
try
    % Read JSON file
    if xASL_exist(pathASL4D,'file')
        val = spm_jsonread(pathASL4D);
		valXASL = xASL_bids_parms2BIDS([], val, 0, 1);
        % Set vendor to manufacturer and remove slice readout time for 3D cases
        if ~isfield(valXASL,'Vendor'), data.x.Vendor = valXASL.Vendor; end
		
		if isfield(valXASL,'readout_dim')
			if strcmp(valXASL.readout_dim,"2D")
				data.x.Q.SliceReadoutTime = x_temporary.SliceReadoutTime;
			end
		end
    end
catch
    warning('Something went wrong...');
end

% Write data to JSON file
spm_jsonwrite(DataParFile,data);





