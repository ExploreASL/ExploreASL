function xASL_io_Docker(root)
%xASL_io_Docker Function for the docker version of ExploreASL.
%
% FORMAT:       xASL_io_Docker(root);
%
% INPUT:        root - Path to DICOM dataset with BIDS(-like) structure.
%
% OUTPUT:       n/a
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function serves as a wrapper for the ExploreASL
%               workflow.
%
%               1. A predefined data structure is converted into NIFTI/JSON format.
%               The function expects an 'incoming' folder, containing a
%               'raw' folder. The next levels are a 'subject' and a 'visist'
%               folder. Within the 'visit' folder are optimally four
%               subfolders: ASL, T1, FLAIR and M0.
%               2. The output of the ExploreASL_Import function is
%               restructured to enable the following step.
%               3. The last step is the execution of ExploreASL_Master on
%               the generated NIFTI/JSON data structure.
%
% Detailed description of the incoming data structure:
% 
%               - /incoming/raw/sub-###/visit-###/ASL
%               - /incoming/raw/sub-###/visit-###/T1
%               - /incoming/raw/sub-###/visit-###/FLAIR
%               - /incoming/raw/sub-###/visit-###/M0
%
% Structure we want to have for ExploreASL_Master:
%
%               - ...
%
% EXAMPLE:      xASL_io_Docker('/opt/incoming/');
%
%               xASL_io_Docker('M:\SoftwareDevelopment\MATLAB\m.stritt\tmp_data_siemens\incoming');
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Workflow

% Define path to DataParFile
DataParFile = fullfile(root,'data','DataParFile.json');

% (1) Read in DICOM data and convert the data to the NIFTI/JSON structure
ExploreASL_Import(ExploreASL_ImportConfig(root),false,true);

% (2) Make the structure readable for ExploreASL_Master

% Rename folder
movefile(fullfile(root,'analysis'),fullfile(root,'data'));

% Generate DataParFile
fID = fopen(DataParFile,'w');

% Parameters
data.x.name = "incoming";
data.x.subject_regexp = "^sub$";
data.x.Quality = 1;
data.x.bNativeSpaceAnalysis = 1;

% Write data to JSON file
JSONstr = jsonencode(data);
if fID == -1, error('Cannot create JSON file...'); end
fwrite(fID, JSONstr, 'char');
fclose(fID);

% (3) Run ExploreASL_Master
ExploreASL_Master(DataParFile, true);



