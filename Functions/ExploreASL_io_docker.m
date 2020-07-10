function ExploreASL_io_docker(root, pathDataParFile)
%ExploreASL_docker Function for the docker version of ExploreASL.
%
% FORMAT:       ExploreASL_io_docker(root, pathDataParFile);
%
% INPUT:        root                - Path to DICOM dataset with BIDS(-like) structure.
%               pathDataParFile     - Path of the DataParFile within the docker container.
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
% EXAMPLE:      ExploreASL_io_docker('/opt/incoming/', '/opt/incoming/analysis/DataParFile.json');
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Workflow

% (1) Read in DICOM data and convert the data to the NIFTI/JSON structure
ExploreASL_Import(ExploreASL_ImportConfig(root),false,true);


% (2) Make the structure readable for ExploreASL_Master

% ... WORK IN PROGRESS ...

% ... xASL copy, rename etc. ...


% (3) Run ExploreASL_Master
ExploreASL_Master(pathDataParFile, true);



