function xASL_dcm_Import(root)
%xASL_dcm_Import Function for the docker version of ExploreASL.
%
% FORMAT:       xASL_dcm_Import(root);
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
%               - /incoming/raw/sub-###/visit-###/ASL
%               - /incoming/raw/sub-###/visit-###/T1
%               - /incoming/raw/sub-###/visit-###/FLAIR
%               - /incoming/raw/sub-###/visit-###/M0
%
% EXAMPLE:      xASL_dcm_Import('/opt/incoming/');
%
%               xASL_dcm_Import('M:\SoftwareDevelopment\MATLAB\m.stritt\tmp_data_siemens\incoming');
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL

%% Workflow

% Initialize ExploreASL
x = ExploreASL_Initialize('',0);

% Define path to DataParFile
DataParFile = fullfile(root,'data','DataParFile.json');

% Read in DICOM data and convert the data to the NIFTI/JSON structure
ExploreASL_Import(ExploreASL_ImportConfig(root),false,true);

% Make the structure readable for ExploreASL_Master

% Rename folder
movefile(fullfile(root,'analysis'),fullfile(root,'data'));

% Create log file
diary(fullfile(root,'data','xASL_module_Import.log'))
try
    % Open file
    fIDcsv = fopen(fullfile(root,'data','import_summary.csv'),'r');
    % Get individual text lines
    while ~feof(fIDcsv)
        tline = fgetl(fIDcsv);
        lineElements = strsplit(tline,',');
        % Get subjects
        if strcmp('"sub"',char(lineElements(1)))
            fprintf('Subject:\t%s\n', lineElements{2})
            fprintf('Visit:\t\t%s\n', lineElements{4})
            fprintf('NIFTI:\t\t%s\n', lineElements{5})
        end
        % disp(tline)
    end
    % Close file
    fclose(fIDcsv);

catch
    
end
% Turn off diary
diary OFF

% Generate DataParFile
fID = fopen(DataParFile,'w');

% Parameters
data.x.name = "incoming";
data.x.subject_regexp = "^sub$";
data.x.Quality = 1;
data.x.bNativeSpaceAnalysis = 1;
data.x.DELETETEMP = 1;
data.x.readout_dim = "2D";

% Write data to JSON file
JSONstr = jsonencode(data);
if fID == -1, error('Cannot create JSON file...'); end
fwrite(fID, JSONstr, 'char');
fclose(fID);

% Compare ASL4D and M0 JSON files with list of parameters (from TestDataSet)
pathASL4D = fullfile(root,'data','sub','ASL_1','ASL4D.json');
xASL_par_Fix(DataParFile,pathASL4D);














