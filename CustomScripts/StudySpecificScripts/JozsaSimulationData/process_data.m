function process_data(start_num,end_num,data_dir,xASL_path,patient_info)
% This is a function created by Tamas Jozsa to process EPAD data
%
% start_num and end_num - integers specifying the range of patients
% data_dir - string determining the work directory
% xASL_path - a string with the path to the ExploreASL repository
% patient info - string of the .mat file containing patients' features
%
% example: process_data(1,20,'./EPAD_stats/','./ExploreASL-1.1.2/','patient_info.mat')
clc;

%% 1. define files and folders

stat_dir = data_dir;

outfile = 'patient_stats.csv';


%% 2. initialise xASL
CurrentPath = pwd;
cd(xASL_path);
ExploreASL_Master('',0); % initialize ExploreASL, without processing data
cd(CurrentPath); % back to the initial path

%% 3. compute statistical data for patients

% list of folders with ASL data
load(fullfile(stat_dir,patient_info));
ASL_data_folder = 'ASL_1/';



for ii = start_num:end_num
    if pass_check(ii)
        disp( [patient_folders{ii}, ' will be processed'] );
        cd(xASL_path);
        ExploreASL_Master([fullfile(CurrentPath,stat_dir,patient_folders{ii}),'/DATA_PAR_',patient_folders{ii},'.json'], true, true);
        cd(CurrentPath);
    else
        disp( [patient_folders{ii}, ' cannot be processed because of missing files'] );
    end
    
end
