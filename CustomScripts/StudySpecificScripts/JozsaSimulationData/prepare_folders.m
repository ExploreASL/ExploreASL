function prepare_folders(EPAD_data_dir,copied_data_dir,data_par_template_file,quality,xASL_path,outfile)
% This is a script created by Tamas Jozsa to copy EPAD data
%
% EPAD_data_dir,copied_data_dir - strings determining the source and target locations
% data_par_template_file - a template of the parameter file used by ExploreASL
% quality - 0 for testing, 1 for HiFi processing
% xASL_path - a string with the path to the ExploreASL repository
% outfile - name of a .mat file as string, will include patients' features (age, sex, etc.) in copied_data_dir
%
% example:
% prepare_folders('/home/henk/ExploreASL/ASL/VirtualBrain/analysis_EPAD/','./EPAD_stats/','DATA_PAR_template.json',1,'./ExploreASL-1.1.2/','patient_info.mat')
clc;

%% 1. define files, folders, and quality
wdir = EPAD_data_dir; % note that its better to copy this outside of ExploreASL, to separate code and data
stat_dir = copied_data_dir;

data_par_template = regexp(fileread(data_par_template_file), '\r?\n', 'split');
data_par_template{10} = ['"Quality":', num2str(quality), ','];
data_par_template(end) = [];

%% 2. initialise xASL
CurrentPath = pwd;
cd(xASL_path);
ExploreASL_Master('',0); % initialize ExploreASL, without processing data
cd(CurrentPath); % back to the initial path

%% 3. compute statistical data for patients

% list of folders with ASL data
mkdir( fullfile(stat_dir) );
copyfile( fullfile(wdir,'participants.tsv'), ...
                  fullfile(stat_dir,'participants.tsv') );
[participants] = xASL_tsvRead(fullfile(stat_dir,'participants.tsv'));

patient_folders = participants(2:end,1);
Age = participants(2:end,4);
Sex = participants(2:end,6);
n_patients = length(patient_folders);

ASL_data_folder = 'ASL_1/';

parfile=cell(n_patients,1);
pass_check = false(n_patients,1);

for ii = 1:n_patients
    data_par_template{2} = ['"subject_regexp":"^',patient_folders{ii},'$",'];
    parfile{ii} = data_par_template;
end

for ii = 1:n_patients
    try
        % copy data to EPAD_stats
        mkdir( fullfile(stat_dir,patient_folders{ii}) );
        mkdir( fullfile(stat_dir,patient_folders{ii},patient_folders{ii}) );
        copyfile( fullfile(wdir,patient_folders{ii},'/FLAIR.nii.gz'), ...
                  fullfile(stat_dir,patient_folders{ii},patient_folders{ii},'/FLAIR.nii.gz') );
        copyfile( fullfile(wdir,patient_folders{ii},'/T1.nii.gz'), ...
                  fullfile(stat_dir,patient_folders{ii},patient_folders{ii},'/T1.nii.gz') );

        mkdir( fullfile(stat_dir,patient_folders{ii},patient_folders{ii},'/',ASL_data_folder) );
        copyfile( fullfile(wdir,patient_folders{ii},'/',ASL_data_folder,'/ASL4D.json'), ...
                  fullfile(stat_dir,patient_folders{ii},patient_folders{ii},'/',ASL_data_folder,'ASL4D.json') );
        copyfile( fullfile(wdir,patient_folders{ii},'/',ASL_data_folder,'/ASL4D.nii.gz'), ...
                  fullfile(stat_dir,patient_folders{ii},patient_folders{ii},'/',ASL_data_folder,'ASL4D.nii.gz') );
        copyfile( fullfile(wdir,patient_folders{ii},'/',ASL_data_folder,'/M0.json'), ...
                  fullfile(stat_dir,patient_folders{ii},patient_folders{ii},'/',ASL_data_folder,'M0.json') );
        copyfile( fullfile(wdir,patient_folders{ii},'/',ASL_data_folder,'/M0.nii.gz'), ...
                  fullfile(stat_dir,patient_folders{ii},patient_folders{ii},'/',ASL_data_folder,'M0.nii.gz') );

        % create setting file for xASL
        filePh = fopen([fullfile(stat_dir,patient_folders{ii}),'/DATA_PAR_',patient_folders{ii},'.json'],'w');
        fprintf(filePh,'%s\n',parfile{ii}{:});
        fclose(filePh);
        
        pass_check(ii)=true;
    catch
        disp( [patient_folders{ii}, ' cannot be processed because of missing files'] );
    end
end


for ii = 1:n_patients
    if pass_check(ii)
        disp( [patient_folders{ii}, ' will be processed'] );
    end

save(fullfile(stat_dir,outfile),'pass_check','parfile','patient_folders','n_patients','Age','Sex')
    
end
