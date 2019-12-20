%% Admin
clear
x.MYPATH   = 'c:\ASL_pipeline_HJ';
AdditionalToolboxDir    = 'C:\ASL_pipeline_HJ_toolboxes'; % provide here ROOT directory of other toolboxes used by this pipeline, such as dip_image & SPM12
if ~isdeployed
    addpath(x.MYPATH);

    subfolders_to_add = { 'ANALYZE_module_scripts', 'ASL_module_scripts', fullfile('Development','dicomtools'), fullfile('Development','Filter_Scripts_JanCheck'), 'MASTER_scripts', 'spm_jobs','spmwrapperlib' };
    for ii=1:length(subfolders_to_add)
        addpath(fullfile(x.MYPATH,subfolders_to_add{ii}));
    end
end

addpath(fullfile(AdditionalToolboxDir,'DIP','common','dipimage'));

[x.SPMDIR, x.SPMVERSION] = xASL_adm_CheckSPM('FMRI',fullfile(AdditionalToolboxDir,'spm12') );
addpath( fullfile(AdditionalToolboxDir,'spm12','compat') );

if isempty(which('dip_initialise'))
    fprintf('%s\n','CAVE: Please install dip_image toolbox!!!');
else dip_initialise
end

%% 

ROOT        = 'C:\Backup\ASL\PatriciaMood\GIFMI_STUMOOD\analysis';
SUBJECTS    = xASL_adm_GetFsList(ROOT, '^SUB\d{3}$',1);

for iS=1:length(SUBJECTS)
    Sex{iS,1}       = SUBJECTS{iS};
    NUMBERsex       = str2num(SUBJECTS{iS}(4:end));
    if      NUMBERsex<25
            Sex{iS,2}   = 2; % female
    else    Sex{iS,2}   = 1; % male
    end
end

save( fullfile(ROOT,'Sex.mat'),'Sex');
