addpath(genpath('c:\ASL_pipeline_HJ'))

ROOT    = 'E:\Backup\ASL_E\KCL_stats\INOX\analysis\dartel\StatsCompare\SVMtrial';

load('E:\Backup\ASL_E\KCL_stats\INOX\analysis\dartel\StatsCompare\SVMtrial\Expt_def.mat');

for iSub=1:15
    Expt_def.groups{1}.subjects{iSub}.data_files{1}   = ['E:\Backup\ASL_E\KCL_stats\INOX\analysis\dartel\StatsCompare\SVMtrial\CBF1_' num2str(iSub) '.nii'];
    Expt_def.groups{2}.subjects{iSub}.data_files{1}   = ['E:\Backup\ASL_E\KCL_stats\INOX\analysis\dartel\StatsCompare\SVMtrial\CBF2_' num2str(iSub) '.nii'];
end

save('E:\Backup\ASL_E\KCL_stats\INOX\analysis\dartel\StatsCompare\SVMtrial\Expt_def.mat','Expt_def');





% Expt_def.groups{1,1}.subjects{1,1} = rmfield(Expt_def.groups{1,1}.subjects{1,1},'data_files');

% Create mask

for iSub=1:15
    tnii                = xASL_io_ReadNifti( fullfile( ROOT, ['CBF1_' num2str(iSub) '.nii'] ) );
    tnii2               = xASL_io_ReadNifti( fullfile( ROOT, ['CBF2_' num2str(iSub) '.nii'] ) );
    TOTnii(:,:,:,iSub)  = tnii.dat(:,:,:);
    TOTnii(:,:,:,iSub*2)= tnii2.dat(:,:,:);
end

TOTnii  = xASL_stat_MeanNan(TOTnii,4)>0;
xASL_io_SaveNifti( fullfile( ROOT, ['CBF1_1.nii']), fullfile( ROOT, ['MaskSVM.nii']), TOTnii ) ;

piet =   xASL_io_ReadNifti( fullfile( ROOT, ['CBF1_1.nii']) );
piet    = piet.dat(:,:,:);
piet    = piet(TOTnii);
piet    = sort(piet);
klaas   = sort(class_n);
whos piet
whos TOTnii
