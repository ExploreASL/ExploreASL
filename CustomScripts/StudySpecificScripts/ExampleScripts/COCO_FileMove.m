%% FileMove

addpath(genpath('c:\ASL_pipeline_HJ'));

Root    = 'E:\Backup\ASL_E\COCO\nifti_files_HJ';
List    = xASL_adm_GetFileList(Root, '^ASL4D\.(nii|nii\.gz)$', 'FPListRec');
for iList =2:length(List);
    [path file ext]     = fileparts(List{iList});
    mkdir(fullfile(path,'ASL_1'));
    NewFile             = fullfile(path,'ASL_1',[file ext]);
    xASL_Move( List{iList}, NewFile);
    clear path file ext
end

List    = xASL_adm_GetFileList(Root, '^ASL4D_parms\.mat$', 'FPListRec');
for iList =1:length(List);
    [path file ext]     = fileparts(List{iList});
    mkdir(fullfile(path,'ASL_1'));
    NewFile             = fullfile(path,'ASL_1',[file ext]);
    xASL_Move( List{iList}, NewFile);
    clear path file ext
end


%% Unzip & rename
addpath(genpath('c:\ASL_pipeline_HJ'));

Root    = 'E:\Backup\ASL_E\COCO\segm\segm';
List    = xASL_adm_GetFileList(Root, '^.*\.nii\.gz$', 'FPListRec');

for iList=1:length(List) 
    %gunzip(List{iList});
    %delete(List{iList});
    xASL_adm_UnzipNifti(List{iList});
end

List    = xASL_adm_GetFileList(Root, '^lesprob.*\.(nii|nii\.gz)$', 'FPListRec');
for iList=1:length(List)
    xASL_Rename( List{iList}, 'WMH_SEGM.nii', 1);
end
List    = xASL_adm_GetFileList(Root, '^wnu.*FLAIR\.(nii|nii\.gz)$', 'FPListRec');
for iList=1:length(List)
    xASL_Rename( List{iList}, 'FLAIR.nii', 1);
end
