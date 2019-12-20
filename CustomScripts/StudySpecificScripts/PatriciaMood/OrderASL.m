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

%% BEFORE IMPORT

ROOT    = 'C:\Backup\ASL\PatriciaMood\GIFMI_STUMOOD\raw';
Dlist   = xASL_adm_GetFsList(ROOT,'^SUB\d{3}$',1);

for iD=1:length(Dlist)
    % label
    Dlist2  = xASL_adm_GetFsList( fullfile(ROOT, Dlist{iD},Dlist{iD}),'^NS_TE00_TI1700_\d{4}$',1);
    for iD2=1:length(Dlist2)
        if  iD2<5
            xASL_Rename( fullfile(ROOT, Dlist{iD},Dlist{iD},Dlist2{iD2}),[Dlist2{iD2}(4:end-4) '0' num2str(iD2*2)]);
        else
            xASL_Rename( fullfile(ROOT, Dlist{iD},Dlist{iD},Dlist2{iD2}),[Dlist2{iD2}(4:end-4) num2str(iD2*2)]);
        end
    end

    % control
    Dlist2  = xASL_adm_GetFsList( fullfile(ROOT, Dlist{iD},Dlist{iD}),'^SS_TE00_TI1700_\d{4}$',1);
    for iD2=1:length(Dlist2)
        if  iD2<6
            xASL_Rename( fullfile(ROOT, Dlist{iD},Dlist{iD},Dlist2{iD2}),[Dlist2{iD2}(4:end-4) '0' num2str(iD2*2-1)]);
        else
            xASL_Rename( fullfile(ROOT, Dlist{iD},Dlist{iD},Dlist2{iD2}),[Dlist2{iD2}(4:end-4) num2str(iD2*2-1)]);
        end
    end
end

%% AFTER IMPORT
clear
ROOT    = 'C:\Backup\ASL\PatriciaMood\GIFMI_STUMOOD\analysis';
Dlist   = xASL_adm_GetFsList(ROOT,'^SUB\d{3}$',1);

for iD=1:length(Dlist)

    for ii=1:19
        clear ASLname1 ASLname2 ASLmat1 ASLmat2 tempNII1 tempNII2 parms1 parms2 IM

        ASLname1     = fullfile( ROOT, Dlist{iD}, ['ASL_' num2str(ii*2-1)],'ASL4D.nii');
        ASLmat1      = fullfile( ROOT, Dlist{iD}, ['ASL_' num2str(ii*2-1)],'ASL4D_parms.mat');
        ASLname2     = fullfile( ROOT, Dlist{iD}, ['ASL_' num2str(ii*2)],'ASL4D.nii');
        ASLmat2      = fullfile( ROOT, Dlist{iD}, ['ASL_' num2str(ii*2)],'ASL4D_parms.mat');

        NewDir       = fullfile(ROOT, Dlist{iD}, ['ASL_' num2str(ii)]);

        tempNII1     = xASL_io_ReadNifti(ASLname1);
        tempNII2     = xASL_io_ReadNifti(ASLname2);

        parms1       = load(ASLmat1);
        parms2       = load(ASLmat2);

        if      parms1.parms.RepetitionTime~=parms2.parms.RepetitionTime || parms1.parms.EchoTime~=parms2.parms.EchoTime
                error('Parameters');
        elseif  parms1.parms.AcquisitionTime>parms2.parms.AcquisitionTime
                error('Parameters');
        end

        if      min(size(tempNII1.dat)~=size(tempNII2.dat))
                error('Nifti matrix size');
        end

        if      size(tempNII1.dat,4)>4 || size(tempNII2.dat,4)>4
                error('SizeNii');
        end

        IM(:,:,:,1)     = tempNII1.dat(:,:,:);
        IM(:,:,:,2)     = tempNII2.dat(:,:,:);

        delete(ASLname1);
        delete(ASLmat1);
        xASL_adm_CreateDir( NewDir );
        xASL_io_SaveNifti( ASLname2, fullfile(NewDir,'ASL4D.nii'), IM, size(IM,4) );
        xASL_Move(ASLmat2,fullfile(NewDir,'ASL4D_parms.mat'));

        delete(ASLname2);
        rmdir(fullfile( ROOT, Dlist{iD}, ['ASL_' num2str(ii*2)]));
        if  ii>1
            rmdir(fullfile( ROOT, Dlist{iD},['ASL_' num2str(ii*2-1)]));
            xASL_Copy(fullfile( ROOT, Dlist{iD},'ASL_1','M0.nii'),fullfile(NewDir,'M0.nii') );
            xASL_Copy(fullfile( ROOT, Dlist{iD},'ASL_1','M0_parms.mat'),fullfile(NewDir,'M0_parms.mat') );
        end
    end
end
