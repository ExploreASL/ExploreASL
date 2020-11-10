ExploreASL_Master('',0);

%% Prisma 3D pCASL

% mDir    = 'C:\Backup\ASL\BioFinder\Prisma_3D_pCASL';
% Flist   = xASL_adm_GetFileList(mDir,'t1_mprage', 'FPListRec');
% for iL=1:length(Flist)
%     [Fpath Ffile Fext]  = xASL_fileparts(Flist{iL});
%     NewName             = fullfile(Fpath,'T1.nii.gz');
%     xASL_Move(Flist{iL},NewName);
% end
%
% Flist   = xASL_adm_GetFileList(mDir,'t2_spc_da', 'FPListRec');
% for iL=1:length(Flist)
%     [Fpath Ffile Fext]  = xASL_fileparts(Flist{iL});
%     NewName             = fullfile(Fpath,'FLAIR.nii.gz');
%     xASL_Move(Flist{iL},NewName);
% end
%
% Flist   = xASL_adm_GetFileList(mDir,'tgse_pcasl', 'FPListRec');
% for iL=1:length(Flist)
%     [Fpath Ffile Fext]  = xASL_fileparts(Flist{iL});
%     ASLdir              = fullfile(Fpath,'ASL_1');
%     xASL_adm_CreateDir(ASLdir);
%     NewName             = fullfile(ASLdir,'ASL4D.nii.gz');
%     xASL_Move(Flist{iL},NewName);
% end
%
% Flist   = xASL_adm_GetFileList(mDir,'ASL4D\.(nii|nii\.gz)', 'FPListRec');
% for iL=1:length(Flist)
%     xASL_io_SplitASL_M0(Flist{iL},1);
% end

% Flist   = xASL_adm_GetFileList(mDir,'M0\.(nii|nii\.gz)', 'FPListRec');
% for iL=1:length(Flist)
%     clear parms
%     [Fpath Ffile Fext]      = xASL_fileparts(Flist{iL});
%     ParmsMat                = fullfile(Fpath,[Ffile '_parms.mat']);
%     parms.RepetitionTime    = 4600;
%     save(ParmsMat,'parms');
% end


%% Skyra FAIR/2D PCASL
mDir    = 'C:\Backup\ASL\BioFinder\Skyra_FAIR_and_2D_pCASL';

Flist   = xASL_adm_GetFileList(mDir,'t1_mprage', 'FPListRec');
for iL=1:length(Flist)
    [Fpath Ffile Fext]  = xASL_fileparts(Flist{iL});
    NewName             = fullfile(Fpath,'T1.nii.gz');
    xASL_Move(Flist{iL},NewName);
end

Flist   = xASL_adm_GetFileList(mDir,'ASL_3D_tra_iso', 'FPListRec');
for iL=1:length(Flist)
    [Fpath Ffile Fext]  = xASL_fileparts(Flist{iL});
    ASLdir              = fullfile(Fpath,'ASL_1');
    xASL_adm_CreateDir(ASLdir);
    NewName             = fullfile(ASLdir,'ASL4D.nii.gz');
    xASL_Move(Flist{iL},NewName);
end

Flist   = xASL_adm_GetFileList(mDir,'ep2d_.*asl_.*ASL', 'FPListRec');
for iL=1:length(Flist)
    [Fpath Ffile Fext]  = xASL_fileparts(Flist{iL});
    ASLdir              = fullfile(Fpath,'ASL_2');
    xASL_adm_CreateDir(ASLdir);
    NewName             = fullfile(ASLdir,'ASL4D.nii.gz');
    xASL_Move(Flist{iL},NewName);
end

Flist   = xASL_adm_GetFileList(mDir,'ASL4D\.(nii|nii\.gz)', 'FPListRec');
for iL=1:length(Flist)
    xASL_TrackProgress(iL,length(Flist));
    [Fpath Ffile Fext]  = xASL_fileparts(Flist{iL});
    ParmsMat            = fullfile(Fpath,'ASL4D_parms.mat');

    if     ~isempty(strfind(Flist{iL},'ASL_1')) % 3D GRASE
            clear parms
            parms.M0                = 1050;
            parms.readout_dim       = '3D';
            parms.EchoTime          = 20; % Dummy parameter!!!!!!!!!!!!!!!!!!!!!!!
            parms.BackgroundSuppressionNumberPulses = 2;
            parms.RepetitionTime    = 5000;
            save(ParmsMat, 'parms');
    elseif ~isempty(strfind(Flist{iL},'ASL_2')) % 2D EPI
            clear parms
            parms.M0                = 'no_background_suppression';
            parms.readout_dim       = '2D';
            parms.RepetitionTime    = 5400;
            parms.SliceReadoutTime  = 30;
            save(ParmsMat, 'parms');

%             xASL_io_SplitASL_M0(Flist{iL},1);
    else
        error('Unknown session');
    end
end


%% Trio
mDir    = 'C:\Backup\ASL\BioFinder\Trio';

Flist   = xASL_adm_GetFileList(mDir,'t1_MP-RAGE', 'FPListRec');
for iL=1:length(Flist)
    [Fpath Ffile Fext]  = xASL_fileparts(Flist{iL});
    NewName             = fullfile(Fpath,'T1.nii.gz');
    xASL_Move(Flist{iL},NewName);
end

Dlist   = xASL_adm_GetFsList(mDir,'\d{3}',1);
for iD=1:length(Dlist)
    cDir        = fullfile(mDir,Dlist{iD});
    FlistSS     = xASL_adm_GetFileList(cDir,'ss');
    FlistNS     = xASL_adm_GetFileList(cDir,'ns');

    if  length(FlistSS)>0 && length(FlistNS)>0
        clear IM IM1 IM2
        IM1     = xASL_io_Nifti2Im(FlistSS{1});
        IM2     = xASL_io_Nifti2Im(FlistNS{1});

        % Check
        if ~min(size(IM1)==size(IM2))
            error('NS & SS nifti files weren''t same dimensions');
        end

        Dim4                    = size(IM1,4)*2;
        IM(:,:,:,1:2:Dim4-1)    = IM1;
        IM(:,:,:,2:2:Dim4-0)    = IM2;

        ASLdir          = fullfile(cDir,'ASL_1');
        xASL_adm_CreateDir(ASLdir);
        ASLpath         = fullfile(ASLdir,'ASL4D.nii');
        xASL_io_SaveNifti(FlistSS{1},ASLpath,IM,[],0);

    end
end

Flist   = xASL_adm_GetFileList(mDir,'ASL4D\.(nii|nii\.gz)', 'FPListRec');
for iL=1:length(Flist)
    xASL_TrackProgress(iL,length(Flist));
    clear parms
    [Fpath Ffile Fext]      = xASL_fileparts(Flist{iL});
    ParmsMat                = fullfile(Fpath,'ASL4D_parms.mat');
    parms.RepetitionTime    = 4200;
    parms.SliceReadoutTime  = 30;
    save(ParmsMat, 'parms');
    ParmsMat                = fullfile(Fpath,'M0_parms.mat');
    save(ParmsMat, 'parms');
end
