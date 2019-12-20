%% DARTEL parallel

% Later implement this as normal DARTEL module, running it parallel simply
% by selecting different subjects for each Matlab run

% KernelSize                      = [8;6.25;4.5;2.75;1];
% No smoothing required, images are almost good enough

% Combine pGM and pWM templates into single 4D NIfTI

Nrs         = [0 1 2 3 4 5 6];
rparam1     = [2   1    0.5   0.25   0.25  0.125];
rparam2     = [1   0.5  0.25  0.125  0.125 0.0625];
K           = [0   2    4     6      8     16];
DummyVar    = '';


LockDir     = fullfile(x.D.ROOT,'lock','DARTEL_T1','DARTEL_module');

%% From SubjectList select only the first time points
clear ToWarpList
nextN   = 1;


for iS=1:x.nSubjects
    x.P.SubjectID = x.SUBJECTS{iS};
    [~, ~, ~, IsSubject] = xASL_init_LongitudinalRegistration( x );
    % To check whether or not we will run longitudinal registration

    if  IsSubject==1 % only run for first volumes
        ToWarpList{nextN}      = x.SUBJECTS{iS};
        nextN                  = nextN+1;
    end
end

%% DivideIn4
% length(ToWarpList)/4 = 25

ToWarpList  = ToWarpList([1:2]);  % PARALLEL 1
% ToWarpList  = ToWarpList(26:50); % PARALLEL 2
% ToWarpList  = ToWarpList(51:75); % PARALLEL 3
% ToWarpList  = ToWarpList(76:end);% PARALLEL 4

% for ii=1:6
ii=1

    LockFileSave    = fullfile(LockDir,['001_create_templates_' num2str(ii) '.status']);

%     if ~exist(LockFileSave,'file')

        %% Save average image
        List{1} = xASL_adm_GetFileList(x.D.PopDir, '^rc1T1.*\.(nii|nii\.gz)$', 'FPList', [0 Inf]);
        List{2} = xASL_adm_GetFileList(x.D.PopDir, '^rc2T1.*\.(nii|nii\.gz)$', 'FPList', [0 Inf]);
        clear IM
        for iT=1:2
            for iL=1:length(List{iT})
                IM(:,:,:,iL) = xASL_io_Nifti2Im(List{iT}{iL});
            end
            MeanIM(:,:,:,iT) = xASL_stat_MeanNan(IM, 4);
            MadiM(:,:,:,iT) = xASL_stat_MadNan(IM, [], 4);
            MeanIM(:,:,:,iT) = (MeanIM(:,:,:,iT) + MadiM(:,:,:,iT))./2;

            SmoothAv = [3 2 1 0 0 0 0 0];
            if SmoothAv(ii)~=0
                MeanIM(:,:,:,iT) = xASL_im_ndnanfilter(MeanIM(:,:,:,iT), 'gauss', [SmoothAv SmoothAv SmoothAv]);
            end
        end
        TemplateFile= fullfile(x.D.PopDir,['T1_template_' num2str(Nrs(ii)) '.nii']);
        xASL_io_SaveNifti(x.D.ResliceRef,TemplateFile,MeanIM, [], 0);
        clear IM


        %% Create flowfields
        LockDirII   = fullfile(LockDir,['Run_' num2str(ii)]);
        xASL_adm_CreateDir(LockDirII);

        for  iS=1:length(ToWarpList)
            % Create flowfield
            IMpath1     = fullfile(x.D.PopDir,['rc1T1_' ToWarpList{iS} '.nii']);
            IMpath2     = fullfile(x.D.PopDir,['rc2T1_' ToWarpList{iS} '.nii']);
            LockFile    = fullfile(LockDirII,[ToWarpList{iS} '.status']);
            u_file1     = fullfile(x.D.PopDir,['u_rc1T1_' ToWarpList{iS} '.nii']);
            u_file2     = ['u_rc1T1_' ToWarpList{iS} '_T1_template.nii'];
            u_Path2     = fullfile(x.D.PopDir, u_file2);

            if  exist(IMpath1,'file') && exist(IMpath2,'file') && ~exist(LockFile,'file')
                DARTELjobWrap(IMpath1,IMpath2,TemplateFile,rparam1(ii),rparam2(ii),K(ii));

                xASL_Rename(u_file1,u_file2,1);

                INPUTname1      = fullfile(x.D.ROOT, ToWarpList{iS}, 'c1T1.nii');
                INPUTname2      = fullfile(x.D.ROOT, ToWarpList{iS}, 'c2T1.nii');

                % Create new image
                xASL_spm_deformations(x,{IMpath1;IMpath2},{IMpath1;IMpath2}, 2, [], [], u_Path2);
                save(LockFile,DummyVar);
            end
        end
        toc
        save(LockFileSave,DummyVar);
%     end



% end
