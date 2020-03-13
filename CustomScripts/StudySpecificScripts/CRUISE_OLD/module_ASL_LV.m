function result = module_ASL(job, x)
% $Rev:: 475                   $:  Revision of last commit
% $Author:: hjmutsaerts        $:  Author of last commit
% $Date:: 2017-03-21 13:12:48 #$:  Date of last commit

%% ASL module of ExploreASL Paul, Matthan, Henk-Jan


%% 0. INPUT

[path   x.P.SessionID ext]     = fileparts(x.SESSIONDIR);
[dummy x.P.SubjectID ext]      = fileparts(path);

xASL_adm_CreateDir(x.D.ASLCheckDir);
xASL_adm_CreateDir(x.MotionDir);   xASL_adm_CreateDir(x.ExclusionDir);
xASL_adm_CreateDir(x.SNRdir);      xASL_adm_CreateDir(x.RawEPIdir);

if  strcmp(x.M0,'separate_scan')
    xASL_adm_CreateDir(x.D.M0CheckDir);  xASL_adm_CreateDir(x.D.M0regASLdir);
end
if  strcmp(x.M0,'no_background_suppression')
    xASL_adm_CreateDir(x.D.RawDir);
end

x     = GenericMutexModules( x, 'ASL' ); % starts mutex locking process to ensure that everything will run only once
result      = false;

%% 0.9 change working directory to make sure that unspecified output will go there...
oldFolder = cd(x.SESSIONDIR);

UnZipnativeNii4Processing(x.SUBJECTDIR); % in case the script searches for *.nii

%% Check file existence
% get input path
asl4D_raw_nii = xASL_adm_GetFsList(x.SESSIONDIR, x.ASL4DFILE);
nFiles = length(asl4D_raw_nii);
if nFiles==0
    error('AslPipeline:missingFile', 'No file found that matches "%s" in "%s"', x.ASL4DFILE, x.SESSIONDIR);
elseif nFiles>1
    error('AslPipeline:incorrectNrOfFiles', 'Expected exactly one nifti file "%s", but found %d', x.ASL4DFILE, nFiles);
end
x.P.ASL4D = fullfile(x.SESSIONDIR, asl4D_raw_nii{1}); % remove cell array from this single element

if ~isfield(x,'motion_correction')
    x.motion_correction   = 0;  % for turbo-quasar cbf maps
    disp('Motion_correction will not be applied');   %THIS IS A CHECK TO MAKE SURE I CHANGE THIS BACK TO 1
end

%% Check existence structural reference files,
%% of which this module is dependent
pGM         = fullfile(x.SUBJECTDIR,['c1' x.P.STRUCT '.nii']);
pWM         = fullfile(x.SUBJECTDIR,['c2' x.P.STRUCT '.nii']);
FlowField   = fullfile(x.SUBJECTDIR,['y_' x.P.STRUCT '.nii']);

if  ~xASL_exist(pGM,'file') || ~xASL_exist(pWM,'file') || ~xASL_exist(FlowField,'file')
    error('Module skipped because structural files missing, please rerun structural module');
end


%% 1    Motion correction ASL

if ~x.motion_correction
    fprintf('%s\n','Motion correction was disabled, skipping');
else


    if ~x.mutex.HasState('002_realign_ASL')

        % First, solve dimensionality
        TempNii     = xASL_io_ReadNifti( x.P.ASL4D);
        TempNii.dat = squeeze(TempNii.dat);
        create(TempNii);

        % First, check matrix size: throw error if 2D data with 3 dimensions only

        if  size(TempNii.dat,4)>1 && size(TempNii.dat,4)/2~=round(size(TempNii.dat,4)/2)
            error('Uneven number of control-label frames, either incomplete pairs or M0 image in time-series!');
        end

        % 1    Estimate motion from mean head position using SPM_realign_asl (saved in matrix)
        % 2    Calculate and plot position and motion parameters

        if  size(TempNii.dat,4)>2

            % First create mean PWI, tSD and tSNR image before MoCo, for
            % comparison

            fprintf('%s\n','Computing mean, SD and SNR ASL images before motion correction');
            fprintf('%s\n','This enables to investigate the effect of motion correction');
            [Co La]         = Check_control_label_order(TempNii.dat(:,:,:,:));
            PWI             = Co-La;

            IM{1}           = mean(Co,4);
            IM{2}           = std(Co,[],4);
            IM{3}           = IM{1}./IM{2};

            IM{4}           = mean(PWI,4);
            IM{5}           = std(PWI,[],4);
            IM{6}           = IM{4}./IM{5};

            Fname           = {'mean_control_beforeMoCo' 'SD_control_beforeMoCo' 'SNR_control_beforeMoCo'};
            Fname(4:6)      = {'mean_PWI_beforeMoCo' 'SD_PWI_beforeMoCo' 'SNR_PWI_beforeMoCo'};

            for ii=1:length(Fname)
                NativePath{ii,1}  = fullfile(x.SESSIONDIR,[Fname{ii} '.nii']);
                DARTELpath{ii,1}  = fullfile(x.D.PopDir, [Fname{ii} '_' x.P.SubjectID '_' x.P.SessionID '.nii']);
                NativeMat{ii,1}   = fullfile(x.SESSIONDIR,[Fname{ii} '.mat']);
                xASL_io_SaveNifti(x.P.ASL4D,NativePath{ii},IM{ii},16,0);
                if  xASL_exist(NativeMat{ii},'file');xASL_delete(NativeMat{ii});end
            end

            xASL_spm_deformations(x,x.P.SubjectID,NativePath,DARTELpath);

            for ii=1:length(NativePath)
                if  xASL_exist(NativePath{ii},'file')
                    xASL_delete(NativePath{ii});
                end
            end

        else fprintf('%s\n',['Skipping creation pre-MoCo QC images for ' x.P.SubjectID '_' x.P.SessionID ' because it had only ' num2str(size(TempNii.dat,4)) ' 3D frames.']);

        end
        if  size(TempNii.dat,4)>1
            % Run Motion Correction
            realign_ASL(x);
        else
            fprintf('%s\n',['Skipping motion correction for ' x.P.SubjectID '_' x.P.SessionID ' because it had only ' num2str(size(TempNii.dat,4)) ' 3D frames.']);
        end

        x.mutex.AddState('002_realign_ASL');
        x.mutex.DelState('0025_register_ASL');
    else    fprintf('%s\n','002_realign_ASL session has already been performed, skipping...');
    end

end

% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
[path file ext]             = fileparts(x.P.ASL4D);
x.despiked_raw_asl    = fullfile(path,['despiked_' file '.nii']);
despiked_raw_asl2           = fullfile(path,['rtemp_despiked_' file '.nii']);
if ~xASL_exist( x.despiked_raw_asl ,'file') && ~xASL_exist( despiked_raw_asl2 ,'file')
    x.despiked_raw_asl = x.P.ASL4D;
end

%% 2    Registration ASL sessions to T1w

if ~x.mutex.HasState('0025_register_ASL')

    register_ASL(x);

    x.mutex.AddState('0025_register_ASL');
    x.mutex.DelState('003_reslice_ASL');
%     x.mutex.DelState('004_realign_reslice_M0'); % also dependent on ASL registration
else    fprintf('%s\n','0025_register_ASL       has already been performed, skipping...');
end



%% 3    Reslice ASL data to high resolution standard space
if ~x.mutex.HasState('003_reslice_ASL') && x.mutex.HasState('0025_register_ASL')

    % 2    Create slice gradient image for quantification reference, in case of 2D ASL
    % 3    Reslice ASL time series to MNI space (currently 1.5 mm^3)
    % 4    Create mean control image, masking to 20% of max value if used as M0 (no background suppression)
    % 5    Smart smoothing mean_control if used as M0

    reslice_ASL(x);

    x.mutex.AddState('003_reslice_ASL');
%     x.mutex.DelState('004_realign_reslice_M0'); % not required to
%     repeat if already performed
%     x.mutex.DelState('0035_PreparePV');
    x.mutex.DelState('005_quantification');
        x.mutex.DelState('006_visualize');
else    fprintf('%s\n','003_reslice_average_ASL       has already been performed, skipping...');
end


%% 4    Resolution estimation & prepare partial volume maps

pGM_PV      = fullfile(x.D.PopDir,['PV_pGM_' x.P.SubjectID '.nii']);
pWM_PV      = fullfile(x.D.PopDir,['PV_pWM_' x.P.SubjectID '.nii']);

%% We only have to do this once per subject, we do not anticipate resolution differences between ASL sessions
if  ~(xASL_exist(pGM_PV,'file') && xASL_exist(pWM_PV,'file')) || strcmp(x.P.SessionID,x.SESSIONS{1})
    % Or if this is the first session, redo this (if ExploreASL is rerun)
    if ~x.mutex.HasState('0035_PreparePV') && x.mutex.HasState('003_reslice_ASL')

        % 1 Upsample ASL image to 1.5x1.5x1.5 mm
        % 2 Reslice pGM & pWM native images to this orientation (will rotate them
        % into ASL space)
        % 3 Estimate ASL sequence smoothness
        % 4 Apply same smoothness to tissue priors
        % 5 Transform tissue priors to MNI space
        PreparePV_ASL(x);

        x.mutex.AddState('0035_PreparePV');

    else    fprintf('%s\n','0035_PreparePV       has already been performed, skipping...');
    end
else fprintf('%s\n','pGM_PV & pWM_PV already existed, skipping...');

end

%% 5    Process M0
if ~x.mutex.HasState('004_realign_reslice_M0') && x.mutex.HasState('0025_register_ASL')

%     Any M0 will be processed here. Even if part of the subjects doesn't
%     have an M0, since this can be later imputed, or an average population
%     M0 image could be used. Also, without background suppression and
%     without an M0, the MeanControl image is before saved as M0, and will
%     be processed here as well.
        M0_raw_nii = xASL_adm_GetFileList( x.SESSIONDIR, ['^' x.P.M0 '\.(nii|nii\.gz)$'], 'FPList', [0 Inf]);
        if ~isempty(M0_raw_nii)

            % 1)    Motion correction if there are multiple frames
            % 2)    Registration M0 -> mean_control_ASL image
            % 3)    Smart smoothing
            % 4)    Reslice M0 to MNI space (currently 1.5 mm^3)
            % 5)    Averaging if multiple frames
            % 6)    Masking
            % 7)    Correction for scale slopes & incomplete T1 recovery

            realign_reslice_M0(x);

            x.mutex.AddState('004_realign_reslice_M0');
            x.mutex.DelState('005_quantification');

        else    fprintf('%s\n','M0 post-processing skipped, because no M0 available');
        end

else    fprintf('%s\n','004_realign_reslice_M0                 has already been performed, skipping...');
end


%% 6    Quantification
% Quantification is performed here according to ASL consensus paper (Alsop, MRM 2016)
if ~x.mutex.HasState('005_quantification') && x.mutex.HasState('003_reslice_ASL')

    % 1    Prepare M0 image
    % 2    Prepare PWI image
    % 3    Load slice gradient if 2D
    % 4    CBF quantification
    % 5    Division by M0 & scale slopes
    % 6    Vendor-specific quantification correction
    % 7    Repair negative artifacts (perfusion fluctuation)
    % 8    Compression high macro-vascular signal
    % 9    Save files
    % 10   Computation FEAST transit times
    rtemp_ASL4D = xASL_adm_GetFileList( x.SESSIONDIR, ['^rtemp_' x.P.ASL4D '\.nii$'], 'FPList', [0 Inf]);

    % this is to allow pre-prepared rtemp_ASL4D.nii files from the 140
    % dynamics 4D timeseries pipeline that were subsequently split into
    % pre and post because those had better filtering than the 35 dynamics
    if ~isempty(rtemp_ASL4D);
        InputPath = rtemp_ASL4D;
    else
        InputPath = fullfile(x.D.PopDir, ['PWI_' x.P.SubjectID '_' x.P.SessionID '.nii']);
    end

    fprintf('%s\n','Quantifying ASL');
    average_quantify_ASL( x, InputPath);

    x.mutex.AddState('005_quantification');
    x.mutex.DelState('006_visualize');
else    fprintf('%s\n','005_quantification                     has already been performed, skipping...');
end


%% 7    Visual check
if ~x.mutex.HasState('006_visualize') && x.mutex.HasState('005_quantification')

    fprintf('%s\n','print visual quality assurance checks');

    close all % close all Figures to avoid capturing & saving the wrong Figure

    % print visual quality assurance checks
    % PWI (scaled for GM)

    xASL_im_CreateVisualFig(       x.D.PopDir, ['^qCBF_' x.P.SubjectID '_' x.P.SessionID '\.(nii|nii\.gz)$'], x.D.ASLCheckDir      , ['qCBF_' x.P.SubjectID '_' x.P.SessionID], x.CLIP_ZERO ); % zero clipping
    xASL_im_CreateVisualFig(       x.D.PopDir, ['^qCBF_untreated_' x.P.SubjectID '_' x.P.SessionID '\.(nii|nii\.gz)$'], x.D.ASLCheckDir      , ['qCBF_untreated_' x.P.SubjectID '_' x.P.SessionID], x.CLIP_ZERO ); % zero clipping

    % Check registration ASL with T1
    visual_registration_check_MNI_1_5(x.D.PopDir, ['^PWI_' x.P.SubjectID '_' x.P.SessionID '\.(nii|nii\.gz)$'] ,['^rc2' x.P.STRUCT '_' x.P.SubjectID '\.(nii|nii\.gz)$'], x.D.T1_ASLREGDIR,0.5,1);

    % Check PV_pGM & PV_pWM
    xASL_im_CreateVisualFig(       x.D.PopDir, ['^PV_pGM_' x.P.SubjectID '\.(nii|nii\.gz)$'], x.D.ASLCheckDir      , ['PV_pGM_' x.P.SubjectID], x.CLIP_ZERO ); % zero clipping

    TempNii     = xASL_io_ReadNifti( x.P.ASL4D );

    if  size(TempNii.dat,4)>1 % if there are individual frames
        % Print large overview of all raw frames & SD
%         visual_nii_check_frames( x.P.ASL4D, x.RawEPIdir, [x.P.SubjectID '_' x.P.SessionID],2,1);
%         visual_nii_check_SD_timeseries( x.P.ASL4D, x.RawEPIdir, [x.P.SubjectID '_' x.P.SessionID],0);

        % mean control
        xASL_im_CreateVisualFig(       x.D.PopDir, ['^mean_control_' x.P.SubjectID '_' x.P.SessionID '\.(nii|nii\.gz)$'],x.D.RawDir      , ['mean_control_' x.P.SubjectID '_' x.P.SessionID], x.CLIP_ZERO ); % zero clipping
    end

    if  size(TempNii.dat,4)>2 % if there are more than two frames
        % Print overview of SD of paired subtractions
%         visual_nii_check_SD_timeseries( x.P.ASL4D, x.RawEPIdir, ['subtr_' x.P.SubjectID '_' x.P.SessionID],1);

        % SD & SNR maps
        xASL_im_CreateVisualFig(       x.D.PopDir, ['^SD_'  x.P.SubjectID '_' x.P.SessionID '\.(nii|nii\.gz)$'], x.SNRdir      , ['SD_map_'  x.P.SubjectID '_' x.P.SessionID], x.CLIP_ZERO );
        xASL_im_CreateVisualFig(       x.D.PopDir, ['^SNR_' x.P.SubjectID '_' x.P.SessionID '\.(nii|nii\.gz)$'], x.SNRdir      , ['SNR_map_' x.P.SubjectID '_' x.P.SessionID], x.CLIP_ZERO );
    end

        % slice gradient
        xASL_adm_CreateDir(x.D.SliceCheckDir);
        visual_registration_check_MNI_1_5_cor_sag(x.D.PopDir, ['^PWI_' x.P.SubjectID '_' x.P.SessionID '\.(nii|nii\.gz)$'],x.D.PopDir, ['^slice_gradient_' x.P.SubjectID '_' x.P.SessionID '\.(nii|nii\.gz)$'],x.D.SliceCheckDir,['slice_gradient_' x.P.SubjectID '_' x.P.SessionID],x.jet256);

    % Check registration M0 with ASL, if instead separate M0 scan (if mean_control was used in case of 'no background suppression', than already perfect registration)
    if  strcmp(x.M0, 'separate_scan')
        xASL_im_CreateVisualFig( x.D.PopDir, ['^' x.P.M0 '_' x.P.SubjectID '_' x.P.SessionID '\.nii'], x.D.M0CheckDir, [x.P.M0 '_' x.P.SubjectID '_' x.P.SessionID] , 1); % zero clipping
        xASL_im_CreateVisualFig( x.D.PopDir, ['^' x.P.M0 '_' x.P.SubjectID '_' x.P.SessionID '\.nii'], x.D.M0CheckDir, [x.P.M0 '_' x.P.SubjectID '_' x.P.SessionID] , 1,'x.jet256' ); % zero clipping
    end

    if  x.nSessions>1
        if strcmp(x.session.options{1},'crushed') && strcmp(x.session.options{2},'non-crushed')
            if  strcmp(x.SESSIONDIR(length(x.SUBJECTDIR)+2:end),'ASL_2')
                xASL_im_CreateVisualFig(   x.D.PopDir, ['^TT_' x.P.SubjectID '\.(nii|nii\.gz)$'], x.D.TTCheckDir, ['TT_map_'  x.P.SubjectID], x.CLIP_ZERO );
                xASL_im_CreateVisualFig(   x.D.PopDir, ['^TT_' x.P.SubjectID '\.(nii|nii\.gz)$'], x.D.TTCheckDir, ['TT_map_'  x.P.SubjectID], x.CLIP_ZERO,x.jet256 );
                visual_registration_check_MNI_1_5(x.D.PopDir, ['^TT_' x.P.SubjectID  '\.(nii|nii\.gz)$'] ,['^' 'rc2' x.P.STRUCT '_' x.P.SubjectID '\.nii'], x.D.TTCheckDir,1,1,[],x);
            end
        end
    end
        x.mutex.AddState('006_visualize');
else    fprintf('%s\n','006_visualize                          has already been performed, skipping...');
end

%% 999 Ready
x.mutex.AddState('999_ready');
cd(oldFolder);

x.mutex.Unlock();
x.result  = true;
result          = true;

end
