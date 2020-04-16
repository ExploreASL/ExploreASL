function [result, x] = xASL_module_dwi(x)
% $Rev:: 475                   $:  Revision of last commit
% $Author:: hjmutsaerts        $:  Author of last commit
% $Date:: 2017-03-21 13:12:48 #$:  Date of last commit

%% DWI module of ExploreASL Paul, Matthan, Henk-Jan
% scripts taken from Chris Vriends script (run_topup_DWI.sh)
% Some extra options for Eddy current correction:
% --repol does some additional correction of motion-related artifacts
% --verbose to see it running

%% 0. INPUT

x = xASL_init_GenericMutexModules( x, 'dwi' ); % starts mutex locking process to ensure that everything will run only once
result = false;

%% 0.9 change working directory to make sure that unspecified output will go there...
oldFolder = cd(x.SESSIONDIR);

Fpath = fileparts(x.SESSIONDIR);
[~, Ffile] = fileparts(Fpath);
iSubject = find(strcmp(x.SUBJECTS,Ffile));

if ~isfield(x,'DoWADQCDC')
    x.DoWADQCDC = false;
end

%% Check whether we already ran this (to know whether we should give verbose feedback)
if ~x.mutex.HasState('999_ready')
    bO = true; % generate output, some processing has and some has not been yet done
else
    bO = false; % skip output, as all processing has been performed
end

x.Sequence      = '2D_EPI';
x.readout_dim   = '2D';


%% Skip if existing
if x.mutex.HasState('010_TopUp_dwi') && x.mutex.HasState('030_RegistrationDWI2T1w')
    if x.mutex.HasState('020_EddyCurrent') && x.mutex.HasState('040_dtiFit')
        if x.mutex.HasState('050_resliceDWI') && x.mutex.HasState('060_visualize')
            if x.mutex.HasState('999_ready')
                result = true;
                return;
            end
        end
    end
end

%% Define paths
FileDWI = xASL_adm_GetFileList(x.SESSIONDIR, '^dwi(|_run-\d)(|_dwi)(?!.*ADC)\.nii$', 'List', [0 Inf]);
if isempty(FileDWI)
    if bO; fprintf('%s\n','No DWI found, skipping'); end
    result = true;
    return; % skip if no DWI data found
elseif length(FileDWI)>1
    warning('More DWI files found than 1');
end
[~, FileDWI] = xASL_fileparts(FileDWI{1});

x = xASL_init_FileSystem(x);

PathX = fullfile(x.SUBJECTDIR, 'x.mat');

x.P.Path_dwi = fullfile(x.SESSIONDIR, [FileDWI '.nii']);
x.P.Path_dwi_ORI = fullfile(x.SESSIONDIR, [FileDWI '_ORI.nii']);
x.P.Path_dwiMat = fullfile(x.SESSIONDIR, [FileDWI '.mat']);
x.P.Path_dwi_mean = fullfile(x.SESSIONDIR, [FileDWI '_mean.nii']);
x.P.Path_dwi_mask = fullfile(x.SESSIONDIR, [FileDWI '_mean_mask.nii']);
% x.P.Path_wdwi = fullfile(x.SESSIONDIR, ['w' FileDWI '.nii']);
% x.P.Path_wdwiMat = fullfile(x.SESSIONDIR, ['w' FileDWI '.mat']);

% Derivatives
PathTopUp = 'TopUp'; % base name of TopUp output files (spline coefficient & movparms)
PathEddy = fullfile(x.SESSIONDIR, 'Eddy'); % base name for output
PathEddyNii = fullfile(x.SESSIONDIR, 'Eddy.nii'); % NIfTI name
PathEddyMat = fullfile(x.SESSIONDIR ,'Eddy.mat');
PathWithoutOutlierRemoval = fullfile(x.SESSIONDIR, 'Eddy.eddy_outlier_free_data.nii');
PathFit = fullfile(x.SESSIONDIR, 'DTIfit'); % base name
PathAcqp = fullfile(x.SESSIONDIR, 'TopUpAcqParms_dwi.txt');
PathBvec = fullfile(x.SESSIONDIR, [FileDWI '.bvec']);
PathBval = fullfile(x.SESSIONDIR, [FileDWI '.bval']);
IndexPath = fullfile(x.SESSIONDIR ,'Index.txt');
PathB0 = fullfile(x.SESSIONDIR ,'B0.nii');
PathB0Mat = fullfile(x.SESSIONDIR ,'B0.mat');
PathUnwarped = fullfile(x.SESSIONDIR ,'Unwarped.nii');
PathUnwarpedMat = fullfile(x.SESSIONDIR ,'Unwarped.mat');
PathFA = fullfile(x.SESSIONDIR ,'DTIfit_FA.nii');
PathSSE = fullfile(x.SESSIONDIR ,'DTIfit_sse.nii');
PathField = fullfile(x.SESSIONDIR ,'Field.nii');

PathPopB0 = fullfile(x.D.PopDir, ['rdwi_B0_' x.SUBJECTS{iSubject} '.nii']);
PathPopUnwarped = fullfile(x.D.PopDir, ['rdwi_Unwarped_' x.SUBJECTS{iSubject} '.nii']);
PathPopFA = fullfile(x.D.PopDir, ['rdwi_FA_' x.SUBJECTS{iSubject} '.nii']);
PathPopSSE = fullfile(x.D.PopDir, ['rdwi_SSE_' x.SUBJECTS{iSubject} '.nii']);
PathPopMean = fullfile(x.D.PopDir, ['rdwi_mean_' x.SUBJECTS{iSubject} '.nii']);


%% Check existence structural reference files, of which this module is dependent
if  (~xASL_exist(x.P.Path_c1T1,'file') || ~xASL_exist(x.P.Path_c2T1,'file') || ~xASL_exist(x.P.Path_y_T1,'file')) && ~x.mutex.HasState('999_ready')
    xASL_adm_UnzipNativeFiles4Processing(x.SUBJECTDIR, x); % in case the script searches for *.nii
    warning('Structural files missing');
end

%% Delete derivatives from a previous run
if ~x.mutex.HasState('030_RegistrationDWI2T1w') && ~x.mutex.HasState('020_EddyCurrent')
    xASL_adm_DeleteFileList(x.SESSIONDIR, '^(B0|dwi.*mean|dwi.*mask|Eddy|wdwi|DTIfit|Field|Index\.txt|Unwarped|y_DWI|TopUp|xASL_qc|qcdc).*', false, [0 Inf]);
    
    [Fpath1, Ffile1] = xASL_fileparts(x.P.Path_dwi_ORI);
    [Fpath2, Ffile2] = xASL_fileparts(x.P.Path_dwi);
    ExtOri = {'.nii' '.json' '.mat' '_sn.mat' '_parms.mat'};
    for iExt=1:length(ExtOri)
        PathOrig = fullfile(Fpath1, [Ffile1 ExtOri{iExt}]);
        PathDest = fullfile(Fpath2, [Ffile2 ExtOri{iExt}]);
        if xASL_exist(PathOrig,'file')
            xASL_Move(PathOrig, PathDest, true);
        end
    end    
end


%% -----------------------------------------------------------------------------
%% 1    TopUp
if ~x.mutex.HasState('010_TopUp_dwi') || ~xASL_exist(PathAcqp,'file') || ~xASL_exist(PathField,'file')
    xASL_adm_DeleteFileList(x.SESSIONDIR, '^(B0|Field|Unwarped|TopUp).*', false, [0 Inf]);
    bSuccess = xASL_fsl_TopUp(x.SESSIONDIR, 'dwi', x);

    if bSuccess
        x.mutex.AddState('010_TopUp_dwi');
        x.mutex.DelState('020_EddyCurrent');
    else
        warning('TopUp didnt run correctly, rest of the DTI results may not be trusted');
    end
elseif   bO; fprintf('%s\n','010_TopUp_dwi has already been performed, skipping...');
end


%% -----------------------------------------------------------------------------
%% 2    Eddy current (distortion correct, motion correction etc)
if ~x.mutex.HasState('020_EddyCurrent') || ~xASL_exist(PathEddyNii,'file')

    xASL_adm_DeleteFileList(x.SESSIONDIR, '^Eddy.*', false, [0 Inf]);

    % Create index (this is a reference from the image volumes to acqparameters & TopUp,
    %  with respect to the phase encoding (PE) direction
    %  We assume here that all image volumes have the same PE direction
    Dim4 = size(xASL_io_Nifti2Im(x.P.Path_dwi),4);

    fclose all;
    FID = fopen(IndexPath,'wt');

    for iD=1:Dim4
        fprintf(FID, '1 ');
    end
    fclose(FID);

    % Determine CUDA or openMP (CUDA = multi-core GP, openMP = multi-core CPU)
    fprintf('Readying for eddy current-correction...\n');
    [FSLdir, x, RootWSLdir] = xASL_fsl_SetFSLdir(x, x.bAutomaticallyDetectFSL);

    % Compute brain mask
    fprintf('Computing mean DWI image...\n');
    xASL_fsl_RunFSL(['/bin/fslmaths ' xASL_adm_UnixPath(x.P.Path_dwi) ' -Tmean '...
        xASL_adm_UnixPath(x.P.Path_dwi_mean)], x);
    fprintf('Running Bet for brainmasking...\n');
    xASL_fsl_RunFSL(['/bin/bet ' xASL_adm_UnixPath(x.P.Path_dwi_mean) ' '...
        xASL_adm_UnixPath(x.P.Path_dwi_mean) ' -R -m'], x);    
    
    % skip determining this for now, simply use CPU version, GPU not tested
    % yet
%     CudaLocations{1} = 'ls -d /usr/local/cuda-9.1*;';
%     CudaLocations{2} = 'ls -d /usr/local/cuda-8.0*;';
%     CudaLocations{3} = 'ls -d /usr/local/cuda*;';
%     CudaLocations{4} = ['ls -f ' FSLdir '/bin/eddy_cuda9.1;'];
%     CudaLocations{5} = ['ls -f ' FSLdir '/bin/eddy_cuda8.0;'];
%     CudaLocations{6} = ['ls -f ' FSLdir '/bin/eddy_cuda;'];
%
%     for iCuda=1:length(CudaLocations)
%         if ispc
%             CudaLocations{iCuda} = ['wsl ' CudaLocations{iCuda}];
%         end
%         StatusCheck(iCuda) = system(CudaLocations{iCuda});
%     end
%
%     CudaVersionUnknown = StatusCheck(1)~=0 && StatusCheck(2)~=0 && StatusCheck(3)==0;
%
%     if (StatusCheck(1)==0 || CudaVersionUnknown) && StatusCheck(4)==0
%         EddyCommand = '/bin/eddy_cuda9.1';
%         fprintf('Found CUDA (9.1) installed, trying if this works\n');
%     elseif (StatusCheck(2)==0 || CudaVersionUnknown) && StatusCheck(5)==0
%         EddyCommand = '/bin/eddy_cuda8.0';
%         fprintf('Found CUDA (8.0) installed, trying if this works\n');
%     elseif CudaVersionUnknown && StatusCheck(6)==0
%         EddyCommand = '/bin/eddy_cuda';
%         fprintf('Found CUDA installed, trying if this works\n');
%     elseif CudaVersionUnknown && StatusCheck(4)==0
%         EddyCommand = '/bin/eddy_cuda9.1';
%         fprintf('Found CUDA installed, trying if this works\n');
%     else
%         EddyCommand = '/bin/eddy_openmp';
%         fprintf('No CUDA found, running CPU only\n');
%     end

    if exist(fullfile(RootWSLdir,'bin','eddy_openmp'), 'file')
        EddyCommand = '/bin/eddy_openmp';
    elseif exist(fullfile(RootWSLdir,'bin','eddy'), 'file')
        EddyCommand = '/bin/eddy';
    else
        error('Cannot find correct eddy function!');
    end
        

    % Run eddy current
    fprintf('%s\n', ['Computing eddy current correction using ' EddyCommand]);

    EddyCommand = [EddyCommand ' --imain=' xASL_adm_UnixPath(x.P.Path_dwi) ' --mask='...
        xASL_adm_UnixPath(x.P.Path_dwi_mask) ' --acqp=' xASL_adm_UnixPath(PathAcqp)...
        ' --index=' xASL_adm_UnixPath(IndexPath) ' --bvecs='...
        xASL_adm_UnixPath(PathBvec) ' --bvals=' xASL_adm_UnixPath(PathBval) ' --out='...
        xASL_adm_UnixPath(PathEddy) ' --topup=' xASL_adm_UnixPath(PathTopUp) ' --repol --verbose'];
    % --repol does some additional correction of motion-related artifacts
    % --verbose to see it running
    
    if x.Quality
        EddyCommand = [EddyCommand ' --niter=6 --flm=quadratic --interp=spline'];
    else
        EddyCommand = [EddyCommand ' --niter=1 --flm=linear --interp=trilinear'];
    end

    xASL_fsl_RunFSL(EddyCommand, x);
    xASL_delete(PathWithoutOutlierRemoval); % this is Eddy without outlier removal, don't need it
    xASL_delete(x.P.Path_dwi_mean);
    xASL_delete(x.P.Path_dwi_mask);
    
    x.mutex.AddState('020_EddyCurrent');
    x.mutex.DelState('030_RegistrationDWI2T1w');
elseif   bO; fprintf('%s\n','020_EddyCurrent has already been performed, skipping...');
end

%% -----------------------------------------------------------------------------
%% 3    Registration to T1w
if ~x.mutex.HasState('030_RegistrationDWI2T1w')

    % Registration to T1w
    OtherList = xASL_adm_GetFileList(x.SESSIONDIR, '^(B0|dwi|Field|TopUp|Unwarped).*\.nii', 'FPList', [0 Inf]);

    xASL_im_CenterOfMass(PathEddyNii, OtherList);
    % Create mutual information reference image for best registration reference
    % Disabled for now, simply using T1w as everyone
%     NewIM = xASL_io_Nifti2Im(x.P.Path_c1T1)+xASL_io_Nifti2Im(x.P.Path_c2T1).*2+xASL_io_Nifti2Im(x.P.Path_c3T1).*3;
%     xASL_io_SaveNifti(x.P.Path_c1T1, PathMIref, NewIM, [], 0);
    xASL_spm_coreg(x.P.Path_T1, PathEddyNii, OtherList, x, [], true);

    % Remove the mat-files (that have equal parameters and won't be used by
    % FSL, only by SPM), so this would incorrectly not apply the FSL
    % registration parameters in the NIfTI header on the other volumes
    xASL_delete(PathB0Mat);
    xASL_delete(x.P.Path_dwiMat);
    xASL_delete(PathEddyMat);
    xASL_delete(PathUnwarpedMat);

    x.mutex.AddState('030_RegistrationDWI2T1w');
    x.mutex.DelState('040_dtiFit');
elseif   bO; fprintf('%s\n','030_RegistrationDWI2T1w has already been performed, skipping...');
end


%% -----------------------------------------------------------------------------
%% 4    DTI fit to create SSE parameter
if ~x.mutex.HasState('040_dtiFit')

    xASL_adm_DeleteFileList(x.SESSIONDIR, '^DTIfit.*', false, [0 Inf]);
    xASL_delete(x.P.Path_dwi_mean);
    
    % Compute brain mask
    fprintf('Computing mean DWI image...\n');
    xASL_fsl_RunFSL(['/bin/fslmaths ' xASL_adm_UnixPath(PathEddyNii) ' -Tmean '...
        xASL_adm_UnixPath(x.P.Path_dwi_mean)], x);
    fprintf('Running Bet for brainmasking...\n');
    xASL_fsl_RunFSL(['/bin/bet ' xASL_adm_UnixPath(x.P.Path_dwi_mean) ' '...
        xASL_adm_UnixPath(x.P.Path_dwi_mean) ' -R -m'], x);     
    
    fprintf('Computing DTI fit...');
    xASL_fsl_RunFSL(['/bin/dtifit --data=' xASL_adm_UnixPath(PathEddy) ' --out='...
        xASL_adm_UnixPath(PathFit) ' --mask=' xASL_adm_UnixPath(x.P.Path_dwi_mask) ' --bvecs='...
        xASL_adm_UnixPath(PathBvec) ' --bvals=' xASL_adm_UnixPath(PathBval) ' --sse'], x, [], [], false);

% output files:
% <basename>_V1 - 1st eigenvector
% <basename>_V2 - 2nd eigenvector
% <basename>_V3 - 3rd eigenvector
% <basename>_L1 - 1st eigenvalue
% <basename>_L2 - 2nd eigenvalue
% <basename>_L3 - 3rd eigenvalue
% <basename>_MD - mean diffusivity
% <basename>_FA - fractional anisotropy
% <basename>_MO - mode of the anisotropy (oblate ~ -1; isotropic ~ 0; prolate ~ 1)
% <basename>_S0 - raw T2 signal with no diffusion weighting
% optional output
% <basename>_sse - Sum of squared error


    x.mutex.AddState('040_dtiFit');
    x.mutex.DelState('050_resliceDWI');
elseif   bO; fprintf('%s\n','040_dtiFit has already been performed, skipping...');
end


%% -----------------------------------------------------------------------------
%% 5    Reslice DTI data
if ~x.mutex.HasState('050_resliceDWI')
    x.P.Path_y_ASL = fullfile(x.SESSIONDIR, 'y_DWI.nii');

    % 1    First create dwi-standard space flow field
    %       If no T1 flow field exists, create an identity flowfield.
    %       We now assume the structural module hasn't run, and we simply want to run the ASL module to quickly check how the images look like
    %       So we only run the automatic Center of Mass ACPC alignment
    if ~xASL_exist(x.P.Path_y_T1,'file') && ~x.Quality
        IDmatrixPath = fullfile(x.D.MapsDir,'Identity_Deformation_y_T1.nii');
        xASL_Copy(IDmatrixPath, x.P.Path_y_ASL, true);

    elseif ~xASL_exist(x.P.Path_y_T1,'file')
        error('Structural module did not run correctly yet');
    else
        xASL_wrp_CreateASLDeformationField(x, true, [2.5 2.5 2.5], x.P.Path_dwi_mean); % assume bit lower resolution than the 2 mm that is usually used
    end

    % output QC parameters: average SSE? over mask?
    % SSE image?
    % FA image
    % first image (B0 image)
    % images in multiple orientations

    InputPaths = {PathB0, PathUnwarped, PathFA, PathSSE, x.P.Path_dwi_mean};
    OutputPaths = {PathPopB0, PathPopUnwarped, PathPopFA, PathPopSSE, PathPopMean};

    xASL_spm_deformations(x, InputPaths, OutputPaths, [], [], [], x.P.Path_y_ASL);

    x.mutex.AddState('050_resliceDWI');
    x.mutex.DelState('060_visualize');
elseif  bO; fprintf('%s\n','050_resliceDWI has already been performed, skipping...');
end


%% -----------------------------------------------------------------------------
%% 6    Visual check
if ~x.mutex.HasState('060_visualize')

    %% Initialization
    fprintf('%s\n','print visual quality assurance checks');
    Parms.ModuleName = 'dwi';
    close all % close all Figures to avoid capturing & saving the wrong Figure

    x = xASL_adm_LoadX(x, PathX, true); % assume x.mat is newer than x

    % Clear any previous QC images or parameters
    if isfield(x,'Output_im') && isfield(x.Output_im,'dwi')
       x.Output_im = rmfield(x.Output_im,'dwi');
    end   
    if isfield(x,'Output') && isfield(x.Output,'dwi')
       x.Output = rmfield(x.Output,'dwi');
    end
    
    x.D.dwiCheckDir = fullfile(x.D.PopDir, 'dwiCheck');
    xASL_adm_CreateDir(x.D.dwiCheckDir);

    if ~isfield(x,'Output')
        x.Output = struct;
    end
    if ~isfield(x.Output,'dwi')
        x.Output.dwi = struct;
    end

    %% Visual QC of BVEC values
    DTI_VisualQC_BVEC(PathBvec, x);

    %% Visual QC of TopUp
    [Output1, Output2] = xASL_im_VisualQC_TopUp(PathPopB0, PathPopUnwarped, x, iSubject, x.D.dwiCheckDir);
    
    x.Output.dwi.MeanAI_PreTopUp_Perc = Output1;
    x.Output.dwi.MeanAI_PostTopUp_Perc = Output2;
    %% Get visualization settings
    % Parameters for creating visual QC Figures:

    % Clip SSE image
    tIM = xASL_io_Nifti2Im(PathPopSSE);
    tIM(tIM<0) = 0;
    tIM(tIM>2) = 2;
    xASL_io_SaveNifti(PathPopSSE, PathPopSSE, tIM, [], false);

    ImIn = {{PathPopB0} {PathPopB0 x.P.Pop_Path_rc2T1}};
    ImIn(3) = {{PathPopFA}};
    ImIn(4) = {{PathPopSSE}};

    DirOut = {x.D.dwiCheckDir x.D.dwiCheckDir x.D.dwiCheckDir x.D.dwiCheckDir x.D.dwiCheckDir};
    x.V.IntScale = {1 [1 0.5] 1 1 1};
    x.V.ColorMapIs = {[] [] x.S.jet256 x.S.jet256};
    x.V.NameExt = {[] 'Reg_pWM_' [] []};


    %%  Perform the visualization
    fprintf('%s','Printing images...  ');
    for iM=1:length(ImIn)
        xASL_TrackProgress(iM, length(ImIn));
        % Default parameters
        Pars    = {'ClipZero' 'IntScale' 'NameExt' 'ColorMapIs'};
        for iP=1:length(Pars)
            if ~isfield(x.V,Pars{iP})
                x.V.(Pars{iP}){iM} = [];
            elseif length(x.V.(Pars{iP}))<iM
                x.V.(Pars{iP}){iM} = [];
            end
        end
        % visualize
        Parms.ModuleName = 'dwi';
        Parms.IM = xASL_im_CreateVisualFig(x, ImIn{iM}, DirOut{iM}, x.V.IntScale{iM}, x.V.NameExt{iM}, x.V.ColorMapIs{iM});
        % add single slice to QC collection
        if  sum(~isnan(Parms.IM(:)))>0 % if image is not empty
            x = xASL_im_AddIM2QC(x, Parms);
        end
    end

    %% Compute overlap (intersection)/DICE coefficient of DWI with T1w brainmask
    x.Output.dwi.CoveragePerc = xASL_qc_ComputeFoVCoverage(PathB0, x);

    %% Orientation check
    x.Output.dwi = xASL_qc_ComputeNiftiOrientation(x, x.P.Path_dwi, x.Output.dwi);
    
    %% Calculate DWI stats within the WM
	% Presmooth the image
	xASL_im_PreSmooth(PathEddyNii, x.P.Path_c2T1, x.P.Path_rc2T1, [], [], [], 1);
	
	% Run the transformation
    xASL_spm_reslice(PathEddyNii, x.P.Path_rc2T1, [], 1, x.Quality, x.P.Path_rc2T1, 1);
	
    MaskIM = xASL_io_Nifti2Im(x.P.Path_rc2T1)>0.9;

    DWIdataPaths = {PathSSE, PathFA, PathEddyNii};
    DWIParms = {'SSE' 'FA' 'B0'};

    for iDWI=1:length(DWIdataPaths)
        % Define defaults
        x.Output.dwi.([DWIParms{iDWI} '_WM_Mean']) = NaN;
        x.Output.dwi.([DWIParms{iDWI} '_WM_SD']) = NaN;
        x.Output.dwi.([DWIParms{iDWI} '_WM_Min']) = NaN;
        x.Output.dwi.([DWIParms{iDWI} '_WM_Max']) = NaN;        
        
        DataIM = xASL_io_Nifti2Im(DWIdataPaths{iDWI});
        DataIM = DataIM(:,:,:,1);
        DataIM = DataIM(MaskIM);

        if sum(MaskIM(:))==0
            warning(['Empty maskIM, check ' x.P.Path_rc2T1]);
        elseif numel(DataIM(:))==0
            warning(['Empty image, check ' DWIdataPaths{iDWI}]);
        else
            MeanParm(iDWI) = xASL_stat_MeanNan(DataIM);
            SDevParm(iDWI) = xASL_stat_StdNan(DataIM);
            MinParm(iDWI) = min(DataIM);
            MaxParm(iDWI) = max(DataIM);

            x.Output.dwi.([DWIParms{iDWI} '_WM_Mean']) = MeanParm(iDWI);
            x.Output.dwi.([DWIParms{iDWI} '_WM_SD']) = SDevParm(iDWI);
            x.Output.dwi.([DWIParms{iDWI} '_WM_Min']) = MinParm(iDWI);
            x.Output.dwi.([DWIParms{iDWI} '_WM_Max']) = MaxParm(iDWI);
        end
    end

    xASL_delete(x.P.Path_rc2T1);
    xASL_delete(PathX);

    xASL_qc_PrintOrientation(x.SESSIONDIR, x.P.Path_dwi, x.SESSIONDIR, 'RigidRegdwi');
    % This function summarizes the dwi orientation. Especially check the determinant, for left-right flips
    PathOrientationResults = fullfile(x.SESSIONDIR,'xASL_qc_PrintOrientation_RigidRegdwi.tsv');
    x.Output.dwi = xASL_im_DetermineFlip(x, iSubject, PathOrientationResults, x.Output.dwi);    
    
    save(PathX, 'x'); % future: do this in each xWrapper

    x.mutex.AddState('060_visualize');
elseif  bO; fprintf('%s\n','060_visualize has already been performed, skipping...');
end    
   
%% -----------------------------------------------------------------------------
%% 7    Store QC
if ~x.mutex.HasState('070_QC')
    
    x = xASL_adm_LoadX(x, PathX, true); % assume x.mat is newer than x
    x = xASL_qc_CollectParameters(x, iSubject, 'dwi');
    save(PathX, 'x'); % future: do this in each xWrapper 

    x.mutex.AddState('070_QC');
elseif  bO; fprintf('%s\n','070_QC has already been performed, skipping...');
end    

%% -----------------------------------------------------------------------------
%% 7    WAD-QC
if ~x.mutex.HasState('070_WADQC') && x.DoWADQCDC
    xASL_qc_WADQCDC(x, iSubject, 'dwi');
    x.mutex.AddState('070_WADQC');
elseif x.mutex.HasState('070_WADQC') && bO
    fprintf('%s\n', '070_WADQC has already been performed, skipping...');
end


%% Householding
if x.DELETETEMP
    xASL_delete(PathB0);
    xASL_delete(PathUnwarped);
end



%% -----------------------------------------------------------------------------
%% 999 Ready
x.mutex.AddState('999_ready');
cd(oldFolder);

x.mutex.Unlock();
x.result  = true;
result    = true;

end




%% ========================================================================
% =========================================================================
function DTI_VisualQC_BVEC(PathBvec, x)
%DTI_VisualQC_BVEC Plots distribution of BVECs

Path_QC_Image = fullfile(x.D.dwiCheckDir, ['VisualQC_Bvec_' x.P.SubjectID '.png']);

if ~exist(PathBvec, 'file')
    fprintf('No DTI BVEC file found, skipping Bvec Visual QC');
    return;
else
    fprintf('Running Bvec Visual QC');
    bvecs = load(PathBvec);
end

fig  = figure('Visible','off', 'position',[100 100 500 500]);
subplot(2,2,1);
plot(bvecs(1,:),bvecs(2,:),'*r');
axis([-1 1 -1 1]);
xlabel('Bvec x');
ylabel('Bvec y');
title('Distribution Bvectors');

subplot(2,2,2);
plot(bvecs(2,:),bvecs(3,:),'*r');
axis([-1 1 -1 1]);
xlabel('Bvec y');
ylabel('Bvec z');
title('Distribution Bvectors');

subplot(2,2,3);
plot(bvecs(1,:),bvecs(3,:),'*r');
axis([-1 1 -1 1]);
xlabel('Bvec x');
ylabel('Bvec z');
title('Distribution Bvectors');

subplot(2,2,4);
plot3(bvecs(1,:),bvecs(2,:),bvecs(3,:),'*r');
axis([-1 1 -1 1 -1 1]);
xlabel('Bvec x');
ylabel('Bvec y');
zlabel('Bvec z');
title('Distribution Bvectors');

saveas(fig, Path_QC_Image);

fprintf('\n');

end
