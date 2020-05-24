function [result, x] = xASL_module_func(x)
%xASL_module_func ExploreASL module for fMRI/func pre-processing
%
% FORMAT: [result, x] = xASL_module_func(x)
%
% INPUT:
%   x  - x structure containing all input parameters (REQUIRED)
%   x.SUBJECTDIR  -  anatomical directory, containing the derivatives of anatomical images (REQUIRED)
%   x.SESSIONDIR  -  ASL directory, containing the derivatives of perfusion images (REQUIRED)
%
%
% OUTPUT:
%   result  - true for successful run of this module, false for insuccessful run
%   x       - x structure containing all output parameters (REQUIRED)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL module processes fMRI
% images, this module is created 'quick and dirty' from the ExploreASL ASL
% module. It works and has tested with EPAD data, after initial import, but
% use for own risk in other circumstances.
%
% This module has the following submodules/wrappers:
%
% 1    TopUp
% 2    Motion correction fMRI
% 3    Registration fMRI sessions to T1w
% 4    Reslice fMRI data
% 5    Visual check
% 6    QC
% 7    WAD-QC
%
% EXAMPLE: [~, x] = xASL_module_func(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL




%% -----------------------------------------------------------------------------
%% 0    Admin
x = xASL_init_InitializeMutex(x, 'func'); % starts mutex locking process to ensure that everything will run only once
result = false;

if ~isfield(x,'SavefMRI')
    x.SavefMRI=0;
end

%% 0.9 change working directory to make sure that unspecified output will go there...
oldFolder = cd(x.SESSIONDIR);

[Fpath, Ffile] = fileparts(x.SESSIONDIR);
iSess = find(strcmp(x.SESSIONS,Ffile));
[~, Ffile] = fileparts(Fpath);

x.Sequence      = '2D_EPI';
x.readout_dim   = '2D';

%% Only continue if fMRI exists
x.P.Path_func_bold = fullfile(x.SESSIONDIR, 'func_run-1_bold.nii');
x.P.Path_func_bold_ORI = fullfile(x.SESSIONDIR, 'func_run-1_bold_ORI.nii');
if ~xASL_exist(x.P.Path_func_bold, 'file') && ~xASL_exist(x.P.Path_func_bold_ORI, 'file')
    fprintf('%s\n','No fMRI found, skipping');
    result = true;
    return;
elseif ~xASL_exist(x.P.Path_func_bold, 'file') && xASL_exist(x.P.Path_func_bold_ORI, 'file')
    xASL_Move(x.P.Path_func_bold_ORI, x.P.Path_func_bold);
end

x = xASL_init_FileSystem(x); % this reinitiates x.P

x.P.Path_func_bold = fullfile(x.SESSIONDIR, 'func_run-1_bold.nii'); % needs to be here for re-initiation
x.P.Path_func_bold_mat = fullfile(x.SESSIONDIR, 'func_run-1_bold.mat');
x.P.Path_func_bold_ORI = fullfile(x.SESSIONDIR, 'func_run-1_bold_ORI.nii'); % needs to be here for re-initiation
x.P.Path_func_NormPE = fullfile(x.SESSIONDIR, 'func_run-1_NormPE.nii');
x.P.Path_func_RevPE = fullfile(x.SESSIONDIR, 'func_run-1_RevPE.nii');
x.P.Path_despiked_func_bold = fullfile(x.SESSIONDIR, 'despiked_func_run-1_bold.nii');
x.P.Path_despiked_func_bold_mat = fullfile(x.SESSIONDIR, 'despiked_func_run-1_bold.mat');
x.P.Path_func_bold_parms_mat = fullfile(x.SESSIONDIR, 'func_run-1_bold_parms.mat');

PathField = fullfile(x.SESSIONDIR ,'Field.nii');
PathFieldCoef = fullfile(x.SESSIONDIR ,'TopUp_fieldcoef.nii');
PathB0 = fullfile(x.SESSIONDIR ,'B0.nii');
PathUnwarped = fullfile(x.SESSIONDIR ,'Unwarped.nii');
PathPopB0 = fullfile(x.D.PopDir, ['rFunc_B0_' x.SUBJECTS{x.iSubject} '.nii']);
PathPopUnwarped = fullfile(x.D.PopDir, ['rFunc_Unwarped_' x.SUBJECTS{x.iSubject} '.nii']);

PathX = fullfile(x.SUBJECTDIR,'x.mat'); % later optimize this, do this in each xWrapper


%% Delete old files
if ~x.mutex.HasState('030_register_func') && ~x.mutex.HasState('020_realign_func')
    DelList = {x.P.Path_FoV x.P.Path_mean_control x.P.Path_mean_PWI_Clipped x.P.Path_PWI x.P.Path_mean_PWI_Clipped_sn_mat x.P.Path_SliceGradient x.P.Path_rrM0};
    for iD=1:length(DelList)
        xASL_delete(DelList{iD});
    end
    xASL_adm_DeleteFileList(x.SESSIONDIR,'^(B0|Field|rp|SD|SNR|TopUp|Unwarped|qcdc|xASL_qc|y_ASL|wfunc).*$'); % backwards compatibility
    
    [Fpath1, Ffile1] = xASL_fileparts(x.P.Path_func_bold_ORI);
    [Fpath2, Ffile2] = xASL_fileparts(x.P.Path_func_bold);
    ExtOri = {'.nii' '.json' '.mat' '_sn.mat' '_parms.mat'};
    for iExt=1:length(ExtOri)
        PathOrig = fullfile(Fpath1, [Ffile1 ExtOri{iExt}]);
        PathDest = fullfile(Fpath2, [Ffile2 ExtOri{iExt}]);
        if xASL_exist(PathOrig,'file')
            xASL_Move(PathOrig, PathDest, true);
        end
    end
end

if ~isfield(x,'motion_correction')
    x.motion_correction   = 1;
end



if ~x.mutex.HasState('999_ready')
    bO = true; % generate output, some processing has and some has not been yet done
    % If fMRI_4D_parms.mat (fMRI parameters) file exist, load and overwrite existing parameters in the xASL file (inheritance principle)
    [~, x] = xASL_adm_LoadParms(x.P.Path_func_bold_parms_mat, x);
else
    bO = false; % skip output, as all processing has been performed
end

%% Check existence structural reference files, of which this module is dependent
if (~xASL_exist(x.P.Path_c1T1,'file') || ~xASL_exist(x.P.Path_c2T1,'file') || ~xASL_exist(x.P.Path_y_T1,'file')) && bO
    warning('Structural files missing');
end

[Fpath, Ffile, Fext] = xASL_fileparts(x.P.Path_func_bold);

if ~isfield(x,'DoWADQCDC')
    x.DoWADQCDC = false;
end

%% -----------------------------------------------------------------------------
%% 1    TopUp
if ~x.mutex.HasState('010_TopUp_func') || ~xASL_exist(PathFieldCoef, 'file')
    
    bSuccess = xASL_fsl_TopUp(x.SESSIONDIR, 'func', x, x.P.Path_func_bold);

    if bSuccess
        InputPaths = {PathB0, PathUnwarped};
        OutputPaths = {PathPopB0, PathPopUnwarped};
        xASL_spm_deformations(x, InputPaths, OutputPaths); % that they are there for visual QC at the end

        x.mutex.AddState('010_TopUp_func');
        x.mutex.DelState('020_realign_fMRI');
    end
elseif   bO; fprintf('%s\n','001_TopUp has already been performed, skipping...');
end


%% -----------------------------------------------------------------------------
%% 2    Motion correction fMRI

if ~x.motion_correction
    if bO; fprintf('%s\n','Motion correction was disabled, skipping'); end
else
    if ~x.mutex.HasState('020_realign_func')

        % Remove previous files
        DelList = {x.P.Path_func_bold_mat x.P.Path_despiked_func_bold x.P.Path_despiked_func_bold_mat 'rp_func_bold.txt'};
        for iD=1:length(DelList)
            xASL_delete(fullfile(x.SESSIONDIR, DelList{iD}));
        end

%         % First, solve dimensionality (in case there are empty dims, that need restructuring)
        nVol = size(xASL_io_Nifti2Im(x.P.Path_func_bold),4);

        % Then, check matrix size: throw error if 2D data with 3 dimensions only


        if  nVol>1
            % Before motion correction, we align the images with ACPC
            OtherList = {x.P.Path_func_bold_ORI, PathB0, PathUnwarped, x.P.Path_func_NormPE, x.P.Path_func_RevPE, PathField, PathFieldCoef};
            xASL_im_CenterOfMass(x.P.Path_func_bold, OtherList, 10);

            % Run motion Correction
            xASL_wrp_RealignASL(x, false); % same as ASL, but no subtraction/pairwise data
        else
            warning(['Skipping motion correction, fMRI ' x.P.SubjectID '_' x.P.SessionID ' only had 1 volume!!!!!']);
        end

        x.mutex.AddState('020_realign_func');
        x.mutex.DelState('030_register_func');
    elseif   bO; fprintf('%s\n','020_realign_func session has already been performed, skipping...');
    end
end


% no spike removal performed for fMRI
nVol = size(xASL_io_Nifti2Im(x.P.Path_func_bold),4);

%% -----------------------------------------------------------------------------
%% 3    Registration fMRI sessions to T1w

if ~x.mutex.HasState('030_register_func')

    xASL_wrp_Register_func(x);

    x.mutex.AddState('030_register_func');
    x.mutex.DelState('040_reslice_func');
elseif  bO; fprintf('%s\n','030_register_func has already been performed, skipping...');
end




%% -----------------------------------------------------------------------------
%% 4    Reslice fMRI data
if ~x.mutex.HasState('040_reslice_func') && x.mutex.HasState('030_register_func')

    xASL_wrp_Reslice_func(x);

    x.mutex.AddState('040_reslice_func');
    x.mutex.DelState('050_visualize');
    x.mutex.DelState('060_QC');
elseif  bO; fprintf('%s\n','040_reslice_func has already been performed, skipping...');
end




%% -----------------------------------------------------------------------------
%% 5    Visual check
if ~x.mutex.HasState('050_visualize')

    fprintf('%s\n','print visual quality assurance checks');
    Parms.ModuleName = 'func';
    close all % close all Figures to avoid capturing & saving the wrong Figure

    x = xASL_adm_LoadX(x, fullfile(x.SUBJECTDIR,'x.mat'), true); % assume x.mat is newer than x

    % Clear any previous QC images & parameters
    if isfield(x,'Output_im') && isfield(x.Output_im,'func')
       x.Output_im = rmfield(x.Output_im,'func');
    end    
    if isfield(x,'Output') && isfield(x.Output,'func')
       x.Output = rmfield(x.Output,'func');
    end        
    
    %% Get visualization settings
    % Parameters for creating visual QC Figures:
    % CBF, CBF with overlay c2T1, CBF with overlay c2T1
    % MeanControl SD
    % M0 NoSmoothM0 NoSmoothM0 with overlay c1T1
    % TT TT with overlay c2T1

    ImIn = {{x.P.Pop_Path_mean_control}  {x.P.Pop_Path_mean_control x.P.Pop_Path_rc2T1}};
    ImIn(3:4) = {{x.P.Pop_Path_SD} {x.P.Pop_Path_SNR}};

    x.D.FuncCheckDir = fullfile(x.D.PopDir, 'FuncCheck');

    DirOut          = {x.D.FuncCheckDir x.D.FuncCheckDir x.D.FuncCheckDir x.D.FuncCheckDir x.D.FuncCheckDir};

    x.V.IntScale(2)   = {[0.5 0.5]};

    x.V.NameExt         = {[] 'Reg_pWM_'};

    %%  Perform the visualization
    fprintf('%s','Printing images...  ');
    for iM=1:length(ImIn)
        xASL_TrackProgress(iM, length(ImIn));
        % Default parameters
        Pars    = {'ClipZero' 'IntScale' 'NameExt' 'ColorMapIs'};
        for iP=1:length(Pars)
            if     ~isfield(x.V,Pars{iP})
                    x.V.(Pars{iP}){iM} = [];
            elseif  length(x.V.(Pars{iP}))<iM
                    x.V.(Pars{iP}){iM} = [];
            end
        end
        % visualize
        Parms.ModuleName = 'func';
        Parms.IM      = xASL_im_CreateVisualFig( x, ImIn{iM}, DirOut{iM}, x.V.IntScale{iM}, x.V.NameExt{iM}, x.V.ColorMapIs{iM});
        % add single slice to QC collection
        if  sum(~isnan(Parms.IM(:)))>0 % if image is not empty
            x   = xASL_im_AddIM2QC(x,Parms);
        end
    end

    if xASL_exist(PathPopB0,'file') && xASL_exist(PathPopUnwarped,'file')% if we have TopUp results
        [Output1, Output2] = xASL_im_VisualQC_TopUp(PathPopB0, PathPopUnwarped, x, x.iSubject, x.D.FuncCheckDir);
        x.Output.func(x.iSubjectSession).MeanAI_PreTopUp_Perc = Output1;
        x.Output.func(x.iSubjectSession).MeanAI_PostTopUp_Perc = Output2;
        xASL_delete(PathPopB0);
        xASL_delete(PathPopUnwarped);
    end

    save(PathX,'x'); % future: do this in each xWrapper

        x.mutex.AddState('050_visualize');
        x.mutex.DelState('060_QC');        
elseif  bO; fprintf('%s\n','050_visualize has already been performed, skipping...');
end

%% -----------------------------------------------------------------------------
%% 6    QC

if ~x.mutex.HasState('060_QC')
    x = xASL_adm_LoadX(x, PathX, true); % assume x.mat is newer than x

    x.Output.func.CoveragePerc = xASL_qc_ComputeFoVCoverage(x.P.Path_func_bold, x);

    %  Run this part only if we have > 10 time points
    if nVol>10
        if  xASL_exist(x.P.Path_c1T1,'file') && xASL_exist(x.P.Path_c2T1,'file')
			% Presmooth 
			xASL_im_PreSmooth(x.P.Path_func_bold,x.P.Path_c1T1,x.P.Path_rc1T1,[],[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);
			xASL_im_PreSmooth(x.P.Path_func_bold,x.P.Path_c2T1,x.P.Path_rc2T1,[],[],x.P.Path_mean_PWI_Clipped_sn_mat, 1);

			% first reslice tissue maps to ASL space
			xASL_spm_reslice(x.P.Path_func_bold, x.P.Path_rc1T1, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality,x.P.Path_rc1T1);
			xASL_spm_reslice(x.P.Path_func_bold, x.P.Path_rc2T1, x.P.Path_mean_PWI_Clipped_sn_mat, 1, x.Quality,x.P.Path_rc2T1);
			% Save time series, required for ASL temporal SNR
			
			TempFunc = xASL_qc_temporalSNR(x.P.Path_func_bold, {x.P.Path_rc1T1 x.P.Path_rc2T1});
			FuncFields = fields(TempFunc);
			for iA=1:length(FuncFields)
				x.Output.func.(FuncFields{iA}) = TempFunc.(FuncFields{iA});
			end
			
		else
			warning('Skipping SPM UP QC because structural files didnt exist');
        end
    end

    xASL_qc_PrintOrientation(x.SESSIONDIR, x.P.Path_func_bold, x.SESSIONDIR, 'RigidRegfunc');
    % This function summarizes the func orientation. Especially check the determinant, for left-right flips    

    x = xASL_qc_CollectParameters(x, x.iSubject, 'func');
    
    xASL_delete(PathX);
    save(PathX,'x');

    x.mutex.AddState('060_QC');
    x.mutex.DelState('070_WADQC');
elseif  bO; fprintf('%s\n','060_QC has already been performed, skipping...');
end

%% -----------------------------------------------------------------------------
%% 7    WAD-QC
if ~x.mutex.HasState('070_WADQC') && x.DoWADQCDC
    xASL_qc_WADQCDC(x, x.iSubj, 'func');
    x.mutex.AddState('070_WADQC');
elseif x.mutex.HasState('070_WADQC') && bO
    fprintf('%s\n', '070_WADQC has already been performed, skipping...');
end



%% -----------------------------------------------------------------------------
%% 999 Ready
x.mutex.AddState('999_ready');
cd(oldFolder);

x.mutex.Unlock();
x.result  = true;
result    = true;

end
