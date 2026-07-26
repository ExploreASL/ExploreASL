function [result, x] = xASL_module_dce(x)
%xASL_module_dce ExploreASL module for DCE pre-processing
%
% FORMAT: [result, x] = xASL_module_dce(x)
%
% INPUT:
%   x  - x structure containing all input parameters (REQUIRED)
%
% OUTPUT:
%   result  - true for successful run of this module, false for insuccessful run
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL module processes  DCE
% images, this module is created 'quick and dirty' from the ExploreASL ASL
% module. It has been tested with TheriniBio data only, use for own risk in other circumstances.
%
% Scripts used (and based on) OSIPY
%
% This module has the following submodules/wrappers:
%
% 0     000_CopyRawdata         Copy from rawdata
% 1     010_MotionEstimation    Estimate motion
% 2     020_RegisterDCE2T1      Registration DCE to T1
% 3     030_MotionCorrection    Motion correction
% 4     040_PVmapsMasks         Create PV maps & masks
% 5     050_AIF                 Get AIF
% 6     060_T1mapping           T1 mapping
% 7     070_OSIPY               OSIPY signal->concentration->Ktrans fitting
% 8     080_QC                  QC
% 9     090_Visualization       QC visualization



%
% EXAMPLE: [~, x] = xASL_module_dce(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2026-... ExploreASL





% x.modules.dce.Model = 'patlak';
% x.modules.dce.TemporalResolution
% x.modules.dce.bMotionCorrection = true;
% x.modules.dce.bUseMeasuredAIF = true;
% x.modules.dce.bRunQuantification = true;


%% -----------------------------------------------------------------------------
%% 0    Admin
x = xASL_init_InitializeMutex(x, 'dce' ); % starts mutex locking process to ensure that everything will run only once
result = false;


% Check whether we already ran this (to know whether we should give verbose feedback)
if ~x.mutex.HasState('999_ready')
    bO = true; % generate output, some processing has and some has not been yet done
else
    bO = false; % skip output, as all processing has been performed
end

if ~isfield(x,'DoWADQCDC')
    x.DoWADQCDC = false;
end

% if ~isfield(x,'bAutomaticallyDetectPython')
%     x.bAutomaticallyDetectPython = false;
% end

% Skip if existing
if x.mutex.HasState('999_ready')
    result = true;
    x.mutex.Unlock();
    return;
end


%% Parameters
x.modules.dce.InjectionFrame = 4;

% Python
x.modules.dce.PathPython = '/Users/hjmutsaerts/venvs/py312/bin/python3.12';
x.modules.dce.PythonCodePath = fullfile(x.opts.MyPath, 'WorkInProgress', 'DCE');


%% Paths
x = xASL_init_FileSystem(x);
PathX = fullfile(x.dir.SUBJECTDIR, 'x.mat');

x.SESSIONS = {'DCE_1'};
x.SESSION = 'DCE_1';
x.dir.SESSIONDIR = fullfile(x.dir.SUBJECTDIR, x.SESSION);

% destination paths
% dir_dceDest -> x.SESSIONDIR
% nii_dceDest = x.P.Path_DCE4D
% json_dceDest = x.P.Path_DCE4D_json

x.P.Path_DCE4D = fullfile(x.dir.SESSIONDIR, 'DCE4D.nii');
x.P.Path_DCE4D_json = fullfile(x.dir.SESSIONDIR, 'DCE4D.json');
x.P.Path_rDCE4D = fullfile(x.dir.SESSIONDIR, 'rDCE4D.nii');

x.P.Path_despotIR = fullfile(x.dir.SESSIONDIR, 'despotIR.nii');
x.P.Path_despotFA1 = fullfile(x.dir.SESSIONDIR, 'despotFA_flip-01.nii');
x.P.Path_despotFA2 = fullfile(x.dir.SESSIONDIR, 'despotFA_flip-02.nii');
x.P.Path_despotIR_json = fullfile(x.dir.SESSIONDIR, 'despotIR.json');
x.P.Path_despotFA1_json = fullfile(x.dir.SESSIONDIR, 'despotFA_flip-01.json');
x.P.Path_despotFA2_json = fullfile(x.dir.SESSIONDIR, 'despotFA_flip-02.json');

x.P.Path_VFA = fullfile(x.dir.SESSIONDIR, 'despot_vfa.nii'); % variable flip angle
x.P.Path_rVFA = fullfile(x.dir.SESSIONDIR, 'rdespot_vfa.nii'); % resampled

% Derivatives
x.P.Path_DCE_mean = fullfile(x.dir.SESSIONDIR, 'DCE_mean.nii');
x.P.Path_DCE_SD = fullfile(x.dir.SESSIONDIR, 'DCE_SD.nii');
x.P.Path_DCE_SNR = fullfile(x.dir.SESSIONDIR, 'DCE_SNR.nii');
x.P.Path_DCE_CoV = fullfile(x.dir.SESSIONDIR, 'DCE_CoV.nii');

path_RealignParameters = fullfile(x.dir.SESSIONDIR, 'rp_DCE4D.txt');

x.P.Path_PVgm = fullfile(x.dir.SESSIONDIR, 'PVgm.nii');
x.P.Path_PVwm = fullfile(x.dir.SESSIONDIR, 'PVwm.nii');
x.P.Path_PVcsf = fullfile(x.dir.SESSIONDIR, 'PVcsf.nii');
x.P.Path_PV_T1w = fullfile(x.dir.SESSIONDIR, 'PVt1w.nii');

x.P.Path_y_DCE = fullfile(x.dir.SESSIONDIR, 'y_DCE.nii');
x.P.Path_AIF = fullfile(x.dir.SESSIONDIR, 'AIF.mat');

x.P.Path_DCE4D_Concentration = fullfile(x.dir.SESSIONDIR, 'DCE4D_Concentration.nii');
x.P.Path_Ktrans = fullfile(x.dir.SESSIONDIR, 'Ktrans.nii');
x.P.Path_Vp = fullfile(x.dir.SESSIONDIR, 'Vp.nii');

x.P.Path_DCE_RSquared = fullfile(x.dir.SESSIONDIR, 'DCE_RSquared.nii');
x.P.Path_DCE_FitMask = fullfile(x.dir.SESSIONDIR, 'DCE_FitMask.nii');

x.D.DCECheckDir = fullfile(x.D.PopDir, 'DCECheck');

% Standard space
x.P.Pop_Path_DCE_mean = fullfile(x.D.PopDir, ['DCE_mean_' x.SUBJECT '.nii']);
x.P.Pop_Path_DCE_CoV = fullfile(x.D.PopDir, ['DCE_CoV_' x.SUBJECT '.nii']);
x.P.Pop_Path_DCE_Ktrans = fullfile(x.D.PopDir, ['DCE_Ktrans_' x.SUBJECT '.nii']);
x.P.Pop_Path_DCE_Vp = fullfile(x.D.PopDir, ['DCE_Vp_' x.SUBJECT '.nii']);


%% Original paths
SubjVisit = x.SUBJECT; % change this for the subject/session we're processing
visitIndex = strfind(SubjVisit, '_');

if numel(visitIndex)~=1
    error(['Too many underscores found in subject_visit for ' SubjVisit]);
else
    subject = SubjVisit(1:visitIndex-1);
    visit = SubjVisit(visitIndex+1:end);
end

dir_dceOrig = fullfile(x.dir.RawData, subject, ['ses-' visit], 'dce');
nii_dceOrig = xASL_adm_GetFileList(dir_dceOrig, ['^' subject '_ses-' visit '_dce.*\.nii'], 'FPList');
json_dceOrig = xASL_adm_GetFileList(dir_dceOrig, ['^' subject '_ses-' visit '_dce.*\.json'], 'FPList');

% Despot T1 mapping
dir_despotOrig = fullfile(x.dir.RawData, subject, ['ses-' visit], 'anat');
nii_despotOrig_despotIR = xASL_adm_GetFileList(dir_despotOrig, ['^' subject '_ses-' visit '_acq-despotIR.*\.nii'], 'FPList');
nii_despotOrig_despotFA1 = xASL_adm_GetFileList(dir_despotOrig, ['^' subject '_ses-' visit '_acq-despotFA_flip-01.*\.nii'], 'FPList');
nii_despotOrig_despotFA2 = xASL_adm_GetFileList(dir_despotOrig, ['^' subject '_ses-' visit '_acq-despotFA_flip-02.*\.nii'], 'FPList');


%% -----------------------------------------------------------------------------
%% 0.   Copy rawdata->derivatives
if ~x.mutex.HasState('000_CopyRawdata')

    fprintf('\n\n\n%s\n', '=========================================================================');
    fprintf('%s\n', '0. Copying rawdata->derivatives');
    fprintf('%s\n\n\n\n', '=========================================================================');

    xASL_delete(x.dir.SESSIONDIR, true);
    xASL_adm_CreateDir(x.dir.SESSIONDIR);

    xASL_Copy(nii_dceOrig{1}, x.P.Path_DCE4D);
    xASL_Copy(json_dceOrig{1}, x.P.Path_DCE4D_json);
    
    % Despot T1 mapping
    xASL_Copy(nii_despotOrig_despotIR{1}, x.P.Path_despotIR);
    xASL_Copy(nii_despotOrig_despotFA1{1}, x.P.Path_despotFA1);
    xASL_Copy(nii_despotOrig_despotFA2{1}, x.P.Path_despotFA2);
    xASL_Copy([nii_despotOrig_despotIR{1}(1:end-4) '.json'], x.P.Path_despotIR_json);
    xASL_Copy([nii_despotOrig_despotFA1{1}(1:end-4) '.json'], x.P.Path_despotFA1_json);
    xASL_Copy([nii_despotOrig_despotFA2{1}(1:end-4) '.json'], x.P.Path_despotFA2_json);
    
    if  (~xASL_exist(x.P.Path_c1T1, 'file') || ~xASL_exist(x.P.Path_c2T1, 'file') || ~xASL_exist(x.P.Path_y_T1, 'file')) && ~x.mutex.HasState('999_ready')
        warning('Structural files missing');
    end

    % Unzipping
    xASL_adm_UnzipNifti(x.P.Path_DCE4D);
    xASL_adm_UnzipNifti(x.P.Path_T1);

    % Convert to single precision
    xASL_io_SaveNifti(x.P.Path_DCE4D, x.P.Path_DCE4D, xASL_io_Nifti2Im(x.P.Path_DCE4D), 32);
    xASL_io_SaveNifti(x.P.Path_despotIR, x.P.Path_despotIR, xASL_io_Nifti2Im(x.P.Path_despotIR), 32);
    xASL_io_SaveNifti(x.P.Path_despotFA1, x.P.Path_despotFA1, xASL_io_Nifti2Im(x.P.Path_despotFA1), 32);
    xASL_io_SaveNifti(x.P.Path_despotFA2, x.P.Path_despotFA2, xASL_io_Nifti2Im(x.P.Path_despotFA2), 32);

    x.mutex.AddState('000_CopyRawdata');
    x.mutex.DelState('010_MotionEstimation');
end





%% Check existence structural reference files, of which this module is dependent


% % % % % % % % %% Delete derivatives from a previous run
% % % % % % % % if ~x.mutex.HasState('030_RegistrationDWI2T1w') && ~x.mutex.HasState('020_EddyCurrent')
% % % % % % % %     xASL_adm_DeleteFileList(x.SESSIONDIR, '^(B0|dwi.*mean|dwi.*mask|Eddy|wdwi|DTIfit|Field|Index\.txt|Unwarped|y_DWI|TopUp|xASL_qc|qcdc).*', false, [0 Inf]);
% % % % % % % % 
% % % % % % % %     [Fpath1, Ffile1] = xASL_fileparts(x.P.Path_dwi_ORI);
% % % % % % % %     [Fpath2, Ffile2] = xASL_fileparts(x.P.Path_dwi);
% % % % % % % %     ExtOri = {'.nii' '.json' '.mat' '_sn.mat' '_parms.mat'};
% % % % % % % %     for iExt=1:length(ExtOri)
% % % % % % % %         PathOrig = fullfile(Fpath1, [Ffile1 ExtOri{iExt}]);
% % % % % % % %         PathDest = fullfile(Fpath2, [Ffile2 ExtOri{iExt}]);
% % % % % % % %         if xASL_exist(PathOrig,'file')
% % % % % % % %             xASL_Move(PathOrig, PathDest, true);
% % % % % % % %         end
% % % % % % % %     end    
% % % % % % % % end



nii = xASL_io_ReadNifti(x.P.Path_DCE4D);
nVolumes = nii.hdr.dim(5);
optimFWHM_Res_mm = nii.hdr.pixdim(2:4);



%% -----------------------------------------------------------------------------
%% 1    Motion estimation
if ~x.mutex.HasState('010_MotionEstimation')

    fprintf('\n\n\n%s\n', '=========================================================================');
    fprintf('%s\n', '1. Motion estimation');
    fprintf('%s\n\n\n\n', '=========================================================================');

    flags.quality = 1;
    
    flags.sep = double(min(nii.hdr.pixdim(2:4)));
    flags.rtm = 1; % realign to mean
    flags.interp = 1;
    flags.graphics = 0;
    
    spm_realign(spm_vol(x.P.Path_DCE4D), flags, 0);


    x.mutex.AddState('010_MotionEstimation');
    x.mutex.DelState('020_RegisterDCE2T1');
elseif   bO; fprintf('%s\n','010_MotionEstimation has already been performed, skipping...');
end



%% -----------------------------------------------------------------------------
%% 2    Registration to T1w
if ~x.mutex.HasState('020_RegisterDCE2T1')

    fprintf('\n\n\n%s\n', '=========================================================================');
    fprintf('%s\n', '2. Registration DCE to T1w');
    fprintf('%s\n\n\n\n', '=========================================================================');

    % Create temporary DCE_mean
    temp_meanDCEpath = fullfile(x.dir.SESSIONDIR, 'temp_DCE_mean.nii');
    temp_meanDCEim = xASL_stat_MeanNan(xASL_io_Nifti2Im(x.P.Path_DCE4D), 4);
    xASL_io_SaveNifti(x.P.Path_DCE4D, temp_meanDCEpath, temp_meanDCEim);

    % Registration to T1w
    OtherList = xASL_adm_GetFileList(x.dir.SESSIONDIR, '^(DCE4D|despot).*\.nii', 'FPList');

    xASL_spm_coreg(x.P.Path_T1, temp_meanDCEpath, OtherList, x);

    xASL_delete(temp_meanDCEpath);

    x.mutex.AddState('020_RegisterDCE2T1');
    x.mutex.DelState('030_MotionCorrection');
elseif   bO; fprintf('%s\n','020_RegisterDCE2T1 has already been performed, skipping...');
end


%% -----------------------------------------------------------------------------
%% 3    Motion correction & resampling in native space
if ~x.mutex.HasState('030_MotionCorrection')

    fprintf('\n\n\n%s\n', '=========================================================================');
    fprintf('%s\n', '3. Motion correction');
    fprintf('%s\n\n\n\n', '=========================================================================');

    for iVol=1:nVolumes
	    matlabbatch{1}.spm.spatial.realign.write.data{iVol, 1} = [x.P.Path_DCE4D ',' num2str(iVol)];
    end
    
    matlabbatch{1}.spm.spatial.realign.write.roptions.which     = [2 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.interp    = 4;
    matlabbatch{1}.spm.spatial.realign.write.roptions.wrap      = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.mask      = 1;
    matlabbatch{1}.spm.spatial.realign.write.roptions.prefix    = 'r';
    
    spm_jobman('run',matlabbatch);
    matlabbatch = cell(0);

    x.mutex.AddState('030_MotionCorrection');
    x.mutex.DelState('040_PVmaps');
elseif   bO; fprintf('%s\n','030_MotionCorrection has already been performed, skipping...');
end



%% -----------------------------------------------------------------------------
%% 4    Create PV maps & y_DCE
if ~x.mutex.HasState('040_PVmaps')

    fprintf('\n\n\n%s\n', '=========================================================================');
    fprintf('%s\n', '4. Create PV maps & y_DCE');
    fprintf('%s\n\n\n\n', '=========================================================================');

    xASL_im_PreSmooth(x.P.Path_rDCE4D, x.P.Path_c1T1, x.P.Path_PVgm, optimFWHM_Res_mm);
    xASL_im_PreSmooth(x.P.Path_rDCE4D, x.P.Path_c2T1, x.P.Path_PVwm, optimFWHM_Res_mm);
    xASL_im_PreSmooth(x.P.Path_rDCE4D, x.P.Path_c3T1, x.P.Path_PVcsf, optimFWHM_Res_mm);
    xASL_im_PreSmooth(x.P.Path_rDCE4D, x.P.Path_T1, x.P.Path_PV_T1w, optimFWHM_Res_mm); % for computing TC

    xASL_spm_reslice(x.P.Path_rDCE4D, x.P.Path_PVgm, [], [], x.settings.Quality, x.P.Path_PVgm);
    xASL_spm_reslice(x.P.Path_rDCE4D, x.P.Path_PVwm, [], [], x.settings.Quality, x.P.Path_PVwm);
    xASL_spm_reslice(x.P.Path_rDCE4D, x.P.Path_PVcsf, [], [], x.settings.Quality, x.P.Path_PVcsf);
    xASL_spm_reslice(x.P.Path_rDCE4D, x.P.Path_PV_T1w, [], [], x.settings.Quality, x.P.Path_PV_T1w);

    % Create y_DCE transformation field
    niiT1w = xASL_io_ReadNifti(x.P.Path_T1);
    resSrc = niiT1w.hdr.pixdim(2:4);
    xASL_im_PreSmooth(x.P.Path_rDCE4D, x.P.Path_y_T1, x.P.Path_y_DCE, [], resSrc);

    x.mutex.AddState('040_PVmaps');
    x.mutex.DelState('050_AIF');
elseif   bO; fprintf('%s\n','040_PVmaps has already been performed, skipping...');
end


%% -----------------------------------------------------------------------------
%% 5    Get AIF
if ~x.mutex.HasState('050_AIF')

    fprintf('\n\n\n%s\n', '=========================================================================');
    fprintf('%s\n', '5. Get AIF');
    fprintf('%s\n\n\n\n', '=========================================================================');

    % Create CSF mask
    pvGM = xASL_io_Nifti2Im(x.P.Path_PVgm);
    pvWM = xASL_io_Nifti2Im(x.P.Path_PVwm);
    pvCSF = xASL_io_Nifti2Im(x.P.Path_PVcsf);
    maskCSF = pvCSF>(pvGM+pvWM) & (pvGM+pvWM+pvCSF)>0.1;

    % Define ROI for venous output function 
    % Warp dilated MNI confluence sinus mask to DCE native space
    pathConfluenceMNI = fullfile(x.D.MapsSPMmodifiedDir, 'ConfluenceSinusMaskDilated.nii');
    pathConfluenceNative = fullfile(x.dir.SESSIONDIR, 'maskConfluenceSinus.nii');
    xASL_spm_deformations(x, pathConfluenceMNI, pathConfluenceNative, 1, x.P.Path_rDCE4D, [], x.P.Path_y_DCE);
    
    imConfluenceMask = xASL_io_Nifti2Im(pathConfluenceNative);
    imConfluenceMask = imConfluenceMask>0.5;
    
    imConfluenceMaskLocal = imConfluenceMask & maskCSF;

    % Get all time curves
    imDCE4D = xASL_io_Nifti2Im(x.P.Path_rDCE4D);
    data2D = reshape(imDCE4D, [], size(imDCE4D, 4));
    voxelTimeCurves = data2D(imConfluenceMaskLocal(:), :);

    nii = xASL_io_ReadNifti(x.P.Path_rDCE4D);
    voxelVolume_mL = prod(nii.hdr.pixdim(2:4))/1000;
    targetVolume_mL = 0.3;
    nVoxels = ceil(targetVolume_mL/voxelVolume_mL);

    % Properties AIF: 1. largest signal difference
    deltaVoxels = max(voxelTimeCurves, [], 2) - min(voxelTimeCurves, [], 2); % difference
    deltaVoxels(:,2) = 1:length(deltaVoxels); % voxel indices
    deltaVoxels = sortrows(deltaVoxels, 1, 'descend'); % largest differences first
    deltaVoxels(:,3) = 1:length(deltaVoxels); % large difference indices (smallest is largest)
    deltaVoxels = sortrows(deltaVoxels, 2, 'ascend'); % sort back to voxel indices

    % Properties AIF: 2. early enhancement (after injection, first 1/3 of frames should have higher intensity than last 1/3 of frames)
    nFramesAfterInjection = nVolumes-x.modules.dce.InjectionFrame-2;
    nFrames_Third = floor(nFramesAfterInjection/3);
    earlyFrames = [x.modules.dce.InjectionFrame+2:x.modules.dce.InjectionFrame+2+nFrames_Third];
    lateFrames = [nVolumes-nFrames_Third:nVolumes];
    sumEarlySignal = sum(voxelTimeCurves(:, earlyFrames), 2);
    sumLateSignal = sum(voxelTimeCurves(:, lateFrames), 2);
    
    diffEarlyLate = sumEarlySignal - sumLateSignal; % difference
    diffEarlyLate(:,2) = 1:length(diffEarlyLate); % voxel indices
    diffEarlyLate = sortrows(diffEarlyLate, 1, 'descend'); % largest differences first
    diffEarlyLate(:,3) = 1:length(diffEarlyLate); % large difference indices (smallest is largest)
    diffEarlyLate = sortrows(diffEarlyLate, 2, 'ascend'); % sort back to voxel indices

    % 3. Voxels that fullfill both criteria, are low on the ranks (:,3) of both
    voxelCandidate = deltaVoxels(:,3).*diffEarlyLate(:,3); % multiply both ranks
    voxelCandidate(:,2) = 1:length(voxelCandidate); % voxel indices
    voxelCandidate = sortrows(voxelCandidate, 1, 'ascend'); % sort based on being low on both criteria (i.e., scoring high on both criteria)

    % % % Select by properties AIF: precontrast stability
    % % sdPreContrast = std(voxelTimeCurves(:,1:injectionFrame-1), [], 2);
    % % sdPreContrast(:,2) = 1:length(sdPreContrast);
    % % sdPreContrast = sortrows(sdPreContrast, 1);
    
    % Outwards in; restrict ROI? -> /Users/hjmutsaerts/ExploreASL/ExploreASL/External/SPMmodified/MapsAdded/brainCentralityMap.nii
    
    AIF = mean(voxelTimeCurves(voxelCandidate(1:nVoxels, 2), :), 1);
    save(x.P.Path_AIF, 'AIF');

    x.mutex.AddState('050_AIF');
elseif   bO; fprintf('%s\n','050_AIF has already been performed, skipping...');
end



%% -----------------------------------------------------------------------------
%% 6    T1 mapping
if ~x.mutex.HasState('060_T1mapping')

    fprintf('\n\n\n%s\n', '=========================================================================');
    fprintf('%s\n', '6. T1 mapping');
    fprintf('%s\n\n\n\n', '=========================================================================');

    % Dummy T1 map
    
    % Typical tissue T1 values at 3T:
    % 
    % White matter: ~832 ms (Wansapura 1999, PMID 10232510)
    % Gray matter: ~1331 ms (Wansapura 1999, PMID 10232510)
    % Blood: ~1664 ms arterial at Hct 0.42 (Lu 2004, PMID 15334591); ISMRM ASL consensus uses 1650 ms (Alsop 2015, PMC4190138)
    % CSF: ~4000+ ms (Rooney et al., MRM 2007;57:308, PMID 17260370)
    % If your values are outside these ranges, check your flip angles and TR.
    % 0.475*1331+0.475*832+0.05*1664 -> 1100 ms (5% CBV with blood T1, rest equally divided between WM and GM)
    
        % % T1gm = 1331;
        % % T1wm = 832;
        % % T1csf = 4000;
        % % 
        % % T1im = T1gm.*pvGM + T1wm.*pvWM + T1csf.*pvCSF;
        % % 
        % % pathT1map = fullfile(dir_dceDest, 'T1map.nii');
        % % xASL_io_SaveNifti(Path_pvGM, pathT1map, T1im);
    
    
    % Try-out T1-mapping
    % According to https://osipi.github.io/osipy/tutorials/dce-analysis/#background
    
    % combinedIM_VFA(:,:,:,1) = xASL_io_Nifti2Im(x.P.Path_despotIR);
    combinedIM_VFA(:,:,:,1) = xASL_io_Nifti2Im(x.P.Path_despotFA1);
    combinedIM_VFA(:,:,:,2) = xASL_io_Nifti2Im(x.P.Path_despotFA2);
    
    xASL_io_SaveNifti(x.P.Path_despotFA1, x.P.Path_VFA, combinedIM_VFA);
    xASL_spm_reslice(x.P.Path_rDCE4D, x.P.Path_VFA);

    x.mutex.AddState('060_T1mapping');
elseif   bO; fprintf('%s\n','060_T1mapping has already been performed, skipping...');
end
    

%% -----------------------------------------------------------------------------
%% 7    OSIPY signal->concentration->Ktrans fitting
if ~x.mutex.HasState('070_OSIPY')

    fprintf('\n\n\n%s\n', '=========================================================================');
    fprintf('%s\n', '7. OSIPY signal->concentration->Ktrans fitting');
    fprintf('%s\n\n\n\n', '=========================================================================');
    
    xASL_adm_UnzipNifti(x.P.Path_rDCE4D);
    xASL_adm_UnzipNifti(x.P.Path_rVFA);
    
    % Start Python within the venv:
    pyenv('Version', x.modules.dce.PathPython);
    insert(py.sys.path,int32(0), x.modules.dce.PythonCodePath);
    
    % -> Convert rDCE_4D NIfTI to concentration
    mod = py.importlib.import_module('Step4_Signal2Concentration');
    % mod = py.importlib.reload(mod);
    mod.nifti_signal_to_concentration(x.P.Path_rDCE4D, x.P.Path_rVFA, x.P.Path_rVFA, x.P.Path_DCE4D_Concentration); % PathT1map is a dummy, not used yet
    
    % -> Convert AIF to concentration
    load(x.P.Path_AIF);
    concentrationAIF = mod.signal_to_concentration(AIF, x.P.Path_rVFA, x.P.Path_rVFA);

    % -> Convert concentration 2 Ktrans & Vp
    mod = py.importlib.import_module('Step5_Concentration2Ktrans');
    % mod = py.importlib.reload(mod);
    mod.concentration_to_ktrans(x.P.Path_DCE4D_Concentration, x.P.Path_Ktrans, x.P.Path_Vp, concentrationAIF, x.P.Path_DCE_RSquared, x.P.Path_DCE_FitMask);

    x.mutex.AddState('070_OSIPY');
elseif   bO; fprintf('%s\n','070_OSIPY has already been performed, skipping...');
end


%% -----------------------------------------------------------------------------
%% 8    QC statistics
if ~x.mutex.HasState('080_QC')

    fprintf('\n\n\n%s\n', '=========================================================================');
    fprintf('%s\n', '8. QC statistics');
    fprintf('%s\n\n\n\n', '=========================================================================');

    x = xASL_adm_LoadX(x, PathX, true); % assume x.mat is newer than x

    % Clear any previously stored QC parameters
    if isfield(x,'Output') && isfield(x.Output,'DCE')
       x.Output = rmfield(x.Output,'DCE');
    end

    % Create masks
    pvGM = xASL_io_Nifti2Im(x.P.Path_PVgm);
    pvWM = xASL_io_Nifti2Im(x.P.Path_PVwm);
    pvCSF = xASL_io_Nifti2Im(x.P.Path_PVcsf);

    maskGM = pvGM>(pvWM+pvCSF) & (pvGM+pvWM+pvCSF)>0.1;
    maskWM = pvWM>(pvGM+pvCSF) & (pvGM+pvWM+pvCSF)>0.1;
    maskWB = (pvGM+pvWM)>pvCSF & (pvGM+pvWM+pvCSF)>0.1;

    % A. Mean, SD, SNR of DCE
    IM_rDCE4D = xASL_io_Nifti2Im(x.P.Path_rDCE4D);
    meanIM = xASL_stat_MeanNan(IM_rDCE4D, 4);
    sdIM = xASL_stat_StdNan(IM_rDCE4D, [], 4);
    % snrIM = meanIM./sdIM;
    covIM = sdIM./meanIM;
    
    xASL_io_SaveNifti(x.P.Path_rDCE4D, x.P.Path_DCE_mean, meanIM, [], 0);
    xASL_io_SaveNifti(x.P.Path_rDCE4D, x.P.Path_DCE_SD, sdIM, [], 0);
    % xASL_io_SaveNifti(x.P.Path_rDCE4D, x.P.Path_DCE_SNR, snrIM, [], 0);
    xASL_io_SaveNifti(x.P.Path_rDCE4D, x.P.Path_DCE_CoV, covIM, [], 0);


    % B. WB SNR
    meanIM = xASL_stat_MeanNan(IM_rDCE4D, 4);
    sdIM = xASL_stat_StdNan(IM_rDCE4D, [], 4);
    snrIM = meanIM./sdIM;
    
    x.Output.DCE.DCE_SNR_WB_ratio = xASL_stat_MeanNan(snrIM(maskWB));

    
    % C. GM-WM CNR
    SD = xASL_stat_StdNan(sdIM(maskWB));
    meanGM = xASL_stat_MeanNan(meanIM(maskGM));
    meanWM = xASL_stat_MeanNan(meanIM(maskWM));
    diffGM_WM = abs(meanGM-meanWM);
    
    x.Output.DCE.DCE_GMWM_CNR_ratio = diffGM_WM/SD;


    %% D. Compute Motion parameters
    % rp = realign parameters
    % FD = framewise displacement

    path_RealignParameters = fullfile(x.dir.SESSIONDIR, 'rp_DCE4D.txt');
    rp = load(path_RealignParameters, '-ascii'); % load the 3 translation and 3 rotation values
    MeanRadius = 50; % typical distance center head to cerebral cortex (Power et al., NeuroImage 2012)
    
    tx = rp(:,1); ty = rp(:,2); tz  = rp(:,3); % translations
    rx = rp(:,4); ry = rp(:,5); rz  = rp(:,6); % rotations (pitch, roll, yaw)
    
    PartTranslation = tx.^2 + ty.^2 + tz.^2;
    PartRotation = 0.2*MeanRadius^2* ((cos(rx)-1).^2 + (sin(rx)).^2 + (cos(ry)-1).^2 + (sin(ry)).^2 + (cos(rz)-1).^2 + (sin(rz)).^2);
    
    NDV = sqrt(PartTranslation + PartRotation); % net displacement vector
    NDV_FD = abs(NDV(2:end) - NDV(1:end-1));
    
    x.Output.DCE.DCE_meanFD_mm = mean(NDV_FD); % mean framewise displacement (mm)
    x.Output.DCE.DCE_maxFD_mm = max(NDV_FD); % max framewise displacement (mm)


    %% E. Compute alignment with T1w (tanimoto coefficient)
    imT1w = xASL_io_Nifti2Im(x.P.Path_PV_T1w);
    meanIM_T1 = xASL_io_Nifti2Im(x.P.Path_DCE_mean);
    
    % Normalize the image intensities
    minT1w = min(imT1w(maskWB));
    minDCE = min(meanIM_T1(maskWB));
    
    imT1w = imT1w-minT1w;
    meanIM_T1 = meanIM_T1-minDCE;
    
    imT1w = imT1w./mean(imT1w(maskWB));
    meanIM_T1 = meanIM_T1./mean(meanIM_T1(maskWB));
    
    % imshow3D([imT1w meanIM])
    
    x.Output.DCE.DCE2T1w_TC_Perc = xASL_qc_TanimotoCoeff(imT1w, meanIM_T1, maskWB, 3);


    %% G. AIF QC
    % 1. timeToPeak
    load(x.P.Path_AIF);
    [~, maxIndex] = max(AIF);
    x.Output.DCE.AIF_timeToPeak_nFrames = maxIndex - x.modules.dce.InjectionFrame;
    
    % 2. Relative peak signal
    x.Output.DCE.AIF_relativePeak_Ratio = (max(AIF) - min(AIF)) / std(AIF); % how much larger is max-min signal than SD signal


    %% H. Ktrans

    % Use R^2 as goodness-of-fit parameter
    Rsquared = xASL_io_Nifti2Im(x.P.Path_DCE_RSquared);
    x.Output.DCE.DCE_Ktrans_GM_R_squaredMean_Perc = xASL_stat_MeanNan(Rsquared(maskGM));
    x.Output.DCE.DCE_Ktrans_WM_R_squaredMean_Perc = xASL_stat_MeanNan(Rsquared(maskWM));

    Rsquared(Rsquared<0) = 0;

    imKtrans = xASL_io_Nifti2Im(x.P.Path_Ktrans);
    
    robustGMMask = maskGM & isfinite(imKtrans);
    kTransValues = imKtrans(robustGMMask);
    kTransWeights = Rsquared(robustGMMask);

    x.Output.DCE.DCE_Ktrans_GM_robustMean = sum(kTransValues.*kTransWeights)./sum(kTransWeights(:));
    x.Output.DCE.DCE_Ktrans_GM_mean = xASL_stat_MeanNan(imKtrans(maskGM));
    x.Output.DCE.DCE_Ktrans_GM_SD = xASL_stat_StdNan(imKtrans(maskGM));
    x.Output.DCE.DCE_Ktrans_GM_min = min(imKtrans(maskGM));
    x.Output.DCE.DCE_Ktrans_GM_max = max(imKtrans(maskGM));
    
    robustWMMask = maskWM & isfinite(imKtrans);
    kTransValues = imKtrans(robustWMMask);
    kTransWeights = Rsquared(robustWMMask);

    x.Output.DCE.DCE_Ktrans_WM_robustMean = sum(kTransValues.*kTransWeights)./sum(kTransWeights(:));
    x.Output.DCE.DCE_Ktrans_WM_mean = xASL_stat_MeanNan(imKtrans(maskWM));
    x.Output.DCE.DCE_Ktrans_WM_SD = xASL_stat_StdNan(imKtrans(maskWM));
    x.Output.DCE.DCE_Ktrans_WM_min = min(imKtrans(maskWM));
    x.Output.DCE.DCE_Ktrans_WM_max = max(imKtrans(maskWM));

    robustParenchymaMask = maskWB & isfinite(imKtrans);
    kTransValues = imKtrans(robustParenchymaMask);
    kTransWeights = Rsquared(robustParenchymaMask);

    x.Output.DCE.DCE_Ktrans_Parenchyma_robustMean = sum(kTransValues.*kTransWeights)./sum(kTransWeights(:));
    x.Output.DCE.DCE_Ktrans_Parenchyma_mean = xASL_stat_MeanNan(imKtrans(maskWB));
    x.Output.DCE.DCE_Ktrans_Parenchyma_SD = xASL_stat_StdNan(imKtrans(maskWB));
    x.Output.DCE.DCE_Ktrans_Parenchyma_min = min(imKtrans(maskWB));
    x.Output.DCE.DCE_Ktrans_Parenchyma_max = max(imKtrans(maskWB));

    
    %% I. Vp 
    imVp = xASL_io_Nifti2Im(x.P.Path_Vp);
    
    robustGMMask = maskGM & isfinite(imVp);
    VpValues = imVp(robustGMMask);
    VpWeights = Rsquared(robustGMMask);

    x.Output.DCE.DCE_Vp_GM_robustMean = sum(VpValues.*VpWeights)./sum(VpWeights(:));
    x.Output.DCE.DCE_Vp_GM_mean = xASL_stat_MeanNan(imVp(maskGM));
    x.Output.DCE.DCE_Vp_GM_SD = xASL_stat_StdNan(imVp(maskGM));
    x.Output.DCE.DCE_Vp_GM_min = min(imVp(maskGM));
    x.Output.DCE.DCE_Vp_GM_max = max(imVp(maskGM));
    
    robustWMMask = maskWM & isfinite(imVp);
    VpValues = imVp(robustWMMask);
    VpWeights = Rsquared(robustWMMask);

    x.Output.DCE.DCE_Vp_WM_robustMean = sum(VpValues.*VpWeights)./sum(VpWeights(:));
    x.Output.DCE.DCE_Vp_WM_mean = xASL_stat_MeanNan(imVp(maskWM));
    x.Output.DCE.DCE_Vp_WM_SD = xASL_stat_StdNan(imVp(maskWM));
    x.Output.DCE.DCE_Vp_WM_min = min(imVp(maskWM));
    x.Output.DCE.DCE_Vp_WM_max = max(imVp(maskWM));

    robustParenchymaMask = maskWB & isfinite(imVp);
    VpValues = imVp(robustParenchymaMask);
    VpWeights = Rsquared(robustParenchymaMask);

    x.Output.DCE.DCE_Vp_Parenchyma_robustMean = sum(VpValues.*VpWeights)./sum(VpWeights(:));
    x.Output.DCE.DCE_Vp_Parenchyma_mean = xASL_stat_MeanNan(imVp(maskWB));
    x.Output.DCE.DCE_Vp_Parenchyma_SD = xASL_stat_StdNan(imVp(maskWB));
    x.Output.DCE.DCE_Vp_Parenchyma_min = min(imVp(maskWB));
    x.Output.DCE.DCE_Vp_Parenchyma_max = max(imVp(maskWB));

    x.Output.DCE.DCE_nVolumes = nVolumes;

    save(PathX, 'x'); % future: do this in each xWrapper

    x.mutex.AddState('080_QC');
elseif   bO; fprintf('%s\n','080_QC has already been performed, skipping...');
end    



%% -----------------------------------------------------------------------------
%% 9    Visualization
if ~x.mutex.HasState('090_Visualization')

    fprintf('\n\n\n%s\n', '=========================================================================');
    fprintf('%s\n', '9. Visualization');
    fprintf('%s\n\n\n\n', '=========================================================================');

    close all % close all Figures to avoid capturing & saving the wrong Figure


    %% AIF
    if usejava('jvm')
        fig = figure('Visible','off');
        plot(AIF, 'r-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 30, 'MarkerEdgeColor', 'b');
    
        xlabel('DCE frame (n)');
        ylabel('DCE value (a.u.)');
        title('Arterial Input Function');
    
        pathAIFcheck = fullfile(x.D.DCECheckDir, ['AIF_' SubjVisit '.jpg']);
        fprintf('Saving AIF plot to %s\n', pathAIFcheck);
        xASL_adm_CreateDir(x.D.DCECheckDir);
        saveas(fig, pathAIFcheck, 'jpg');
        close all;
    else
        fprintf('Skipping motion vs exclusion overview, missing JVM\n');
    end

    % visualization debugging
    % visIM = IM(:,:,:,1);
    % visIM = visIM./max(visIM(:));
    % 
    % imshow3D([visIM+maskCSF+imConfluenceMask visIM+imConfluenceMaskLocal])


    %% Resample images to standard space
    InputPaths = {x.P.Path_DCE_mean, x.P.Path_DCE_CoV, x.P.Path_Ktrans, x.P.Path_Vp};
    OutputPaths = {x.P.Pop_Path_DCE_mean, x.P.Pop_Path_DCE_CoV, x.P.Pop_Path_DCE_Ktrans, x.P.Pop_Path_DCE_Vp};
    
    xASL_spm_deformations(x, InputPaths, OutputPaths, [], [], [], x.P.Path_y_DCE);

    %% Visualization
    xASL_vis_CreateVisualFig(x, {x.P.Pop_Path_DCE_mean}, x.D.DCECheckDir, [], 'DCE_mean_');
    xASL_vis_CreateVisualFig(x, {x.P.Pop_Path_DCE_CoV}, x.D.DCECheckDir, [], 'DCE_CoV_', [], [], {x.S.masks.skull});
    xASL_vis_CreateVisualFig(x, {x.P.Pop_Path_DCE_Ktrans}, x.D.DCECheckDir, [], 'DCE_Ktrans_', [], [], {x.S.masks.skull});
    xASL_vis_CreateVisualFig(x, {x.P.Pop_Path_DCE_Vp}, x.D.DCECheckDir, [], 'DCE_Vp_', [], [], {x.S.masks.skull});

    x.mutex.AddState('090_Visualization');
elseif   bO; fprintf('%s\n','090_Visualization has already been performed, skipping...');
end


%% Householding
% if x.settings.DELETETEMP
% %    xASL_delete(.....);
% %    xASL_delete(.....);
% end



%% -----------------------------------------------------------------------------
%% 999 Ready
x.mutex.AddState('999_ready');

x.mutex.Unlock();
x.result  = true;
result    = true;

end