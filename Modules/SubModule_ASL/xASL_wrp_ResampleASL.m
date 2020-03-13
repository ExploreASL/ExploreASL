function xASL_wrp_ResampleASL(x)
%xASL_wrp_ResampleASL Submodule of ExploreASL ASL Module, that reslices native
%space images to standard space
%
% FORMAT: xASL_wrp_ResampleASL(x)
%
% INPUT:
%   x  - structure containing fields with all information required to run this submodule (REQUIRED)
%
% OUTPUT:
% OUTPUT FILES: NIfTIs in standard space
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule resamples native space NIfTIs to standard space, using the deformation fields computed in the structural module
%              after smoothing these transformation fields to the ASL resolution.
%              The applied interpolation is a combination of all transformations (e.g. motion correction, registration to
%              T1w, and transformation of T1w to standard space. This submodule performs the following steps:
%
%   1    Create slice gradient image for quantification reference, in case of 2D ASL
%   2    Reslice ASL time series to MNI space (currently 1.5 mm^3)
%   3    Create mean control image, if available
%         This also applies a bilateral filter, if requested
%         If x.M0 is set as UseControlAsM0, this mean control NIfTI will be
%         copied to an M0 NIfTI (and processed in the M0 submodule)
%   4    Perform pair-wise subtraction (to create PWI.nii), in native space
%   5    Same in standard space
%   6    Save PWI NIfTI & time-series-related maps (SD, SNR)
%
% EXAMPLE: xASL_wrp_ResampleASL(x);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL



%% ------------------------------------------------------------------------------------------
%% 0) Administration
if strcmp(x.M0,'no_background_suppression')
    x.M0 = 'UseControlAsM0'; % backward compatibility
end

% Use either original or motion estimated ASL4D
% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
if ~xASL_exist(x.P.Path_despiked_ASL4D,'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end
tempnii = xASL_io_ReadNifti(x.P.Path_despiked_ASL4D);
nVolumes = double(tempnii.hdr.dim(5));

if ~isfield(x,'SavePWI4D')
    x.SavePWI4D   = 0;
end

xASL_wrp_CreateASLDeformationField(x); % make sure we have the deformation field in ASL resolution

%% TopUp files
PathB0 = fullfile(x.SESSIONDIR ,'B0.nii');
PathUnwarped = fullfile(x.SESSIONDIR ,'Unwarped.nii');
PathPopB0 = fullfile(x.D.PopDir, ['rASL_B0_' x.P.SubjectID '_' x.P.SessionID '.nii']);
PathPopUnwarped = fullfile(x.D.PopDir, ['rASL_Unwarped_' x.P.SubjectID '_' x.P.SessionID '.nii']);

InputPaths = {PathB0, PathUnwarped};
OutputPaths = {PathPopB0, PathPopUnwarped};
xASL_spm_deformations(x, InputPaths, OutputPaths, [], [], [], x.P.Path_y_ASL);

%% ------------------------------------------------------------------------------------------
%% 1    Create slice gradient image for quantification reference, in case of 2D ASL
xASL_wrp_CreateSliceGradient(x);


%% ------------------------------------------------------------------------------------------
%% 2    Reslice ASL time series to MNI space (currently 1.5 mm^3)
fprintf('%s\n','Convert ASL time series to single precision');
xASL_io_SaveNifti(x.P.Path_despiked_ASL4D, x.P.Path_temp_despiked_ASL4D, tempnii.dat(:,:,:,:),32,0);

% Convertion to single precision for precious data, to avoid
% digitization artifacts in the spatial processing
% Plus this creates a temporary copy to not touch the original ASL file

xASL_spm_deformations(x, x.P.Path_temp_despiked_ASL4D,x.P.Path_rtemp_despiked_ASL4D,4,[],x.P.Path_mean_PWI_Clipped_sn_mat,x.P.Path_y_ASL);


%% ------------------------------------------------------------------------------------------
%% 3    Create mean control image, if available
%       This also applies a bilateral filter, if requested
%       If x.M0 is set as UseControlAsM0, this mean control NIfTI will be
%       copied to an M0 NIfTI (and processed in the M0 submodule)
if  nVolumes>1 % this is when a mean control image can be created

    %% First create SD and SNR images for control series, to use for QC

    for iS=1:nVolumes
        matlabbatch{1}.spm.spatial.realign.write.data{iS,1} = [x.P.Path_despiked_ASL4D ',' num2str(iS)];
    end

    matlabbatch{1}.spm.spatial.realign.write.roptions.which     = [2 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.interp    = 2;
    matlabbatch{1}.spm.spatial.realign.write.roptions.wrap      = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.mask      = 1;
    matlabbatch{1}.spm.spatial.realign.write.roptions.prefix    = 'r';

    spm_jobman('run',matlabbatch); % this applies the motion correction in native space

    if ~xASL_exist(x.P.Path_rdespiked_ASL4D,'file')
        [Fpath, Ffile, Fext]        = xASL_fileparts(x.P.Path_despiked_ASL4D);
        x.P.Path_rdespiked_ASL4D    = fullfile(Fpath,['r' Ffile Fext]);
    end


    %% ------------------------------------------------------------------------------------------
    %% Bilateral filter to remove smooth artifacts
    if ~isdeployed % skip this part for compilation, to avoid DIP image issues
        if  x.BILAT_FILTER>0 && nVolumes>9
            volIM                       = xASL_io_ReadNifti(x.P.Path_rtemp_despiked_ASL4D); % load resliced time-series
            VoxelSize                   = double(volIM.hdr.pixdim(2:4));
            volIM                       = single(volIM.dat(:,:,:,:,:,:,:,:));

            mask                        = x.skull; % get standard space mask
            mask(isnan(mean(volIM,4)))  = 0; % remove outside FoV voxels

            ovol                        = xASL_wrp_Filter(volIM, mask, VoxelSize, x); % run filter

            xASL_io_SaveNifti( x.P.Path_rtemp_despiked_ASL4D, x.P.Path_rtemp_despiked_ASL4D, ovol,32,0 ); % store in file
        end
    end

    ControlIm = xASL_quant_GetControlLabelOrder(xASL_io_Nifti2Im(x.P.Path_rdespiked_ASL4D));
    IM_mean = xASL_stat_MeanNan(ControlIm,4);
    xASL_io_SaveNifti(x.P.Path_rdespiked_ASL4D, x.P.Path_mean_control, IM_mean, [], 0);

    InputFiles  = {x.P.Path_mean_control};
    OutputFiles = {x.P.Pop_Path_mean_control};
    xASL_spm_deformations(x, InputFiles, OutputFiles, [], [], x.P.Path_mean_PWI_Clipped_sn_mat, x.P.Path_y_ASL);

    xASL_adm_DeleteFilePair(x.P.Path_SD_control, 'mat');
    xASL_adm_DeleteFilePair(x.P.Path_SNR_control, 'mat');
    xASL_adm_DeleteFilePair(x.P.Path_rdespiked_ASL4D, 'mat');

    % Visual check of M0-pGM registration for masking
    xASL_im_CreateVisualFig(x, {x.P.Pop_Path_mean_control x.P.Pop_Path_rc1T1}, x.D.M0regASLdir,[0.5 0.2]);

    if strcmp(x.M0,'UseControlAsM0') && nVolumes==1
        warning('Couldnt create mean control image to be used as M0, timeseries missing');
    end

    if  strcmp(x.M0,'UseControlAsM0')
        % if there is no background suppression, we use the mean control
        % image as M0 image, which has perfect registration

        if  xASL_exist(x.P.Path_M0,'file') && ~xASL_exist(x.P.Path_M0_backup,'file')
            % Backup M0, if not already performed in previous run
            xASL_Copy(x.P.Path_M0,x.P.Path_M0_backup);
        end
        if  exist(x.P.Path_M0_parms_mat,'file') && ~exist(x.P.Path_M0_backup_parms_mat,'file')
            % Same for the parms file
            xASL_Copy(x.P.Path_M0_parms_mat,x.P.Path_M0_backup_parms_mat);
        end

        xASL_Copy(x.P.Path_mean_control, x.P.Path_M0,1); % overwrite M0 in case of previous processing
        if exist(x.P.Path_ASL4D_parms_mat,'file')
            xASL_Copy(x.P.Path_ASL4D_parms_mat, x.P.Path_M0_parms_mat,1);
        end
        if exist(x.P.Path_ASL4D_json,'file')
            xASL_Copy(x.P.Path_ASL4D_json, x.P.Path_M0_json,1);
        end
        fprintf('%s\n','Copying mean control image, for use as M0');
    end

else
    fprintf('%s\n','Creation mean control image was skipped, because there was only 1 frame');
end


%% ------------------------------------------------------------------------------------------
%% 4    Pair-wise subtraction native space
% Load ASL time series (after being pre-processed)
ASL_im = xASL_io_Nifti2Im(x.P.Path_temp_despiked_ASL4D); % Load time-series nifti
dim4 = size(ASL_im,4);

if dim4>1 && round(dim4/2)==dim4/2 && dim4==nVolumes
    % if has an even time series, same length as before (not changed by filter)
    % do a paired subtraction
    [ControlIm, LabelIm] = xASL_quant_GetControlLabelOrder(ASL_im);
    ASL_im = ControlIm - LabelIm;
    xASL_io_SaveNifti(x.P.Path_temp_despiked_ASL4D, x.P.Path_PWI, xASL_stat_MeanNan(ASL_im, 4), [], false);

    if x.SavePWI4D % option to store subtracted time-series (PWI4D)
        xASL_io_SaveNifti(x.P.Path_temp_despiked_ASL4D, x.P.Path_PWI4D, ASL_im, 32, false);
        xASL_spm_reslice(x.P.Path_PWI, x.P.Path_PWI4D, [], [], x.Quality, x.P.Path_PWI4D);
    end
elseif dim4 == 1
    xASL_io_SaveNifti(x.P.Path_temp_despiked_ASL4D, x.P.Path_PWI, ASL_im, [], false);
end


%% ------------------------------------------------------------------------------------------
%% 5    Pair-wise subtraction standard space
% Load ASL time series (after being pre-processed)
ASL_im = xASL_io_Nifti2Im(x.P.Path_rtemp_despiked_ASL4D); % Load time-series nifti
dim4 = size(ASL_im,4);

if dim4>1 && round(dim4/2)==dim4/2 && dim4==nVolumes
    % if has an even time series, same length as before (not changed by filter)
    % do a paired subtraction
    [ControlIm, LabelIm] = xASL_quant_GetControlLabelOrder(ASL_im);
    ASL_im = ControlIm - LabelIm;

    if x.SavePWI4D % option to store subtracted time-series (PWI4D)
        xASL_io_SaveNifti(x.P.Path_rtemp_despiked_ASL4D, x.P.Pop_Path_PWI4D, ASL_im,32, false);
    end
end


%% ------------------------------------------------------------------------------------------
% 6    Save PWI NIfTI & time-series-related maps (SD, SNR)
PWI = xASL_stat_MeanNan(ASL_im,4);
MaskIM = xASL_io_Nifti2Im(fullfile(x.D.MapsSPMmodifiedDir, 'rgrey.nii'));
MaskIM = MaskIM>(0.7*max(MaskIM(:)));

CoV0 = xASL_stat_ComputeSpatialCoV(PWI,MaskIM)*100;

if  size(ASL_im,4)>1
    % Calculate two halves, for reproducibility (wsCV) CBF & spatial CoV
    HalfVol                         = round(size(ASL_im,4)/2);
    Part1                           = xASL_stat_MeanNan(ASL_im(:,:,:,1:HalfVol),4);
    Part2                           = xASL_stat_MeanNan(ASL_im(:,:,:,HalfVol+1:end),4);

    MaskIM                          = MaskIM & isfinite(Part1) & isfinite(Part2);
    mean1                           = mean(Part1(MaskIM));
    std1                            = std(Part1(MaskIM));
    CoV1                            = std1/mean1;
    mean2                           = mean(Part2(MaskIM));
    std2                            = std(Part2(MaskIM));
    CoV2                            = std2/mean2;

    wsCV_mean                       = 100*abs(mean1-mean2)/(0.5*(mean1+mean2));
    wsCV_CoV                        = 100*abs(CoV1-CoV2)/(0.5*(CoV1+CoV2));

    fprintf('%s\n',['pGM>0.7: wsCV mean CBF = ' num2str(wsCV_mean,3) '%, wsCV spatial CoV = ' num2str(wsCV_CoV,3) '%'])

    new_SD                          = xASL_stat_StdNan( ASL_im,0,4);
    new_SNR                         = PWI./new_SD;
    new_SNR(new_SNR<0)              = 0; % clip @ zero

    xASL_io_SaveNifti(x.P.Path_rtemp_despiked_ASL4D,fullfile(x.D.PopDir,['PWI_part1_'  x.P.SubjectID '_' x.P.SessionID '.nii']),Part1,32);
    xASL_io_SaveNifti(x.P.Path_rtemp_despiked_ASL4D,fullfile(x.D.PopDir,['PWI_part2_'  x.P.SubjectID '_' x.P.SessionID '.nii']),Part2,32);

    xASL_io_SaveNifti(x.P.Path_rtemp_despiked_ASL4D,x.P.Pop_Path_SD ,new_SD ,32,0);
    xASL_io_SaveNifti(x.P.Path_rtemp_despiked_ASL4D,x.P.Pop_Path_SNR,new_SNR,32,0);
    fprintf('%s\n','Standard space timeseries-related images saved (SD & SNR), & part1 & part2 for reproducibility');
else
    fprintf('%s\n',['Standard space SD & SNR maps were not created because of only ' num2str(size(ASL_im,4)) ' frame(s)']);
end

xASL_io_SaveNifti(x.P.Path_rtemp_despiked_ASL4D, x.P.Pop_Path_PWI ,PWI,32, 0); % single precision for better precision

if  x.DELETETEMP
    xASL_delete(x.P.Path_rtemp_despiked_ASL4D);
    xASL_adm_DeleteFilePair( x.P.Path_temp_despiked_ASL4D, 'mat');
else
    xASL_io_SaveNifti(x.P.Path_rtemp_despiked_ASL4D, x.P.Path_rtemp_despiked_ASL4D, ASL_im, 32);
end


fprintf('%s\n','Standard space PWI maps saved');
fprintf('%s\n',['Spatial CoV pGM>0.7 = ' num2str(CoV0,3) '%']);

end
