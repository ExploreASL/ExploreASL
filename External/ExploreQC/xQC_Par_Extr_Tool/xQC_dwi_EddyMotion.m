function [EddyQC] = xQC_dwi_EddyMotion (dwi_directory, subject, GMpath, WMpath, CSFpath, TopUp_Path )
% xQC_dwi_Eddy compute QC parameters derived by Eddy Correction in Diffusion weighted images
%
% FORMAT: [EddyQC] = xASL_qc_dwi_EddyMotion (dwi_directory)
%
% INPUT:
%   dwi_directory                            - Path to the subject DWI directory
%   subject                                   - subject ID
% OUTPUT:
% EddyQC.motion.Average_x_translation_mm      - Average Volume to Volume motion on the X axis (in mm)
% EddyQC.motion.Average_y_translation_mm      - Average Volume to Volume motion on the Y axis (in mm)
% EddyQC.motion.Average_z_translation_mm      - Average Volume to Volume motion on the Z axis (in mm)
% EddyQC.motion.Average_x_rotation_deg        - Average Volume to Volume % rotation on the X axis (in degrees)
% EddyQC.motion.Average_y_rotation_deg        - Average Volume to Volume % rotation on the Y axis (in degrees)
% EddyQC.motion.Average_z_rotation_deg        - Average Volume to Volume % rotation on the X axis (in degrees)
% EddyQC.motion.Average_abs_motion_mm         - Average motion of each volume compared to the first
% EddyQC.motion.Average_rel_motion_mm         - Average motion of each volume compared to the previous
% EddyQC.Outlier_Perc                         - Percentage of outlier slices as computed by FSL Eddy Current
% EddyQC.Induced_Distortion                   - Standard deviation of EddyCurrent induced distortion
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Description : This function computes QC parameters from the FSL Eddy
%               Current correction parameters
%
% EXAMPLE: EddyQC = xASL_qc_dwi_Eddy ('/examplepath/DWIfolder/');

% Set DWI directory

%% Set Paths to NIfTIs and load images
GMim = xASL_io_Nifti2Im(GMpath)>0.5;
WMim = xASL_io_Nifti2Im(WMpath)>0.5;
CSFim = xASL_io_Nifti2Im(CSFpath)>0.5;
TopUp = xASL_io_Nifti2Im(TopUp_Path);

% Create Brain Mask
BrainMask = (GMim+WMim+CSFim);
xASL_io_CreateNifti(fullfile(dwi_directory, 'tmp_brainmask.nii'), BrainMask);

%Set Path to  Eddy Current output files
Motion_Pars_file = fullfile(dwi_directory, 'Eddy.eddy_parameters');
Absolute_Relative_Motion_file = fullfile(dwi_directory, 'Eddy.eddy_movement_rms');
Outlier_file = fullfile(dwi_directory, 'Eddy.eddy_outlier_map');

if ~exist(Motion_Pars_file) || ~exist(Absolute_Relative_Motion_file) || ~exist(Outlier_file)
    
    fprintf(['Eddy output not found for ' subject 'skipping' ])
    EddyQC = struct();
    
    EddyQC.motion.Average_x_translation_mm = NaN;
    EddyQC.motion.Average_y_translation_mm = NaN;
    EddyQC.motion.Average_z_translation_mm = NaN;
    
    EddyQC.motion.Average_x_rotation_deg = NaN;
    EddyQC.motion.Average_y_rotation_deg = NaN;
    EddyQC.motion.Average_z_rotation_deg = NaN;
    
    EddyQC.motion.Average_abs_motion_mm = NaN;
    EddyQC.motion.Average_rel_motion_mm = NaN;
    
    EddyQC.Outlier_Perc = NaN;
    
    EddyQC.Induced_Distortion = NaN;
    EddyQC.Susc_Induced_Distortion = NaN;
else
    
    %we do not have this for now, we need to run eddy current with -mporder
    %This can only be done when running EddyCurrent with GPU acceleration,
    %because of the computation intensivity of this scheme.
    %Motion_over_time = fullfile(dwi_directory,'Eddy.eddy_movement_over_time' )
    
    % Reslice Brain Mask into dwi space and extract topup for Eddy
    % Current Susceptibility induced distortion parameters for the whole brain
    xASL_spm_reslice(TopUp_Path, fullfile(dwi_directory, 'tmp_brainmask.nii'),[], [], [], fullfile(dwi_directory, 'tmp_brainmask_resliced'), 0);
    BrainMask = xASL_io_Nifti2Im(fullfile(dwi_directory, 'tmp_brainmask_resliced.nii'));
    BrainMask = BrainMask>0.5;
    
    xASL_delete(fullfile(dwi_directory, 'tmp_brainmask.nii'));
    xASL_delete(fullfile(dwi_directory, 'tmp_brainmask_resliced.nii'));
    
    TopUpBrain = TopUp(BrainMask);
    
    % Import data from Eddy Current outputs
    Total_Motion = importdata(Absolute_Relative_Motion_file);
    Motion_xyz = importdata(Motion_Pars_file);
    Induced_distortion_coeff = Motion_xyz(:, 7:9);
    Outlier_count = importdata(Outlier_file);
    Outlier_count = str2double(Outlier_count.textdata(2:end,:));
    
    % Compute QC measures
    EddyQC = struct();
    
    EddyQC.motion.Average_x_translation_mm = xASL_stat_MeanNan(Motion_xyz(2:end, 1));
    EddyQC.motion.Average_y_translation_mm = xASL_stat_MeanNan(Motion_xyz(2:end, 2));
    EddyQC.motion.Average_z_translation_mm = xASL_stat_MeanNan(Motion_xyz(2:end, 3));
    
    EddyQC.motion.Average_x_rotation_deg = xASL_stat_MeanNan(Motion_xyz(2:end, 4));
    EddyQC.motion.Average_y_rotation_deg = xASL_stat_MeanNan(Motion_xyz(2:end, 5));
    EddyQC.motion.Average_z_rotation_deg = xASL_stat_MeanNan(Motion_xyz(2:end, 6));
    
    EddyQC.motion.Average_abs_motion_mm = mean(Total_Motion(2:end, 1));
    EddyQC.motion.Average_rel_motion_mm = mean(Total_Motion(2:end, 2));
    
    EddyQC.Outlier_Perc = (xASL_stat_SumNan(Outlier_count(:))/numel(Outlier_count))*100;
    
    EddyQC.Induced_Distortion = xASL_stat_MeanNan(xASL_stat_StdNan(Induced_distortion_coeff));
    EddyQC.Susc_Induced_Distortion = xASL_stat_StdNan(TopUpBrain);
    
    
end