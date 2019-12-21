function [EddyQC] = xASL_qc_dwi_EddyMotion (subj_directory)
% xASL_qc_dwi_Eddy compute QC parameters derived by Eddy Correction in Diffusion weighted images 
%
% FORMAT: [EddyQC] = xASL_qc_dwi_EddyMotion (dwi_directory)
%
% INPUT:
%   subjdir  	     - Path to the subject directory 
%                          
%  
%
% OUTPUT:
% EddyQC.motion.Average_x_translation_mm      - Average Volume to Volume motion on the X asxis (in mm)                                                   
% EddyQC.motion.Average_y_translation_mm      - Average Volume to Volume motion on the Y asxis (in mm)
% EddyQC.motion.Average_z_translation_mm      - Average Volume to Volume motion on the Z asxis (in mm)
% EddyQC.motion.Average_x_rotation_deg        - Average Volume to Volume % rotation on the X axis (in degree)
% EddyQC.motion.Average_y_rotation_deg        - Average Volume to Volume % rotation on the Y axis (in degree)
% EddyQC.motion.Average_z_rotation_deg        - Average Volume to Volume % rotation on the X axis (in degree)
% EddyQC.motion.Average_abs_motion_mm         - Average motion of each volume compared to the first 
% EddyQC.motion.Average_rel_motion_mm         - Average motion of each volume compared to the previous
% EddyQC.Outlier_Perc                         - Percentage of outliers slices as computed by FSL Eddy Current
% EddyQC.Induced_Distortion                   - Standard deviation of EddyCurrent induced distortion 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Description : This function computes  the QC parameters derived from the
% FSL Eddy Current correction
%
% 
% EXAMPLE: EddyQC = xASL_qc_dwi_Eddy ('/examplepath/DWIfolder/');

%set DWI directory
dwi_directory = fullfile(subj_directory, 'dwi');

%%Set Paths to NIfTIs and load images
C1Path= fullfile (subj_directory, 'c1T1.nii');
C2Path= fullfile (subj_directory, 'c2T1.nii');
C3Path= fullfile (subj_directory, 'c3T1.nii');
TopUp_Path = fullfile(dwi_directory, 'TopUp_fieldcoef.nii');


GMim = (xASL_io_Nifti2Im(C1Path)>0.5);
WMim = (xASL_io_Nifti2Im(C2Path)>0.5);
CSFim = (xASL_io_Nifti2Im(C3Path)>0.5);
TopUp = xASL_io_Nifti2Im(TopUp_Path);

% Create Brain Mask
BrainMask = (GMim+WMim+CSFim);
xASL_io_CreateNifti(fullfile(dwi_directory, 'tmp_brainmask.nii'), BrainMask);

%Set Path to  Eddy Current output files
Motion_Pars_file = fullfile(dwi_directory, 'Eddy.eddy_parameters');
Absolute_Relative_Motion_file = fullfile(dwi_directory, 'Eddy.eddy_movement_rms');
Outlier_file = fullfile(dwi_directory, 'Eddy.eddy_outlier_map');
%we do not have this for now, we need to run eddy current with -mporder 
%This can only be done when running EddyCurrent with GPU acceleration, 
%because of the computation intensivity of this scheme.
%Motion_over_time = fullfile(dwi_directory,'Eddy.eddy_movement_over_time' )

% Reslice Brain Mask into dwi space and extract brain from topup for Eddy
% Current Susceptibility induced distortion computation 
xASL_spm_reslice(TopUp_Path, fullfile(dwi_directory, 'tmp_brainmask.nii'),[], [], [], fullfile(dwi_directory, 'tmp_brainmask_resliced'), 0)
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





