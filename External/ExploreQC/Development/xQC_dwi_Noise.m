function [Noise_mean , Noise_SD] = xQC_dwi_Noise(SubjDir, filename, bSaveDenoised , bSaveNoise)
% Computes Noise on Diffusion Image
%
% FORMAT: [Noise_mean , Noise_SD] = xQC_dwiNoise(SubjDir, filename, bSaveDenoised , bSaveNoise)
% INPUT:
%   SubjDir         - Path to subject directory (REQUIRED)
%   Filename        - Name of diffusion nifti file (REQUIRED)
%   bSaveDenoised   - 1 will save the denoised image 
%   bSaveNoise      - 1 will save the noise image 
% 
% OUTPUT:
%   Noise_mean      - Mean Noise
%   Noise_SD        - Standard Deviation of noise map
% EXAMPLE:
% 
% [Noise_mean , Noise_SD] = xQC_dwiNoise('Path\to\Subj', 'dwi_run-1.nii', 1, 0)

% DESCRIPTION:
% This funtion uses the MPdenoising toolbox (Veraart, J. et al., 2016) to estimate a noise map 
% of diffusion mri and extract mean and standard eviation value
% REFERENCES: 
%      Veraart, J.; Fieremans, E. & Novikov, D.S. Diffusion MRI noise mapping
%      using random matrix theory Magn. Res. Med., 2016, early view, doi:
%      10.1002/mrm.26059


dwiPath = fullfile(SubjDir, 'dwi', filename)

dwiiM= xASL_io_Nifti2Im(dwiPath);


% Make a mask with the C2T1 image, this makes noise estimation faster

PathC2 = fullfile(SubjDir, 'c2T1.nii');

if ~exist(PathC2)
    warning (['No C2 file available in ' SubjDir ' calculating noise on the whole brain'])
    
    [Signal, NoiseMap] = MPdenoising(dwiiM);
else
    
    %Reslice C2 into dwi space
    xASL_spm_reslice(dwiPath, PathC2, [], [], [], fullfile(SubjDir, 'WM_tmp_reslicedmask'), 0);
    
    C2im = xASL_io_Nifti2Im(fullfile(SubjDir, 'WM_tmp_reslicedmask.nii'));
    C2im = C2im>0.5;
    
    xASL_delete(fullfile(SubjDir, 'WM_tmp_reslicedmask.nii'))
    
    [Signal, NoiseMap] = MPdenoising(dwiiM, C2im);
    
    Noise_mean = xASL_stat_MeanNan(NoiseMap(:))
    Noise_SD = xASL_stat_StdNan(NoiseMap(:))
end

if bSaveDenoised
    DenoisedPath = fullfile(SubjDir, 'dwi', ['denoised_' filename ])
    xASL_io_SaveNifti(dwiPath,DenoisedPath, Signal)
end

if bSaveNoise
    NoisePath = fullfile(SubjDir, 'dwi', 'Noise_dwi.nii');
    xASL_io_SaveNifti(dwiPath,NoisePath, NoiseMap)
end



end
