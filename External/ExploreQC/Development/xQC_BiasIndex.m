function BiasIndex = xQC_BiasIndex(SubjDir, SPMdir, ScanType)
% Computes BiasIndex on FLAIR or T1w images
%
% FORMAT: BiasIndex = xQC_BiasIndex(SubjDir, SPMdir, ScanType)
% INPUT:
%   SubjDir     - Path to subject directory (REQUIRED)
%   SPMdir      - Path to directory where SPM is stored (REQUIRED)
%   ScanType    - Either 'FLAIR' or 'T1w' (REQUIRED)
% 
% OUTPUT:
%   BiasIndex   - Mean Bias Index 
% 
% EXAMPLE:
% 
% BiasIndex = xQC_BiasIndex('Path\to\Subj' , 'path\to\SPM', 'FLAIR')
% 
% DESCRIPTION:
%  This function computes the Bias Index (BI) as the mean signal
%  inhomogeneity, i.e. the mean difference between the bias field corrected
%  image and the original one. (Peltonen, J. et al., 2018) 



if strcmp('FLAIR', ScanType)
    PathIn = fullfile(SubjDir, 'FLAIR.nii');
    PathOut = fullfile(SubjDir,'FLAIRcorr.nii');
elseif strcmp('T1', ScanType)
    PathIn = fullfile(SubjDir, 'T1.nii');
    PathOut = fullfile(SubjDir,'T1corr.nii');
else 
    error('Please Specify ScanType')
end 


xASL_spm_BiasfieldCorrection(PathIn, SPMdir, [], [], PathOut);

BiasfieldCorrectedNIfTI = xASL_io_Nifti2Im(PathOut);
UncorrectedNIfTI = xASL_io_Nifti2Im(PathIn);


SubtractedIm = BiasfieldCorrectedNIfTI-UncorrectedNIfTI;

BiasIndex = xASL_stat_MeanNan(SubtractedIm(:));

end 
