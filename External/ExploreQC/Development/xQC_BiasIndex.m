function [BI_mean, BI_SD] = xQC_BiasIndex(SubjDir, SPMdir, ScanType)
% Computes BiasIndex on FLAIR or T1w images
%
% FORMAT: BiasIndex = xQC_BiasIndex(SubjDir, SPMdir, ScanType)
% INPUT:
%   SubjDir     - Path to subject directory (REQUIRED)
%   SPMdir      - Path to directory where SPM is stored (REQUIRED)
%   ScanType    - Either 'FLAIR' or 'T1w' (REQUIRED)
% 
% OUTPUT:
%   BI_mean   - Mean Bias Index 
%   BI_SD     - Bias Index Standard Deviation
% EXAMPLE:
% 
% BiasIndex = xQC_BiasIndex('Path\to\Subj' , 'path\to\SPM', 'FLAIR')
% 
% DESCRIPTION:
%  This function computes the Bias Index (BI) as the mean signal
%  inhomogeneity, i.e. the mean difference between the bias field corrected
%  image and the original one. (Peltonen, J. et al., 2018) 



if strcmp('FLAIR', ScanType)
    NiftiPath = fullfile(SubjDir, 'FLAIR.nii');
elseif strcmp('T1', ScanType)
    NiftiPath = fullfile(SubjDir, 'T1.nii');
else 
    error('Please Specify ScanType')
end 


matlabbatch = spm_segm_BI(NiftiPath, SPMdir);

spm_jobman('run', matlabbatch);

BIpath = fullfile(SubjDir, ['BiasField_' ScanType '.nii']);
BiasField = xASL_io_Nifti2Im(BIpath);

BI_mean= xASL_stat_MeanNan(BiasField(:));
BI_SD= xASL_stat_StdNan(BiasField(:));

end 




function matlabbatch = spm_segm_BI(NiftiPath, SPMdir)

matlabbatch{1}.spm.spatial.preproc.channel.vols = {[NiftiPath ',1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(SPMdir, 'tpm', 'TPM.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(SPMdir, 'tpm', 'TPM.nii,2')};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(SPMdir, 'tpm', 'TPM.nii,3')};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(SPMdir, 'tpm', 'TPM.nii,4')};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(SPMdir, 'tpm', 'TPM.nii,5')};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(SPMdir, 'tpm', 'TPM.nii,6')};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 0;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 0;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = '';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 9;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];


end 
