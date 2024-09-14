% Copyright 2015-2024 ExploreASL (Works In Progress code)
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

function [x] = xASL_spm_GLMmodel(x)
%xASL_spm_GLMmodel ExploreASL wrapper for SPM GLM - statistical model estimation
fprintf('%s\n','Running GLM model')
% This part crashes if the implicitly created mask contains no voxels 
%% Admin
load(x.S.SPMmat1);
save(x.S.SPMmat,'SPM');
%% Delete previous files
xASL_adm_DeleteFileList(x.S.SPMdir,'^beta_.*\.(nii|nii\.gz)$');
xASL_adm_DeleteFileList(x.S.SPMdir,'^mask\.(nii|nii\.gz)$');
xASL_adm_DeleteFileList(x.S.SPMdir,'^ResMS\.(nii|nii\.gz)$');
xASL_adm_DeleteFileList(x.S.SPMdir,'^RPV\.(nii|nii\.gz)$');   
%% Run model estimation
matlabbatch{1}.spm.stats.fmri_est.spmmat            = { x.S.SPMmat };
matlabbatch{1}.spm.stats.fmri_est.write_residuals   = 0; % don't write residual maps
matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1; % restrict maximum likelihood
spm_jobman('run',matlabbatch);
% Save SPM2
load(x.S.SPMmat);
save(x.S.SPMmat2,'SPM');
end
