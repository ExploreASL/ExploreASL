function xASL_spm_affine(srcPath, refPath, fwhmSrc, fwhmRef)
% ExploreASL wrapper for SPM affine registration function (a.k.a. 'old normalize' but without DCT)
%
% FORMAT: xASL_spm_affine(srcPath, refPath, fwhmSrc, fwhmRef)
% 
% INPUT:
%   refPath - path to reference space (NifTI image) you want to register the source image to (REQUIRED)
%   srcPath - path to source image (NifTI image) you want to register (REQUIRED)
%   fwhmSrc - Gaussian smoothing to be applied to the source image before estimating the affine registration, in FWHM (mm) (REQUIRED)
%   fwhmRef - Gaussian smoothing to be applied to the reference image before estimating the affine registration, in FWHM (mm) (REQUIRED)
%
% OUTPUT: n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This SPM wrapper runs SPM's old normalize-estimate function, which calculates the affine transformation (i.e. linear + zooming and shearing) that is required to
% align the source image with the reference image. Note that the old normalize SPM function by default also estimates a low-degree Discrete
% Cosine Transform (DCT) but this is disabled in this wrapper to have a simple affine transformation. Also note that this affine transformation uses a
% correlation cost function, hence it requires the source and reference images to have similar contrasts.
% As opposed to other SPM estimation functions (e.g. xASL_spm_coreg) this
% function does not change the orientation header, not does it allow to change those of others (e.g. the OtherList in xASL_spm_coreg). Instead,
% it stores its estimated affine transformation as orientation difference matrix in a file with the same path but _sn.mat extension.
% For the provided smoothing FWHM, note that smoothnesses combine with Pythagoras' rule (i.e. square summing)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_spm_affine('/MyStudy/Subject1/T1.nii.gz', '/MyStudy/Subject1/mean_control.nii', 5, 5);
% __________________________________
% Copyright 2015-2019 ExploreASL
% 2015-2019 HJ
   
% Check parameters
if nargin < 4
	error('xASL_spm_affine: Requires 4 input parameters.');
end

if ~xASL_exist(srcPath) || ~xASL_exist(refPath)
	error('xASL_spm_affine: Cannot find input images.');
end

% Unzip the input files 
srcPath = xASL_adm_UnzipNifti(srcPath);
refPath = xASL_adm_UnzipNifti(refPath);

% Old Normalize settings
matlabbatch{1}.spm.tools.oldnorm.est.subj.source            = xASL_spm_admin(srcPath);
matlabbatch{1}.spm.tools.oldnorm.est.subj.wtsrc             = '';
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.template      = xASL_spm_admin(refPath);
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.weight        = '';
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.smosrc        = fwhmSrc;
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.smoref        = fwhmRef;
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.regtype       = 'subj';
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.cutoff        = 25; % biasfield correction
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.nits          = 0; % 16 default (includes DCT)
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.reg           = 1;

%  ------------------------------------------------------------------------------------------
%  Print & run

fprintf('\n%s\n','------------------------------------------------------------------------------------------');
fprintf('%s\n',['Affine registering ' srcPath ' to ' refPath]);   

spm_jobman('run',matlabbatch);

end
