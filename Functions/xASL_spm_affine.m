function xASL_spm_affine(srcPath, refPath, fwhmSrc, fwhmRef, otherList, bDCT, bQuality)
% ExploreASL wrapper for SPM affine registration function (a.k.a. 'old normalize'). On default run without DCT.
%
% FORMAT: xASL_spm_affine(srcPath, refPath, fwhmSrc, fwhmRef[,otherList, bDCT, bQuality])
% 
% INPUT:
%   refPath   - path to reference space (NifTI image) you want to register the source image to (REQUIRED)
%   srcPath   - path to source image (NifTI image) you want to register (REQUIRED)
%   fwhmSrc   - Gaussian smoothing to be applied to the source image before estimating the affine registration, in FWHM (mm) (REQUIRED)
%   fwhmRef   - Gaussian smoothing to be applied to the reference image before estimating the affine registration, in FWHM (mm) (REQUIRED)
%   otherList - a list of NIFTIs to which should this registration be applied (OPTIONAL, default EMPTY)
%   bDCT      - boolean specifying to perform the low-degree Discrete Cosine Transform (DCT) (OTIONAL, default FALSE)
%   bQuality  - boolean for quality mode (TRUE = high, FALSE = low) - decreases the number of DCT coefficients for DCT (OPTIONAL, default TRUE)
%
% OUTPUT: n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This SPM wrapper runs SPM's old normalize-estimate function, which calculates the affine transformation (i.e. linear + zooming and shearing) that is required to
% align the source image with the reference image. By default it does not estimate the low-degree Discrete Cosine Transform (DCT) to have a simple affine transformation 
% but this can be enabled in this wrapper. Also note that this affine transformation uses a correlation cost function, hence it requires the source and reference images 
% to have similar contrasts and resolution - or provide the resolution estimates so the smoothing can be done.
% This function does not change the orientation header by default, but it allows to change those of others through the otherList. If bDCT is used or no otherList given,
% it stores its estimated affine transformation as orientation difference matrix in a file with the same path but _sn.mat extension.
% For the provided smoothing FWHM, note that smoothnesses combine with Pythagoras' rule (i.e. square summing).
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_spm_affine('/MyStudy/Subject1/T1.nii.gz', '/MyStudy/Subject1/mean_control.nii', 5, 5);
% __________________________________
% Copyright 2015-2020 ExploreASL

   
% Check parameters
if nargin<5
    otherList = '';
else
    otherList = xASL_adm_OtherListSPM(otherList);
end
if nargin < 4
	error('Requires 4 input parameters');
end
if ~xASL_exist(srcPath) || ~xASL_exist(refPath)
	error('Cannot find input images');
end

if nargin < 6 || isempty(bDCT)
	bDCT = false;
end

if ~isempty(otherList) && bDCT
	warning('otherList not empty and bDCT==1. DCT produces a _sn.mat file with transformation parameters, it cannot apply the transformation to the files in the otherList');
end

if nargin < 7 || isempty(bQuality)
	bQuality = true;
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
if bDCT
	% Includes also the Direct Cosine Transform
	if bQuality
		matlabbatch{1}.spm.tools.oldnorm.est.eoptions.nits      = 16; % 16 is the default for SPM
	else
		matlabbatch{1}.spm.tools.oldnorm.est.eoptions.nits      = 4; % low quality mode
	end
else
	% Excludes DCT
	matlabbatch{1}.spm.tools.oldnorm.est.eoptions.nits      = 0;
end
matlabbatch{1}.spm.tools.oldnorm.est.eoptions.reg           = 1;

%  ------------------------------------------------------------------------------------------
%  Print & run

fprintf('\n%s\n','------------------------------------------------------------------------------------------');
fprintf('%s\n',['Affine registering ' srcPath ' to ' refPath]);   

spm_jobman('run',matlabbatch);

%% Apply to NIfTIs
% If we apply this to NIfTIs, we remove the transformation,
% to keep this transparant. This is done by providing otherList,
% even if we want to apply it to the srcPath
% If we don't provide an otherList, this estimated affine transformation
% is only saved as _sn.mat, not applied

if ~isempty(otherList) && ~bDCT
    [Fpath, Ffile] = xASL_fileparts(srcPath);
    PathMat = fullfile(Fpath, [Ffile '_sn.mat']);
    SnMat = load(PathMat);
    Affine = (SnMat.VG.mat/SnMat.Affine)/SnMat.VF.mat;

    for iList=1:length(otherList)
        nii = [];
        
        % find comma
        StartComma = regexp(otherList{iList},',');
        if ~isempty(StartComma)
            Index4D = xASL_str2num(otherList{iList}(StartComma+1:end));
            otherList{iList} = otherList{iList}(1:StartComma-1);
            if Index4D==1
                % apply affine transformation to NIfTI file
                nii = xASL_io_ReadNifti(otherList{iList});
                nii.mat = Affine*nii.mat;
                create(nii);
            else
                % apply affine transformation to motion orientation sidecar
                [Fpath, Ffile] = xASL_fileparts(otherList{iList});
                PathMotion = fullfile(Fpath, [Ffile '.mat']);
                Mat = load(PathMotion,'-mat');
                mat = Mat.mat;
                mat(:,:,Index4D) = Affine*mat(:,:,Index4D);
                save(PathMotion,'mat');
            end
        else % apply affine transformation to NIfTI file
            nii = xASL_io_ReadNifti(otherList{iList});
            nii.mat = Affine*nii.mat;
            create(nii);
        end
    end
    xASL_delete(PathMat);
end


end
