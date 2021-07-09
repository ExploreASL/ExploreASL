function xASL_spm_coreg(refPath, srcPath, OtherList, x, sep, FastReg)
%xASL_spm_coreg ExploreASL wrapper for SPM coregistration function (a.k.a. 'coreg')
%
% FORMAT: xASL_spm_coreg(refPath, srcPath[, OtherList, x, sep, FastReg])
%
% INPUT:
%   refPath             - path to reference space (NifTI image) you want to register the source image to (REQUIRED)
%   srcPath             - path to source image (NifTI image) you want to register (REQUIRED)
%   OtherList           - cell list with paths to NIfTI images that you want to apply the calculated transformation to (OPTIONAL)
%                         the srcPath NIfTI orientation header is always adapted with the calculated transformation
%   x                   - ExploreASL structure containing fields with global information about the pipeline environment
%                         and settings (e.g. x.settings.Quality), useful when you want this script to copy the options of an ExploreASL pipeline run (OPTIONAL)
%   sep                 - separation setting used by SPM coregistration, that defines the average distance between sampled points (in mm)
%                         can also have a vector of multiple values to allow a coarse registration with increasingly finer ones (OPTIONAL, DEFAULT=[4 2])
%   FastReg             - New experimental feature to speed up registration of low resolution images.
%                         When low resolution images are registered to high resolution images, the high resolution matrix and
%                         sep values can unnecessarily increase the precision and duration of the coregistration (e.g. from 10 to 55 seconds)
%                         This parameter is WIP, so please leave it at default setting for now(OPTIONAL, default=false)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This SPM wrapper runs SPMs coregister-estimate function, which calculates the 6 parameter rigid-body transformation (a.k.a. linear) that is required to
% align the source image with the reference image. This 6 parameter transformation (i.e. 3 XYZ translations and 3 rotations) is applied to
% the orientation header of the source NIfTI image, and also to the images provided in OtherList (optional).
% Note that this SPM registration function uses a normalized mutual information (NMI) by default, enabling registration between two images
% with different contrast.
% Note that this algorithm will use the first volume of a multi-volume NIfTI
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_spm_coreg('/MyStudy/Subject1/T1.nii.gz', '/MyStudy/Subject1/mean_control.nii', {'/MyStudy/Subject1/M0.nii'});
% __________________________________
% Copyright 2015-2019 ExploreASL

%%  ------------------------------------------------------------------------------------------
%%   Admin, manage input, unzipping if needed
[~, refFile, refExt] = xASL_fileparts(refPath);
[~, srcFile, srcExt] = xASL_fileparts(srcPath);

refFile = [refFile refExt];
srcFile = [srcFile srcExt];

if nargin<4 || isempty(x) || ~isfield(x,'Quality') || isempty(x.settings.Quality)
    x.settings.Quality = true; % default quality is high
    x.DELETETEMP = true;
end
if nargin<3 || isempty(OtherList)
    OtherList = {''};
end
if nargin<5 || isempty(sep)
    if x.settings.Quality
        sep = [4 2];
    else
        sep = 6;
    end
end
if nargin<6 || isempty(FastReg)
    FastReg = false; % fast registration by resampling to resolution of source
    % use e.g. when registering 2D EPI to high resolution T1w,
    % knowing that NMI registration will suffice, which is case e.g. with
    % DTI to c1+c2+c3
end

if FastReg
    [Fpath, Ffile] = xASL_fileparts(refPath);
    tempResliced = fullfile(Fpath, [Ffile '_tempResliced.nii']);
    niiRef = xASL_io_ReadNifti(refPath);
    niiSrc = xASL_io_ReadNifti(srcPath);
    VoxelSizeRef = niiRef.hdr.pixdim(2:4);
    VoxelSizeSrc = niiSrc.hdr.pixdim(2:4);
    VoxelSizeRefMin = min(VoxelSizeRef);
    VoxelSizeSrcMin = min(VoxelSizeSrc);
    
    if VoxelSizeRef>2*VoxelSizeSrc
        ResliceIs = 2;
    elseif VoxelSizeSrc>2*VoxelSizeRef
        % this is the case e.g. when registering fMRI/ASL image to T1w
        ResliceIs = 1;
    else
        ResliceIs = 0;
    end
    
    if ResliceIs==1
        % reslice ref
        xASL_Copy(refPath, tempResliced, true);
        NewVoxelSize = repmat(VoxelSizeSrcMin,[1,3]);
        xASL_spm_smooth(tempResliced, NewVoxelSize, tempResliced);
        xASL_im_Upsample(tempResliced, tempResliced, NewVoxelSize);

        matlabbatch{1}.spm.spatial.coreg.estimate.ref = xASL_spm_admin(tempResliced);
        matlabbatch{1}.spm.spatial.coreg.estimate.source = xASL_spm_admin(srcPath);
    elseif ResliceIs==2
        % reslice source
        xASL_Copy(srcPath, tempResliced, true);
        NewVoxelSize = repmat(VoxelSizeRefMin,[1,3]);
        xASL_spm_smooth(tempResliced, NewVoxelSize, tempResliced);
        xASL_im_Upsample(tempResliced, tempResliced, NewVoxelSize);

        matlabbatch{1}.spm.spatial.coreg.estimate.ref = xASL_spm_admin(refPath);
        matlabbatch{1}.spm.spatial.coreg.estimate.source = xASL_spm_admin(tempResliced);
    else
        % reslice none
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = xASL_spm_admin(refPath);
        matlabbatch{1}.spm.spatial.coreg.estimate.source = xASL_spm_admin(srcPath);
    end
else
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = xASL_spm_admin(refPath);
    matlabbatch{1}.spm.spatial.coreg.estimate.source = xASL_spm_admin(srcPath);    
end


%%  ------------------------------------------------------------------------------------------
%%  Default SPM settings

matlabbatch{1}.spm.spatial.coreg.estimate.other = xASL_adm_OtherListSPM(OtherList);

matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = sep;


%%  ------------------------------------------------------------------------------------------
%%  Print & run
fprintf('\n%s\n\n','------------------------------------------------------------------------------------------');
fprintf('%s\n',['Rigid-body registering ' srcFile ' to ' refFile]);

spm_jobman('run',matlabbatch);
close all

if x.DELETETEMP
    srcPath = xASL_spm_admin(srcPath);
    xASL_adm_DeleteFileList(fileparts(srcPath{1}), '^spm_.*\.ps$', false, [0 Inf]);

    if FastReg
        xASL_delete(tempResliced);
    end
end

end
