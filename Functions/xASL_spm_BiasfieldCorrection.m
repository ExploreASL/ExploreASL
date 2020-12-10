function xASL_spm_BiasfieldCorrection(PathIn, SPMdir, Quality, PathMask, PathOut)
%xASL_spm_BiasfieldCorrection Remove biasfield from image
%
% FORMAT: xASL_spm_BiasfieldCorrection(PathIn, SPMdir, Quality, MaskName, PathOut)
%
% INPUT:
%   PathIn      - path to the NIfTI from which to compute the biasfield
%   SPMdir      - path to the SPM installation (REQUIRED, in ExploreASL use x.SPMdir)
%   Quality     - boolean specifying smoothness of the biasfield, with 1 for high and 0 for low quality (OPTIONAL, DEFAULT=high quality)
%   MaskName    - legacy (unused) parameter
%   PathOut     - path to output NIfTI (OPTIONAL, DEFAULT = PathIn prefixed by 'm'
%
% OUTPUT:
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function is a wrapper around the SPM "old segment"
% function, for biasfield removal. It is tested for M0 and mean control
% images. It conducts the following steps:
%
% 1. Create implicit mask
% 2. Define SPM 'old segmentation' settings
% 3. Run SPM 'old segmentation'
% 4. Delete temporary files
% 5. Rename temporary SPM file into output file
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_spm_BiasfieldCorrection('/MyStudy/sub-001/T1.nii', x.SPMdir);
% __________________________________
% Copyright 2015-2020 ExploreASL



%% ------------------------------------------------------------------------------------------------------
%% Admin
[Fpath, Ffile, Fext] = xASL_fileparts(PathIn);
PathIn = xASL_spm_admin(PathIn);
PathTemp = fullfile(Fpath, ['m' Ffile Fext]);

if nargin<5 || isempty(PathOut)
    PathOut = PathTemp;
end
if nargin<4 || isempty(PathMask)
    NoMaskInput = true;
else
    NoMaskInput = false;
end
if nargin<3 || isempty(Quality)
    Quality = true; % default to high quality
end
if nargin<2 || isempty(SPMdir)
    % try to redeem this path
    SPMdir = which('spm');
    if ~isempty(SPMdir)
        SPMdir = fileparts(SPMdir);
    else
        error('Cannot find SPM, path needed for tissue priors');
    end
end

% Delete any existing transformations
snMat{1} = fullfile(Fpath,[Ffile '_seg_sn.mat']);
snMat{2} = fullfile(Fpath,[Ffile '_seg_inv_sn.mat']);

for ii=1:2
    xASL_delete(snMat{ii});
end

%% 1) Create implicit mask
if NoMaskInput
    % create an implicit mask from NaNs or zeros, to speed up this
    % segmentation, but SPM probably does this itself (implicit mask)?

    PathMask = fullfile(Fpath, ['ImplicitMask_' Ffile Fext]);
    tempIM = xASL_io_Nifti2Im(PathIn{1}(1:end-2));
    MaskIM = ones(size(tempIM));
    MaskIM(tempIM==0 | isfinite(tempIM)) = 0;
    xASL_io_SaveNifti(PathIn{1}(1:end-2), PathMask, logical(MaskIM), [], 0);
end




%% ------------------------------------------------------------------------------------------------------
%% 2) Define SPM 'old segmentation' settings
matlabbatch{1}.spm.tools.oldseg.data            = PathIn;
matlabbatch{1}.spm.tools.oldseg.output.GM       = [0 0 0];
matlabbatch{1}.spm.tools.oldseg.output.WM       = [0 0 0];
matlabbatch{1}.spm.tools.oldseg.output.CSF      = [0 0 0];
matlabbatch{1}.spm.tools.oldseg.output.biascor  = 1;
matlabbatch{1}.spm.tools.oldseg.output.cleanup  = 0;
matlabbatch{1}.spm.tools.oldseg.opts.tpm = {
                                            fullfile(SPMdir,'toolbox','OldSeg','grey.nii');
                                            fullfile(SPMdir,'toolbox','OldSeg','white.nii');
                                            fullfile(SPMdir,'toolbox','OldSeg','csf.nii');
                                            };
matlabbatch{1}.spm.tools.oldseg.opts.ngaus      = [2;2;2;4];
matlabbatch{1}.spm.tools.oldseg.opts.regtype    = 'mni';
matlabbatch{1}.spm.tools.oldseg.opts.warpreg    = 1;
matlabbatch{1}.spm.tools.oldseg.opts.warpco     = Inf; % 25; % disables the registration, saves time
matlabbatch{1}.spm.tools.oldseg.opts.biasreg    = 0.0001;
matlabbatch{1}.spm.tools.oldseg.opts.biasfwhm   = 60;
matlabbatch{1}.spm.tools.oldseg.opts.msk        = {PathMask};

if Quality
    matlabbatch{1}.spm.tools.oldseg.opts.samp = 9;
    % PM: can make this dependent on the image resolution
else
    matlabbatch{1}.spm.tools.oldseg.opts.samp = 32; % sufficient for smooth biasfields, saves time
end

%% ------------------------------------------------------------------------------------------------------
%% 3) Run SPM 'old segmentation'
spm_jobman('run',matlabbatch);



%% ------------------------------------------------------------------------------------------------------
%% 4) Delete temporary files
for ii=1:2
    xASL_delete(snMat{ii});
end
if NoMaskInput
    xASL_delete(PathMask); % this is the created mask only
end


%% ------------------------------------------------------------------------------------------------------
%% 5) Rename temporary SPM file into output file
if ~strcmp(PathTemp, PathOut)
    xASL_Move(PathTemp, PathOut, true);
end


end