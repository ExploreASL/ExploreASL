function xASL_spm_BiasfieldCorrection(PathIn, SPMdir, Quality, MaskName, PathOut)
%xASL_spm_BiasfieldCorrection ExploreASL wrapper around SPM "old segment" function,
% for biasfield removal only
% tested for M0 & mean_control images
% PathIn = path to input NIfTI
% SPMdir = folder where SPM is installed [optional]
% Quality = option to allow low quality (0), default is higher quality (1)
% MaskName = option for explicit masking [optional, default = no masking]
% PathOut = path to output NIfTI [optional]

%% ------------------------------------------------------------------------------------------------------
%% Admin
[Fpath, Ffile, Fext] = xASL_fileparts(PathIn);
PathIn = xASL_spm_admin(PathIn);
PathTemp = fullfile(Fpath, ['m' Ffile Fext]);

if nargin<5 || isempty(PathOut)
    PathOut = PathTemp;
end
if nargin<4 || isempty(MaskName)
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

if NoMaskInput
    % create an implicit mask from NaNs or zeros, to speed up this
    % segmentation, but SPM probably does this itself (implicit mask)?

    ImplicitMaskName            = fullfile(Fpath, ['ImplicitMask_' Ffile Fext]);
    MaskName                    = ImplicitMaskName;
    tIM                         = xASL_io_Nifti2Im(PathIn{1}(1:end-2));
    MaskIM                      = ones(size(tIM));
    MaskIM(tIM==0 | tIM==NaN)   = 0;
    xASL_io_SaveNifti(PathIn{1}(1:end-2),ImplicitMaskName,logical(MaskIM),16,0);
end




%% ------------------------------------------------------------------------------------------------------
%% Segment including creation biasfield
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
matlabbatch{1}.spm.tools.oldseg.opts.msk        = {MaskName};

if Quality
    matlabbatch{1}.spm.tools.oldseg.opts.samp = 9; % sufficient for smooth biasfields, saves time
    % PM: can make this dependent on the image resolution
else
    matlabbatch{1}.spm.tools.oldseg.opts.samp = 32; % sufficient for smooth biasfields, saves time
end

spm_jobman('run',matlabbatch);






%% ------------------------------------------------------------------------------------------------------
% %Delete temporary files
for ii=1:2
    xASL_delete(snMat{ii});
end
xASL_delete(ImplicitMaskName); % this is the created mask only


%% ------------------------------------------------------------------------------------------------------
%% Rename temporary SPM file into output file

if ~strcmp(PathTemp,PathOut)
    xASL_Move(PathTemp,PathOut, true);
end

end