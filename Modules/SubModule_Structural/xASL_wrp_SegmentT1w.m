function [x] = xASL_wrp_SegmentT1w(x, SegmentSPM12)
%xASL_wrp_SegmentT1w Submodule of ExploreASL Structural Module, that segments 3D T1 (or T2) scan
%
% FORMAT: [x] = xASL_wrp_SegmentT1w(x, SegmentSPM12)
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P     - paths with NIfTIs for which this function should be applied to (REQUIRED)
%   SegmentSPM12 - Whether to run SPM12 (true) or CAT12 (false) (OPTIONAL, DEFAULT = false)
%   x.bFixResolution - resample to a resolution that CAT12 accepts (OPTIONAL, DEFAULT=false)
%   x.Pediatric_Template - boolean specifying if we use a pediatric
%             template instead of adult one (OPTIONAL, DEFAULT = false)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule segments high resolution structural/anatomical scans into GM/WM/CSF/soft tissue/bone/air tissue classes.
% It will save GM/WM/CSF in native space, and the transformation from native to standard space.
% This transformation includes Geodesic Shooting/DARTEL for CAT12.
%
% This submodule contains the following steps:
% 1) Administration
% 2) Extra segmentation options by Jan Petr
% 3) Segmentation using CAT12
%    -> If CAT12 fails, it will be repeated with higher contrast, higher strength affine preprocessing & less biasfield regularization
%    -> If CAT12 fails twice, it will be skipped & SPM12 will be run
% 4) Segmentation using SPM12
% 5) File management CAT12
% 6) File management lesions
% 7) Resample lesions to standard space
%    -> for the lesion masking. MORE EXPLANATION NEEDED BY JAN
% 8) Manage flowfields
%    -> smooth combination non-linear flowfield outside the lesion & uniform flowfield within the lesion
% 9) File management
%
%
% EXAMPLE: xASL_wrp_SegmentT1w(x);
%
% REFERENCE:
% Gaser, C., 2009. Partial volume segmentation with adaptive maximum a posteriori (MAP) approach. NeuroImage 47, S121.
% Mendrik AM, Vincken KL, Kuijf HJ, et al. MRBrainS Challenge: Online Evaluation Framework for Brain Image Segmentation in 3T MRI Scans. Comput.Intell.Neurosci. 2015. p. 813696-mrbrains13.isi.uu.nl/results.php.
% __________________________________
% Copyright 2015-2019 ExploreASL






if nargin<2 || isempty(SegmentSPM12)
    SegmentSPM12 = false; % by default use CAT12, not SPM12 to segment
end
if ~isfield(x,'bFixResolution') || isempty(x.bFixResolution)
    x.bFixResolution = false;
end
if ~isfield(x,'Pediatric_Template') || isempty(x.Pediatric_Template)
    x.Pediatric_Template = false;
end

%% --------------------------------------------------------------------------------
%% 1) Administration
if x.Quality~=0 && x.Quality~=1
    error('Wrong quality definition');
end

if ~isfield(x,'Seg')
	x.Seg = {};
end

% Check whether we should do normal or strong biasfield correction
x.T1BiasFieldRegularization = true; % default
if isfield(x,'Vendor') && ~isempty(regexp(x.Vendor,'GE'))
    x.T1BiasFieldRegularization = false; % SPM12
    % GE has wider bore scanners, resulting in a wide biasfield
end



% Study file management
xASL_io_ReadNifti(x.P.Path_T1); % unzip if necessary
xASL_adm_DeleteFileList(x.SUBJECTDIR, ['^c[1-3]' x.P.STRUCT '\.(nii|nii\.gz)$']); % make sure no old files are left when job fails
xASL_adm_DeleteFileList(x.SUBJECTDIR, ['^' x.P.STRUCT '_seg8\.mat$']);
xASL_adm_DeleteFileList(x.SUBJECTDIR, ['^y_' x.P.STRUCT '\.nii$']);
xASL_adm_DeleteFileList(x.SUBJECTDIR, '^catreport.*\.pdf$');

xASL_adm_RemoveTempFilesCAT12(x, true); % remove previous temporary CAT12 files

%% --------------------------------------------------------------------------------
%% 1.5 Fix resolution for CAT12 % if we would activate this, we should also
%     fix the FLAIR resolution, otherwise this crashes when cleaning up the
%     FLAIR WMH segmentation
if x.bFixResolution % if we use CAT12 to segment
    % we force sufficient isotropy, otherwise CAT12 will crash
    tNii = xASL_io_ReadNifti(x.P.Path_T1);
    CurrentVoxelSize = tNii.hdr.pixdim(2:4);
    if max(CurrentVoxelSize./[1.5 1.5 1.5]>[1.1 1.1 1.1]) % if any dimension is smaller than 1.5 mm resolution:
        NewVoxelSize = min(CurrentVoxelSize,[1.5 1.5 1.5]);
    
        % first backup T1, if backup doest exist yet
        if ~xASL_exist(x.P.Path_T1_ORI)
            xASL_Copy(x.P.Path_T1, x.P.Path_T1_ORI);
        end
        
        xASL_im_Upsample(x.P.Path_T1, x.P.Path_T1, NewVoxelSize, [], [], 'spline');
    end
end


%% --------------------------------------------------------------------------------
%% 2) Extra segmentation options by Jan Petr (NEED MERGING WITH ABOVE)
% NOVICE study
% x.Seg.DisableDARTEL = 1;

% PICTURE study
% x.Seg.Method = 'GS'; % Or 'DARTEL'
% x.Seg.SaveSPMFlowField = 1;
% x.Seg.SaveOriginalFlowField = 1;
% x.Seg.SaveMixedFlowField = 1;
% x.Seg.SaveIntermedFlowField = 1;


% Disable the DARTEL/CAT12 registration completely and run only DCT+Affine registration (if TRUE)
if ~isfield(x.Seg,'DisableDARTEL')
	x.Seg.DisableDARTEL = false;
end

% Sets the basic method for the MNI space registration after CAT12 segmentation
% 'default' = Automatic choosing between GS and DARTEL.
% 'GS' = does only Geodesic shooting
% 'DARTEL' = does only DARTEL
if ~isfield(x.Seg,'Method')
	x.Seg.Method = 'default';
elseif ~strcmp(x.Seg.Method, 'default') && ~strcmp(x.Seg.Method, 'DARTEL') && ~strcmp(x.Seg.Method, 'GS')
    warning(['Wrong x.Seg.Method: ' xASL_num2str(x.Seg.Method) ', using default setting instead']);
    x.Seg.Method = 'default';
end

% If TRUE, then additionally saves the flow field of the affine+DCT normalization (low degree of freedom non-linear)
% that is executed before DARTEL or GS for pre-registration
% Filename y_T1_SPM.nii
if ~isfield(x.Seg,'SaveSPMFlowField')
	x.Seg.SaveSPMFlowField = false;
end

% In the case of lesion masking and mixing of the AFFINE+DARTEL flow field, this option saves additionally the
% non-modified DARTEL/GS flow-field
% Filename y_T1_orig.nii
if ~isfield(x.Seg,'SaveOriginalFlowField')
	x.Seg.SaveOriginalFlowField = false; % Additionally save the original flowfield
end

% In the case of lesion masking it saves additionally the  mix of the AFFINE+DARTEL flow field and the
% DARTEL/GS flow-field
% Filename y_T1_mixed.nii
if ~isfield(x.Seg,'SaveMixedFlowField')
	x.Seg.SaveMixedFlowField = false; % Additionally saves an improved mixing that takes care to be invertible
end

% Saves the intermediate steps for DARTEL and GS y_T1_1.nii, y_T1_2.nii,... y_T1_6.nii
if ~isfield(x.Seg,'SaveIntermedFlowField')
	x.Seg.SaveIntermedFlowField = false;
end

% If the DARTEL/GS is disabled, then it does not make sense to activate any of these extra options
if x.Seg.DisableDARTEL
	x.Seg.Method = 'default';
	x.Seg.SaveSPMFlowField = false;
	x.Seg.SaveOriginalFlowField = false;
	x.Seg.SaveMixedFlowField = false;
	x.Seg.SaveIntermedFlowField = false;
end





%% -------------------------------------------------------------------------------------------
%% 3) Segmentation using CAT12
%  This runs by default (default = x.SegmentSPM == 0)
%  When it fails, it will pass x.SegmentSPM == 1
if ~SegmentSPM12
    SegmentSPM12 = xASL_wrp_CAT12Segmentation(x);
end






%% -------------------------------------------------------------------------------------------
%% 4) Segmentation using SPM12
%  This usually gives poorer results than CAT12, but can be chosen if CAT12 doesnt work
if SegmentSPM12
    xASL_wrp_SPM12Segmentation(x);

    if ~x.Quality % With low quality, registration was performed on lower resolution,
        % & transformation needs to be upsampled to correct resolution
        xASL_spm_reslice(x.D.ResliceRef, x.P.Path_y_T1, [], [], x.Quality, x.P.Path_y_T1, 1);
    end
    return; % rest of the function is housekeeping for CAT12
end





%% ----------------------------------------------------------------------------------------
%% 5) File management CAT12 names
InFile{1} = fullfile(x.SUBJECTDIR,'mri',['p1' x.P.STRUCT '.nii']);
InFile{2} = fullfile(x.SUBJECTDIR,'mri',['p2' x.P.STRUCT '.nii']);
InFile{3} = fullfile(x.SUBJECTDIR,'mri',['p3' x.P.STRUCT '.nii']);
InFile{4} = fullfile(x.SUBJECTDIR,'mri',['y_' x.P.STRUCT '.nii']);
InFile{5} = fullfile(x.SUBJECTDIR,'label',['catROI_' x.P.STRUCT '.mat']);
InFile{6} = fullfile(x.SUBJECTDIR,'report',['catreport_' x.P.STRUCT '.pdf']);
InFile{7} = fullfile(x.SUBJECTDIR,'report',['cat_' x.P.STRUCT '.mat']);

OutFile{1} = x.P.Path_c1T1; % GM segmentation
OutFile{2} = x.P.Path_c2T1; % WM segmentation
OutFile{3} = x.P.Path_c3T1; % CSF segmentation
OutFile{4} = x.P.Path_y_T1; % deformation field to common space
OutFile{5} = fullfile(x.D.TissueVolumeDir,['catROI_' x.P.STRUCT '_' x.P.SubjectID '.mat']); % contains ROI volume values from several atlases
OutFile{6} = fullfile(x.SUBJECTDIR,['catreport_' x.P.STRUCT '.pdf']);
OutFile{7} = fullfile(x.D.TissueVolumeDir,['cat_' x.P.STRUCT '_' x.P.SubjectID '.mat']);

% If saving separate files - copy the intermediate steps
if x.Seg.SaveIntermedFlowField
    switch x.Seg.Method
        case 'DARTEL'
            nSteps = 6;
        case 'GS'
            nSteps = 6;
    end
    for iii = 1:nSteps
        if xASL_exist(fullfile(x.SUBJECTDIR, 'mri', ['y_' x.P.STRUCT '_' num2str(iii) '.nii']), 'file')
            xASL_Move(fullfile(x.SUBJECTDIR,'mri',['y_' x.P.STRUCT '_' num2str(iii) '.nii']),...
            fullfile(x.SUBJECTDIR, ['y_' x.P.STRUCT '_' num2str(iii) '.nii']), true, false);
        end
    end
end

for ii=1:length(InFile)
    if xASL_exist(InFile{ii},'file')
        xASL_Move(InFile{ii}, OutFile{ii}, true, false); % move file, overwrite old file, no verbose
    end
end


%% ----------------------------------------------------------------------------------------
%% 6) File management lesion segmentations
%  Obtain paths of lesion files if they exist, to correct flowfields for
Lesion_T1_list     = xASL_adm_GetFileList(x.SUBJECTDIR,['^Lesion_' x.P.STRUCT '_\d*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);
Lesion_FLAIR_list  = xASL_adm_GetFileList(x.SUBJECTDIR,['^Lesion_' x.P.FLAIR '_\d*\.(nii|nii\.gz)$'],'FPList',[0 Inf]);
Path_Transf_SPM      = fullfile(x.SUBJECTDIR,'mri',['y_' x.P.STRUCT '_withoutDARTEL.nii']);
%         Path_Transf_SPM_subj = fullfile(x.SUBJECTDIR,['y_' x.P.STRUCT '_withoutDARTEL.nii']);
Path_Transf_DARTEL   = fullfile(x.SUBJECTDIR,'mri',['y_' x.P.STRUCT '_withDARTEL.nii']);

% 1) Resample all lesion masks to standard space (avoids differences T1 & FLAIR spaces)
% 2) Combine all lesions into single Total_Lesion.nii, save this, remove the others
% 3) correct flow fields with the Total_Lesion.nii

% Replace flow fields, to use the SPM flow field for resampling of lesions to standard space

if xASL_exist(x.P.Path_y_T1,'file')
    xASL_Move(x.P.Path_y_T1, Path_Transf_DARTEL, true);
end
if xASL_exist(Path_Transf_SPM,'file')
    xASL_Move(Path_Transf_SPM, x.P.Path_y_T1, true);
end

% Delete existing zipped lesion masks
xASL_adm_DeleteFileList(x.D.PopDir, ['^rLesion.*' x.P.SubjectID '\.nii\.gz'], 0, [0 Inf]);



%% ----------------------------------------------------------------------------------------
%% 7) Resample all lesion masks to standard space (avoid differences between input spaces)
for iL=1:length(Lesion_FLAIR_list)
    rLesion_FLAIR_list{iL} = fullfile(x.D.PopDir, ['rLesion_' x.P.FLAIR '_' num2str(iL) '_' x.P.SubjectID '.nii']);
    xASL_spm_deformations(x,Lesion_FLAIR_list{iL}, rLesion_FLAIR_list{iL}, 2);
end
for iL=1:length(Lesion_T1_list)
    rLesion_T1_list{iL} = fullfile(x.D.PopDir, ['rLesion_' x.P.STRUCT '_' num2str(iL) '_' x.P.SubjectID '.nii']);
    xASL_spm_deformations(x, Lesion_T1_list{iL}, rLesion_T1_list{iL}, 2);
end
xASL_Move(x.P.Path_y_T1, Path_Transf_SPM, true); % move the flow field back

% Combine all lesions
rLesionList = xASL_adm_GetFileList(x.D.PopDir, ['^rLesion.*' x.P.SubjectID '\.nii'], 'FPList', [0 Inf]);
LesionIM = zeros([121 145 121]); % assuming 1.5 mm MNI
for iL=1:length(rLesionList)
    LesionIM(xASL_im_ConvertMap2Mask(xASL_io_Nifti2Im(rLesionList{iL}))) = 1;
    xASL_delete(rLesionList{iL});
end
% Lesion saving in standard space happens later, with the improved flow fields





%% ----------------------------------------------------------------------------------------
%% 8) Manage flowfields
if xASL_stat_SumNan(LesionIM(:))>0
    % Mix the DARTEL (non-linear) and SPM (uniform) flowfields based on Euclidian distance matrix
    Fields2Save = xASL_wrp_CombineFlowFields(x, Path_Transf_SPM, Path_Transf_DARTEL, LesionIM);
    xASL_Move(Path_Transf_DARTEL, x.P.Path_y_T1, true);
else % if no lesion existed, keep the non-linear flowfield (if it exists, otherwise take the SPM flowfield)
    if xASL_exist(Path_Transf_DARTEL, 'file')
        xASL_Move(Path_Transf_DARTEL, x.P.Path_y_T1, true);
    elseif xASL_exist(Path_Transf_SPM, 'file')
           xASL_Move(Path_Transf_SPM, x.P.Path_y_T1, true);
    end

    Fields2Save = xASL_io_Nifti2Im(x.P.Path_y_T1);
end

% Extrapolate over NaNs for smooth edges of flow fields
for iE=1:size(Fields2Save, 5)
    Fields2Save(:,:,:,1,iE) = xASL_im_ExtrapolateOverNaNs(squeeze(Fields2Save(:,:,:,1,iE)));
end


xASL_io_SaveNifti(x.P.Path_y_T1, x.P.Path_y_T1, Fields2Save, [], false);



%% ----------------------------------------------------------------------------------------
%% 9) File management
xASL_adm_RemoveTempFilesCAT12(x);

if ~x.Quality % With low quality, registration was performed on lower resolution,
              % & transformation needs to be upsampled to correct resolution
    xASL_spm_reslice(x.D.ResliceRef, x.P.Path_y_T1, [], [], x.Quality, x.P.Path_y_T1, 1);
end



end















%% ===================================================================================================================
%% ===================================================================================================================
function xASL_adm_RemoveTempFilesCAT12(x, bForce)
%xASL_adm_RemoveTempFilesCAT12 Removes residual/temporal files in mri directory

if nargin<2 || isempty(bForce)
    if ~x.DELETETEMP
        bForce = false;
    else
        bForce = true; % default
    end
end

% delete previous error files
% Remove residual files in mri directory
if bForce
    DirList{1} = fullfile(x.SUBJECTDIR,'mri');
    DirList{2} = fullfile(x.SUBJECTDIR,'label');
    DirList{3} = fullfile(x.SUBJECTDIR,'report');
    ErrDir     = fullfile(x.SUBJECTDIR,'err');
    if exist(ErrDir, 'dir')
        TL = xASL_adm_GetFileList(ErrDir,'^cat.*$','List',[0 Inf], true); % find dirs
        if ~isempty(TL)
            DirList{4} = fullfile(ErrDir,TL{1});
        end
        DirList{end+1} = ErrDir;
    end

    for iL=1:length(DirList)
        DelList = xASL_adm_GetFileList(DirList{iL},'^.*\.(pdf|nii|xml|mat|ini|db|DS_Store)$','FPlist',[0 Inf]);
        for iD=1:length(DelList)
            xASL_delete(DelList{iD});
        end
        if exist(DirList{iL}, 'dir')
           try
              rmdir(DirList{iL}, 's');
           catch % do nothing
           end
        end
    end
end


end



%% ===================================================================================================================
%% ===================================================================================================================
function xASL_wrp_SPM12Segmentation(x)
%xASL_wrp_SPM12Segmentation Run the SPM12 segmentation


    fprintf('%s\n','Segmenting structural image using SPM12');
    warning('off', 'MATLAB:nearlySingularMatrix');

    matlabbatch{1}.spm.spatial.preproc.channel.vols = {x.P.Path_T1};

    if x.T1BiasFieldRegularization
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0;
    else
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.01;
    end

    %% Fill the templates used for segmentation
	if x.Pediatric_Template
		for iDim=1:6
			matlabbatch{1}.spm.spatial.preproc.tissue(iDim).tpm = {fullfile(x.SPMDIR, 'toolbox', 'cat12', 'templates_pediatric', ['infant-1yr-TPM.nii,' num2str(iDim)])};
		end
	else
		for iDim=1:6
			matlabbatch{1}.spm.spatial.preproc.tissue(iDim).tpm = {fullfile(x.SPMDIR, 'tpm', ['TPM.nii,' num2str(iDim)])};
		end
	end
	

    %% Tissue-class specific options
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus      = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native     = [1 0]; % store this
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped     = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus      = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native     = [1 0]; % store this
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped     = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus      = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native     = [1 0]; % store this
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped     = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus      = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native     = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped     = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus      = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native     = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped     = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus      = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native     = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped     = [0 0];

    %% General options
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm     = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write        = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf             = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup         = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg             = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg          = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm            = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.write           = [0 1];

    if x.Quality
        matlabbatch{1}.spm.spatial.preproc.warp.samp        = 3;
    else
        matlabbatch{1}.spm.spatial.preproc.warp.samp        = 9;
    end

    %% Run the segmentation
    spm_jobman('run', matlabbatch);
    warning('on', 'MATLAB:nearlySingularMatrix');

end






%% ===================================================================================================================
%% ===================================================================================================================
function [SegmentSPM12] = xASL_wrp_CAT12Segmentation(x)
%xASL_wrp_CAT12Segmentation Run the CAT12 segmentation

SegmentSPM12 = true; % by default, run SPM12 when CAT12 crashes


%% --------------------------------------------------------------------
%% CAT12 template & registration settings
SPMTemplateNII    = fullfile(x.SPMDIR, 'tpm', 'TPM.nii');
DartelTemplateNII = fullfile(x.SPMDIR, 'toolbox', 'cat12', 'templates_1.50mm', 'Template_1_IXI555_MNI152.nii');
GSTemplateNII     = fullfile(x.SPMDIR, 'toolbox', 'cat12', 'templates_1.50mm', 'Template_0_IXI555_MNI152_GS.nii');

if x.Pediatric_Template
	DartelTemplateNII = fullfile(x.SPMDIR, 'toolbox', 'cat12', 'templates_pediatric', 'Template_1_infant-1yr_DARTEL.nii');
	GSTemplateNII     = fullfile(x.SPMDIR, 'toolbox', 'cat12', 'templates_pediatric', 'Template_0_infant-1yr_CAT.nii');
	x.Seg.Method = 'DARTEL';
end

if ~xASL_exist(SPMTemplateNII, 'file')
    error('SPM tissue priors missing, is SPM12 installed in the correct folder?');
elseif ~xASL_exist(DartelTemplateNII,'file') || ~xASL_exist(GSTemplateNII,'file')
    error('CAT12 DARTEL/GS template missing, is CAT12 installed in the correct folder? It should be installed in SPM12/toolbox/cat12');
end

matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm                         = {SPMTemplateNII};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.xasl_savesteps           = x.Seg.SaveIntermedFlowField;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.xasl_quality             = x.Quality;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.xasl_disabledartel       = x.Seg.DisableDARTEL;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.xasl_lesion              = {xASL_im_Lesion2CAT(x.P.Path_T1)};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm   = {DartelTemplateNII}; % Runs DARTEL to this n=555 subjects template
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm = {GSTemplateNII}; % Runs Geodesic Shooting to this n=555 subjects template
switch x.Seg.Method
    case 'default'
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr = 0.5; % chooses between 0 DARTEL & 4 Geodesic Shooting, this is the optimized 0.5 Geodesic Shooting
    case 'GS'
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr = 4;
    case 'DARTEL'
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr = 0;
end


% matlabbatch{1}.spm.tools.cat.estwrite.extopts.cleanupstr = 0.5; % default line in cat12 segmentation batch script
% Here we don't define CleanUpStr, which gives an error. For some reason it doesn"t want to be loaded through
% the SPM conf initialization script. Nevertheless, cleanupstr = 0.5 by default in cat_defaults.m
% so we can remove it here and ignore it

%% --------------------------------------------------------------------
%% CAT12 segmentation quality settings
matlabbatch{1}.spm.tools.cat.estwrite.extopts.SLC = 0;
if x.Quality
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP           = 1070; % light cleanup
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr        = 0.5; % strength local adaptive segmentation
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr       = 0; % using GCUT instead of SPM approach (more robust)
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox           = 1.5; % voxelsize on which registration is run (1.5 == default)
    matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr          = 0.5; % SPM bias-correction strenght
    matlabbatch{1}.spm.tools.cat.estwrite.opts.samp             = 3;   % spm sampling distance
    

else
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP           = 0; % light cleanup
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr        = 0; % strength local adaptive segmentation
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr       = 0; % default SPM approach
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox           = 3; % voxelsize on which registration is run (1.5 == default)
    matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr          = eps;
    matlabbatch{1}.spm.tools.cat.estwrite.opts.samp             = 9;   % spm sampling distance
end

if x.T1BiasFieldRegularization
    matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5; % CAT12
else
    % disable biasfield regularization for large biasfields (e.g. GE wide bore scanner)
    matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.75; % CAT12
end

%% --------------------------------------------------------------------
%% Residual CAT12 segmentation settings
matlabbatch{1}.spm.tools.cat.estwrite.data                  = {x.P.Path_T1}; % T1.nii
matlabbatch{1}.spm.tools.cat.estwrite.nproc                 = 0; % don't split the segmentation in multiple processes.
% This will run multi-threaded, but not start up other Matlab instances

matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg           = 'mni'; % regularize affine registration for MNI European brains
matlabbatch{1}.spm.tools.cat.estwrite.output.surface        = 0;   % don't do surface modeling
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.noROI  = struct([]); % don't do ROI estimations
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native      = 1;   % save c1T1 in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod         = 0;   % don't save modulation
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel      = 0;   % don't save DARTEL space c1T1, this happens below in the reslice part
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native      = 1;   % save c2T1 in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod         = 0;   % don't save modulation
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel      = 0;   % don't save DARTEL space c2T1, this happens below in the reslice part
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native      = 1;   % save c3T1 in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod         = 0;   % don't save modulation
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel      = 0;   % don't save DARTEL space c2T1, this happens below in the reslice part
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobian.warped= 0;   % don't save Jacobians
matlabbatch{1}.spm.tools.cat.estwrite.output.warps          = [1 0]; % save warp to MNI
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped    = 0;   % don't save bias-corrected T1.nii
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native      = 1;   % save c3T1 in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod         = 0;   % don't save modulation
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel      = 0;

if ~x.bFixResolution
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed= [1 0.1]; % process everything on 1 mm fixed resolution (default)
else
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.best = [0.5 0.1]; % process everything on best resolution (is probably lower as we forced this to be 1.5 mm)
end

 % PM: we can use this image, which is skull-stripped and in common space, to use as common space mask
 % to subtract c1 & c2 from to obtain c3 (CSF)

%% --------------------------------------------------------------------
%% Run CAT12 segmentation
% 1) We try the above settings
% 2) If this fails, we try to improve the contrast & repeat CAT12
% 3) If this fails, we exit this functions and pass an argument to run the SPM12 segmentation instead

try % 1) First attempt CAT12
    spm_jobman('run',matlabbatch); % Run CAT12
    SegmentSPM12 = false;
catch
    if ~x.Quality
        warning('CAT12 failed with x.Quality==0, try x.Quality==1 instead!');
    end

    try % 2) Second attempt CAT12
        %    Increase strength affine preprocessing (APP)
        %    improve contrast
        %    increase biasfield correction
        %    (This emperically outperformed repetition with SPM12)

        if ~xASL_exist(x.P.Path_T1_ORI, 'file') && xASL_exist(x.P.Path_T1, 'file')
            xASL_Copy(x.P.Path_T1, x.P.Path_T1_ORI);
        end

        xASL_io_SaveNifti(x.P.Path_T1,x.P.Path_T1, xASL_io_Nifti2Im(x.P.Path_T1).^(2^0.5), [], false); % increase contrast
        matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.75; % increase biasfield correction
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1; % Increase strenght affine preprocessing

        xASL_adm_RemoveTempFilesCAT12(x); % Delete previous CAT12 derivatives
        spm_jobman('run', matlabbatch); % Run CAT12
        SegmentSPM12 = false;
    catch % 3) Return to run SPM12 segmentation
        xASL_adm_RemoveTempFilesCAT12(x); % Delete previous CAT12 derivatives
        SegmentSPM12 = true;
    end
end


end










%% ===================================================================================================================
%% ===================================================================================================================
function [Fields2Save] = xASL_wrp_CombineFlowFields(x, Path_Transf_SPM, Path_Transf_DARTEL, LesionIM)
%xASL_wrp_CombineFlowFields Create Euclidian distance matrix, to mix the DARTEL (non-linear) and SPM (uniform) flow


%% Admin
IM_spm = xASL_io_Nifti2Im(Path_Transf_SPM);
IM_dartel = xASL_io_Nifti2Im(Path_Transf_DARTEL);
distFromLesion = xASL_im_DistanceTransform(logical(LesionIM));
distInLesion = xASL_im_DistanceTransform(logical(1-LesionIM));

%% Simple mixing by Henk
dist          = distFromLesion;
dist(dist==1) = 0;
dist(dist<0)  = 0;
dist          = 1-(dist./max(dist(:)));
dist          = dist.^8;
dist          = repmat(dist,[1 1 1 1 3]);

%% Save the SPM8 flow-field
if x.Seg.SaveSPMFlowField
    xASL_io_SaveNifti(Path_Transf_DARTEL, [x.P.Path_y_T1(1:(end-4)) '_SPM.nii'], IM_spm, [], false);
end

%% Save the original unaltered GS or DARTEL flow-field
if x.Seg.SaveOriginalFlowField
    xASL_io_SaveNifti(Path_Transf_DARTEL, [x.P.Path_y_T1(1:(end-4)) '_orig.nii'], IM_dartel, [], false);
end

%% Save the improved mixing
if x.Seg.SaveMixedFlowField && x.Quality
    % Looks at border points of the lesion
    pointsWithin = (distInLesion == 1);
    % Calculate the difference in transformation movement between the SPM and DARTEL/GS
    borderLineTrans = IM_dartel(repmat(pointsWithin,[1 1 1 1 3])) - IM_spm(repmat(pointsWithin,[1 1 1 1 3]));
    borderLineTrans = reshape(borderLineTrans,[],3);
    borderLineTrans = sqrt(sum(borderLineTrans.^2,2));

    % Finds the maximum discrepancy = MAXdist
    distMax = round(max(borderLineTrans)/1.5);

    % Finds the maximum thickness of the lesion = MAXlesion
    distMaxLesion = max(distInLesion(:));
    distMaxLesion = floor(distMaxLesion/4);

    % Sets the fixed borders as the MAX and 25% of MAXlesion within the lesion
    distMap = distFromLesion;
    distMap(distMap > distMax) = distMax;
    distMap = distMap + distMaxLesion;
    distMap(distInLesion>distMaxLesion) = 0;
    indx = (distInLesion<=distMaxLesion) .* (distInLesion>0);
    distMap(indx > 0) = distMaxLesion-distInLesion(indx > 0)+1;

    % Sets this as 0 and 1 and interpolate the registration transformations quadratically within
    distMap = distMap/(distMax+distMaxLesion);
    distMap(distMap>1) = 1;
    distMap(distMap<0) = 0;
    distMap = distMap.^2;
    distMap = repmat(distMap,[1 1 1 1 3]);
    IM_mixed = (1-distMap).*IM_spm + distMap.*IM_dartel;
    xASL_io_SaveNifti(Path_Transf_DARTEL,[x.P.Path_y_T1(1:(end-4)) '_mixed.nii'], IM_mixed, [], false);
end

%% mix the two flow fields, in non-linear proportions
% Deprecated function. In case of a large deformation - the IM_spm field will be completely off, pointing to
% a wrong location, because it will not manage to find the tumor. It might be more nicely regularized, but
% it will point to a wrong location, because the deformation will not be captured.
if x.Quality
    IM_dartel = dist.*IM_spm + (1-dist).*IM_dartel;
end
%IM_dartel_old       = dist.*IM_spm + (1-dist).*IM_dartel;
%xASL_io_SaveNifti(Path_Transf_DARTEL,[x.P.Path_y_T1(1:(end-4)) '_old.nii'],IM_dartel_old,[],0);

% Instead, we might assume linear interpolation of the border shift - linear deformation of the tumor according to the
% position of the border of it
%xASL_io_SaveNifti(Path_Transf_DARTEL,[x.P.Path_y_T1(1:(end-4)) '_orig.nii'],IM_dartel,[],0);
%IM_dartel(repmat(logical(LesionIM),[1 1 1 1 3])) = NaN;

%% Fix the edges of the flowfields, extrapolate over the NaNs
Fields2Save = IM_dartel;


end
