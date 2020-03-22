function [result, x] = xASL_module_Structural(x)
%xASL_module_Structural ExploreASL module for structural processing
%
% FORMAT: [result, x] = xASL_module_Structural(x)
%
% INPUT:
%   x  - x structure containing all input parameters (REQUIRED)
%   x.SUBJECTDIR  -  anatomical directory, containing the derivatives of anatomical images (REQUIRED)
%
%
% OUTPUT:
%   result  - true for successful run of this module, false for insuccessful run
%   x       - x structure containing all output parameters (REQUIRED)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This first ExploreASL module processes the structural
% images, i.e. high-resolution T1w and FLAIR (if present), on an individual (i.e. subject-to-subject) basis.
% If a FLAIR is present, this is processed first to obtain a WMH mask to fill the hypointense lesions on the T1w,
% before segmenting the T1w. For the T1w segmentation this module uses CAT12
% by default but if this fails it falls back to SPM after trying to
% optimize the image contrast. This module has the following steps/submodules/wrappers:
%
% 010_LinearReg_T1w2MNI         - Ensure the alignment of subjects' anterior commissure (AC) with the AC in MNI & apply this to all images
% 020_LinearReg_FLAIR2T1w       - Align the FLAIR (if present) with T1w
% 030_FLAIR_BiasfieldCorrection - Perform a biasfield correction (if not performed  by LST in following steps)
% 040_LST_Segment_FLAIR_WMH     - Segment WMH lesions on FLAIR (if present)
% 050_LST_T1w_LesionFilling_WMH - Use WMH segmentation to fill lesions on T1w
% 060_Segment_T1w               - Tissue segmentation on T1w
% 070_CleanUpWMH_SEGM           - Extra WMH cleanup of some over- and under-segmentation
% 080_Resample2StandardSpace    - Clone all images to standard space
% 090_GetVolumetrics            - Obtain whole-brain volumes of GM, WM, CSF, WMH
% 100_VisualQC                  - Obtain QC parameters & save QC Figures
% 110_DoWADQCDC                 - QC for WAD-QC DICOM server (OPTIONAL)
%
% EXAMPLE: [~, x] = xASL_module_Structural(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2019 ExploreASL



%% -------------------------------------------------------------------------------------------
%% Administration

x = xASL_init_GenericMutexModules(x, 'T1' ); % starts mutex locking process to ensure that everything will run only once
x = xASL_init_FileSystem(x); % initialize FileSystem, quick & dirty
oldFolder = cd(x.SUBJECTDIR); % make sure that unspecified output will go here

if x.mutex.HasState('999_ready')
    bO = false; % no Output, as everything has been done already
else
    bO = true; % yes Output, not completely done, so we want to know what has or has not been done
end

if ~isfield(x, 'DoWADQCDC')
    x.DoWADQCDC = false;
end

% Start with throwing a warning if the current subject doesn"t have a structural image
if ~xASL_exist(x.P.Path_T1,'file') && ~xASL_exist(x.P.Path_T1_ORI,'file')
    % First try to find one with a more BIDS-compatible name & rename it (QUICK & DIRTY FIX)
    Flist = xASL_adm_GetFileList(x.SUBJECTDIR, '(.*T1.*run.*|.*run.*T1.*)\.(nii|nii\.gz)', 'FPList', [0 Inf]);
    if ~isempty(Flist)
        xASL_Move(Flist{1}, x.P.Path_T1);
    else
        warning([x.SUBJECTDIR ' didnt contain a T1w structural image, skipping...']);
        result = true;
        return;
    end
end

if ~xASL_exist(x.P.Path_FLAIR,'file') && ~xASL_exist(x.P.Path_FLAIR_ORI,'file')
    % First try to find one with a more BIDS-compatible name & rename it (QUICK & DIRTY FIX)
    Flist = xASL_adm_GetFileList(x.SUBJECTDIR, '(.*FLAIR.*run.*|.*run.*FLAIR.*)\.(nii|nii\.gz)', 'FPList', [0 Inf]);
    if ~isempty(Flist)
        xASL_Move(Flist{1}, x.P.Path_FLAIR);
    end
end

%% Which WMH segmentation algorithm to use
% Rules are here:
% 1) for existing WMH_SEGM we force a dummy LPA run, which is
%    faster, and will be overwritten later by the pre-existing WMH_SEGM
% 2) default is 'LPA' unless LGA was chosen
if xASL_exist(x.P.Path_WMH_SEGM, 'file') && bO
    x.WMHsegmAlg = 'LPA';
	% force fast WMH segmentation, because already WMH_SEGM.nii existing, which we assume to be high quality
    fprintf('%s\n', 'Warning: existing WMH_SEGM.nii detected, so using this WMH_SEGM.nii instead of segmenting ourselves');
    fprintf('%s\n', 'However, if this WMH_SEGM.nii was created in a previous run, consider removing it first!');
elseif xASL_exist(x.P.Path_WMH_SEGM, 'file') && ~bO
    % suppress output when this subject is ready
    x.WMHsegmAlg = 'LPA';
elseif ~isfield(x,'WMHsegmAlg') || isempty(x.WMHsegmAlg)
    x.WMHsegmAlg = 'LPA'; % default WMH segmentation algorithm, seems more robust & faster than LGA
elseif strcmp(x.WMHsegmAlg, 'LGA')
    % keep LGA
elseif xASL_exist(x.P.Path_FLAIR, 'file')
    warning('Unknown WMH segmentation option, selecting LPA');
    x.WMHsegmAlg = 'LPA';
end
switch x.WMHsegmAlg
    case 'LPA'
        rWMHPath = fullfile(x.SUBJECTDIR, ['ples_lpa_mr' x.P.FLAIR '.nii']);
    case 'LGA'
        rWMHPath = fullfile(x.SUBJECTDIR, ['ples_lga_0.3_rmr' x.P.FLAIR '.nii']);
end


%% -----------------------------------------------------------------------------
%% Check whether we require certain files to be present, or otherwise to skip this subject
% By default we don't check this (because a population might not contain a
% FLAIR or M0, or this differs between subjects). This is useful when the
% data is incomplete, but one wants to start image processing nevertheless

if ~isfield(x,'SkipIfNoFlair'); x.SkipIfNoFlair = 0; end
if ~isfield(x,'SkipIfNoASL');   x.SkipIfNoASL   = 0; end
if ~isfield(x,'SkipIfNoM0');    x.SkipIfNoM0    = 0; end

bContinue = true;
if x.SkipIfNoFlair && ~xASL_exist(x.P.Path_FLAIR,'file')
    bContinue = false;
end
for iSess=1:x.nSessions % need to define various ASL4D_session files still
    if ~xASL_exist(x.P.Path_ASL4D,'file') && x.SkipIfNoASL
       bContinue = false;
    elseif ~xASL_exist(x.P.Path_M0,'file') && x.SkipIfNoM0
       bContinue = false;
    end
end


if ~bContinue
    fprintf('%s\n','Skipping this subject, because no FLAIR, ASL or M0 present whereas this was requested');
    x.mutex.Unlock();
    cd(oldFolder);
    result = false;
    return;
end


StateName{ 1} = '010_LinearReg_T1w2MNI';
StateName{ 2} = '020_LinearReg_FLAIR2T1w';
StateName{ 3} = '030_FLAIR_BiasfieldCorrection';
StateName{ 4} = '040_LST_Segment_FLAIR_WMH';
StateName{ 5} = '050_LST_T1w_LesionFilling_WMH';
StateName{ 6} = '060_Segment_T1w';
StateName{ 7} = '070_CleanUpWMH_SEGM';
StateName{ 8} = '080_Resample2StandardSpace';
StateName{ 9} = '090_GetVolumetrics';
StateName{10} = '100_VisualQC_Structural';
StateName{11} = '110_DoWADQCDC';


%% -----------------------------------------------------------------------------
%% Check if we need to reset the original T1w/FLAIR,
%  This is the case if any of the T1w .status before segmentation is missing

if ~x.mutex.HasState(StateName{1}) || ~x.mutex.HasState(StateName{5})

    % if we rerun the structural module, then clean all modules (as they
    % are dependent of the structural module for segmentation/registration)
    % don't remove any previous WMH_SEGM
    xASL_adm_CleanUpBeforeRerun(x.D.ROOT, [1 2], false, false, x.P.SubjectID);
    
    if xASL_exist(x.P.Path_T1, 'file') % Fix multiple T1w iterations
        IM = xASL_io_Nifti2Im(x.P.Path_T1);
        if size(IM,4)>1 || size(IM,5)>1 || size(IM,6)>1 || size(IM,7)>1
            warning(['Too many dims, using first: ' x.P.Path_T1]);
            xASL_Copy(x.P.Path_T1, x.P.Path_T1_ORI, true);
            xASL_io_SaveNifti(x.P.Path_T1, x.P.Path_T1, IM(:,:,:,1,1,1,1), [], false);
        end
    end
    if xASL_exist(x.P.Path_FLAIR, 'file') % Fix multiple T1w iterations
        IM = xASL_io_Nifti2Im(x.P.Path_FLAIR);
        if size(IM,4)>1 || size(IM,5)>1 || size(IM,6)>1 || size(IM,7)>1
            warning(['Too many dims, using first: ' x.P.Path_T1]);
            xASL_Move(x.P.Path_FLAIR, x.P.Path_FLAIR_ORI, true);
            xASL_io_SaveNifti(x.P.Path_FLAIR, x.P.Path_FLAIR, IM(:,:,:,1,1,1,1), [], false);
        end
    end
end



%% -----------------------------------------------------------------------------
%% 1    Register T1w -> MNI
iState = 1;
if ~x.mutex.HasState(StateName{iState}) % tracks progress through lock/*.status files, & locks current run

    xASL_wrp_LinearReg_T1w2MNI(x, x.bAutoACPC);

    x.mutex.AddState(StateName{iState});
    xASL_adm_CompareDataSets([], [], x); % unit testing
    x.mutex.DelState(StateName{iState+1});
else
	xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	if  bO
		fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);
	end
end



%% -----------------------------------------------------------------------------
%% 2    Register FLAIR -> T1w (if FLAIR.nii exists)
iState = 2;
if ~x.mutex.HasState(StateName{iState}) % tracks progress through lock/*.status files, & locks current run
    if xASL_exist(x.P.Path_FLAIR, 'file')

        xASL_wrp_LinearReg_FLAIR2T1w(x, x.bAutoACPC);

        x.mutex.AddState(StateName{iState});
        xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState(StateName{iState+1});
	else
		if bO; fprintf('%s\n',[StateName{iState} ': No FLAIR data found, skipping...']);end
    end

else
	if xASL_exist(x.P.Path_FLAIR, 'file')
		xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	end
	if bO
		fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);
	end
end



%% -----------------------------------------------------------------------------
%% 3  T1w-based FLAIR biasfield correction (temporary disabled, deprecated?)
iState = 3;
% Here we acquire a biasfield from the T1 and apply it to the FLAIR, as with large lesions, the biasfield
% cannot be reliably estimated on the FLAIR. By default = disabled, enable this in case of large lesions.

if xASL_exist(x.P.Path_WMH_SEGM,'file') || ~xASL_exist(x.P.Path_FLAIR,'file')
    if bO; fprintf('%s\n','WMH segmentation already existing, or FLAIR.nii missing, skipping FLAIR biasfield correction'); end % We skip now
else
    if ~x.mutex.HasState(StateName{iState})
		if xASL_exist(x.P.Path_FLAIR,'file')
			if  isfield(x,'T1wBasedFLAIRBiasfieldCorrection') && x.T1wBasedFLAIRBiasfieldCorrection==1
				xASL_wrp_FLAIR_BiasfieldCorrection(x);
			else
				fprintf('%s\n',[StateName{iState} ' not requested, skipping...']); % this is off by default
			end

			x.mutex.AddState(StateName{iState}); % tracks progress through lock/*.status files, & locks current run
			xASL_adm_CompareDataSets([], [], x); % unit testing
			x.mutex.DelState(StateName{iState+1});
		else
			if bO; fprintf('%s\n',[StateName{iState} ': No FLAIR data found, skipping...']);end
		end
	else
		if xASL_exist(x.P.Path_FLAIR,'file')
			xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
		end
		if  bO; fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);end
    end
end





%% -----------------------------------------------------------------------------
%% 4    Lesion Segmentation Toolbox: segment FLAIR.nii if available
iState = 4;
if ~x.mutex.HasState(StateName{iState})  % tracks progress through lock/*.status files, & locks current run
    if xASL_exist(x.P.Path_FLAIR,'file')

        xASL_wrp_LST_Segment_FLAIR_WMH(x, rWMHPath, x.WMHsegmAlg);

        x.mutex.AddState(StateName{iState});
        xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState(StateName{iState+1});

    elseif bO; fprintf('%s\n',[StateName{iState} ': No FLAIR data found, skipping...']);
    end
else
	if xASL_exist(x.P.Path_FLAIR,'file')
		xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	end
	if  bO; fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);end
end





%% -----------------------------------------------------------------------------
%% 5    Lesion filling
iState = 5;
if ~x.mutex.HasState(StateName{iState}) && x.bLesionFilling  % tracks progress through lock/*.status files, & locks current run

    if xASL_exist(rWMHPath,'file')
        % check if the FLAIR segmentation exists from LST
        xASL_wrp_LST_T1w_LesionFilling_WMH(x, rWMHPath);

        x.mutex.AddState(StateName{iState});
        xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState(StateName{iState+1});
	else
		if bO; fprintf('%s\n','050_LesionFilling: No FLAIR data found, skipping...'); end
    end

else
	if xASL_exist(rWMHPath,'file')
		xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	end
	if  bO; fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);end
end







%% -----------------------------------------------------------------------------
%% 6    CAT12 segmentation
iState = 6;
% CAT12 outperforms SPM12. Therefore, always run CAT12, unless this crashes, then we try SPM12

if ~isfield(x,'Segment_SPM12')
    x.Segment_SPM12 = false; % by default, use CAT12, not SPM12 for segmentation
end
if ~isfield(x,'bFixResolution')
    x.bFixResolution = false; % by default, keep the original resolution
end

% Now check if the segmentation results exist
catVolFile = fullfile(x.D.TissueVolumeDir,['cat_' x.P.STRUCT '_' x.P.SubjectID '.mat']);
MatFile   = fullfile(x.SUBJECTDIR, [x.P.STRUCT '_seg8.mat']);
VolumetricResultsMissing = ~xASL_exist(catVolFile, 'file') && ~xASL_exist(MatFile, 'file');
Reprocessing = ~xASL_exist(x.P.Path_c1T1, 'file') || ~xASL_exist(x.P.Path_c2T1, 'file') || VolumetricResultsMissing;

if x.mutex.HasState(StateName{iState}) && Reprocessing
    warning('Segmentation results were missing, repeating');
end
if ~x.mutex.HasState(StateName{iState}) || Reprocessing

    x = xASL_wrp_SegmentT1w(x, x.Segment_SPM12);

    x.mutex.AddState(StateName{iState});  % tracks progress through lock/*.status files, & locks current run
    xASL_adm_CompareDataSets([], [], x); % unit testing
    x.mutex.DelState(StateName{iState+1});
else
	xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	if bO; fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);end
end



%% -----------------------------------------------------------------------------
%% 7    CleanUp WMH segmentation
iState = 7;
if ~x.mutex.HasState(StateName{iState})
    if xASL_exist(x.P.Path_WMH_SEGM, 'file')

        UniqueN = unique(xASL_io_Nifti2Im(x.P.Path_WMH_SEGM));
        
        if length(UniqueN)==2 && UniqueN(1)==0 && UniqueN(2)==1
            warning('Skipping WMH cleanup, because WMH_SEGM seems imported as binary mask');
        else
            xASL_wrp_CleanUpWMH_SEGM(x);
        end

        x.mutex.AddState(StateName{iState});
        xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState(StateName{iState+1});
    end
else
	if xASL_exist(x.P.Path_WMH_SEGM, 'file')
		xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	end
end

%% -----------------------------------------------------------------------------
%% 8    Resample to standard space
iState = 8;
FileMissing = false; % Check availability NIfTIs in standard space, else redo warping
List2Check = {x.P.Pop_Path_rT1 x.P.Pop_Path_rc1T1 x.P.Pop_Path_rc2T1};
for iCheck=1:length(List2Check)
    if ~xASL_exist(List2Check{iCheck},'file')
        FileMissing = true;
    end
end

if ~x.mutex.HasState(StateName{iState}) || FileMissing
    if ~x.mutex.HasState(StateName{6})
        fprintf('%s\n',[StateName{6} ' not performed, ' StateName{iState} ' skipped']);
    else % Run only if hasn't been run, if segmentation has been run & rigid-body registrations (for all volumes) have been run

        xASL_wrp_Resample2StandardSpace(x);

        x.mutex.AddState(StateName{iState});
        xASL_adm_CompareDataSets([], [], x); % unit testing
        x.mutex.DelState(StateName{iState+2}); % resampling has no effect on the tissue volumetrics, those are copied from CAT12
    end

else
	xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	if  bO; fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);end
end



%% -----------------------------------------------------------------------------
%% 9    Get WMH & tissue volumes
iState = 9;
Reprocessing = false;
WMHvolFileMissing = isempty(xASL_adm_GetFileList(x.D.TissueVolumeDir,['^WMH_LST_(LGA|LPA)_' x.P.SubjectID '\.tsv$'],'FPList',[0 Inf]));
if xASL_exist(x.P.Path_FLAIR,'file') && WMHvolFileMissing && x.mutex.HasState(StateName{iState})
    Reprocessing = true;
    warning(['WMH/LST volumetric results were missing, rerunning ' StateName{iState}]);
end
T1volFileMissing = ~exist(fullfile(x.D.TissueVolumeDir,['TissueVolume_' x.P.SubjectID '.tsv']), 'file');
if xASL_exist(x.P.Path_T1,'file') && T1volFileMissing && x.mutex.HasState(StateName{iState})
    Reprocessing = true;
    warning(['T1/CAT12 volumetric results were missing, rerunning ' StateName{iState}]);
end
if ~x.mutex.HasState(StateName{iState}) || Reprocessing % tracks progress through lock/*.status files, & locks current run

    % PM: this should be corrected for any
    % Lesion_(T1|FLAIR)_(1|2|3|..)\.nii files

    if ~x.mutex.HasState(StateName{6})
        fprintf('%s\n', [StateName{6} ' not performed, ' StateName{iState} ' skipped']);
    else % run only if segmentation has been run

        xASL_wrp_GetVolumetrics(x);
        x.mutex.AddState(StateName{iState});
    end
else
	xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	if  bO; fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);end
end





%% -----------------------------------------------------------------------------
%% 10    Visual QC
iState = 10;
if ~x.mutex.HasState(StateName{iState}) % tracks progress through lock/*.status files, & locks current run
    if ~x.mutex.HasState(StateName{6}) || ~x.mutex.HasState(StateName{8})
        fprintf('%s\n', [StateName{6} ' or ' StateName{8} ' not performed, skipping ' StateName{iState}]);
    else

        xASL_wrp_VisualQC_Structural(x);

        x.mutex.AddState(StateName{iState});
        x.mutex.DelState(StateName{iState+1});
    end

else
	xASL_adm_CompareDataSets([], [], x,2,StateName{iState}); % unit testing - only evaluation
	if  bO; fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);end
end

%% -----------------------------------------------------------------------------
%% 11    WAD-QC
iState = 11;

if ~x.mutex.HasState(StateName{iState}) && x.DoWADQCDC
    xASL_qc_WADQCDC(x, x.iSubject, 'Structural');
    x.mutex.AddState(StateName{iState});
elseif x.mutex.HasState(StateName{iState}) && bO
    fprintf('%s\n', [StateName{iState} ' has already been performed, skipping...']);
end




%% -----------------------------------------------------------------------------
%% 999 Ready
x.mutex.AddState('999_ready');
x.mutex.Unlock();
cd(oldFolder);
result = true;


end
