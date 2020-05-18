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

x = xASL_init_InitializeMutex(x, 'T1' ); % starts mutex locking process to ensure that everything will run only once
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

%% -----------------------------------------------------------------------------
%% Check whether we require certain files to be present, or otherwise to skip this subject
% By default we don't check this (because a population might not contain a
% FLAIR or M0, or this differs between subjects). This is useful when the
% data is incomplete, but one wants to start image processing nevertheless



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
%% 8    Resample to standard space
iState = 8;

    if ~x.mutex.HasState(StateName{iState}) && xASL_exist(x.P.Path_FLAIR,' file')

        xASL_wrp_Resample2StandardSpace(x);

        x.mutex.AddState(StateName{iState});
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
%% 999 Ready
x.mutex.AddState('999_ready');
x.mutex.Unlock();
cd(oldFolder);
result = true;


end
