function [result, x] = xASL_module_Population(x)
%xASL_module_Population ExploreASL module for population-based/group-based processing
%
% FORMAT: [result, x] = xASL_module_Population(x)
%
% INPUT:
%   x       - x structure containing all input parameters (REQUIRED)
%
% OUTPUT:
%   result  - true for successful run of this module, false for insuccessful run
%   x       - x structure containing all output parameters
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL module processes all available images on the
% group level. It assumes that all images were adequately processed in the
% previous modules. It will perform the following group-wise processing and
% checks:
%
% - 010_CreatePopulationTemplates - Create population average images, to compare scanners, cohorts etc without physiological variance
% - 020_CreateAnalysisMask        - Generate a group-level mask by combining individuals masks, for ROI-based analysis & VBA
% - 030_CreateBiasfield           - When there are multiple scanners, create scanner-specific biasfields (uses Site.mat for this)
% - 040_GetDICOMStatistics        - Create TSV file with overview of DICOM parameters
% - 050_GetVolumeStatistics       - Create TSV file with overview of volumetric parameters
% - 060_GetMotionStatistics       - Create TSV file with overview of motion parameters
% - 065_GetRegistrationStatistics - Create TSV file with overview of the registration statistics
% - 070_GetROIstatistics          - Create TSV file with overview of regional values (e.g. qCBF, mean control, pGM etc)
% - 080_SortBySpatialCoV          - Sort ASL_Check QC images by their spatial CoV in quality bins
% - 090_DeleteAndZip              - Delete temporary files and gzip all NIfTIs
%
% EXAMPLE: [~, x] = xASL_module_Population(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2020 ExploreASL



%% ------------------------------------------------------------------------------------------------------------
%% Admin

if x.iWorker>1 % run population module only once when ExploreASL is called multiple times in parallel
    warning(['I am worker ' xASL_num2str(x.iWorker) ', population module should not run in parallel, skipping']);
    result = true;
end

if ~isfield(x,'bNativeSpaceAnalysis') || isempty(x.bNativeSpaceAnalysis)
    x.bNativeSpaceAnalysis = 0;
end

% Check if we have ASL or not, to know if we need to run ASL-specific stuff/warnings
bHasASL = ~isempty(xASL_adm_GetFileList(x.D.PopDir, '^.*ASL_\d\.nii$'));
if ~bHasASL
    warning('Detected no ASL scans, skipping ASL-specific parts of the Population module');
end

x = xASL_init_InitializeMutex(x, 'QC'); % starts mutex locking process to ensure that everything will run only once
x = xASL_init_FileSystem(x);
xASL_adm_CreateDir(x.S.StatsDir);

% Obtain ASL sequence
x = xASL_adm_DefineASLSequence(x);

StateName{1}  = '010_CreatePopulationTemplates';
StateName{2}  = '020_CreateAnalysisMask';
StateName{3}  = '030_CreateBiasfield';
StateName{4}  = '040_GetDICOMStatistics';
StateName{5}  = '050_GetVolumeStatistics';
StateName{6}  = '060_GetMotionStatistics';
StateName{7}  = '065_GetRegistrationStatistics';
StateName{8}  = '070_GetROIstatistics';
StateName{9}  = '080_SortBySpatialCoV';
StateName{10} = '090_DeleteAndZip';



%% ------------------------------------------------------------------------------------------------------------
%% 1    Create template images
if ~x.mutex.HasState(StateName{1})
    xASL_wrp_CreatePopulationTemplates(x);  % this doesn't work nicely yet with sessions, should be changed after new BIDS is implemented

    % Save FoV mask as susceptibility mask for 3D spiral
    % as 3D spiral doesnt have a susceptibility artifact (or negligible)

    FoVPath = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^FoV_n' xASL_num2str(x.nSubjectsSessions) '_bs-mean_Unmasked\.nii$'], 'FPList');
    SusceptPath = fullfile(x.D.TemplatesStudyDir,['MaskSusceptibility_n' xASL_num2str(x.nSubjectsSessions) '_bs-mean.nii']);

    if strcmpi(x.Sequence,'3d_spiral') && ~isempty(FoVPath)
        xASL_io_SaveNifti(FoVPath{1}, SusceptPath, xASL_io_Nifti2Im(FoVPath{1}), [], false);
    end

    x.mutex.AddState(StateName{1});
    fprintf('%s\n',[StateName{1} ' was performed']);
else
    fprintf('%s\n',[StateName{1} ' has already been performed, skipping...']);
end


%% General settings
x = xASL_adm_CreateFileReport(x);
% xASL_wrp_PVC_HiRes( x ); % PVEc correction in standard space high resolution, using B-splines

%% ------------------------------------------------------------------------------------------------------------
%% 2    Create population-based analysis mask for ROI-based analysis & VBA
if ~x.mutex.HasState(StateName{2}) && bHasASL
    x = xASL_im_CreateAnalysisMask(x);
    x.mutex.AddState(StateName{2});
    fprintf('%s\n',[StateName{2} ' was performed']);
elseif bHasASL
    fprintf('%s\n',[StateName{2} ' has already been performed, skipping...']);
end



%% -----------------------------------------------------------------------------
%% 3    Multi-sequence equalization
if ~x.mutex.HasState(StateName{3}) && bHasASL
    xASL_wrp_CreateBiasfield(x); % later to include: smoothness equalization, geometric distortion correction etc
    x.mutex.AddState(StateName{3});
    fprintf('%s\n',[StateName{3} ' was performed']);
elseif bHasASL
    fprintf('%s\n',[StateName{3} ' has already been performed, skipping...']);
end




%% -----------------------------------------------------------------------------
%% 4    Print DICOM header parameters & check whether there are outliers
if ~x.mutex.HasState(StateName{4})
    ScanType = {'ASL4D' 'M0'};
    HasSessions = {1 1};

    for iType=1:length(ScanType)
        xASL_stat_GetDICOMStatistics(x, ScanType{iType}, HasSessions{iType});
    end

    xASL_stat_GetAcquisitionTime(x); % This provides an overview of Acquisition times

    x.mutex.AddState(StateName{4});
    fprintf('%s\n',[StateName{4} ' was performed']);
else
    fprintf('%s\n',[StateName{4} ' has already been performed, skipping...']);
end


%% -----------------------------------------------------------------------------
%% 5    Summarize volume statistics (uses native space)
if ~x.mutex.HasState(StateName{5})

    xASL_stat_GetVolumeStatistics(x);

    x.mutex.AddState(StateName{5});
    fprintf('%s\n',[StateName{5} ' was performed']);
else
    fprintf('%s\n',[StateName{5} ' has already been performed, skipping...']);
end



%% -----------------------------------------------------------------------------
%% 6    Summarize motion statistics (using generated net displacement vector (NDV) motion results from ASL-realign module)
if ~x.mutex.HasState(StateName{6}) && bHasASL
    try
        xASL_stat_GetMotionStatistics(x);
        x.mutex.AddState(StateName{6});
        fprintf('%s\n',[StateName{6} ' was performed']);
    catch ME
        warning('Motion summarizing failed:');
        fprintf('%s\n',ME.message);
    end
elseif bHasASL
    fprintf('%s\n',[StateName{6} ' has already been performed, skipping...']);
end

%% -----------------------------------------------------------------------------
%% 6.5  Summarize registration statistics (using the Tanimoto coefficients calculated in the ASL and Structural submodules)
if ~x.mutex.HasState(StateName{7})
    try
        xASL_stat_GetRegistrationStatistics(x);
        x.mutex.AddState(StateName{7});
        fprintf('%s\n',[StateName{7} ' was performed']);
    catch ME
        warning('Registration summarizing failed:');
        fprintf('%s\n',ME.message);
    end
else
    fprintf('%s\n',[StateName{7} ' has already been performed, skipping...']);
end


%% -----------------------------------------------------------------------------
%% 7    ROI statistics
if ~x.mutex.HasState(StateName{8})
    
    % Default & Fallback for Atlases
    if ~isfield(x.S,'Atlases')
        x.S.Atlases = {'TotalGM','DeepWM'}; % Default
    end

    x = xASL_init_LoadMetadata(x); % Add statistical variables, if there are new ones
    % if exist('ASL','var')
    %     xASL_vis_OverlapT1_ASL(x, ASL.Data.data); % Overlap T1 GM probability map & CBF, Create image showing spatial/visual agreement between T1 GM segmentation & ASL
    % end

    xASL_stat_ComputeWsCV(x); % This computes wsCV & bsCV to compute power   
    
    % Iterate over atlases in cell structure
    
    % ROI statistics
    % x.S.SubjectWiseVisualization =1; % set this on to visualize the subject-wise masks
    % over CBF maps (takes lot of extra time though)
    
    % x.S.InputDataStr              = 'SD'; % 'SD' 'TT' 'M0' 'R1' 'ASL_HctCohort' 'ASL_HctCorrInd'
    % x.S.InputAtlasPath            = fullfile(x.D.MapsSPMmodifiedDir,'WholeBrain.nii');
    % % x.S.InputAtlasPath            = fullfile(x.D.AtlasDir,'HOcort_CONN.nii');
    % xASL_wrp_GetROIstatistics( x);
    % x.S.InputDataStr              = 'SNR'; % 'SD' 'TT' 'M0' 'R1' 'ASL_HctCohort' 'ASL_HctCorrInd'
    % x.S.InputAtlasPath            = fullfile(x.D.MapsSPMmodifiedDir,'WholeBrain.nii');
    % xASL_wrp_GetROIstatistics( x);

    x.S.InputDataStr = 'qCBF'; % 'SD' 'TT' 'M0' 'R1' 'ASL_HctCohort' 'ASL_HctCorrInd'
	x.S.InputDataStrNative = 'CBF'; % 'SD' 'TT' 'M0' 'R1' 'ASL_HctCohort' 'ASL_HctCorrInd'
    
    % Iterate over atlases
    for iAtlas=1:length(x.S.Atlases)        
        
        % Check if atlas name is in path list
        if isfield(x.P.Atlas,x.S.Atlases{iAtlas})
            x.S.InputAtlasPath = x.P.Atlas.(x.S.Atlases{iAtlas});
        end
        
        % ROI statistics (default: standard space)
        x.S.InputNativeSpace = 0;
        xASL_wrp_GetROIstatistics(x);
        % ROI statistics (optional: native space)
        if x.bNativeSpaceAnalysis
            x.S.InputNativeSpace = 1;
            x.S.InputAtlasNativeName = x.S.Atlases{iAtlas};
            xASL_wrp_GetROIstatistics(x);
        end
        
        
    end
    
    % Check if we should do the same for Lesion or ROI masks (i.e.
    % individual "atlases") -> NB not yet developed/tested in native space
    LesionROIList = xASL_adm_GetFileList(x.D.PopDir, 'r(Lesion|ROI).*\.nii$', 'FPList', [0 Inf]);
    x.S.InputNativeSpace = 0;
    for iAtlas = 1:length(LesionROIList)
        x.S.InputAtlasPath = LesionROIList{iAtlas};
        xASL_wrp_GetROIstatistics(x);
    end
    
    % % Do the same for volumetrics
    % x.S.IsASL = false;
    % x.S.IsVolume = true;
    % x.S.InputDataStr = 'mrc1T1'; % GM volume
    % x.S.InputAtlasPath = fullfile(x.D.AtlasDir,'HOcort_CONN.nii');
    % xASL_wrp_GetROIstatistics(x);
    % x.S.InputAtlasPath = fullfile(x.D.AtlasDir,'HOsub_CONN.nii');
    % xASL_wrp_GetROIstatistics(x);
    % x.S.InputAtlasPath = fullfile(x.D.AtlasDir,'Hammers.nii');
    % xASL_wrp_GetROIstatistics(x);

    x.mutex.AddState(StateName{8});
    fprintf('%s\n',[StateName{8} ' was performed']);
else
    fprintf('%s\n',[StateName{8} ' has already been performed, skipping...']);
end

%% -----------------------------------------------------------------------------
%% 8    QC categorization based on spatial CoV:
if ~x.mutex.HasState(StateName{9}) && bHasASL
    xASL_qc_SortBySpatialCoV(x);

    % When this has been visually corrected, following function will obtain the QC categories
    % xASL_qc_ObtainQCCategoriesFromJPG(x);
    x.mutex.AddState(StateName{9});
    fprintf('%s\n',[StateName{9} ' was performed']);
elseif bHasASL
    fprintf('%s\n',[StateName{9} ' has already been performed, skipping...']);
end


%% -----------------------------------------------------------------------------
%% 9    Reduce data size
if ~x.mutex.HasState(StateName{10})
    if ~x.bReproTesting && x.DELETETEMP
        xASL_adm_DeleteManyTempFiles(x);
    end

    % Zip temporary files
    % This way, temporary files that take up a lot of data, will take up have
    % the disc space. These files may be re-used in the event of a
    % re-processing, which is why we zip them instead of deleting them.
    xASL_adm_GzipAllFiles(x.D.ROOT);
    x.mutex.AddState(StateName{10});
    fprintf('%s\n',[StateName{10} ' was performed']);
else
    fprintf('%s\n',[StateName{10} ' has already been performed, skipping...']);
end


%% -----------------------------------------------------------------------------
%% 999 Ready
x.mutex.AddState('999_ready');
x.mutex.Unlock();
result = true;
close all;

end
