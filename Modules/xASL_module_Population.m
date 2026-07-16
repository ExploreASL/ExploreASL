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
% - `010_CreatePopulationTemplates` - Create population average images, to compare scanners, cohorts etc without physiological variance
% - `020_CreateAnalysisMask`        - Generate a group-level mask by combining individuals masks, for ROI-based analysis & VBA
% - `030_CreateBiasfield`           - When there are multiple scanners, create scanner-specific biasfields (uses Site.mat for this)
% - `040_GetDICOMStatistics`        - Create TSV file with overview of DICOM parameters
% - `050_GetVolumeStatistics`       - Create TSV file with overview of volumetric parameters
% - `060_GetMotionStatistics`       - Create TSV file with overview of motion parameters
% - `065_GetRegistrationStatistics` - Create TSV file with overview of the registration statistics
% - `070_GetROIstatistics`          - Create TSV file with overview of regional values (e.g. qCBF, mean control, pGM etc)
%                                   7.a Perform statistics for normal atlases
%                                   7.b Perform statistics for Lesion and ROI files
%                                   7.c Parse TSVs & add to participants.tsv
%                                   7.d Generate the participants.json sidecar of participants.tsv
% - `080_SortBySpatialCoV`          - Sort ASL_Check QC images by their spatial CoV in quality bins
% - `090_DeleteTempFiles`           - Delete temporary files
% - `100_GZipAllFiles`              - Zip files to reduce disc space usage of temporary and non-temporay NIfTI files
%
% EXAMPLE: [~, x] = xASL_module_Population(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________


%% ------------------------------------------------------------------------------------------------------------
%% Admin
x = xASL_init_SubStructs(x);

% Input check
if x.opts.nWorkers>1 % don't run population module when ExploreASL is parallelized
    warning('Population module should not run in parallel, skipping...');
    fprintf('%s\n', 'Best to run this module after you have ran ExploreASL in parallel, by restarting ExploreASL non-parallel');
    result = true;
    return;
end

% Check again for Atlases (main checking is done when loading dataPar)
if ~isfield(x.S,'Atlases') && ~isfield(x.S, 'TissueMasking')
	% If no values are provided, then we provide the defaults
	% GM WM GM tissues with the Total & DeepWM & Tatu_ACA_MCA_PCA
	
	if isfield(x.S, 'TissueThreshold')
		% If TissueThresholds is set, then add Total-GM, Total-WM, and vascular territory atlas
		x.S.Atlases = {'Total', 'DeepWM', 'Tatu_ACA_MCA_PCA'}; 
		x.S.TissueMasking = {'GM', 'WM', 'GM'};
	else
		x.S.Atlases = {'Total', 'Total', 'DeepWM', 'Tatu_ACA_MCA_PCA', 'Tatu_ACA_MCA_PCA'}; 
		x.S.TissueMasking = {'GM', 'GM', 'WM', 'GM', 'GM'};
		x.S.TissueThreshold = [0.7, 0.5, 0.5, 0.7, 0.5];
	end

    % Note that Atlases and TissueMasking  should be in the same order as the atlases/ROIs
    % A mismatch (e.g. TissueMasking=GM for Atlases=deepWM) would result in an empty ROI, producing a NaN in the .tsv table
    % You can also use CSF or combinations like GM+WM or GM+CSF
elseif ~isfield(x.S, 'Atlases') || ~isfield(x.S, 'TissueMasking') || length(x.S.Atlases)~=length(x.S.TissueMasking)
	% Incorrect values are provided
	error('You need to provide x.S.Atlases and x.S.TissueMasking with the same length');
end

if ~isfield(x.S, 'TissueThreshold')
	% The default threshold is 0.7
	x.S.TissueThreshold = ones(1,length(x.S.TissueMasking)) * 0.7;
elseif length(x.S.TissueMasking) ~= length(x.S.TissueThreshold)
	error('x.S.TissueThreshold has to match the length of x.S.TissueMasking');
elseif max(x.S.TissueThreshold) > 1
	error('Maximum value for x.S.TissueThreshold is 1');
elseif min(x.S.TissueThreshold) < 0
	error('Minimum value for x.S.TissueThreshold is 0');
end

if ~isfield(x.S, 'LesionROIThreshold')
	% The default LesionROIThreshold is 0.5
	x.S.LesionROIThreshold = 0.5;
elseif max(x.S.LesionROIThreshold) > 1
	error('Maximum value for x.S.LesionROIThreshold is 1');
elseif min(x.S.LesionROIThreshold) < 0
	error('Minimum value for x.S.LesionROIThreshold is 0');
end

% Print the used atlases	
fprintf('\nThe following atlases have been selected with the following tissue masks and tissue thresholds:\n')
for iAtlas=1:length(x.S.Atlases)
    fprintf('%s\n', [x.S.TissueMasking{iAtlas} ' > ' num2str(x.S.TissueThreshold(iAtlas)) ' within ' x.S.Atlases{iAtlas} ' atlas']);
end
fprintf('\n');

% Default datatypes
if ~isfield(x.S,'DataTypes') || isempty(x.S.DataTypes)
	x.S.DataTypes = {'qCBF'}; % Default
	% Alternatives: 'Tex' 'ATT' 'SD' 'M0' 'ABV' 'ITT'
	% These can be added in the dataPar manually
end

% Create population directory
xASL_adm_CreateDir(x.D.PopDir);

if ~isfield(x.modules.population,'bNativeSpaceAnalysis') || isempty(x.modules.population.bNativeSpaceAnalysis)
    x.modules.population.bNativeSpaceAnalysis = 0;
end

% Check if we have ASL or not, to know if we need to run ASL-specific stuff/warnings
bHasASL = ~isempty(xASL_adm_GetFileList(x.D.PopDir, '^.*ASL_\d+\.nii$'));
if ~bHasASL
    warning('Detected no ASL scans, skipping ASL-specific parts of the Population module');
end

x = xASL_init_InitializeMutex(x, 'Population'); % starts mutex locking process to ensure that everything will run only once

if x.mutex.bAnyModuleLocked
    % If any module is locked, we skip this module
    result = true;
    return;
end


x = xASL_init_FileSystem(x);


StateName{1}  = '010_CreatePopulationTemplates';
StateName{2}  = '020_CreateAnalysisMask';
StateName{3}  = '030_CreateBiasfield';
StateName{4}  = '040_GetDICOMStatistics';
StateName{5}  = '050_GetVolumeStatistics';
StateName{6}  = '060_GetMotionStatistics';
StateName{7}  = '065_GetRegistrationStatistics';
StateName{8}  = '070_GetROIstatistics';
StateName{9}  = '080_SortBySpatialCoV';
StateName{10} = '090_DeleteTempFiles';
StateName{11} = '100_GZipAllFiles';


x.S.TemplateNumberName = ['_n' xASL_num2str(x.dataset.nSubjects)];


%% ------------------------------------------------------------------------------------------------------------
%% 1.   Create template images
if ~x.mutex.HasState(StateName{1})
    % This generates templates per session/run
    xASL_wrp_CreatePopulationTemplates(x);

    % Save FoV mask as susceptibility mask for 3D spiral
    % as 3D spiral doesnt have a susceptibility artifact (or negligible)

    FoVPath = xASL_adm_GetFileList(x.D.TemplatesStudyDir, ['^FoV' x.S.TemplateNumberName '_bs-mean_Unmasked\.nii$'], 'FPList');
    PathTemplateSusceptibilityMask = fullfile(x.D.TemplatesStudyDir,['MaskSusceptibility' x.S.TemplateNumberName '_bs-mean.nii']);

    if ~xASL_exist(PathTemplateSusceptibilityMask, 'file')
        warning('Susceptibility mask template was missing...');
        
        if ~isempty(FoVPath)
            xASL_io_SaveNifti(FoVPath{1}, PathTemplateSusceptibilityMask, xASL_io_Nifti2Im(FoVPath{1}), [], false);
            fprintf('and was replaced by FoV mask...\n');
        end
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
%% 2.   Create population-based analysis mask for ROI-based analysis & VBA
if ~x.mutex.HasState(StateName{2}) && bHasASL
    x = xASL_im_CreateGroupAnalysisMask(x);
    x.mutex.AddState(StateName{2});
    fprintf('%s\n',[StateName{2} ' was performed']);
elseif bHasASL
    fprintf('%s\n',[StateName{2} ' has already been performed, skipping...']);
end


%% -----------------------------------------------------------------------------
%% 3.   Multi-sequence equalization
if ~x.mutex.HasState(StateName{3}) && bHasASL
    xASL_wrp_CreateBiasfield(x); % later to include: smoothness equalization, geometric distortion correction etc
    x.mutex.AddState(StateName{3});
    fprintf('%s\n',[StateName{3} ' was performed']);
elseif bHasASL
    fprintf('%s\n',[StateName{3} ' has already been performed, skipping...']);
end


%% -----------------------------------------------------------------------------
%% 4.   Print DICOM header parameters & check whether there are outliers
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
%% 5.   Summarize volume statistics (uses native space)
if ~x.mutex.HasState(StateName{5})

    xASL_stat_GetVolumeStatistics(x);

    x.mutex.AddState(StateName{5});
    fprintf('%s\n',[StateName{5} ' was performed']);
else
    fprintf('%s\n',[StateName{5} ' has already been performed, skipping...']);
end


%% -----------------------------------------------------------------------------
%% 6.   Summarize motion statistics (using generated net displacement vector (NDV) motion results from ASL-realign module)
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
%% 6.5   Summarize registration statistics (using the Tanimoto coefficients calculated in the ASL and Structural submodules)
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
%% 7.    ROI statistics
%% 7.a   Perform statistics for normal atlases
if ~x.mutex.HasState(StateName{8})
    
    x = xASL_init_LoadMetadata(x); % Add statistical variables, if there are new ones
    % if exist('ASL','var')
    %     xASL_vis_OverlapT1_ASL(x, ASL.Data.data); % Overlap T1 GM probability map & CBF, Create image showing spatial/visual agreement between T1 GM segmentation & ASL
    % end

    xASL_stat_ComputeWsCV(x); % This computes wsCV & bsCV to compute power   

    % ROI statistics
    % x.S.SubjectWiseVisualization =1; % set this on to visualize the subject-wise masks
    % over CBF maps (takes lot of extra time though)
    
    % Iterate over DataTypes
    for iDataType=1:length(x.S.DataTypes)
        x.S.InputDataStr = x.S.DataTypes{iDataType};
    
        % Iterate over atlases
        x.dir.dirAtlas = fullfile(x.opts.MyPath, 'external', 'Atlases');

        for iAtlas=1:length(x.S.Atlases)
            % Note that the number of ROI atlases here should be the same the number of tissue masking chosen
            % If needed, an atlas or tissue type can be provided multiple times in different combinations

            % We use the specified tissue type
			% 'GM' = gray matter
            % 'WM' = white matter
			% 'CSF' = cerebrospinal fluid
            % 'GM+WM' = whole brain parenchyma GM+WM, previously was 'WB', alternatively can be defined as 'WM+GM'
			% 'GM+CSF' = GM+CSF combination, 'CSF+GM' does the same
			% 'WM+CSF' = WM+CSF combination
			% 'GM+WM+CSF' = GM+WM+CSF combination

            x.S.TissueMaskingLocal = x.S.TissueMasking{iAtlas};
			x.S.TissueThresholdLocal = x.S.TissueThreshold(iAtlas);
            
            % Find the path of the atlas
            pathAtlas = fullfile(x.dir.dirAtlas, [x.S.Atlases{iAtlas} '.nii']);
            
            % Check if atlas name is in path list
            if isfield(x.D.Atlas, x.S.Atlases{iAtlas})
                % Atlas is found in the default atlas list
                x.S.InputAtlasPath = x.D.Atlas.(x.S.Atlases{iAtlas});
            elseif xASL_exist(pathAtlas, 'file')
                    % try to find the atlas in the default folder
                    x.D.Atlas.(x.S.Atlases{iAtlas}) = pathAtlas;
                    x.S.InputAtlasPath = x.D.Atlas.(x.S.Atlases{iAtlas});
            else
                warning(['Unknown atlas: ' x.S.Atlases{iAtlas} ', skipping']);
            end

            % ROI statistics (default: standard space)
            x.S.InputNativeSpace = 0;
			x.S.bSubjectSpecificROI = false; % lesion/ROIs designated per subject (e.g., Lesion_T1_2.nii)
            % x.S.SubjectWiseVisualization = true; defaulted to false,
            % set this to true for visualization ROIs

            xASL_wrp_GetROIstatistics(x);
            % ROI statistics (optional: native space)
            if x.modules.population.bNativeSpaceAnalysis
                x.S.InputNativeSpace = 1;
                x.S.InputAtlasNativeName = [x.S.Atlases{iAtlas} '_Atlas'];
                xASL_wrp_GetROIstatistics(x);
            end
        end
        
		%% -----------------------------------------------------------------------------
		%% 7.b Perform statistics for Lesion and ROI files
		% Read the names of the lesion files
		LesionROIList = xASL_adm_GetFileList(x.D.PopDir, '(?i)^r(Lesion|ROI)_(T1|FLAIR|T2)_\d*_.*\.nii', 'List', [0 Inf]);
		% Go through the lesions and remove the subject names
		for iROI = 1:length(LesionROIList)
			[~, iEnd] = regexpi(LesionROIList{iROI}, '^r(Lesion|ROI)_(T1|FLAIR|T2)_\d*_');
			if isempty(iEnd)
				LesionROIList{iROI} = '';
			else
				LesionROIList{iROI} = LesionROIList{iROI}(1:iEnd);
			end
		end
		
		% Obtain a unique list of lesion names without the subject name
		LesionUniqueROIList = unique(LesionROIList);

		% Standard space analysis in a specific ROI with no tissue restriction
        x.S.InputNativeSpace = 0;
		x.S.bSubjectSpecificROI = true;
		x.S.TissueMaskingLocal = 'GM+WM+CSF';
		x.S.TissueThresholdLocal = x.S.LesionROIThreshold;
		for iROI = 1:length(LesionUniqueROIList)
            x.S.InputAtlasPath = fullfile(x.D.PopDir, LesionUniqueROIList{iROI});
            xASL_wrp_GetROIstatistics(x);
		end

		% Lesion/ROI statistics in native space
		if x.modules.population.bNativeSpaceAnalysis
			x.S.InputNativeSpace = 1;
			x.S.bSubjectSpecificROI = true;
			x.S.TissueMaskingLocal = 'GM+WM+CSF';
			x.S.TissueThresholdLocal = x.S.LesionROIThreshold;
			for iROI = 1:length(LesionUniqueROIList)
				x.S.InputAtlasPath = fullfile(x.D.PopDir, LesionUniqueROIList{iROI});
				% Remove 'r' at the start
				x.S.InputAtlasNativeName = LesionUniqueROIList{iROI}(2:end-1);
				xASL_wrp_GetROIstatistics(x);
			end
        end
    end

    %% -----------------------------------------------------------------------------
    %% 7.c Parse TSVs & add to participants.tsv
    regExp_Type = {'qCBF' 'ATT' 'M0' 'meanControl' 'ITT' 'Tex'};
    key_Type = {'cbf' 'att' 'm0' 'control' 'itt' 'tex'};
    regExp_Stats = {'mean' 'CoV'};
    regExp_Atlas = {'Total' 'DeepWM' 'Tatu_ACA_MCA_PCA'};
    regExp_Tissue = {'GM' 'WM' 'GM'};
    regExp_PVC = {'PVC0' 'PVC2'};

    for iType=1:length(regExp_Type)
        for iStat=1:length(regExp_Stats)
            for iAtlas=1:length(regExp_Atlas)
                for iPVC=1:length(regExp_PVC)
                    regExp = [regExp_Stats{iStat} '_' regExp_Type{iType} '.*StandardSpace_' regExp_Atlas{iAtlas} regExp_Tissue{iAtlas} '_n=' num2str(x.dataset.nSubjects) '_' date '_' regExp_PVC{iPVC} '\.tsv'];
                    fList = xASL_adm_GetFileList(fullfile(x.D.PopDir, 'Stats'), regExp, 'FPList');
                    if length(fList)>1
                        warning('Multiple stats files found to be added to participants.tsv, using the first only');
                    end
                    if ~isempty(fList)
                        % Addition to participants.tsv
                        tableIs = xASL_tsvRead(fList{1});
                        
                        % filter the bilateral values
                        iBilateral = find(cellfun(@(y) strcmp(y(end-1:end),'_B'), tableIs(1,:)));
                        tableROI = tableIs(1,iBilateral);
                        
                        tableSubjRuns = tableIs(3:end,1:2);
                        tableValues = tableIs(3:end,iBilateral);
                        
                        for iKey=1:length(tableROI)
                            dataIn = [tableSubjRuns tableValues(:, iKey)];
                            keyBIDS = [regExp_Stats{iStat} '_' key_Type{iType} '_' tableROI{iKey} '_' regExp_PVC{iPVC}];
                            xASL_bids_Add2ParticipantsTSV(dataIn, keyBIDS, x);
                        end
                    end
                        
                end
            end
        end
    end

    %% -----------------------------------------------------------------------------    
    %% 7.d Generate the participants.json sidecar of participants.tsv
    xASL_bids_GenerateParticipantsJSON(x);

    x.mutex.AddState(StateName{8});
    fprintf('%s\n',[StateName{8} ' was performed']);
else
    fprintf('%s\n',[StateName{8} ' has already been performed, skipping...']);
end


%% -----------------------------------------------------------------------------
%% 8.   QC categorization based on spatial CoV:
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
%% 9.  Reduce data size
if ~x.mutex.HasState(StateName{10})
    if ~x.settings.bReproTesting && x.settings.DELETETEMP
        xASL_adm_DeleteManyTempFiles(x);
    end
    x.mutex.AddState(StateName{10});
    fprintf('%s\n',[StateName{10} ' was performed']);
else
    fprintf('%s\n',[StateName{10} ' has already been performed, skipping...']);
end


%% 10.  xASL_adm_GzipAllFiles
if ~x.mutex.HasState(StateName{11})
    xASL_adm_GzipAllFiles(x.dir.xASLDerivatives,[],[],fullfile(x.opts.MyPath,'External'), false);
    x.mutex.AddState(StateName{11});
    fprintf('%s\n',[StateName{11} ' was performed']);
else
        fprintf('%s\n',[StateName{11} ' has already been performed, skipping...']);
end 


%% -----------------------------------------------------------------------------
%% 999. Ready
x.mutex.AddState('999_ready');
x.mutex.Unlock();
result = true;
close all;


end
