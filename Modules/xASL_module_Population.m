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
% - `080_SortBySpatialCoV`          - Sort ASL_Check QC images by their spatial CoV in quality bins
% - `090_DeleteTempFiles`           - Delete temporary files
% - `100_GZipAllFiles`              - Zip files to reduce disc space usage of temporary and non-temporay NIfTI files
%
% EXAMPLE: [~, x] = xASL_module_Population(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


%% ------------------------------------------------------------------------------------------------------------
%% Admin
[x] = xASL_init_SubStructs(x);

% Input check
if x.opts.nWorkers>1 % don't run population module when ExploreASL is parallelized
    warning('Population module should not run in parallel, skipping...');
    fprintf('%s\n', 'Best to run this module after you have ran ExploreASL in parallel, by restarting ExploreASL non-parallel');
    result = true;
    return;
end

x = xASL_wrp_Population_PrepareAtlas4ROI(x); % Parse x.S.Atlases & x.S.TissueMasking

% Default datatypes
if ~isfield(x.S,'DataTypes') || isempty(x.S.DataTypes)
	x.S.DataTypes = {'qCBF'}; % Default
	% Alternatives: 'Tex' 'ATT' 'SD' 'M0' 'ABV'
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
%% 1    Create template images
if ~x.mutex.HasState(StateName{1})
    xASL_wrp_CreatePopulationTemplates(x);  % this doesn't work nicely yet with sessions, should be changed after new BIDS is implemented

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
%% 2    Create population-based analysis mask for ROI-based analysis & VBA
if ~x.mutex.HasState(StateName{2}) && bHasASL
    x = xASL_im_CreateGroupAnalysisMask(x);
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
        
        if isempty(regexp(x.S.InputDataStr(1), '(q|r)', 'once'))
            % Some DataTypes have a prefix in standard space but not in
            % native space: e.g.,
            % qCBF in standard space is CBF in native space
            % rc1T1 in standard space is c1T1 in native space
            x.S.InputDataStrNative = x.S.InputDataStr;
        else
            x.S.InputDataStrNative = x.S.InputDataStr(2:end);
        end
    
        % Iterate over atlases
        x.dir.dirAtlas = fullfile(x.opts.MyPath, 'external', 'Atlases');

        for iAtlas=1:length(x.S.Atlases)
            % Note that the number of ROI atlases here should be the same the number of tissue masking chosen
            % If needed, an atlas or tissue type can be provided multiple times in different combinations

            % We use the specified tissue type
            % 'GM' = gray matter
            % 'WM' = white matter
            % 'WB' = whole brain (= GM+WM)
            x.S.TissueMaskingLocal = x.S.TissueMasking(iAtlas);
            
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
                % ROI statistics (optional: native space)
                if x.modules.population.bNativeSpaceAnalysis
                    x.S.InputNativeSpace = 1;
                    x.S.InputAtlasNativeName = [x.S.Atlases{iAtlas} '_Atlas'];
                    xASL_wrp_GetROIstatistics(x);
                end
            end
        end
    
        % Check if we should do the same for Lesion or ROI masks (i.e. individual "atlases")
		% PM: not yet developed/tested in native space
        
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

        x.S.InputNativeSpace = 0;
		x.S.bSubjectSpecificROI = true;
        for iROI = 1:length(LesionUniqueROIList)
            x.S.InputAtlasPath = fullfile(x.D.PopDir, LesionUniqueROIList{iROI});
            xASL_wrp_GetROIstatistics(x);
        end
    end

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
    if ~x.settings.bReproTesting && x.settings.DELETETEMP
        xASL_adm_DeleteManyTempFiles(x);
    end
    x.mutex.AddState(StateName{10});
    fprintf('%s\n',[StateName{10} ' was performed']);
else
    fprintf('%s\n',[StateName{10} ' has already been performed, skipping...']);
end

%% 10    xASL_adm_GzipAllFiles
if ~x.mutex.HasState(StateName{11})
    xASL_adm_GzipAllFiles(x.dir.xASLDerivatives,[],[],fullfile(x.opts.MyPath,'External'));
    x.mutex.AddState(StateName{11});
    fprintf('%s\n',[StateName{11} ' was performed']);
else
        fprintf('%s\n',[StateName{11} ' has already been performed, skipping...']);
end
    

%% -----------------------------------------------------------------------------
%% 999 Ready
x.mutex.AddState('999_ready');
x.mutex.Unlock();
result = true;
close all;


end




%% -----------------------------------------------------------------------------
%% -----------------------------------------------------------------------------
function [x] = xASL_wrp_Population_PrepareAtlas4ROI(x)
%xASL_wrp_Population_PrepareAtlas4ROI Parse x.S.Atlases & x.S.TissueMasking

bPrintInstructions = false; % Suboptimal state, print instructions
bAtlasTissueMatch = true; % Necessary condition

if ~isfield(x.S,'Atlases') && ~isfield(x.S, 'TissueMasking')
	% Default atlases/ROIs & tissue masks if nothing is provided
	x.S.Atlases = {'TotalGM','DeepWM'}; % Default
    x.S.TissueMasking = {'GM' 'WM'}; % GM WM, fits with the TotalGM & DeepWM above
    % Note that this should be in the same order as the atlases/ROIs
    % A mismatch (e.g. TissueMasking=GM for ROI=deepWM) would result in an empty ROI, producing a NaN in the .tsv table
elseif ~isfield(x.S, 'Atlases') && isfield(x.S, 'TissueMasking')
	% Missing Atlases, but provided TissueMasking - cannot continue
    warning('Custom tissue-types (x.S.TissueMasking) specified without ROI atlas-selection (x.S.Atlases). Atlases need to be provided. See instructions below:');
    bAtlasTissueMatch = false;
	bPrintInstructions = true;
elseif isfield(x.S, 'Atlases') && isfield(x.S, 'TissueMasking') && length(x.S.Atlases)~=length(x.S.TissueMasking)
	% Non matching lengths, cannot continue
    warning('The same number of ROI atlases as subject-wise tissue-types provided does not match:');
    fprintf('%s\n', ['x.S.Atlases: ' strjoin(x.S.Atlases)]);
    fprintf('%s\n', ['x.S.TissueMasking: ' strjoin(x.S.TissueMasking)]);
    bAtlasTissueMatch = false;
	bPrintInstructions = true;
elseif 	isfield(x.S, 'Atlases') && ~isfield(x.S, 'TissueMasking')
	% TissueMasking not provided, so it has to be extracted from Atlases as previously.
    warning('ROIs provided in x.S.Atlases without the tissue-types for these ROIs in x.S.TissueMasking. This will be fixed automatically, provide x.S.TissueMarking the next time');
    bAtlasTissueMatch = true;
	bPrintInstructions = true;
	% Fill in TissueMasking based on the atlas names and default to GM. The user will see a warning and automatic tissue masks to verify
	for iAtlas = 1:numel(x.S.Atlases)
		if ~isempty(regexpi(x.S.Atlases{iAtlas}, 'WM')) || ~isempty(regexpi(x.S.Atlases{iAtlas}, 'whitematter'))
			x.S.TissueMasking{iAtlas} = 'WM';
		elseif ~isempty(regexpi(x.S.Atlases{iAtlas}, 'WB')) || ~isempty(regexpi(x.S.Atlases{iAtlas}, 'wholebrain'))
			x.S.TissueMasking{iAtlas} = 'GM';
		else
			x.S.TissueMasking{iAtlas} = 'GM';
		end
	end
end
if bPrintInstructions
    fprintf('%s\n', 'When ROI atlases are provided in x.S.Atlases, their tissue types');
    fprintf('%s\n', 'need to be provided as well in dataPar.json, with either option ''GM'', ''WM'', ''WB'' (GM+WM).');
    fprintf('%s\n', 'E.g., when not provided, these default to:');
    fprintf('%s\n', 'x.S.Atlases = [''TotalGM'' ''DeepWM'']');
    fprintf('%s\n\n', 'x.S.TissueMasking = [''GM'' ''WM'']');
end

if ~bAtlasTissueMatch
	% No match means that we have to end it
    error('Not the same number of ROI atlases as subject-wise tissue-types, skipping');
end

fprintf('\nThe following ROI atlases have been selected with the following tissue:\n')
for iAtlas=1:length(x.S.Atlases)
    fprintf('%s\n', [x.S.TissueMasking{iAtlas} ' tissue within ' x.S.Atlases{iAtlas} ' ROIs']);
end
fprintf('\n');

end