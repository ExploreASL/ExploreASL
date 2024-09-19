function [x] = ExploreASL_Process(x)
%ExploreASL_Process Multi-step processing workflow for the STRUCTURAL, ASL and POPULATION module.
%
% FORMAT: [x] = ExploreASL_Process(x)
%
% INPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x             - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Multi-step processing workflow for the STRUCTURAL, ASL and POPULATION module.
%
% xASL_Module_Structure    - Processes structural data, i.e. high-resolution T1w and FLAIR scans, 
%                            type help xASL_module_Structural for more information.
% xASL_Module_ASL          - Processes ASL data, i.e. ASL scans and M0, type help xASL_module_ASL for more information. 
% xASL_Module_Population   - Processes data on population basis, mostly QC, but also template/parametric maps creation, 
%                            etc. Type help xASL_module_Population for more information.
%
%
% DETAILS:
%
% ===== xASL_module_Structural =====
%
%  1. Alignment T1w -> ACPC in MNI
%  2. Coregister FLAIR, T2, T1c -> T1w
%  3. FLAIR biasfield correction
%  4. Segment FLAIR
%  5. Lesion Filling
%  6. Segment T1w
%  7. CleanUpWMH_SEGM
%  8. Resample structural images -> common space
%  9. Get tissue volumes (incl WMH if available)
% 10. Visual & automatic QC
%
%
% ===== xASL_module_ASL =====
%
%  1. TopUp geometric distortion correction
%  2. Motion correction
%  3. Registration ASL sessions to T1w
%  4. Reslice ASL data to high resolution standard space
%  5. Resolution estimation & prepare partial volume maps
%  6. Process M0
%  7. Quantification
%  8. Create analysis mask
%  9. Visual check
% 10. WAD-QC
%
%
% ===== xASL_module_Population =====
%
%  1. CreatePopulationTemplates
%  2. CreateAnalysisMask
%  3. CreateBiasfield
%  4. GetDICOMStatistics
%  5. GetVolumeStatistics
%  6. GetMotionStatistics
%  6. GetRegistrationStatistics
%  7. GetROIstatistics
%  8. SortBySpatialCoV
%  9. DeleteTempFiles
% 10. GzipAllFiles
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright (c) 2015-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________



    %% 0. Workflow for initialization of data loading and processing
    
    x = xASL_init_DataLoading(x); % this initializes all generic data loading parameters
	x = xASL_init_LoadDataPar(x); % Load/create dataPar.json & its settings
    x = xASL_init_SubjectList(x); % create subject list for loading data from rawdata (BIDS2Legacy) or from derivatives (legacy)
    x = xASL_init_Parallelization(x); % choose which subjects this worker processes

    if x.opts.nWorkers==1
        % Remove lock dirs from previous runs that might still exist if the pipeline crashed.
        % This is only performed if ExploreASL is not running in parallel
        % Note that any lock dirs for individual modules/su
        x = xASL_init_RemoveLockDirs(x);
    end

    xASL_init_PrintUserFeedback(x, 1, 0);

    if x.opts.bReadRawdata
        % 0 Run BIDS to Legacy conversion.
        % The rawdata (datasetRoot/rawdata/*) is in BIDS, where the ExploreASL derivatives 
        % (datasetRoot/derivatives/ExploreASL/*) are in ExploreASL's own legacy format.
        % In the future, we will probably move to BIDS derivatives, but for now we keep the legacy format.
        % Therefore, here we copy all BIDS data to the legacy format, and then we run the processing.
        xASL_adm_BreakString('BIDS to ExploreASL LEGACY CONVERSION'); % Print feedback
        [~, x] = xASL_init_Iteration(x, 'xASL_module_BIDS2Legacy');
        x = xASL_init_CreateParticipantsTSV(x);
    end

    % Here, we load all the ExploreASL derivatives data into the Matlab x structure, 
    % such that we can use it for processing.
    % i.e., datasetRoot/derivatives/ExploreASL/*
	
    x = xASL_init_Session_TimePoint_Lists(x);
    % generate session & visit/time point lists
    % this is based on the legacy structure
    x = xASL_init_Process(x); % this initializes all generic processing settings

    
    % -----------------------------------------------------------------------------
    %% 1  xASL_module_Structural
    if x.opts.bProcess(1)==1
        [~, x] = xASL_init_Iteration(x,'xASL_module_Structural');
        % The following DARTEL module is an optional extension of the structural module
        % to create population-specific templates
        if isfield(x.modules.structural,'bSegmentSPM12') && x.modules.structural.bSegmentSPM12 && x.dataset.nSubjects>1
            % in case we used SPM12 instead of CAT12 for segmentation, we have to run DARTEL separately
            [~, x] = xASL_init_Iteration(x,'xASL_module_DARTEL');
        end
        % Now only check the availability of files when not running parallel
        if x.opts.nWorkers==1; xASL_adm_CreateFileReport(x); end        
    end

    %% Optional modules
    % The following are optional extensions of the structural module and not required to run, normally they can be ignored:
    if isfield(x.modules, 'bRunLongReg') && x.modules.bRunLongReg
        % Use this module for longitudinal registration
        [~, x] = xASL_init_Iteration(x,'xASL_module_LongReg');
    end
    if isfield(x.modules, 'bRunDARTEL') && x.modules.bRunDARTEL
        % Use this module for additional additional -subject registration/creating templates
        [~, x] = xASL_init_Iteration(x,'xASL_module_DARTEL');
    end
    
    % -----------------------------------------------------------------------------
    %% 2    xASL_module_ASL  
    if x.opts.bProcess(2)==1
        [~, x] = xASL_init_Iteration(x,'xASL_module_ASL');
        % Now only check the availability of files when not running parallel
        if x.opts.nWorkers==1; xASL_adm_CreateFileReport(x); end    
    end


    % -----------------------------------------------------------------------------
    %% 3    xASL_module_Population
    % Performs all group-level processing & QC
    if x.opts.bProcess(3)==1
        [~, x] = xASL_init_Iteration(x,'xASL_module_Population');
    end
    
    % -----------------------------------------------------------------------------
    
    %%      Zip derivatives
    % xASL_module_Population also zips all, so in that case the zipping
    % here would be skipped

    % Input check
    if x.opts.nWorkers>1 % don't run population module when ExploreASL is parallelized
        fprintf('%s\n', 'Not zipping NIfTIs because running ExploreASL in parallel mode');
	else
		if isfield(x, 'dir') && isfield(x.dir,'xASLDerivatives') && ~isempty(x.dir.xASLDerivatives)
			xASL_adm_GzipAllFiles(x.dir.xASLDerivatives, [], [], fullfile(x.opts.MyPath, 'External'));
		end
    end

end


%% =======================================================================================================================
%% =======================================================================================================================
function [x] = xASL_init_RemoveLockDirs(x)
    %xASL_init_RemoveLockDirs Remove 'lock-dir' if present from aborted previous run, for current subjects only
    %
    % FORMAT: [x] = xASL_init_RemoveLockDirs(x)
    % 
    % INPUT:
    %   x        - ExploreASL x struct (STRUCT, REQUIRED)
    %
    % OUTPUT:
    %   x        - ExploreASL x struct 
    %                         
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % DESCRIPTION:    Remove 'lock-dir' if present from aborted previous run, for current subjects only.
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % EXAMPLE:        [x] = xASL_init_RemoveLockDirs(x);
    % __________________________________
    
    
    %% LockDir within 2 directories (e.g. T1, FLAIR or ASL)
    LockDir = fullfile(x.dir.xASLDerivatives, 'lock');

    if exist(LockDir, 'dir')

        listLockedDirs = xASL_adm_FindByRegExp(fullfile(x.dir.xASLDerivatives, 'lock'), {'(ASL|Structural|LongReg_T1|BIDS2Legacy)', x.dataset.subjectRegexp, '.*module.*','^(locked)$'}, 'Match', 'Directories');
        if ~isempty(listLockedDirs)
            fprintf('\n');
            warning('Locked folders were found, consider removing them before proceeding:');
            for iDir=1:length(listLockedDirs)
                fprintf('%s\n', listLockedDirs{iDir});
            end
            fprintf('\n');
        end

        % LockDir within 2 directories (e.g. DARTEL)
        listLockedDirs = xASL_adm_FindByRegExp(fullfile(x.dir.xASLDerivatives, 'lock'), {'(Population|DARTEL_T1)', '.*module.*','^(locked)$'}, 'Match','Directories');
        if ~isempty(listLockedDirs)
            fprintf('\n');
            warning('Locked folders were found, consider removing them before proceeding:');
            for iDir=1:length(listLockedDirs)
                fprintf('%s\n', listLockedDirs{iDir});
            end
            fprintf('\n');        
        end
    end
    
end