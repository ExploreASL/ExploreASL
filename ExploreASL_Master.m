function [x] = ExploreASL_Master(varargin)
%ExploreASL_Master ExploreASL pipeline master wrapper calling the individual import & pipeline modules
%
% FORMAT: [x] = ExploreASL([DataParPath, ImportModules, ProcessModules, bPause, iWorker, nWorkers])
% 
% INPUT:
%   DataParPath    - Path to data parameter file (OPTIONAL, DEFAULT = prompting user input)
%
%   ImportModules  - [DCM2NII, NII2BIDS, ANONYMIZE, BIDS2LEGACY] (OPTIONAL, BOOLEAN ARRAY)
%                  - DCM2NII = Run the DICOM to NIFTI conversion (BOOLEAN, DEFAULT = 0)
%                  - NII2BIDS = Run the NIFTI to BIDS conversion (BOOLEAN, DEFAULT = 0)
%                  - ANONYMIZE = Run the defacing and full anonymization (BOOLEAN, DEFAULT = 0)
%                  - BIDS2LEGACY = Run the BIDS to LEGACY conversion (BOOLEAN, DEFAULT = 0)
%
%   ProcessModules - [STRUCTURAL, ASL, POPULATION] (OPTIONAL, BOOLEAN ARRAY) 
%                  - STRUCTURAL = Run the Structural Module (BOOLEAN, DEFAULT = 0)
%                  - ASL = Run the ASL Module (BOOLEAN, DEFAULT = 0)
%                  - POPULATION = Run the Population Module (BOOLEAN, DEFAULT = 0)
%
%   bPause         - TRUE = Pause workflow before ExploreASL pipeline (OPTIONAL, DEFAULT = FALSE)
%
%   iWorker        - Allows parallelization when called externally. 
%                    iWorker defines which of the parallel ExploreASL calls we are (OPTIONAL, DEFAULT=1)
%
%   nWorkers       - Allows parallelization when called externally. 
%                    nWorkers defines how many ExploreASL calls are made in parallel (OPTIONAL, DEFAULT=1)
%
% OUTPUT:
%   x              - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    This masterscript starts ExploreASL by first calling ExploreASL_Initialize, 
%                 then running xASL_Module_Structure, xASL_module_ASL and xASL_module_Population.
%                 When ProcessModules is set to 0, the ExploreASL pipeline is not started but only  
%                 initialized for debugging. This pipeline can be run from CLI or using the python GUI.
% 
% ExploreASL_Initialize    - This wrapper initializes ExploreASL: managing paths, deployment, etc.
% ExploreASL_ImportMaster  - Multi-step import workflow for DCM2NII, NII2BIDS & BIDS2LEGACY.
% xASL_Module_Structure    - Processes structural data, i.e. high-resolution T1w and FLAIR scans, 
%                            type help xASL_module_Structural for more information.
% xASL_Module_ASL          - Processes ASL data, i.e. ASL scans and M0, type help xASL_module_ASL for more information. 
% xASL_Module_Population   - Processes data on population basis, mostly QC, but also template/parametric maps creation, 
%                            etc. Type help xASL_module_Population for more information. 
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLES:
%
% [x] = ExploreASL_Master('/MyDisk/MyStudy/DataPar.json', 1, 1, 1, 1, 1);
% [x] = ExploreASL_Master('/MyDisk/MyStudy/DataPar.json', [1 1 1 1], [1 1 1], 1, 1, 1);
% [x] = ExploreASL_Master('/MyDisk/MyStudy/DataPar.json', '[1 1 1 1]', '[1 1 1]', '1', '1', '1');
%
% For additional examples we recommend to check out the tutorial section: https://exploreasl.github.io/Documentation/site/Tutorials/
% __________________________________
% Copyright 2015-2021 ExploreASL

    % -----------------------------------------------------------------------------
    %% 1 Initialization when calling this function

    % NB: *.mat files that contain statistics in data-root folder
    % (e.g.\analysis) should contain the following format:

    % 1st column should be subject-id column
    % 2nd column should be parameter values
    % or 2nd column is session-id (e.g. 'ASL_1'), then
    % 3rd column contains parameter values
    
    % -----------------------------------------------------------------------------
    % Initialization
    x = ExploreASL_Initialize(varargin{:});
    
    % -----------------------------------------------------------------------------
    % Import Master
    x = ExploreASL_ImportMaster(x);
    
    % -----------------------------------------------------------------------------
    % Re-Initialize for potential data loading/processing
    if x.bReinitialize > 0
        x = ExploreASL_Initialize(x.DataParPath, x.ImportModules, x.ProcessModules, x.bPause, x.iWorker, x.nWorkers);
    end
    
    % -----------------------------------------------------------------------------
    % Print user feedback
    if x.bProcessData==0 || x.bProcessData==2
        if x.bProcessData==2 && nargout==0
            warning('Data loading requested but no output structure defined');
            fprintf('%s\n', 'Try adding "x = " to the command to load data into the x structure');
        end
        return; % skip processing
    elseif ~isdeployed && x.bPause % if this is true, we skip the break here
        fprintf('%s\n','Press any key to start processing & analyzing');
        fprintf('Please ensure you have a read-only copy of your original data as they may be overwritten\n');
        fprintf('%s\n','Or press CTRL/command-C to cancel...  ');
        pause;
    end

    
    %% -----------------------------------------------------------------------------
    %% 1  xASL_module_Structural
    %  1) Alignment T1w -> ACPC in MNI
    %  2) Coregister FLAIR, T2, T1c -> T1w
    %  3) FLAIR biasfield correction
    %  4) Segment FLAIR
    %  5) Lesion Filling
    %  6) Segment T1w
    %  7) CleanUpWMH_SEGM
    %  8) Resample structural images -> common space
    %  9) Get tissue volumes (incl WMH if available)
    % 10) Visual & automatic QC

    if x.ProcessModules(1)==1
        [~, x] = xASL_Iteration(x,'xASL_module_Structural');
        % The following DARTEL module is an optional extension of the structural module
        % to create population-specific templates
        if isfield(x,'SegmentSPM12') && x.SegmentSPM12 && x.nSubjects>1
            % in case we used SPM12 instead of CAT12 for segmentation,
            % we have to run DARTEL separately
            [~, x] = xASL_Iteration(x,'xASL_module_DARTEL');
        end
        % Now only check the availability of files when not running parallel
        if x.nWorkers==1; xASL_adm_CreateFileReport(x); end        
    end

    % Optional modules
    % The following are optional extensions of the structural module and not required to run, normally they can be ignored:
    if isfield(x, 'bRunModule_LongReg') && x.bRunModule_LongReg
        [~, x] = xASL_Iteration(x,'xASL_module_LongReg'); % use this module for longitudinal registration
    end
    if isfield(x, 'bRunModule_DARTEL') && x.bRunModule_DARTEL
        [~, x] = xASL_Iteration(x,'xASL_module_DARTEL'); % use this module for additional additional -subject registration/creating templates
    end
    
    % -----------------------------------------------------------------------------
    %% 2    xASL_module_ASL  
    %  1    TopUp geometric distortion correction
    %  2    Motion correction
    %  3    Registration ASL sessions to T1w
    %  4    Reslice ASL data to high resolution standard space
    %  5    Resolution estimation & prepare partial volume maps
    %  6    Process M0
    %  7    Quantification
    %  8    Create analysis mask
    %  9    Visual check
    % 10    WAD-QC
    
    if x.ProcessModules(2)==1
        [~, x] = xASL_Iteration(x,'xASL_module_ASL');
        % Now only check the availability of files when not running parallel
        if x.nWorkers==1; xASL_adm_CreateFileReport(x); end    
    end


    % -----------------------------------------------------------------------------
    %% 3    xASL_module_Population
    % Performs all group-level processing & QC

    if x.ProcessModules(3)==1
        [~, x] = xASL_Iteration(x,'xASL_module_Population');
    end

    % -----------------------------------------------------------------------------    
    %% Finishing touch
    fprintf('Many thanks for using ExploreASL, please don''t forget to cite https://pubmed.ncbi.nlm.nih.gov/32526385/\n');
    fprintf('Note that ExploreASL is a collaborative effort.\n');
    fprintf('Therefore, please don''t hesitate to contribute by feedback, adding code snippets, or clinical experience!\n');
    

end


