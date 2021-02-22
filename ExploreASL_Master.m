function [x] = ExploreASL_Master(DataParPath, ProcessData, SkipPause, iWorker, nWorkers, iModules)
%ExploreASL_Master ExploreASL pipeline master wrapper calling the individual pipeline modules
%
% FORMAT: [x] = ExploreASL([DataParPath, ProcessData, SkipPause, iWorker, nWorkers, iModules])
% 
% INPUT:
%   DataParPath - path to data parameter file (OPTIONAL, required when ProcessData=true, will then be prompted if omitted)
%   ProcessData - 0 = only initialize ExploreASL functionality (e.g. path management etc, REQUIRED, will be prompted if omitted)
%               - 1 = initialize ExploreASL functionality, load data & start processing pipeline, 
%               - 2 = initialize ExploreASL functionality, load data but no processing
%               - (OPTIONAL, default=prompt the user)
%   SkipPause   - true if calling pipeline without pausing for questioning for user prompting (OPTIONAL, DEFAULT=false)
%   iWorker     - allows parallelization when called externally. iWorker defines which of the parallel ExploreASL calls we are (OPTIONAL, DEFAULT=1)
%   nWorkers    - allows parallelization when called externally. nWorkers defines how many ExploreASL calls are made in parallel (OPTIONAL, DEFAULT=1)
%   iModules    - scalar or vector allowing to select modules to run, e.g.
%                  [1] = only structural, [2] = only ASL, [3] = only
%                  Population. [1 2] = Structural & ASL module, etc.
%                  OPTIONAL, DEFAULT = [1 2 3]
% OUTPUT:
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This masterscript starts ExploreASL by first calling ExploreASL_Initialize, 
% then running xASL_Module_Structure, xASL_module_ASL and xASL_module_Population
% When ProcessData is set to false (either as argument or when prompted),
% the ExploreASL pipeline is not started but only initialized for debugging.
% This pipeline can be run from CLI or as GUI (later to be implemented in the SPM GUI)
%
% xASL_Module_Structure  - processes structural data, i.e. high-resolution
%                          T1w and FLAIR scans, type help xASL_module_Structural for more information
% xASL_Module_ASL        - processes ASL data, i.e. ASL scans and M0, type help xASL_module_ASL for more information 
% xASL_Module_Population - processes data on population basis, mostly QC, but also template/parametric maps creation, etc. Type help xASL_module_Population for more information 
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE for GUI: ExploreASL
% EXAMPLE for calling externally: ExploreASL('//MyDisk/MyStudy/DataPar.m', true, true);
% EXAMPLE for calling externally to run the Structural module only: ExploreASL('//MyDisk/MyStudy/DataPar.m', true, true, [], [], 1);
% EXAMPLE for calling externally to run the ASL & Population modules: ExploreASL('//MyDisk/MyStudy/DataPar.m', true, true, [], [], [2 3]);
% EXAMPLE for debugging/initialization only: [x] = ExploreASL('',0);
% __________________________________
% Copyright 2015-2021 ExploreASL

    % -----------------------------------------------------------------------------
    %% 1 Initialization when calling this function
    % -----------------------------------------------------------------------------

    % NB: *.mat files that contain statistics in data-root folder
    % (e.g.\analysis) should contain the following format:

    % 1st column should be subject-id column
    % 2nd column should be parameter values
    % or 2nd column is session-id (e.g. 'ASL_1', then
    % 3rd column contains parameter values
    
    % using exist(var) here as nargin doesnt work when debugging
    DataParPath = char(DataParPath); % convert to char on default
    if exist('ProcessData', 'var') && ~isempty(ProcessData)
        if ischar(ProcessData)
            ProcessData = str2num(ProcessData);
        end
    else
        ProcessData = []; % by default we let the user choose
    end
    if exist('SkipPause', 'var') && ~isempty(SkipPause)
        if ischar(SkipPause)
            SkipPause= str2num(SkipPause);
        end
    else
        SkipPause = false; % by default we don't skip the pause question below
    end
    if exist('iWorker', 'var') && ~isempty(iWorker)
        if ischar(iWorker)
            iWorker= str2num(iWorker);
        end
    else
        iWorker = 1; % by default we don't parallelize ExploreASL instances
    end
    if exist('nWorkers', 'var') && ~isempty(nWorkers)
        if ischar(nWorkers)
            nWorkers = str2num(nWorkers);
        end
    else
        nWorkers = 1; % by default we don't parallelize ExploreASL instances
    end
    if exist('iModules', 'var') && ~isempty(iModules)
        if ischar(iModules)
            iModules = str2num(iModules);
        end
    else
        iModules = [1 2 3]; % by default we run all default modules (1) Structural, 2) ASL, 3) Population)
    end
    
    
    x = ExploreASL_Initialize(DataParPath, ProcessData, iWorker, nWorkers);
    
    if x.ProcessData==0 || x.ProcessData==2
        if x.ProcessData==2 && nargout==0
            warning('Data loading requested but no output structure defined');
            fprintf('%s\n', 'Try adding "x = " to the command to load data into the x structure');
        end
        return; % skip processing
    elseif ~isdeployed && ~SkipPause % if this exists, we skip the break here
        fprintf('%s\n','Press any key to start processing & analyzing');
        fprintf('Please ensure you have a read-only copy of your original data as they may be overwritten\n');
        fprintf('%s\n','Or press CTRL/command-C to cancel...  ');
        pause;
    else % continue
    end

    
    %% -----------------------------------------------------------------------------
    %% 1  xASL_module_Structural
    %  1  Alignment T1w -> ACPC in MNI
    %  2) Coregister FLAIR -> T1w
    %  3) FLAIR biasfield correction
    %  4) Segment FLAIR
    %  5) Lesion Filling
    %  6) Segment T1w
    %  7) CleanUpWMH_SEGM
    %  8) Resample structural images -> common space
    %  9) Get tissue volumes (incl WMH if available)
    % 10) Visual & automatic QC

    if max(iModules==1)
        [~, x] = xASL_Iteration(x,'xASL_module_Structural');
        % The following DARTEL module is an optional extension of the structural module
        % to create population-specific templates
        if isfield(x,'SegmentSPM12') && x.SegmentSPM12 && x.nSubjects>1
            % in case we used SPM12 instead of CAT12 for segmentation,
            % we have to run DARTEL separately
            [~, x] = xASL_Iteration(x,'xASL_module_DARTEL');
        end
        % Now only check the availability of files when not running parallel
        if nWorkers==1; xASL_adm_CreateFileReport(x); end        
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
    
    if max(iModules==2)
        [~, x] = xASL_Iteration(x,'xASL_module_ASL');
        % Now only check the availability of files when not running parallel
        if nWorkers==1; xASL_adm_CreateFileReport(x); end    
    end


    % -----------------------------------------------------------------------------
    %% 3    xASL_module_Population
    % Performs all group-level processing & QC

    if max(iModules==3)
        [~, x] = xASL_Iteration(x,'xASL_module_Population');
    end

    % -----------------------------------------------------------------------------    
    %% Finishing touch
    fprintf('Many thanks for using ExploreASL, please don''t forget to cite https://pubmed.ncbi.nlm.nih.gov/32526385/\n');
    fprintf('Note that ExploreASL is a collaborative effort.\n');
    fprintf('Therefore, please don''t hesitate to contribute by feedback, adding code snippets, or clinical experience!\n');
    

end



