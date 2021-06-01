function [x] = ExploreASL_Master(varargin)
%ExploreASL_Master ExploreASL pipeline master wrapper calling the individual import & pipeline modules
%
% FORMAT: [x] = ExploreASL([DatasetRoot, ImportModules, ProcessModules, bPause, iWorker, nWorkers])
% 
% INPUT:
%   DatasetRoot    - Path to the BIDS dataset root directory (OPTIONAL, DEFAULT = EMPTY)
%
%   ImportModules  - [DCM2NII, NII2BIDS, DEFACE, BIDS2LEGACY] (OPTIONAL, BOOLEAN ARRAY)
%                  - DCM2NII = Run the DICOM to NIFTI conversion (BOOLEAN, DEFAULT = 0)
%                  - NII2BIDS = Run the NIFTI to BIDS conversion (BOOLEAN, DEFAULT = 0)
%                  - DEFACE = Run the defacing (BOOLEAN, DEFAULT = 0)
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
% ExploreASL_ProcessMaster - Multi-step processing workflow for the STRUCTURAL, ASL and POPULATION module.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLES:
%
% [x] = ExploreASL_Master('/MyDisk/MyStudy', 1, 1, 1, 1, 1);
% [x] = ExploreASL_Master('/MyDisk/MyStudy', [1 1 1 1], [1 1 1], 1, 1, 1);
% [x] = ExploreASL_Master('/MyDisk/MyStudy', '[1 1 1 1]', '[1 1 1]', '1', '1', '1');
%
% For additional examples we recommend to check out the tutorial section: https://exploreasl.github.io/Documentation/site/Tutorials/
% __________________________________
% Copyright 2015-2021 ExploreASL

    % -----------------------------------------------------------------------------
    %% Initialization when calling this function

    % NB: *.mat files that contain statistics in data-root folder
    % (e.g.\analysis) should contain the following format:

    % 1st column should be subject-id column
    % 2nd column should be parameter values
    % or 2nd column is session-id (e.g. 'ASL_1'), then
    % 3rd column contains parameter values
    
    % -----------------------------------------------------------------------------
    %% Initialization
    x = ExploreASL_Initialize(varargin{:});
    
    % -----------------------------------------------------------------------------
    %% Import Master
    if x.opts.bImportData
        x = ExploreASL_ImportMaster(x);
    end
    % Store logging information about errors/warnings in backup variable
    if isfield(x,'logging')
    	loggingBackUp = x.logging;
    end
    
    % -----------------------------------------------------------------------------
    %% Re-Initialize for potential data loading/processing
    if x.opts.bReinitialize
        x = ExploreASL_Initialize(x.opts.DatasetRoot, x.opts.ImportModules, x.opts.ProcessModules, x.opts.bPause, x.opts.iWorker, x.opts.nWorkers);
    end
    % Retrieve logging information about errors/warnings from backup variable
    if exist('loggingBackUp', 'var')
    	x.logging = loggingBackUp;
    end
    
    % -----------------------------------------------------------------------------
    %% Print user feedback
    if ~x.opts.bProcessData || x.opts.bOnlyLoad
        if x.opts.bOnlyLoad && nargout==0
            warning('Data loading requested but no output structure defined');
            fprintf('%s\n', 'Try adding "x = " to the command to load data into the x structure');
        end
        return; % skip processing
    elseif ~isdeployed && x.opts.bPause % if this is true, we skip the break here
        fprintf('%s\n','Press any key to start processing & analyzing');
        fprintf('Please ensure you have a read-only copy of your original data as they may be overwritten\n');
        fprintf('%s\n','Or press CTRL/command-C to cancel...  ');
        pause;
    end

    % -----------------------------------------------------------------------------
    %% Processing Master
    [x] = ExploreASL_ProcessMaster(x);

    % -----------------------------------------------------------------------------    
    %% Finishing touch
    fprintf('Many thanks for using <a href="https://github.com/ExploreASL" rel="nofollow">ExploreASL</a>, ');
    fprintf('please don''t forget to cite <a href="https://pubmed.ncbi.nlm.nih.gov/32526385/" rel="nofollow">https://pubmed.ncbi.nlm.nih.gov/32526385/</a>.\n');
    fprintf('Note that ExploreASL is a collaborative effort.\n');
    fprintf('Therefore, please don''t hesitate to contribute by feedback, adding code snippets, or clinical experience!\n');
    

end


