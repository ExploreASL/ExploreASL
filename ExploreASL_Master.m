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


    %% ExploreASL_Master Workflow
    
    % -----------------------------------------------------------------------------
    % Initialization
    x = ExploreASL_Initialize(varargin{:});
    
    % -----------------------------------------------------------------------------
    % Import Master
    x = ExploreASL_ImportMaster(x);
    
    % -----------------------------------------------------------------------------
    % Data loading
    x = xASL_init_DataLoading(x);
    
    % -----------------------------------------------------------------------------
    % Print user feedback
    xASL_init_PrintUserFeedback(x);

    % -----------------------------------------------------------------------------
    % Processing Master
    x = ExploreASL_ProcessMaster(x);

    % -----------------------------------------------------------------------------    
    % Finishing touch
    fprintf('Many thanks for using <a href="https://github.com/ExploreASL" rel="nofollow">ExploreASL</a>, ');
    fprintf('please don''t forget to cite <a href="https://pubmed.ncbi.nlm.nih.gov/32526385/" rel="nofollow">https://pubmed.ncbi.nlm.nih.gov/32526385/</a>.\n');
    fprintf('Note that ExploreASL is a collaborative effort.\n');
    fprintf('Therefore, please don''t hesitate to contribute by feedback, adding code snippets, or clinical experience!\n');
    

end


