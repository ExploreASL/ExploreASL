function [x] = ExploreASL(varargin)
%ExploreASL ExploreASL pipeline master wrapper calling the individual import & pipeline modules
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
% DESCRIPTION:               This master script runs the following six modular workflows:
% 
% I. Initialization:         ExploreASL_Initialize
%                            This wrapper initializes ExploreASL: managing paths, deployment, etc.
% II. Import Master:         ExploreASL_Import
%                            Multi-step import workflow for DCM2NII, NII2BIDS & BIDS2LEGACY.
% III. Processing Master:      ExploreASL_Process
%                            Multi-step processing workflow for the STRUCTURAL, ASL and POPULATION module.
% IV. Print user feedback:   xASL_init_PrintUserFeedback
%                            This wrapper prints a final feedback after the processing pipeline
%
% This pipeline can be run from CLI or using the python GUI.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLES:
%
% [x] = ExploreASL('/MyDisk/MyStudy', 1, 1, 1, 1, 1);
% [x] = ExploreASL('/MyDisk/MyStudy', [1 1 1], [1 1 1 1], 1, 1, 1);
% [x] = ExploreASL('/MyDisk/MyStudy', '[1 1 1]', '[1 1 1 1]', '1', '1', '1');
%
% For additional examples we recommend to check out the tutorial section: 
% https://exploreasl.github.io/Documentation/latest/Tutorials-ASL-BIDS/
% __________________________________
% Copyright (c) 2015-2022 ExploreASL


    %% ExploreASL Workflow
    
    % -----------------------------------------------------------------------------
    % I. Initialization
    x = ExploreASL_Initialize(varargin{:});
    
    % -----------------------------------------------------------------------------
    % II. Import Master
    if x.opts.bImportData
        x = ExploreASL_Import(x);
    end

    % -----------------------------------------------------------------------------
    % III. Processing Master
    if x.opts.bProcessData
        x = ExploreASL_Process(x);
    end

    % -----------------------------------------------------------------------------    
    % IV. Print user feedback (after pipeline)
    xASL_init_PrintUserFeedback(x, nargout, 1);
    

end


