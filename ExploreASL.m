function [x] = ExploreASL(varargin)
%ExploreASL ExploreASL master script which calls the individual import, defacing & processing modules.
%
% FORMAT: [x] = ExploreASL([DatasetRoot, Import, Deface, Process, bPause, iWorker, nWorkers])
% 
% INPUT:
%
%   Definitions:
%
%   DatasetRoot - [DatasetRoot]                              (OPTIONAL, STRING, DEFAULT = [])
%   Import      - [DCM2NII, NII2BIDS]                        (OPTIONAL, BOOLEAN, DEFAULT = 0)
%   Deface      - [DEFACE]                                   (OPTIONAL, BOOLEAN, DEFAULT = 0)
%   Process     - [BIDS2LEGACY, STRUCTURAL, ASL, POPULATION] (OPTIONAL, BOOLEAN, DEFAULT = 0)
%   bPause      - [bPause]                                   (OPTIONAL, BOOLEAN, DEFAULT = 0)
%   iWorker     - [iWorker]                                  (OPTIONAL, INTEGER, DEFAULT = 1)
%   nWorkers    - [nWorkers]                                 (OPTIONAL, INTEGER, DEFAULT = 1)
%
%   Explanations:
%
%   DatasetRoot - Path to the BIDS dataset root directory (this should be the directory containing sourcedata/rawdata/derivatives)
%   Import      - Run the DICOM to BIDS import workflow (scalar or vector input for the import modules)
%   Deface      - Run the defacing (simple true/false boolean to deface the anatomical data)
%   Process     - Run the ExploreASL processing pipeline (scalar or vector input for the processing modules)
%   bPause      - Pause workflow before ExploreASL pipeline (simple true/false boolean to pause before processing)
%   iWorker     - iWorker defines in which of the parallel ExploreASL calls we are (allows parallelization when called externally)
%   nWorkers    - nWorkers defines how many ExploreASL calls are made in parallel (allows parallelization when called externally)              
%
% OUTPUT:
%
%   x              - Struct containing ExploreASL parameters
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:               This master script runs the following five modular workflows:
% 
% I.   Initialization:       ExploreASL_Initialize:       This wrapper initializes ExploreASL.
% II.  Import:               ExploreASL_Import:           Multi-step import workflow for DCM2NII and NII2BIDS.
% III. Deface:               ExploreASL_Deface:           Iterate over subjects and deface the anatomical files.
% IV.  Process:              ExploreASL_Process:          Multi-step processing workflow for the Structural, ASL and Population module.
% V.   Print user feedback:  xASL_init_PrintUserFeedback: This wrapper prints a final feedback after the processing pipeline
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLES: There are different options to run ExploreASL. The most common and easiest way is to use the scalar input shown below.
%
% Scalar input: Here we run all Import, the Deface and all Process modules. We pause before processing and define iWorker and nWorkers.
%
%      DatasetRoot = '/MyDisk/MyStudy';
%      Import = 1;
%      Deface = 1;
%      Process = 1;
%      bPause = 1;
%      iWorker = 1;
%      nWorkers = 1;
%
% Vector input: Alternatively we can also choose to run only specific sub-modules of the import or processing workflow using the vector notation.
%
%      DatasetRoot = '/MyDisk/MyStudy';
%      Import = [1 1];
%      Deface = 1;
%      Process = [1 1 1 1];
%      bPause = 1;
%      iWorker = 1;
%      nWorkers = 1;
%
% Character input: To run ExploreASL in compiled mode from a console it can be useful to have the possibility to use the character input.
%
%      DatasetRoot = '/MyDisk/MyStudy';
%      Import = '[1 1]';
%      Deface = '1';
%      Process = '[1 1 1 1]';
%      bPause = '1';
%      iWorker = '1';
%      nWorkers = '1';
%
% Command: You can use one of the input options shown above with the matlab command shown below.
%
%      [x] = ExploreASL(DatasetRoot, Import, Deface, Process, bPause, iWorker, nWorkers);
%
% For additional examples we recommend to check out the ASL-BIDS or the Basics tutorial sections:
%
% - https://exploreasl.github.io/Documentation/latest/Tutorials-ASL-BIDS/
% - https://exploreasl.github.io/Documentation/latest/Tutorials-Basics/
%
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
    % III. Deface Master
    if x.opts.bDefaceData
        x = ExploreASL_Deface(x);
    end

    % -----------------------------------------------------------------------------
    % IV. Processing Master
    if x.opts.bProcessData
        x = ExploreASL_Process(x);
    end

    % -----------------------------------------------------------------------------    
    % V. Print user feedback (after pipeline)
    xASL_init_PrintUserFeedback(x, nargout, 1);
    

end


