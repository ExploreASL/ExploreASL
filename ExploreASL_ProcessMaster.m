function [x] = ExploreASL_ProcessMaster(x)
%ExploreASL_ProcessMaster Multi-step processing workflow for the STRUCTURAL, ASL and POPULATION module.
%
% FORMAT: [x] = ExploreASL_ProcessMaster(x)
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
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright 2015-2021 ExploreASL


    %% Processing Workflow

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

    if x.opts.ProcessModules(1)==1
        [~, x] = xASL_Iteration(x,'xASL_module_Structural');
        % The following DARTEL module is an optional extension of the structural module
        % to create population-specific templates
        if isfield(x.settings,'SegmentSPM12') && x.settings.SegmentSPM12 && x.nSubjects>1
            % in case we used SPM12 instead of CAT12 for segmentation,
            % we have to run DARTEL separately
            [~, x] = xASL_Iteration(x,'xASL_module_DARTEL');
        end
        % Now only check the availability of files when not running parallel
        if x.opts.nWorkers==1; xASL_adm_CreateFileReport(x); end        
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
    
    if x.opts.ProcessModules(2)==1
        [~, x] = xASL_Iteration(x,'xASL_module_ASL');
        % Now only check the availability of files when not running parallel
        if x.opts.nWorkers==1; xASL_adm_CreateFileReport(x); end    
    end


    % -----------------------------------------------------------------------------
    %% 3    xASL_module_Population
    % Performs all group-level processing & QC

    if x.opts.ProcessModules(3)==1
        [~, x] = xASL_Iteration(x,'xASL_module_Population');
    end

end