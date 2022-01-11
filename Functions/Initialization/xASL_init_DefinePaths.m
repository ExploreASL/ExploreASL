function [x] = xASL_init_DefinePaths(x)
%xASL_init_DefinePaths Define paths used by ExploreASL
%
% FORMAT: [x] = xASL_init_DefinePaths(x)
%
% INPUT:
%   x       - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x       - ExploreASL x structure (STRUCT)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Define paths used by ExploreASL.
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
%
% Copyright (c) 2015-2022 ExploreASL


    %% General

    % Prefixes standard space
    x.D.CBFPreFix_Resliced  = 'qCBF';
    x.D.c_PreFix{1} = 'rc1T1';
    x.D.c_PreFix{2} = 'rc2T1';

    % Atlases and templates
    if isfield(x.opts, 'MyPath')
        x.D.MapsDir             = fullfile(x.opts.MyPath, 'Maps');
        x.D.MapsSPMmodifiedDir  = fullfile(x.opts.MyPath, 'External', 'SPMmodified', 'MapsAdded');
        x.D.ResliceRef          = fullfile(x.opts.MyPath, 'External', 'SPMmodified', 'MapsAdded', 'rgrey.nii');
        x.D.IdentityTransfRef   = fullfile(x.opts.MyPath, 'External', 'SPMmodified', 'MapsAdded', 'Identity_Deformation_y_T1.nii');
        x.D.TemplateDir         = fullfile(x.opts.MyPath, 'Maps', 'Templates');
        x.D.AtlasDir            = fullfile(x.opts.MyPath, 'External', 'Atlases');
    else
        warning('MyPath field not defined...');
    end

    %% Study-specific
    if and(isfield(x.D, 'ROOT'), isfield(x.opts, 'bProcessData'))
        if x.opts.bProcessData || x.opts.bLoadData
            
            % Base directories
            x.D.PopDir              = fullfile(x.D.ROOT,   'Population');
            x.D.perf                = fullfile(x.D.PopDir, 'perf');
            x.D.anat                = fullfile(x.D.PopDir, 'anat');
            x.D.stats               = fullfile(x.D.PopDir, 'stats');

            % Structural module
            x.D.T1w                 = fullfile(x.D.anat,   'T1w');
            x.D.FLAIR               = fullfile(x.D.anat,   'FLAIR');
            x.D.T1CheckDir          = fullfile(x.D.anat,   'T1Check');
            x.D.TissueVolumeDir     = fullfile(x.D.anat,   'TissueVolume');
            x.D.CoregDir            = fullfile(x.D.anat,   'T1wCoregCheck');
            x.D.FLAIR_CheckDir      = fullfile(x.D.anat,   'FLAIRCheck' );
    		x.D.T1c_CheckDir        = fullfile(x.D.anat,   'T1cCheck' );
    		x.D.T2_CheckDir         = fullfile(x.D.anat,   'T2Check' );
            x.D.FLAIR_REGDIR        = fullfile(x.D.anat,   'FLAIRReg'   );
            x.D.FlowFieldCheck      = fullfile(x.D.anat,   'FlowFieldCheck' );
            x.D.LongRegCheckDir     = fullfile(x.D.anat,   'LongRegCheck');
            x.D.LesionCheckDir      = fullfile(x.D.anat,   'LesionCheck');
            x.D.ROICheckDir         = fullfile(x.D.anat,   'ROICheck');

            % ASL module
            x.D.ASL                 = fullfile(x.D.perf,   'ASL');
            x.D.ASLCheckDir         = fullfile(x.D.perf,   'ASLCheck');
            x.D.MotionDir           = fullfile(x.D.perf,   'MotionASL');
            x.D.ExclusionDir        = fullfile(x.D.perf,   'Exclusion');
            x.D.DICOMparameterDir   = fullfile(x.D.perf,   'DICOMparameters');
            x.D.SNRdir              = fullfile(x.D.perf,   'SD_SNR');
            x.D.M0CheckDir          = fullfile(x.D.perf,   'M0Check');
            x.D.M0regASLdir         = fullfile(x.D.perf,   'M0Reg_ASL');
            x.D.SliceCheckDir       = fullfile(x.D.perf,   'SliceGradientCheck');
            x.D.RawDir              = fullfile(x.D.perf,   'RawSourceIMCheck');
            x.D.RawEPIdir           = fullfile(x.D.perf,   'Raw_EPI_Check');
            x.D.T1_ASLREGDIR        = fullfile(x.D.perf,   'T1_ASLReg');
            x.D.TTCheckDir          = fullfile(x.D.perf,   'ATT_Check');
            x.D.TemplatesStudyDir   = fullfile(x.D.perf,   'Templates');

            % POPULATION module
            x.D.SpaghettiDir        = fullfile(x.D.stats,  'SpaghettiPlots');
            x.S.StatsDir            = fullfile(x.D.stats,  'Statistics');
            x.D.HistogramDir        = fullfile(x.D.stats,  'Histograms');
            x.D.StatsMaps           = fullfile(x.D.stats,  'StatsMaps');
            x.D.SliceGradient       = fullfile(x.D.stats,  'SliceGradient');
            
        end
    end


end

