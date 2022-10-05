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
% Copyright 2015-2021 ExploreASL


    %% General

    % Prefixes standard space
    x.D.CBFPreFix_Resliced  = 'qCBF';
    x.D.c_PreFix{1} = 'rc1T1';
    x.D.c_PreFix{2} = 'rc2T1';

    %% Study-specific
    if and(isfield(x.D, 'ROOT'), isfield(x.opts, 'bProcessData'))
        if x.opts.bProcessData || x.opts.bLoadData
            x.D.PopDir = fullfile(x.D.ROOT,'Population');

            % Structural module
            x.D.T1CheckDir          = fullfile(x.D.PopDir, 'T1Check');
            x.D.TissueVolumeDir     = fullfile(x.D.PopDir, 'TissueVolume');
            x.D.CoregDir            = fullfile(x.D.PopDir, 'T1wCoregCheck');
            x.D.FLAIR_CheckDir      = fullfile(x.D.PopDir, 'FLAIRCheck' );
    		x.D.T1c_CheckDir        = fullfile(x.D.PopDir, 'T1cCheck' );
    		x.D.T2_CheckDir         = fullfile(x.D.PopDir, 'T2Check' );
            x.D.FLAIR_REGDIR        = fullfile(x.D.PopDir, 'FLAIRReg'   );
            x.D.FlowFieldCheck      = fullfile(x.D.PopDir, 'FlowFieldCheck' );
            x.D.LongRegCheckDir     = fullfile(x.D.PopDir, 'LongRegCheck');
            x.D.LesionCheckDir      = fullfile(x.D.PopDir, 'LesionCheck');
            x.D.ROICheckDir         = fullfile(x.D.PopDir, 'ROICheck');

            % ASL module
            x.D.ASLCheckDir         = fullfile(x.D.PopDir, 'ASLCheck');
            x.D.MotionDir           = fullfile(x.D.PopDir, 'MotionASL');
            x.D.ExclusionDir        = fullfile(x.D.PopDir, 'Exclusion');
            x.D.DICOMparameterDir   = fullfile(x.D.PopDir, 'DICOMparameters');
            x.D.SNRdir              = fullfile(x.D.PopDir, 'SD_SNR');
            x.D.M0CheckDir          = fullfile(x.D.PopDir, 'M0Check');
            x.D.M0regASLdir         = fullfile(x.D.PopDir, 'M0Reg_ASL');
            x.D.SliceCheckDir       = fullfile(x.D.PopDir, 'SliceGradientCheck');
            x.D.RawDir              = fullfile(x.D.PopDir, 'RawSourceIMCheck');
            x.D.RawEPIdir           = fullfile(x.D.PopDir, 'Raw_EPI_Check');
            x.D.T1_ASLREGDIR        = fullfile(x.D.PopDir, 'T1_ASLReg');
            x.D.TTCheckDir          = fullfile(x.D.PopDir, 'ATT_Check');
			x.D.TExchCheckDir       = fullfile(x.D.PopDir, 'TExch_Check');
            x.D.TemplatesStudyDir   = fullfile(x.D.PopDir, 'Templates');

            % POPULATION module
            x.D.SpaghettiDir        = fullfile(x.D.PopDir, 'SpaghettiPlots');
            x.S.StatsDir            = fullfile(x.D.PopDir, 'Stats');
            x.D.HistogramDir        = fullfile(x.D.PopDir, 'Histograms');
            x.D.StatsMaps           = fullfile(x.D.PopDir, 'StatsMaps');
            
        end
    end


end

