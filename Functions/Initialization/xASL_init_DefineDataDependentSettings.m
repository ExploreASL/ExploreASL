function [x] = xASL_init_DefineDataDependentSettings(x)
%xASL_init_DefineDataDependentSettings Define ExploreASL environment parameters, dependent of loaded data
%
% FORMAT: [x] = xASL_init_DefineDataDependentSettings(x)
%
% INPUT:
%   x       - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x       - ExploreASL x structure (STRUCT)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Define ExploreASL environment parameters, dependent of loaded data.
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
%
% Copyright 2015-2021 ExploreASL


x = xASL_init_DefinePaths(x);
x = xASL_init_Toolboxes(x);


%% --------------------------------------------------------------------------
%% Reproducibility testing
if ~isfield(x,'settings')
    x.settings = struct;
end
if ~isfield(x.settings,'bReproTesting')
    x.settings.bReproTesting = false;
end

% If the reproducibility is on, then delete the old RMS file
if isfield(x.settings, 'bReproTesting')
    if x.settings.bReproTesting
        xASL_delete(fullfile(x.D.ROOT, 'RMS_Reproducibility.mat'))
    end
end

%% --------------------------------------------------------------------------
%% Setting the option for pediatric template (this is normally set only for specific dataset and general xASL initialization does not have it)
% Set if the pediatric template field is set correctly
if ~isfield(x.settings,'Pediatric_Template') || isempty(x.settings.Pediatric_Template)
    x.settings.Pediatric_Template = false;
end

% We need to define which of the pediatric templates to use
if x.settings.Pediatric_Template
	% Set the default
	if ~isfield(x,'Pediatric_Type') || isempty(x.Pediatric_Type)
		x.Pediatric_Type = 'infant-1yr';
	end
	
	switch (x.Pediatric_Type)
		case {'1yr','infant_1yr','Infant_1yr'}
			x.Pediatric_Type = 'infant-1yr';
		case {'2yr','infant_2yr','Infant_2yr'}
			x.Pediatric_Type = 'infant-2yr';
		case {'neo','infant_neo','Infant_neo','neonate'}
			x.Pediatric_Type = 'infant-neo';
	end
end


%% --------------------------------------------------------------------------
%% Atlases and templates in pediatric version
if x.settings.Pediatric_Template
	x.D.MapsDir             = fullfile(x.opts.MyPath,'Maps', x.Pediatric_Type);
	x.D.MapsSPMmodifiedDir  = fullfile(x.opts.MyPath,'External', 'SPMmodified', 'MapsAdded', x.Pediatric_Type);
	x.D.ResliceRef          = fullfile(x.opts.MyPath,'External', 'SPMmodified', 'MapsAdded', x.Pediatric_Type,'rgrey.nii');
	x.D.IdentityTransfRef   = fullfile(x.opts.MyPath,'External', 'SPMmodified', 'MapsAdded', x.Pediatric_Type,'Identity_Deformation_y_T1.nii');
	x.D.TemplateDir         = fullfile(x.opts.MyPath,'Maps', 'Templates', x.Pediatric_Type);
	x.D.AtlasDir            = fullfile(x.opts.MyPath,'External', 'Atlases', x.Pediatric_Type);
end


%% --------------------------------------------------------------------------
%% Manage input parameters ExploreASL course
Fields = {'bLesionFilling' 'bAutoACPC' 'bGetControlLabelOrder'};
Defaults = [true true true];

for iL=1:length(Fields)
    if ~isfield(x.settings,Fields{iL})
        x.settings.(Fields{iL}) = Defaults(iL);
    elseif isempty(x.settings.(Fields{iL}))
        x.settings.(Fields{iL}) = Defaults(iL);
    elseif ~islogical(x.settings.(Fields{iL}))
        if x.settings.(Fields{iL})==1
            x.settings.(Fields{iL}) = true;
        elseif x.settings.(Fields{iL})==0
            x.settings.(Fields{iL}) = false;
        else
            warning([Fields{iL} ' was not true or false, set to default:' num2str(Defaults(iL))]);
            x.settings.(Fields{iL}) = Defaults(iL);
        end
    end
end

% In a future release we could set all defaults based on a JSON here
if ~isfield(x.modules.structural, 'bSegmentSPM12')
	x.modules.structural.bSegmentSPM12 = false;
end
if ~isfield(x.modules.asl, 'M0_conventionalProcessing')
	x.modules.asl.M0_conventionalProcessing = false;
end

%% Atlases and templates (those may have been removed by the clean x function)
x = xASL_init_MapsAndAtlases(x);

%% Visualization (those may have been removed by the clean x function)
x = xASL_init_VisualizationSettings(x);

%% Check deprecated fields
x = xASL_io_CheckDeprecatedFieldsX(x, 0);


end





%% =======================================================================================================================
%% =======================================================================================================================
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
                x.D.ATTCheckDir         = fullfile(x.D.PopDir, 'ATT_Check'); % Putting TT and ATT maps to the same folder
                x.D.TexCheckDir         = fullfile(x.D.PopDir, 'Tex_Check');
                x.D.TemplatesStudyDir   = fullfile(x.D.PopDir, 'Templates');
    
                % POPULATION module
                x.D.SpaghettiDir        = fullfile(x.D.PopDir, 'SpaghettiPlots');
                x.S.StatsDir            = fullfile(x.D.PopDir, 'Stats');
                x.D.HistogramDir        = fullfile(x.D.PopDir, 'Histograms');
                x.D.StatsMaps           = fullfile(x.D.PopDir, 'StatsMaps');
                
            end
        end
    
    
    end
    

    

