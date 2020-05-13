function x  = xASL_init_ExploreASL_directories( x )
%xASL_init_ExploreASL_directories Part of master script ExploreASL, loading definitions & pathnames

%% Setting the option for pediatric template (this is normally set only for specific dataset and general xASL initialization does not have it)
% Set if the pediatric template field is set correctly
if ~isfield(x,'Pediatric_Template') || isempty(x.Pediatric_Template)
    x.Pediatric_Template = false;
end

% We need to define which of the pediatric templates to use
if x.Pediatric_Template
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

%% ExploreASL
% General
x.D.FunctionsDir        = fullfile(x.MyPath,'Functions');
x.D.FunctionsExternDir  = fullfile(x.MyPath,'External','SPMmodified','xASL');

% Atlases and templates in pediatric version
if x.Pediatric_Template
	x.D.MapsDir             = fullfile(x.MyPath,'Maps',x.Pediatric_Type);
	x.D.MapsSPMmodifiedDir  = fullfile(x.MyPath,'External','SPMmodified','MapsAdded',x.Pediatric_Type);
	x.D.ResliceRef          = fullfile(x.MyPath,'External','SPMmodified','MapsAdded',x.Pediatric_Type,'rgrey.nii');
	x.D.IdentityTransfRef   = fullfile(x.MyPath,'External','SPMmodified','MapsAdded',x.Pediatric_Type,'Identity_Deformation_y_T1.nii');
	x.D.TemplateDir         = fullfile(x.MyPath,'Maps','Templates',x.Pediatric_Type);
	x.D.AtlasDir            = fullfile(x.MyPath,'External','AtlasesNonCommercial',x.Pediatric_Type);
else
	% Atlases and templates
	x.D.MapsDir             = fullfile(x.MyPath,'Maps');
	x.D.MapsSPMmodifiedDir  = fullfile(x.MyPath,'External','SPMmodified','MapsAdded');
	x.D.ResliceRef          = fullfile(x.MyPath,'External','SPMmodified','MapsAdded','rgrey.nii');
	x.D.IdentityTransfRef   = fullfile(x.MyPath,'External','SPMmodified','MapsAdded','Identity_Deformation_y_T1.nii');
	x.D.TemplateDir         = fullfile(x.MyPath,'Maps','Templates');
	x.D.AtlasDir            = fullfile(x.MyPath,'External','AtlasesNonCommercial');
end

% Prefixes standard space
x.D.CBFPreFix_Resliced  = 'qCBF';
x.D.c_PreFix{1}         = 'rc1T1';
x.D.c_PreFix{2}         = 'rc2T1';


%% Study-specific
if  x.ProcessData
    x.D.PopDir           = fullfile(x.D.ROOT,'Population');

    % Structural module
    x.D.T1CheckDir          = fullfile( x.D.PopDir, 'T1Check');
    x.D.TissueVolumeDir     = fullfile( x.D.PopDir, 'TissueVolume');
    x.D.CoregDir            = fullfile( x.D.PopDir, 'T1wCoregCheck');
    x.D.FLAIR_CheckDir      = fullfile( x.D.PopDir, 'FLAIRCheck' );
    x.D.FLAIR_REGDIR        = fullfile( x.D.PopDir, 'FLAIRReg'   );
    x.D.FlowFieldCheck      = fullfile( x.D.PopDir, 'FlowFieldCheck' );
    x.D.LongRegCheckDir     = fullfile( x.D.PopDir, 'LongRegCheck');
    x.D.LesionCheckDir      = fullfile( x.D.PopDir, 'LesionCheck');
    x.D.ROICheckDir         = fullfile( x.D.PopDir, 'ROICheck');

    % ASL module
    x.D.ASLCheckDir         = fullfile( x.D.PopDir, 'ASLCheck');
    x.D.MotionDir           = fullfile( x.D.PopDir, 'MotionASL');
    x.D.ExclusionDir        = fullfile( x.D.PopDir, 'Exclusion');
    x.D.DICOMparameterDir   = fullfile( x.D.PopDir, 'DICOMparameters');
    x.D.SNRdir              = fullfile( x.D.PopDir, 'SD_SNR');
    x.D.M0CheckDir          = fullfile( x.D.PopDir, 'M0Check');
    x.D.M0regASLdir         = fullfile( x.D.PopDir, 'M0Reg_ASL');
    x.D.SliceCheckDir       = fullfile( x.D.PopDir, 'SliceGradientCheck');
    x.D.RawDir              = fullfile( x.D.PopDir, 'RawSourceIMCheck');
    x.D.RawEPIdir           = fullfile( x.D.PopDir, 'Raw_EPI_Check');
    x.D.T1_ASLREGDIR        = fullfile( x.D.PopDir, 'T1_ASLReg');
    x.D.TTCheckDir          = fullfile( x.D.PopDir, 'ATT_Check');
    x.D.TemplatesStudyDir   = fullfile( x.D.PopDir, 'Templates');

    % ANALYZE module
    x.SpaghettiDir          = fullfile( x.D.PopDir, 'SpaghettiPlots');
    x.S.StatsDir            = fullfile( x.D.PopDir, 'Stats');
    x.HistogramDir          = fullfile( x.D.PopDir, 'Histograms');
    x.StatsMaps             = fullfile( x.D.PopDir, 'StatsMaps');

    xASL_adm_CreateDir( x.D.PopDir );
end



end
