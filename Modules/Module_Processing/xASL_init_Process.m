function x = xASL_init_Process(x)
%xASL_init_Process Initialization for ExploreASL Processing
%
% FORMAT: x = xASL_init_Process(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Initialization fore ExploreASL processing
% i.e., datasetRoot/derivatives/ExploreASL/*
% This could contain data that was just copied from rawdata (converted to ExploreASL legacy format)
% or (partly) already processed data.
%
% This function is divided into the following (sub)functions:
% 1. Load dataset_description.json
% 2. Initialize environmental variables dependent on loaded data
% 3. xASL_init_DefineStudyStats
% 4. xASL_init_PrintCheckSettings
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright (c) 2015-2023 ExploreASL



    %% 1. Load dataset_description.json
    x = xASL_init_Loaddataset_description(x);

    
    %% 2. Initialize environmental variables dependent on loaded data
    % Define several ExploreASL environment parameters that are dependent on the loaded data,
    % such as the subfolders of the population folder, which maps and templates to use, which visualization settings, and which fields are deprecated.
    x = xASL_init_DefineDataParDependentSettings(x);
    

    %% 3. xASL_init_DefineStudyStats
    % Define study statistical parameters for this pipeline run
    x = xASL_init_DefineStudyStats(x);


    %% 4. xASL_init_PrintCheckSettings
    % Define & print settings (path, iWorker, nWorkers, Quality, DELETETEMP)
    x = xASL_init_PrintCheckSettings(x);
    


end






%% =======================================================================================================================
%% =======================================================================================================================
%% SUBFUNCTIONS start here





%% =======================================================================================================================
%% =======================================================================================================================
function [x] = xASL_init_DefineDataParDependentSettings(x)
    %xASL_init_DefineDataParDependentSettings Define ExploreASL environment parameters, dependent of loaded data
    %
    % FORMAT: [x] = xASL_init_DefineDataParDependentSettings(x)
    %
    % INPUT:
    %   x       - ExploreASL x structure (STRUCT, REQUIRED)
    %
    % OUTPUT:
    %   x       - ExploreASL x structure (STRUCT)
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % DESCRIPTION: Define several ExploreASL environment parameters that are dependent on the loaded data,
    % such as the subfolders of the population folder, which maps and templates to use, which visualization settings, and which fields are deprecated.
    %
    % EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % REFERENCES:  n/a
    %
    % Copyright 2015-2021 ExploreASL
    
    
    x = xASL_init_DefinePaths(x); % here we define the subfolders of the dataset (mostly of the Population folder)
    x = xASL_init_Toolboxes(x); % here we load SPM (PM: should this part go elsewhere? it is not data-dependent)
    
    
    %% --------------------------------------------------------------------------
    %% Reproducibility testing
    if ~isfield(x.settings,'bReproTesting')
        x.settings.bReproTesting = false;
    end
    
    % If the reproducibility is on, then delete any old RMS file
    if isfield(x.settings, 'bReproTesting')
        if x.settings.bReproTesting
            xASL_delete(fullfile(x.dir.xASLDerivatives, 'RMS_Reproducibility.mat'))
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
                x.D.PopDir = fullfile(x.dir.xASLDerivatives,'Population');
    
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
    
                % Population module
                x.D.SpaghettiDir        = fullfile(x.D.PopDir, 'SpaghettiPlots');
                x.S.StatsDir            = fullfile(x.D.PopDir, 'Stats');
                x.D.HistogramDir        = fullfile(x.D.PopDir, 'Histograms');
                x.D.StatsMaps           = fullfile(x.D.PopDir, 'StatsMaps');
                        
            
            
            end
        
    
        
    
        %% =======================================================================================================================
        %% =======================================================================================================================
        function x  = xASL_init_Toolboxes(x)
        %xASL_init_Toolboxes Here we load SPM
        %
        % FORMAT: x  = xASL_init_Toolboxes(x)
        %
        % INPUT:
        %   x       - ExploreASL x structure (STRUCT, REQUIRED)
        %
        % OUTPUT:
        %   x       - ExploreASL x structure (STRUCT)
        %
        % -----------------------------------------------------------------------------------------------------------------------------------------------------
        % DESCRIPTION: Here we load SPM
        %
        % EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
        %
        % -----------------------------------------------------------------------------------------------------------------------------------------------------
        % REFERENCES:  n/a
        %
        % Copyright 2015-2021 ExploreASL
    
    
            x.D.SPMDIR = fullfile(x.opts.MyPath, 'External', 'SPMmodified');
            x.D.SPMpath = x.D.SPMDIR;
            x.external.SPMVERSION = 'SPM12';
    
            if isdeployed
                spm('asciiwelcome');
                spm('defaults','pet');
                spm_jobman('initcfg');
            else
                addpath(fullfile(x.D.SPMpath ,'compat'));
            end
    
            spm_get_defaults('cmd_line',true);
    
        end




%% =======================================================================================================================
%% =======================================================================================================================
function x = xASL_init_PrintCheckSettings(x)
    %xASL_init_PrintCheckSettings Check whether pre-defined settings existed
    %
    % FORMAT: x = xASL_init_PrintCheckSettings(x)
    %
    % INPUT:
    %   x       - ExploreASL x structure (STRUCT, REQUIRED)
    %
    % OUTPUT:
    %   x       - ExploreASL x structure (STRUCT)
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % DESCRIPTION: Check whether pre-defined settings existed
    %
    % Prints these on the screen as the start of the pipeline.
    % Runs following steps:
    %
    % 1. Set default settings if not defined
    % 2. Print data/study specific settings
    % 3. Print warnings
    %
    % EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % REFERENCES:  n/a
    %
    % Copyright 2015-2021 ExploreASL
        
    
        %% -----------------------------------------------------------------------
        %% 1) Set default settings if not defined
        if ~isfield(x.settings,'Quality') || (x.settings.Quality~=0 && x.settings.Quality~=1)
            x.settings.Quality = 1;
            fprintf('%s\n', 'Default Quality=1 used (optimal quality)');
        end
        if ~isfield(x.settings,'DELETETEMP') || (x.settings.DELETETEMP~=0 && x.settings.DELETETEMP~=1)
            x.settings.DELETETEMP = 1;
        %     fprintf('%s\n','Default x.settings.DELETETEMP=1 used (delete files temporarily used for processing)');
        end
    
        %% -----------------------------------------------------------------------
        %% 2) Print data/study specific settings
        xASL_adm_BreakString('Additional Settings',[],[],1);
        
        if x.opts.nWorkers>1
            fprintf(['I am worker ' num2str(x.opts.iWorker) '/' num2str(x.opts.nWorkers) '\n']);
            fprintf('Note that the resulting number of scans mentioned below applies only to this worker\n');
        end
    
        fprintf('%s\n',[num2str(x.dataset.nTotalSubjects) ' scans - ' ...
            num2str(x.dataset.nExcluded) ' exclusions, resulting in ' ...
            num2str(x.dataset.nSubjects) ' scans of: ']);
    
        for iT=1:x.dataset.nTimePointsTotal
            fprintf('%s\n',['Longitudinal timePoint ' num2str(iT) ' = ' ...
                num2str(x.dataset.nTimePointTotalSubjects(iT)) ' scans - ' ...
                num2str(x.dataset.nTimePointExcluded(iT)) ' exclusions = ' ...
                num2str(x.dataset.nTimePointSubjects(iT)) ' scans']);
        end
    
        fprintf('%s\n',['ASL sessions: ' num2str(x.dataset.nSessions)]);
    
        fprintf('\n%s','Ancillary data, sets: ');
        if isfield(x.S,'SetsID')
                fprintf('%s\n',[num2str(size(x.S.SetsID,2)) ' sets are defined for ' ...
                    num2str(size(x.S.SetsID,1)) ' "SubjectsSessions"']);
    
                for iSet=1:size(x.S.SetsID,2)
                    fprintf(['Set ' num2str(iSet) ' = "' x.S.SetsName{iSet} '" options ']);
                    for iOption=1:size(x.S.SetsOptions{iSet},2)
                        fprintf(['"' x.S.SetsOptions{iSet}{iOption} '"']);
                        if iOption~= size(x.S.SetsOptions{iSet},2)
                            fprintf(' & ');
                        end
                    end
                    if      x.S.Sets1_2Sample(iSet)==1
                            fprintf(', codes for paired data');
                    elseif  x.S.Sets1_2Sample(iSet)==2
                            fprintf(', codes for two-sample data');
                    elseif  x.S.Sets1_2Sample(iSet)==3
                            fprintf([', continuous variate (with ' ...
                                num2str(length(unique(x.S.SetsID(:,iSet)))) ' unique values)']);
                    end                    
    
                    fprintf('\n');
                end
        else    
            fprintf('%s\n','No sets are defined');
        end
    
        if ~isfield(x.Q,'M0')
        %     warning('M0 option missing!');
        else
            fprintf('\n%s\n',['M0 option selected is "' num2str(x.Q.M0) '"']);
        end
        
        % Get the ExploreASL text width
        textBreakExploreASL = xASL_adm_BreakString([], [], false, 0, 0);
    
        if length(x.dir.xASLDerivatives)>length(textBreakExploreASL)
            fprintf('x.dir.xASLDerivatives  %s\n', x.dir.xASLDerivatives(1:length(textBreakExploreASL)-10));
            fprintf('          ... %s\n', x.dir.xASLDerivatives(length(textBreakExploreASL)-9:end));
        else
            fprintf('x.dir.xASLDerivatives  %s\n', x.dir.xASLDerivatives);
        end
        fprintf('x.settings.DELETETEMP %s\n',[num2str(x.settings.DELETETEMP) ' (delete temporary files)']);
        fprintf('x.settings.Quality    %s\n',[num2str(x.settings.Quality) ' (0 = fast try-out; 1 = normal high quality)']);
    
        fprintf('\n');
    
        %% -----------------------------------------------------------------------
        %% 3) Print warnings
        xASL_adm_BreakString('');
        fprintf('\n');
        field_symbol = {'subjectRegexp'};
    
        for iField=1:length(field_symbol)
            if ~isfield(x.dataset,field_symbol{iField})
                warning(['x.dataset' field_symbol{iField} ' was not defined!'])
            end
        end
    
        if ~isfield(x,'dir')
            warning('x.dir did not exist');
        else
            field_symbol = {'xASLDerivatives'};
            for iField=1:length(field_symbol)
                if ~isfield(x.dir,field_symbol{iField})
                    warning(['x.dir.' field_symbol{iField} ' was not defined!'])
                end
            end
        end
    
        if ~isfield(x,'dataset')
            warning('x.dataset didn''nt exist');
        else
            if ~isempty(regexp(x.dataset.subjectRegexp, '^(\^|)\.\*(\$|)$', 'once'))
                warning('Subject regexp not specific! Check that no wrong folders are included as subjects');
            end
        end
    
        fprintf('\n');
    
    end