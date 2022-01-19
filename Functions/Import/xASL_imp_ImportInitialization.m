function x = xASL_imp_ImportInitialization(x)
%xASL_imp_ImportInitialization Initialization before xASL_module_Import
%
% FORMAT: x = xASL_imp_ImportInitialization(x)
%
% INPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% OUTPUT:
%   x        - Struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:    Initialization before xASL_module_Import.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:        n/a
% __________________________________
% Copyright (c) 2015-2022 ExploreASL


    %% Initialization

    % Initialize x struct
    x = xASL_init_SubStructs(x);
    
    % We expect x.opts.DatasetRoot to be the study root directory, but if it is not defined,
    % then the user probably used a path to a descriptive JSON file instead
    if isfield(x, 'dir') && isfield(x.dir, 'DatasetRoot') && isempty(x.dir.DatasetRoot)
        x.dir.DatasetRoot = xASL_fileparts(x.opts.DatasetRoot);
    end
    
    % For the xASL_Iterate support we need x.D.ROOT as well
    x.D.ROOT = x.dir.DatasetRoot;
    
    % Currently fixed import settings
    x.modules.import.settings.bCopySingleDicoms = false;
    x.modules.import.settings.bUseDCMTK = true;
    x.modules.import.settings.bCheckPermissions = false;
    
    % Basic parameter checks
    x = xASL_imp_BasicParameterChecks(x);

    % Initialize the import setup
    if x.opts.ImportModules(1) || x.opts.ImportModules(2) || x.opts.Deface || (x.opts.bLoadData && x.opts.bLoadableData)
        x.modules.import.imPar = xASL_imp_Initialize(x.dir.DatasetRoot, x.dir.sourceStructure);
    else
        x.modules.import.imPar = NaN;
    end

    
    %% Determine subject/session/run structure from sourcedata, temp data or rawdata
    if x.opts.bLoadData && x.opts.bLoadableData
        x = xASL_imp_DetermineSubjectStructure(x);
    end
    
    % Create logging directory if it does not exist already
    xASL_adm_CreateDir(fullfile(x.dir.DatasetRoot,'derivatives','ExploreASL','log'));


end



