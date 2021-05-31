%% ==================================================================================
function x  = xASL_init_Toolboxes(x)
%xASL_init_Toolboxes Check & load ancillary toolboxes, versions and paths
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
% DESCRIPTION: Check & load ancillary toolboxes, versions and paths.
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
%
% Copyright 2015-2021 ExploreASL


    x.D.SPMDIR = fullfile(x.MyPath, 'External', 'SPMmodified');
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

