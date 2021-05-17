%% ==================================================================================
function x  = xASL_init_Toolboxes(x)
%xASL_init_Toolboxes Check & load ancillary toolboxes, versions and paths.

    x.SPMDIR = fullfile(x.MyPath, 'External', 'SPMmodified');
    x.SPMpath = x.SPMDIR;
    x.SPMVERSION = 'SPM12';

    if isdeployed
        spm('asciiwelcome');
        spm('defaults','pet');
        spm_jobman('initcfg');
    else
        addpath(fullfile(x.SPMpath ,'compat'));
    end

    spm_get_defaults('cmd_line',true);

end

