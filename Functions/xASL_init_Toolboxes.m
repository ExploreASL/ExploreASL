function x  = xASL_init_Toolboxes(x)
%xASL_init_Toolboxes Check & load ancillary toolboxes, versions and paths.
% Part of ExploreASL_Initialize.m

x.SPMDIR      = fullfile(x.MyPath, 'External', 'SPMmodified');
x.SPMpath     = x.SPMDIR;
x.SPMVERSION  = 'SPM12';

% This is committed since we now provide the SPM toolbox within ExploreASL
% Speeds up start-up of ExploreASL
% if      exist('spm_jobman') % SPM is already loaded
%         [x.SPMDIR, x.SPMVERSION] = xASL_adm_CheckSPM('FMRI' );
%         x.SPMpath     = x.SPMDIR;
% else    [x.SPMDIR, x.SPMVERSION] = xASL_adm_CheckSPM('FMRI', x.SPMpath );
% end
%
% if ~isfield(x, 'SPMpath') || ~isfield(x,'SPMDIR')
%     x.SPMpath  = uigetdir(cd, 'Select the root-directory where SPM is installed ');
%     if x.SPMpath==0
%        return
%     end
% end

if ~isdeployed % Modification for compilation
    addpath(fullfile(x.SPMpath ,'compat') );
end

% Check conflict with SPM fieldtrip
% if ~isempty(which('fieldtrip2spss'))
%     fprintf('%s\n','CAVE: meannan (and t-test?) conflicts between "fieldtrip" (SPM extension) and matlabs original!!!');
%     fprintf('%s\n','Therefore, fieldtrip has been removed from the Matlab pathdef');
% %      rmpath( fullfile(x.SPMDIR, 'external', 'fieldtrip'));
% % @ Paul: waarom heb je dit in comment gezet? Moet dit op andere manier?
% end

if ~isdeployed % Modification for compilation
    if ismac
        %DIPpath = fullfile(x.MyPath, 'External', 'DIP', 'dip_2.9_MacOS');
    elseif isunix
        %DIPpath = fullfile(x.MyPath, 'External', 'DIP', 'dip_2.9_lnx64');
    elseif ispc
        %DIPpath     = fullfile(x.MyPath, 'External', 'DIP', 'dip_2.8.1_win64');
    end
    %addpath(fullfile(DIPpath,'common','dipimage') );

    %if ~isempty(which('dip_initialise'))
    %    try
    %        dip_initialise;
    %    end
    %end
end

if isdeployed
    spm('asciiwelcome');
    spm('defaults','pet');
    spm_jobman('initcfg');
end

spm_get_defaults('cmd_line',true);

end
