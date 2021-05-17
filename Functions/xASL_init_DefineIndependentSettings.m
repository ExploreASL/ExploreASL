%% ==================================================================================
%% ==================================================================================
function [x] = xASL_init_DefineIndependentSettings(x)
%xASL_init_DefineIndependentSettings Define ExploreASL environment
%parameters, independent of loaded data



%% --------------------------------------------------------------------------------------------------------------------
%% Default parameters
x.stopaftererrors = Inf; % set to a high number (or Inf) when running a complete batch overnight
x.dryrun = false; % set to true to skip all real processing & just print all parameters
x.bOverwrite = true;

if ~exist('groot','builtin')
    % before R2012b
    set(0,'DefaultTextInterpreter','none')
else
    % R2012b+
    set(groot,'DefaultTextInterpreter','none')
end

% Get version
if ~isdeployed
    VersionPath = xASL_adm_GetFileList(x.MyPath, '^VERSION.*$', 'FPList', [0 Inf]);
    if isempty(VersionPath)
        warning('Could not obtain ExploreASL version, version file missing');
    else
        [~, Fname, Fext] = fileparts(VersionPath{1});
        x.Version = [Fname(9:end) Fext];
    end
else
    % Output of compiled ExploreASL version
    try
        versionFile = dir(fullfile(x.MyPath, 'VERSION*'));
        versionFile = versionFile.name;
        x.Version = versionFile(9:end);
    catch
        x.Version = '1.0.0 (TMP)';
    end
end

if ~isfield(x,'Q')
    x.Q = struct;
end

x = xASL_init_DefinePaths(x);
x = xASL_init_Toolboxes(x); % Initialize toolboxes
x = xASL_init_VisualizationSettings(x); % visual settings



end


