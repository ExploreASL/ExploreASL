function [x] = xASL_init_DefineIndependentSettings(x)
%xASL_init_DefineIndependentSettings Define ExploreASL environment parameters, independent of loaded data
%
% FORMAT: [x] = xASL_init_DefineIndependentSettings(x)
%
% INPUT:
%   x       - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:
%   x       - ExploreASL x structure (STRUCT)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Define ExploreASL environment parameters, independent of loaded data.
%
% EXAMPLE:     This is part of the initialization workflow. Check out the usage there.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:  n/a
%
% Copyright 2015-2021 ExploreASL


%% --------------------------------------------------------------------------------------------------------------------
%% Default parameters
x.settings.stopAfterErrors = Inf; % set to a high number (or Inf) when running a complete batch overnight
x.settings.dryRun = false; % set to true to skip all real processing & just print all parameters
x.settings.bOverwrite = true;

if ~exist('groot','builtin')
    % before R2012b
    set(0,'DefaultTextInterpreter','none')
else
    % R2012b+
    set(groot,'DefaultTextInterpreter','none')
end

% Get version
try
    versionFile = dir(fullfile(x.opts.MyPath, 'VERSION*'));
    versionFile = versionFile.name;
    x.Version = versionFile(9:end);

    if isempty(x.Version)
        warning('Unknown ExploreASL version number');
        x.Version = 'VERSION_unknown';
    end
catch
    warning('Could not obtain ExploreASL version, version file missing');
    x.Version = 'VERSION_File_Missing';
end


if ~isfield(x,'Q')
    x.Q = struct;
end


%% Atlases and templates
x = xASL_init_MapsAndAtlases(x);

%% Visualization
x = xASL_init_VisualizationSettings(x);


end





