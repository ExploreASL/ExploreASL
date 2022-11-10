function [Config] = xASL_adm_LoadPDFConfig(ConfigPath)
% xASL_adm_LoadPDFConfig Loads parameters from a .json config file
% Essentialy an alias of spm_jsonread.
%
% FORMAT: [config] = xASL_adm_LoadPDFConfig(ConfigPath)
%
% INPUT:
%   ConfigPath  - path to json REQUIRES
%
%   OUTPUT:
%   config      - Struct containing all parameters of the jsonfile
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads a *.json Config file, and is used at 
% xASL_qc_GenerateReport to make a report with defined parameters:
%
% 1. Load JSON file
% 2. Deal with warnings
% 3. return a struct containing the config parameters
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [config] = xASL_adm_LoadPDFConfig('/ExploreASL/Functions/QualityControl/DefaultPDFConfig.json');
% __________________________________
% Copyright 2015-2022 ExploreASL



%% ------------------------------------------------------------------------
% 0. Argument check
Config = struct; % default

if nargin<1 || isempty(ConfigPath)
    warning('ConfigPath was not specified, aborting');
    return
end

if ~isfile(ConfigPath)
    warning(['Could not find ' ConfigPath])
    return 
end

%% ------------------------------------------------------------------------
%% 2. Load JSON file 
if exist(ConfigPath,'file') 
    Config = spm_jsonread(ConfigPath);
end

%% ------------------------------------------------------------------------
%% 3. Deal with warnings
if isempty(fields(Config))
    warning('config struct empty, Check if config JSON is configured correctly.');
end

end