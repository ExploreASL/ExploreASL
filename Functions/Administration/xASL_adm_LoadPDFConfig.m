function [Config] = xASL_adm_LoadPDFConfig(x)
% xASL_adm_LoadPDFConfig Loads parameters from a .json config file
% Essentialy an alias of spm_jsonread.
%
% FORMAT: [config] = xASL_adm_LoadPDFConfig(x)
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
% EXAMPLE: [config] = xASL_adm_LoadPDFConfig(x);
% __________________________________
% Copyright 2015-2022 ExploreASL



%% ------------------------------------------------------------------------
% 0. Argument check
Config = struct; % default

if nargin<1 || isempty(x)
    warning('x was not specified, aborting');
    return
end

ConfigPath = strcat(x.dir.xASLDerivatives, '/ConfigReportPDF.json');

%% ------------------------------------------------------------------------
%% 2. Load JSON file 
if isfile(ConfigPath)
    Config = spm_jsonread(ConfigPath);
else
    fprintf('configReportPDF.json does not exists, using default config \n');
    standardPath = strcat(x.opts.MyPath, '/Functions/QualityControl/templateConfigReportPDF.json');
    Config = spm_jsonread(standardPath);
end

%% ------------------------------------------------------------------------
%% 3. Deal with warnings
if isempty(fields(Config))
    warning('config struct empty, Check if config JSON is configured correctly.');
end

end