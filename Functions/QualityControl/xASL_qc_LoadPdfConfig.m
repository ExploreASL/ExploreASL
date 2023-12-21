function [config] = xASL_qc_LoadPdfConfig(x, configPath)
% xASL_qc_LoadPdfConfig loads parameters from a .json config file to be used in xASL_qc_GenerateReport.
% Using the config file, the user can define the parameters of the report.
% If no config file is given, it will first look in the derivatives folder if configReportPDF.json exists.
% If nothing exists there a default config file is loaded and used. (source at in the ExploreASL/Functions/QualityControl folder)
%
% FORMAT: [config] = xASL_qc_LoadPdfConfig(x[, configPath])
%
% INPUT:
%   x           - Struct containing all ExploreAsl Parameters (REQUIRED)
%   configPath  - Path to the config file (OPTIONAL)
%
%   OUTPUT:
%   config      - Struct containing all parameters of the jsonfile
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function loads a *.json config file, and is used at 
% xASL_qc_GenerateReport to make a report with defined parameters:
%
% 1. Argument check
% 2. Load JSON file
% 3. Deal with warnings and return output
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [config] = xASL_qc_LoadPdfConfig(x);
% __________________________________
% Copyright 2015-2023 ExploreASL



%% ------------------------------------------------------------------------
% 1. Argument check

if nargin<1 || isempty(x)
    error('No x structure provided');
end

if nargin<2 || isempty(configPath)
    configPath = fullfile(x.dir.xASLDerivatives, 'configReportPDF.json');
end


%% ------------------------------------------------------------------------
%% 2. Load JSON file 
if exist(configPath, 'file')
    config = xASL_io_ReadJson(configPath);
else
    fprintf(['Custom PDF definitions not found, using default configuration for PDF generation.\n']);
    xASL_qc_GeneratePdfConfig(x, x.SUBJECT, false);
    config = xASL_io_ReadJson(configPath);
end

%% ------------------------------------------------------------------------
%% 3. Deal with warnings
if isempty(fields(config))
    warning('Check if configuration JSON for the PDF report is correct, no parameters were loaded.');
end

end