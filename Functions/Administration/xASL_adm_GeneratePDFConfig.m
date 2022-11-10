function xASL_adm_GeneratePDFConfig(x, jsonPath, bPrintInvisible)
%xASL_qc_GeneratePDFConfig Print output QC
%
% FORMAT: xASL_qc_CreatePDF(x[, JsonPath, bVerbose])
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   jsonPath    - Path where json Should be placed  (OPTIONAL DEFAULT is .../derivatives/ExploreASL/ConfigReportPDF.json)
%   bPrintInvisible - boolean specifying wether you want the config file to contain invisible parameters as well
%
% OUTPUT: n/a
% OUTPUTFILE:
%   xASL_Report_SubjectName.pdf - Json Structure containing info for the
%   PDF generator
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function checks output struct x and places the values
% of the output in a Json file to be read by 
% 
% EXAMPLE: xASL_qc_GeneratePDFConfig(x);
% __________________________________
% Copyright (C) 2015-2021 ExploreASL

%% -----------------------------------------------------------------------------------------------
%  Admin


if nargin <2 || isempty(jsonPath)
    jsonPath = strcat(x.dir.xASLDerivatives, '/ConfigReportPDF.json');
    fprintf(['jsonPath was not specified, generating config file at ', jsonPath, '\n']);
end

if nargin <3 || isempty(bPrintInvisible)
    bPrintInvisible = false;
end

% Determine x.mat file
PathX = fullfile(x.dir.SUBJECTDIR,'x.mat');

% Check if x.mat exists
if ~exist(PathX, 'file')
    warning([PathX ' didnt exist, skipping xASL_qc_GeneratePDFConfig']);
    return;
end

% Load Output from x.mat
x = xASL_adm_LoadX(x, PathX, false); 
 
try
    %% Load Standard Configuration
    sPath = strcat(x.opts.MyPath, '/Functions/QualityControl/templateConfigReportPDF.json');
    
    if ~isfile(sPath)
        warning('Default PDF configuration not found, check the explore ASL directory.')
        return
    else 
        config = xASL_adm_LoadPDFConfig(sPath);
    end 

    clear OutputFields 
    OutputFields = fieldnames(x.Output);
    OutputStruct = struct([]);
    StandardValue = struct('Visible', false, 'Range', '', 'Unit','');

    for iModule=1:length(OutputFields)
        ModuleFields = fieldnames(x.Output.(OutputFields{iModule}));
        for iField=1:length(ModuleFields)   
           if isfield(config, OutputFields{iModule}) && isfield(config.(OutputFields{iModule}), ModuleFields{iField})
                OutputStruct(1).(OutputFields{iModule}).(ModuleFields{iField}) = config.(OutputFields{iModule}).(ModuleFields{iField});
           elseif bPrintInvisible
                OutputStruct(1).(OutputFields{iModule}).(ModuleFields{iField}) = StandardValue;
           end
        end
    end

    spm_jsonwrite(jsonPath, OutputStruct, struct('indent', ' '));
    
catch ME
    fprintf('%s\n', ['Creation of ' jsonPath ' failed:']);
    warning(ME.message);
end

end