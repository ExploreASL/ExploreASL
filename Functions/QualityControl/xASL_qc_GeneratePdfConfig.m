function [config] = xASL_qc_GeneratePdfConfig(x, subject, bOverwrite)
% xASL_qc_GeneratePdfConfig Generates a JSON based on all values in x.Output
%
% FORMAT: xASL_qc_GeneratePdfConfig(x[, subject, bOverwrite])
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   subject     - subject name (OPTIONAL, default = x.SUBJECT)
%   bOverWrite  - boolean to determine if current configReportPDF.json should be overwritten. (OPTIONAL, default == true)
%
% OUTPUT: 
%   In the ExploreASL Derivatives folder in the subject directory, a configReportPDF.json is generated containing all quality parameters.
%   Quality parameters are taken from x.Output and x.OutputIm
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function dumps all value from x.Output of subject into a configuration file
%               This configuration can be read by xASL_qc_GeneratePdfConfig to create a pdf that contains all quality values.
% 
% EXAMPLE: xASL_qc_GeneratePdfConfig(x, 'sub-001');
% __________________________________
% Copyright (C) 2015-2023 ExploreASL


% check input
if nargin < 1 || isempty(x)
    error('No x structure provided');
end

if nargin < 2 || isempty(subject)
    if isfield(x, 'SUBJECT')
        subject = x.SUBJECT;
    else 
        warning('subject or x.SUBJECT missing, this might go wrong');
    end
else
    x.SUBJECT = subject;
end

if nargin < 3 || isempty(bOverwrite)
   bOverwrite = true;
end

% Fix <SESSION> not existing
if ~isfield(x, 'SESSION')
    if isfield(x, 'SESSIONS') && ~isempty(x.SESSIONS)
        x.SESSION = x.SESSIONS{1};
    else
        warning('x.SESSIONS was missing, this might go wrong');
    end
end

% Fix x.P not existing
if ~isfield(x, 'P') || isempty(x.P) || isempty(fields(x.P))
    if ~isfield(x, 'dir')
        warning('x.dir missing, this might go wrong');
    end
    if ~isfield(x.dir, 'SUBJECTDIR')
        x.dir.SUBJECTDIR = fullfile(x.dir.xASLDerivatives, x.SUBJECT);
    end
    if ~isfield(x.dir, 'SESSIONDIR')
        x.dir.SESSIONDIR = fullfile(x.dir.SUBJECTDIR, x.SESSION);
    end

    x = xASL_init_FileSystem(x);
end

% Determine x.mat file
PathX = fullfile(x.dir.xASLDerivatives, subject, 'x.mat');

% Check if x.mat file exists already
if ~exist(PathX, 'file')
    warning([PathX ' didnt exist, skipping xASL_qc_CreateOutputPDF']);
    return;
end

x = xASL_adm_LoadX(x, PathX, true); % Assume memory x is newer than x.mat

% Make sure that the directory exists
PrintDir = fullfile(x.dir.xASLDerivatives);
xASL_adm_CreateDir(PrintDir);

% Print the title
fprintf('Creating Default Json in:   \n');
fprintf([PrintDir '\n']);

% Parse the entire Json Stack automatically making all the pages.
config = xASL_sub_createDefaultJson(x);

% Write Pdf configuration
xASL_io_WriteJson(fullfile(PrintDir, 'configReportPDF.json'), config, bOverwrite);

end

function [config] = xASL_sub_createDefaultJson(x)
    config = struct();
    if ~isfield(x, 'Output') || isempty(x.Output) || isempty(fields(x.Output))
        warning('x.Output didnt exist, skipping xASL_qc_GenerateJsonTemplate');
        return;
    end
    if ~isfield(x, 'Output_im') || isempty(x.Output_im) || isempty(fields(x.Output_im))
        warning('x.Output_im didnt exist, skipping xASL_qc_GenerateJsonTemplate');
        return;
    end

    config.pages = struct();
    modules = fieldnames(x.Output);
    for module = 1:size(modules)
        config.pages(module).category = 'metadata';
        config.pages(module).type = 'page';
        config.pages(module).pageIdentifier = modules{module};
        config.pages(module).content = xASL_sub_createModuleContent(x, modules{module});
    end

end

function content = xASL_sub_createModuleContent(x, module)
    qc_parameters = fieldnames(x.Output.(module));
    content = cell(size(qc_parameters,1),1);

    for field = 1:size(qc_parameters)
        qc_content = struct();
        qc_content.category = 'content';
        qc_content.type = 'QCValues';
        qc_content.parameter = qc_parameters{field};
        qc_content.module = module;
        content{field} = qc_content;
    end

    % Add the qc_images
    qc_content = struct();
    qc_content.category = 'content';
    qc_content.type = 'QCValues';
    qc_content.parameter = 'qc_images';
    qc_content.position = '[0.4 0.25]';
    qc_content.size = '[0.6 0.6]';
    qc_content.module = module;
    content{end+1} = qc_content;



end