function [x] = xASL_qc_GenerateReport(x, subject, modules, bOverWrite)
% xASL_qc_GenerateReport Generates a PDF report based on a predefined Json configuration file
%
% FORMAT: xASL_qc_GenerateReport(x [,subject, modules, bOverwrite])
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   subject     - subject name (OPTIONAL, default = x.SUBJECT)
%   modules     - structure with name(s) of the modules that need to be added to the PDF report (OPTIONAL, default==all modules in x.output)
%   bOverWrite  - boolean for overwriting configReportPDF.json with the default one (OPTIONAL, default == true)
%
% OUTPUT: 
%   x           - x structure containing fields with all information as well as now the quality parameters loaded in
% OUTPUTFILE:
%   xASL_Report_SubjectName.pdf - printed PDF rapport containing QC images & values
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function generates a PDF report based on a predefined Json configuration file.
%               Using the function xASL_qc_ParsePdfConfig, the Json configuration file is parsed and all the pages are created.
%               The function xASL_qc_ParsePdfConfig is called recursively to parse all the Json objects in the configuration file.
%               Quality control values are obtained from the x structure using the function xASL_adm_LoadX.
%               Scans are generated using the function xASL_vis_CreateVisualFig.
% 
% EXAMPLE: xASL_qc_GenerateReport(x);
%          xASL_qc_GenerateReport(x, 'sub-001', true);
%          xASL_qc_GenerateReport(x, [], false);
% __________________________________
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

%% 0. Admin
if ~usejava('jvm') % only if JVM loaded
    fprintf('Warning: skipping PDF report, JVM missing\n');
    return;
end

% check input
if nargin < 1 || isempty(x)
    error('No x structure provided');
end

% check input
subjectOld = [];
if nargin < 2 || isempty(subject)
    if isfield(x, 'SUBJECT')
        subject = x.SUBJECT;
    else 
        warning('subject or x.SUBJECT missing, this might go wrong');
    end
else
    if ~ischar(subject)
        error('Subject input should be a string (character array)');
    end
    if isfield(x, 'SUBJECT')
        subjectOld = x.SUBJECT;
    end
    x.SUBJECT = subject;
end

if nargin < 3 || isempty(modules)
    modules = []; % This will use all modules in x.Output, in xASL_qc_GeneratePdfConfig>xASL_sub_createDefaultJson
end

if nargin < 4 || isempty(bOverWrite)
    bOverWrite = true; % We need to overwrite because different subjects can have different modules
end

% Fix <SESSION> not existing
bRemoveSESSIONfield = false;
if ~isfield(x, 'SESSION')
    if isfield(x, 'SESSIONS') && ~isempty(x.SESSIONS)
        bRemoveSESSIONfield = true;
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
    warning([PathX ' didnt exist, skipping xASL_qc_GenerateReport']);
    return;
end

x = xASL_adm_LoadX(x, PathX, true); % Assume memory x is newer than x.mat

% Make sure that the directory exists
PrintDir = fullfile(x.dir.xASLDerivatives, subject);
xASL_adm_CreateDir(PrintDir);


%% 1. Load Pdf configuration
config = xASL_qc_LoadPdfConfig(x, [], bOverWrite, [], modules);


%% 2. Create the PDF
% Print the title
fprintf('\nPrinting ExploreASL PDF report(s) in:   \n');
fprintf([PrintDir '\n']);

% Parse the entire Json Stack automatically making all the pages.
xASL_qc_ParsePdfConfig(config, x);

% Householding
if ~isempty(subjectOld)
    % restore x.SUBJECT
    x.SUBJECT = subjectOld;
end
if bRemoveSESSIONfield
    x = rmfield(x, 'SESSION');
end


%% 3. Delete the config file
configPath = fullfile(x.dir.xASLDerivatives, x.SUBJECT, 'configReportPDF.json');
xASL_delete(configPath);


end