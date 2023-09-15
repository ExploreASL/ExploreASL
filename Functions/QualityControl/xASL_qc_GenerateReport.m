function [x] = xASL_qc_GenerateReport(x, subject)
% xASL_qc_GenerateReport Generates a PDF report based on a predifined Json configuration file
%
% FORMAT: xASL_qc_GenerateReport(x)
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   subject     - subject name (OPTIONAL, default = x.SUBJECT)
%
% OUTPUT: 
%
%   x           - structure containing fields with all information required to run this submodule
% OUTPUTFILE:
%   xASL_Report_SubjectName.pdf - printed PDF rapport containing QC images & values
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function generates a PDF report based on a predifined Json configuration file.
%               Using the function xASL_qc_ParsePdfConfig, the Json configuration file is parsed and all the pages are created.
%               The function xASL_qc_ParsePdfConfig is called recursively to parse all the Json objects in the configuration file.
%               Quality control values are obtained from the x structure using the function xASL_adm_LoadX.
%               Scans are generated using the function xASL_vis_CreateVisualFig.
% 
% EXAMPLE: xASL_qc_GenerateReport(x);
% __________________________________
% Copyright (C) 2015-2023 ExploreASL

if ~usejava('jvm') % only if JVM loaded
    fprintf('Warning: skipping PDF report, JVM missing\n');
    config = NaN;
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
    if isfield(x, 'SUBJECT')
        subjectOld = x.SUBJECT;
    end
    x.SUBJECT = subject;
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
    warning([PathX ' didnt exist, skipping xASL_qc_CreateOutputPDF']);
    return;
end

x = xASL_adm_LoadX(x, PathX, true); % Assume memory x is newer than x.mat

% Make sure that the directory exists
PrintDir = fullfile(x.dir.xASLDerivatives, subject);
xASL_adm_CreateDir(PrintDir);

% Delete existing xASL report files
ExistingPrintFiles = xASL_adm_GetFileList(PrintDir, '^xASL_Report_.+$');
if ~isempty(ExistingPrintFiles)
    for iFile = 1:size(ExistingPrintFiles)
        xASL_delete(ExistingPrintFiles{iFile});
    end
end

% Load Pdf configuration
config = xASL_adm_LoadPdfConfig(x);

% Print the title
fprintf('Printing ExploreASL PDF report:   \n');

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
    

end
