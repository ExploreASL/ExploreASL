function [x] = xASL_qc_GenerateReport(x, subject)
% xASL_qc_GenerateReport print output QC
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
% DESCRIPTION: This function iterates over all values in x.Output and all
%              images in x.Output_im, and prints them in a PDF file.
%              x.Output & x.Output_im should contain the QC/result output
%              of all ExploreASL pipeline steps.
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
if nargin < 2 || isempty(subject)
    subject = x.SUBJECT;
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

end
