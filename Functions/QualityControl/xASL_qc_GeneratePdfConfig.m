function [config] = xASL_qc_GeneratePdfConfig(x, subject, bOverwrite, modules)
% xASL_qc_GeneratePdfConfig Generates a JSON based on all values in x.Output
%
% FORMAT: xASL_qc_GeneratePdfConfig(x[, subject, bOverwrite, modules])
%
% INPUT:
%   x           - structure containing ExploreASL fields with all information required to run this function (REQUIRED)
%   subject     - subject name (OPTIONAL, default = x.SUBJECT)
%   bOverWrite  - boolean to determine if current configReportPDF.json should be overwritten. (OPTIONAL, default == true)
%   modules     - structure with name(s) of the modules that need to be added to the PDF report (OPTIONAL, default==all modules in x.output)
%
% OUTPUT: 
%   In the ExploreASL Derivatives folder in the subject directory, a configReportPDF.json is generated containing all quality parameters.
%   Quality parameters are taken from x.Output and x.OutputIm
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function dumps all keys from x.Output of a subject into a configuration file
%               This configuration can be read by xASL_qc_GeneratePdfConfig to create a pdf that contains all quality values.
% 
% EXAMPLE: xASL_qc_GeneratePdfConfig(x, 'sub-001');
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% __________________________________


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
   bOverwrite = true; % We need to always overwrite, because the numbers of modules can differ between subjects
end

if nargin < 4
   modules = []; % By default we use all modules parsed in x.Output, as defined below in xASL_sub_createDefaultJson
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
    warning([PathX ' didnt exist, skipping xASL_qc_GeneratePdfConfig']);
    return;
end

x = xASL_adm_LoadX(x, PathX, true); % Assume memory x is newer than x.mat

% Make sure that the directory exists
PrintDir = fullfile(x.dir.xASLDerivatives, x.SUBJECT);
xASL_adm_CreateDir(PrintDir);
PrintFile = fullfile(PrintDir, 'configReportPDF.json');

% Print the title
% -> we skip this for now, because we need to reprint the configReport for every subject, 
% because the modules can differ between subjects
% fprintf('Creating default PDF configuration file:   \n');
% fprintf([PrintFile '\n']);

% Parse the entire Json Stack automatically making all the pages.
config = xASL_sub_createDefaultJson(x, modules);

% Write Pdf configuration
xASL_io_WriteJson(PrintFile, config, bOverwrite);

end


%% ================================================================================
%% ================================================================================
function [config] = xASL_sub_createDefaultJson(x, modules)
%   x           - structure containing ExploreASL fields with all information required to run this function (REQUIRED)
%   modules     - structure with name(s) of the modules that need to be added to the PDF report (OPTIONAL, default==all modules in x.output)

    config = struct();
    if ~isfield(x, 'Output') || isempty(x.Output) || isempty(fields(x.Output))
        warning('x.Output didnt exist, skipping xASL_qc_GeneratePdfConfig');
        return;
    end
    if ~isfield(x, 'Output_im') || isempty(x.Output_im) || isempty(fields(x.Output_im))
        warning('x.Output_im didnt exist, skipping xASL_qc_GeneratePdfConfig');
        return;
    end

    config.modules = struct();

    if nargin<2 || isempty(modules)
        modules = fieldnames(x.Output);
    end
    
    for module = 1:length(modules)
        if strcmpi(modules{module}, 'Structural') % for this module, we assume a single session/run
            config.modules(module).category = 'metadata';
            config.modules(module).type = 'page';
            config.modules(module).identifier = modules{module};
            config.modules(module).content = xASL_sub_createPageContent(x.Output.(modules{module}), modules{module}, x.SUBJECT);

        elseif strcmpi(modules{module}, 'Population') % PDF reports are created per subject
            error('We cannot create a PDF report for the population module, skipping...');
        elseif strcmpi(modules{module}, 'import') % PDF reports are created per subject
            error('We cannot create a PDF report for the import module, skipping...');

        elseif strcmpi(modules{module}, 'ASL') || strcmpi(modules{module}, 'fMRI') || strcmpi(modules{module}, 'DTI')
            % for all other modules, such as ASL, fMRI, DTI, we allow multiple sessions/runs
            config.modules(module).category = 'metadata';
            config.modules(module).type = 'module'; % here we define a module instead of a page identifier
            config.modules(module).identifier = modules{module};

            if isfield(x, 'SESSION') % if we can define a single session/run
                allSessions = {x.SESSION};
            else % otherwise, use all sessions/runs
                allSessions = fieldnames(x.Output.(modules{module}));
            end            

            for iSession = 1:length(allSessions)
                config.modules(module).content(iSession).category = 'metadata';
                config.modules(module).content(iSession).type = 'page'; % now we define pages
                config.modules(module).content(iSession).identifier = allSessions{iSession};
                config.modules(module).content(iSession).content = xASL_sub_createPageContent(x.Output.(modules{module}).(allSessions{iSession}), modules{module}, x.SUBJECT, allSessions{iSession});
            end
        else
            error('Unknown module name');
        end
    end

end


%% ============================================================================
%% ============================================================================
function content = xASL_sub_createPageContent(module, modulename, subjectname, sessionname)
    
    if nargin < 4 || isempty(sessionname)
        sessionname = '';
    end

    qc_parameters = sort(fieldnames(module));
    content = cell(size(qc_parameters,1),1);

    % Add the qc_images
    qc_content = struct();
    qc_content.category = 'content';
    qc_content.type = 'text';
    
    subjectname = strrep(subjectname, 'sub-', ''); % don't show the prefix
    if strcmpi(modulename, 'structural') % don't show the sessionname for the structural module
        qc_content.text = ['Subject: ' subjectname ' module: ' modulename];
    else
        qc_content.text = ['Subject: ' subjectname ' module: ' sessionname];
    end
    
    qc_content.textSettings = struct();
    qc_content.textSettings.fontSize = '12';
    qc_content.textSettings.fontWeight = 'bold';

    content{1} = qc_content;

    for field = 1:length(qc_parameters)
        qc_content = struct();
        qc_content.category = 'content';
        qc_content.type = 'QCValues';
        qc_content.parameter = qc_parameters{field};
        qc_content.module = modulename;
        qc_content.session = sessionname;
        content{field+1} = qc_content;
    end

    % Add the qc_images
    qc_content = struct();
    qc_content.category = 'content';
    qc_content.type = 'QCValues';
    qc_content.parameter = 'qc_images';
    qc_content.position = '[0.5 0.42]'; % [0.4 0.25]
    qc_content.size = '[0.5 0.5]'; % [0.6 0.6]
    qc_content.module = modulename;
    qc_content.session = sessionname;
    content{end+1} = qc_content;

end