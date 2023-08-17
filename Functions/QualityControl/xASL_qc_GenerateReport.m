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

% Set defaults
settings = xASL_sub_defaultSettings();

% Load Pdf configuration
config = xASL_adm_LoadPDFConfig(x);

% Print the title
fprintf('Printing ExploreASL PDF report:   \n');

% Parse the entire Json Stack automatically making all the pages.
xASL_sub_parseJson(config, x, [], [], settings);

end


%%================================================================================================================================================
%%================================================================================================================================================
function [settings] = xASL_sub_parseJson(json, x, figure, line, settings)
    fields = fieldnames(json);
    
    for iField = 1:length(fields)
        currentField = json.(fields{iField});

        if strcmp(fields{iField}, 'content') || strcmp(fields{iField}, 'pages')
            for iContent = 1:length(currentField)
                if iscell(currentField(iContent))
                    [settings, line] = xASL_sub_parseContent(currentField{iContent}, x, figure, line, settings);
                else
                    [settings, line] = xASL_sub_parseContent(currentField(iContent), x, figure, line, settings);
                end
            end
        else
            [settings, line] = xASL_sub_parseContent(currentField, x, figure, line, settings);
        end 
    end
end


function [settings, line] = xASL_sub_parseContent(currentField, x, figure, line, settings)

    if ~isstruct(currentField)
%        fprintf('xASL_sub_parseContent couldnt find struct of printable content \n');
        return;
    end

    if ~isfield(currentField, 'type')
        error([currentField, 'contains no field specifying the type of content, add field "type": "text" or "type": "image" to this json field for example']);
    end

    switch currentField.type
        case 'text' 
            line = xASL_sub_PrintText(currentField, figure, line, settings);
        case 'qc'
            line = xASL_sub_PrintQC(currentField, x, figure, line, settings);
        case 'image'
            settings = xASL_sub_PrintImage(currentField, x, figure, settings);  
        case 'scan'
            settings = xASL_sub_PrintScan(currentField, x, figure, settings);
        case 'page'
            xASL_sub_printPage(currentField, x, settings);  
        case 'block'
            xASL_sub_printBlock(currentField, x, figure, settings);
        case 'settings'
            settings = xASL_sub_loadSettings(currentField, settings);  
        case 'patients' 
            line = xASL_sub_PrintPatient(x, figure, line, settings);
    end  
end

function  xASL_sub_printPage(json, x, settings)

    clear primaryFig ax  
    primaryFig = figure('visible', 'off', 'Units', 'centimeters', 'Position', [0 0 21 29.7]);

    %% Print the Title
    logoPath = fullfile(x.opts.MyPath, 'Design', 'ExploreASL_logoHeader.png');
    xASL_sub_PrintImage(logoPath, [], primaryFig, settings, [0 0.96 1 0.04]);

    %% Print the Footer
    xASL_sub_PrintText('This report was automatically generated by ExploreASL', primaryFig, [0 0.02 1 0], settings);

    %% Parse JSON and create the file
    xASL_sub_parseJson(json, x, primaryFig, [0 0.93 1 0], settings);

    PrintFile = ['xASL_Report_', json.pageIdentifier];
    PrintPath = fullfile(x.dir.xASLDerivatives, x.SUBJECT, PrintFile); 

    print(primaryFig, PrintPath, '-dpdf', '-bestfit');

    clear primaryFig ax
end

function xASL_sub_printBlock(input, x, pageFig, settings)
    position = xASL_str2num(input.position);
    size = xASL_str2num(input.size);
    [settings.canvas, line] = xASL_sub_createNewCanvas(position, size, settings.canvas);
    xASL_sub_parseJson(input, x, pageFig, line, settings);
end

function line = xASL_sub_PrintPatient(x, figure, line, settings)
    %Check if participant data exists.
    ParticipantsTSV = xASL_adm_GetFileList(x.D.ROOT, 'participants.tsv');
    PatientInfo = {};

    if isempty(ParticipantsTSV)
        return
    end

    structParticipants = xASL_tsvRead(ParticipantsTSV{1});
    nParticipants = size(structParticipants, 1);
    for iPar = 2:nParticipants
        if structParticipants{iPar, 1} == x.SUBJECT
            %Extract Participants table, and participants info
            PatientInfo = [structParticipants(1, :); structParticipants(iPar, :)];
        end
    end

    if isempty(PatientInfo)
        fprintf ('no patient info found to be printed \n');
        return
    end
    
    line = xASL_sub_PrintText('Participant Information', figure, line, settings );
    for iEntry = 1:size(PatientInfo, 2)
        if strcmp(PatientInfo{1, iEntry}, 'participant_id')
            line = xASL_sub_PrintText(['Participant: ', PatientInfo{2, iEntry}], figure, line, settings );
        elseif strcmp(PatientInfo{1, iEntry}, 'age')
            line = xASL_sub_PrintText(['Participant age: ', num2str(PatientInfo{2, iEntry})], figure, line, settings );
        elseif strcmp(PatientInfo{1, iEntry}, 'sex')
            line = xASL_sub_PrintText(['Participant sex: ', PatientInfo{2, iEntry}], figure, line, settings );
        end
    end

end

function line = xASL_sub_PrintQC(json, x, figure, line, settings)
    % Stop if module and field don't exists in output
    if ~isfield(x.Output, (json.module)) && ~isfield(x.Output.(json.module), json.parameter)
        return     
    end

    % Use field specific settings if they exists, otherwise default to monospace for QC values. 
    settings.fontName = 'monospace';
    if isfield(json, 'settings')
        settings = xASL_sub_loadSettings(json.settings, settings);
    end

    TempValue = x.Output.(json.module).(json.parameter);
    
    if ~isfield(json, 'alias') 
        json.alias = json.parameter;
    end

    if ~isfield(json, 'unit') 
        json.unit = '';
    end

    if isfield(json, 'range') 
        [range] = strsplit(json.range, '-');
        if (TempValue < xASL_str2num(range{1})) || (TempValue > xASL_str2num(range{2}))
            settings.color = 'r';
        end
        json.range = ['(' json.range ')'];
    else
        json.range = '';
    end

    if isnumeric(TempValue)
        TempValue = xASL_num2str(TempValue);
    end

    if size(TempValue, 1) == 1
        TextString = sprintf([sprintf('%-20s', [json.alias, ':']), sprintf('%8s', TempValue), sprintf('%-12s', [json.unit, ' ' , json.range]), ' \n']);
    end

    line = xASL_sub_PrintText(TextString, figure, line, settings);

end

function [settings] = xASL_sub_PrintImage(input, x, figure, settings, position)
    switch nargin
    case 4
        position = [xASL_str2num(input.position) xASL_str2num(input.size)];

        if isfield(input, 'xPath') && isfield(x.P, input.name)
            ImagePath = x.P.(input.name);
        elseif isfield(input, 'absolutePath')
            ImagePath = input.absolutePath;
            ImagePath = xASL_sub_WildcardReplace(ImagePath, x);
        elseif isfield(input, 'popPath')
            ImagePath = fullfile(x.dir.xASLDerivatives, 'Population', input.popPath);
            ImagePath = xASL_sub_WildcardReplace(ImagePath, x);
        elseif isfield(input, 'subjPath')
            ImagePath = fullfile(x.dir.xASLDerivatives, x.SUBJECT, input.subjPath);
        else
            warning('xASL_sub_PrintImage didnt have a defined path') ;
            ImagePath = xASL_sub_WildcardReplace(ImagePath, x);
            return
        end

        if isfield(input, 'header')
            header = input.header;
        else
            header = '';
        end
    case 5
        ImagePath = input;
        header = '';
    end

    [canvas] = xASL_sub_createNewCanvas(position(1:2), position(3:4), settings.canvas);
    ax = axes('Position', canvas, 'Visible', settings.axesVisible, 'Parent', figure);
    [img, ~, alphachannel] = imread(ImagePath);
    fg = imshow(img);

    if strcmp('png', ImagePath(end-2:end))
        fg.AlphaData = alphachannel;
    end

    clear fg
    settings.figureCount = xASL_sub_PrintHeader(header, figure, settings, canvas);
end

function [strout] = xASL_sub_WildcardReplace(strin, x)
    strout = strin;
    substring = regexp(strin, '<\w*>', 'match');
    for substringIndex = 1:length(substring)
        if ~isfield(x, substring{substringIndex}(2:end-1))
            warning(['Could not replace ', substring{substringIndex}, ' check if file exists in ExploreASL/Derivatives/Population']);
        else
            strout = strrep(strout, substring{substringIndex}, x.(substring{substringIndex}(2:end-1)));
        end
    end
end

function [settings] = xASL_sub_PrintScan(input, x, figure, settings)
    position = [xASL_str2num(input.position) xASL_str2num(input.size)];
    [canvas] = xASL_sub_createNewCanvas(position(1:2), position(3:4), settings.canvas);

    if ~isfield(x.P, input.name)
        warning (['could not print', input.name, 'check if nifti exists in ExploreASL/Derivatives/Population']);
        return
    end

    ImIn = {x.P.(input.name)};
    header = input.name;
    ax = axes('Position', canvas, 'Visible', settings.axesVisible, 'Parent', figure);
    
    if isfield(input, 'overlay')
        fields = fieldnames(input.overlay);
        for iField = 1:length(fields)
            if isfield(x.P, fields(iField))
                ImIn(end+1) = {x.P.(fields{iField})};
                header = [header ' + ' fields{iField}];
            end
        end
    end

    if isfield(input, 'slice')
        [x.S.TraSlices] = xASL_sub_getSliceFromStruct(input.slice, 'TraSlice');
        [x.S.CorSlices] = xASL_sub_getSliceFromStruct(input.slice, 'CorSlice');
        [x.S.SagSlices] = xASL_sub_getSliceFromStruct(input.slice, 'SagSlice');       
    else
        x.S.TraSlices = [25, 50, 90];
        x.S.SagSlices = [25, 50, 90];
        x.S.CorSlices = [25, 50, 90];
    end

    img = xASL_vis_CreateVisualFig(x, ImIn, [], [], [], [], [], [], [], [], [], [], []);
    fg = imshow(img);
    clear fg
    settings.figureCount = xASL_sub_PrintHeader(header, figure, settings, canvas);
end

function [slice] = xASL_sub_getSliceFromStruct(struct, name)

    if isfield(struct, name) && ~isempty(struct.(name))
        slice = str2num(struct.(name));
    else 
        slice = [];
    end
end

function [newCanvas, line] = xASL_sub_createNewCanvas(position, size, oldCanvas)
    newPosition = [oldCanvas(1) + position(1) * oldCanvas(3) oldCanvas(2) + position(2) * oldCanvas(4)];
    newSize = [oldCanvas(3:4) .* size];
    newCanvas = [newPosition, newSize];
    line = [newPosition(1)  newPosition(2) + newSize(2) newSize(1) 0];
end

function line = xASL_sub_PrintText(input, figure, line, settings)

    switch class(input)
        case {'string', 'char'} 
            String = input;
        case 'struct'
            String = input.text;
            if isfield(input, 'settings')
                settings = xASL_sub_loadSettings(input.settings, settings);
            end
        otherwise
            class(input)
            error('xASL_sub_PrintText couldnt find string of printable text')
    end

    ax = axes('Position', line , 'Visible', settings.axesVisible, 'Parent', figure);
    text(0, 0, String, 'Parent', ax, 'FontSize', settings.fontSize, 'FontWeight', settings.fontWeight, 'Color', settings.color, 'FontName', settings.fontName, 'Interpreter', 'none', 'VerticalAlignment', 'top');
    line = xASL_sub_NewLine(line, settings);
end

function line = xASL_sub_NewLine(line, settings)
    line(2) = line(2) - (settings.lineSpacing) - (settings.fontSize * 0.001);

    if line(2) < 0 
        warning('No space left on page!');
    elseif line(2) < settings.canvas(2) 
        warning('Printing outside canvas, check block settings.');
    end
end

function [figureCount] = xASL_sub_PrintHeader(header, figure, settings, position)

    if isempty(header)
        figureCount = settings.figureCount;
        return
    else
        figureCount = settings.figureCount + 1;
    end

    position(4) = 0;
    settings.HorizontalAlignment = 'center';
    text = ['Fig ' num2str(figureCount) ': ' header];
    xASL_sub_PrintText(text, figure, position, settings);
end

function [settings] = xASL_sub_loadSettings(json, settings)

    fields = fieldnames(json);
    
    for iField = 1:length(fields)
        if strcmp(fields{iField}, 'fontSize') || strcmp(fields{iField}, 'lineSpacing')
            settings.(fields{iField}) = xASL_str2num(json.(fields{iField}));
        else
            settings.(fields{iField}) = json.(fields{iField});
        end 
    end
end

function [settings] = xASL_sub_defaultSettings()
    settings.color = 'k';
    settings.HorizontalAlignment = 'left';
    settings.fontWeight = 'normal';
    settings.axesVisible = 'off';
    settings.fontName = 'default';
    settings.lineSpacing = 0.005;
    settings.figureCount = 0;
    settings.canvas = [0 0 1 1];

    if ispc
        settings.fontSize = 10;
    elseif isunix || ismac
        settings.fontSize = 8;
    end 

end