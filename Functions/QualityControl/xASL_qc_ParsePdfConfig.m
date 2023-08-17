function [pdfSettings] = xASL_qc_ParsePdfConfig(layoutStructure, x, currentFigure, line, pdfSettings)
% xASL_qc_ParsePdfConfig function used by xASL_qc_CreatePDF to parse the configuration file loaded by xASL_adm_LoadPdfConfig.
%   This function will run recursively, and will print all text, images, scans and other content specified in the json file.
%   The json file should contain a structure with fields specified in the manual.
%   
% FORMAT: xASL_qc_ParsePdfConfig(json, x, currentFigure, line, pdfSettings)
%
% INPUT:
%   layoutStructure - json structure containing all information to be printed (REQUIRED)
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   currentFigure      - currentFigure handle to print to (OPTIONAL, defaults to matlab main figure)
%   line        - line to print to (OPTIONAL, defaults to [0 0.93 1 0] when it generates a new page)
%   pdfSettings    - settings used to print (OPTIONAL, will set default settings if not specified)
%
% OUTPUT: 
%   pdfSettings    - settings to print with
%
% OUTPUTFILE:
%   xASL_Report_SubjectName.pdf - printed PDF rapport containing QC images & values
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function iterates over all values in x.Output and all
%              images in x.Output_im, and prints them in a PDF file.
%              x.Output & x.Output_im should contain the QC/result output
%              of all ExploreASL pipeline steps.
% 
% EXAMPLE: xASL_qc_ParsePdfConfig(x);
% __________________________________
% Copyright (C) 2015-2023 ExploreASL

if nargin < 2 
    error('xASL_qc_ParsePdfConfig requires at least 2 input arguments');
end

if nargin < 3 || isempty(currentFigure)
    currentFigure = gcf;
end

if nargin < 4 || isempty(line)
    line = [0 0.93 1 1];
end

if nargin < 5 || isempty(pdfSettings)
    pdfSettings = xASL_qc_ParsePdfConfig_sub_defaultSettings();
end

fields = fieldnames(layoutStructure);

for iField = 1:length(fields)
    currentField = layoutStructure.(fields{iField});

    if strcmp(fields{iField}, 'content') || strcmp(fields{iField}, 'pages')
        for iContent = 1:length(currentField)
            if iscell(currentField(iContent))
                [pdfSettings, line] = xASL_qc_ParsePdfConfig_sub_parseContent(currentField{iContent}, x, currentFigure, line, pdfSettings);
            else
                [pdfSettings, line] = xASL_qc_ParsePdfConfig_sub_parseContent(currentField(iContent), x, currentFigure, line, pdfSettings);
            end
        end
    else
        [pdfSettings, line] = xASL_qc_ParsePdfConfig_sub_parseContent(currentField, x, currentFigure, line, pdfSettings);
    end 
end

end

% ====================================================================================================================================================
%% Subfunctions

function [pdfSettings, line] = xASL_qc_ParsePdfConfig_sub_parseContent(currentField, x, currentFigure, line, pdfSettings)

    if ~isstruct(currentField)
%        fprintf('xASL_qc_ParsePdfConfig_sub_parseContent couldnt find struct of printable content \n');
        return;
    end

    if ~isfield(currentField, 'type')
        error([currentField, 'contains no field specifying the type of content, add field "type": "text" or "type": "image" to this json field for example']);
    end

    switch currentField.type
        case 'text' 
            line = xASL_qc_ParsePdfConfig_sub_PrintText(currentField, currentFigure, line, pdfSettings);
        case 'qc'
            line = xASL_qc_ParsePdfConfig_sub_PrintQC(currentField, x, currentFigure, line, pdfSettings);
        case 'image'
            pdfSettings = xASL_qc_ParsePdfConfig_sub_PrintImage(currentField, x, currentFigure, pdfSettings);  
        case 'scan'
            pdfSettings = xASL_qc_ParsePdfConfig_sub_PrintScan(currentField, x, currentFigure, pdfSettings);
        case 'page'
            xASL_qc_ParsePdfConfig_sub_printPage(currentField, x, pdfSettings);  
        case 'block'
            xASL_qc_ParsePdfConfig_sub_printBlock(currentField, x, currentFigure, pdfSettings);
        case 'settings'
            pdfSettings = xASL_qc_ParsePdfConfig_sub_loadSettings(currentField, pdfSettings);  
        case 'patients' 
            line = xASL_qc_ParsePdfConfig_sub_PrintPatient(x, currentFigure, line, pdfSettings);
    end  
end

function  xASL_qc_ParsePdfConfig_sub_printPage(json, x, pdfSettings)

    %% Create the figure defaults
    primaryFig = figure('visible', 'off', 'Units', 'centimeters', 'Position', [0 0 21 29.7]);
    ax = axes('Position', [0 0 1 1], 'Visible', 'off', 'Parent', primaryFig);

    %% Print the Title
    logoPath = fullfile(x.opts.MyPath, 'Design', 'ExploreASL_logoHeader.png');
    xASL_qc_ParsePdfConfig_sub_PrintImage(logoPath, [], primaryFig, pdfSettings, [0 0.96 1 0.04]);

    %% Print the Footer
    xASL_qc_ParsePdfConfig_sub_PrintText('This report was automatically generated by ExploreASL', primaryFig, [0 0.02 1 0], pdfSettings);

    %% Parse JSON and create the file
    xASL_qc_ParsePdfConfig(json, x, primaryFig, [0 0.93 1 0], pdfSettings);

    PrintFile = ['xASL_Report_', json.pageIdentifier];
    PrintPath = fullfile(x.dir.xASLDerivatives, x.SUBJECT, PrintFile); 

    print(primaryFig, PrintPath, '-dpdf', '-bestfit');

end

function xASL_qc_ParsePdfConfig_sub_printBlock(input, x, pageFig, pdfSettings)
    position = xASL_str2num(input.position);
    size = xASL_str2num(input.size);
    [pdfSettings.canvas, line] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position, size, pdfSettings.canvas);
    xASL_qc_ParsePdfConfig(input, x, pageFig, line, pdfSettings);
end

function line = xASL_qc_ParsePdfConfig_sub_PrintPatient(x, currentFigure, line, pdfSettings)
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
    
    line = xASL_qc_ParsePdfConfig_sub_PrintText('Participant Information', currentFigure, line, pdfSettings );
    for iEntry = 1:size(PatientInfo, 2)
        if strcmp(PatientInfo{1, iEntry}, 'participant_id')
            line = xASL_qc_ParsePdfConfig_sub_PrintText(['Participant: ', PatientInfo{2, iEntry}], currentFigure, line, pdfSettings );
        elseif strcmp(PatientInfo{1, iEntry}, 'age')
            line = xASL_qc_ParsePdfConfig_sub_PrintText(['Participant age: ', num2str(PatientInfo{2, iEntry})], currentFigure, line, pdfSettings );
        elseif strcmp(PatientInfo{1, iEntry}, 'sex')
            line = xASL_qc_ParsePdfConfig_sub_PrintText(['Participant sex: ', PatientInfo{2, iEntry}], currentFigure, line, pdfSettings );
        end
    end

end

function line = xASL_qc_ParsePdfConfig_sub_PrintQC(json, x, currentFigure, line, pdfSettings)
    % Stop if module and field don't exists in output
    if ~isfield(x.Output, (json.module)) && ~isfield(x.Output.(json.module), json.parameter)
        return     
    end

    % Use field specific settings if they exists, otherwise default to monospace for QC values. 
    pdfSettings.fontName = 'monospace';
    if isfield(json, 'settings')
        pdfSettings = xASL_qc_ParsePdfConfig_sub_loadSettings(json.settings, pdfSettings);
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
            pdfSettings.color = 'r';
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

    line = xASL_qc_ParsePdfConfig_sub_PrintText(TextString, currentFigure, line, pdfSettings);

end

function [pdfSettings] = xASL_qc_ParsePdfConfig_sub_PrintImage(input, x, currentFigure, pdfSettings, position)
    switch nargin
    case 4
        position = [xASL_str2num(input.position) xASL_str2num(input.size)];

        if isfield(input, 'xPath') && isfield(x.P, input.name)
            ImagePath = x.P.(input.name);
        elseif isfield(input, 'absolutePath')
            ImagePath = input.absolutePath;
            ImagePath = xASL_qc_ParsePdfConfig_sub_WildcardReplace(ImagePath, x);
        elseif isfield(input, 'popPath')
            ImagePath = fullfile(x.dir.xASLDerivatives, 'Population', input.popPath);
            ImagePath = xASL_qc_ParsePdfConfig_sub_WildcardReplace(ImagePath, x);
        elseif isfield(input, 'subjPath')
            ImagePath = fullfile(x.dir.xASLDerivatives, x.SUBJECT, input.subjPath);
        else
            warning('xASL_qc_ParsePdfConfig_sub_PrintImage didnt have a defined path') ;
            ImagePath = xASL_qc_ParsePdfConfig_sub_WildcardReplace(ImagePath, x);
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

    [canvas] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position(1:2), position(3:4), pdfSettings.canvas);

    cla;
    ax = axes('Position', canvas, 'Visible', pdfSettings.axesVisible, 'Parent', currentFigure);
    [img, ~, alphachannel] = imread(ImagePath);
    fg = imshow(img);

    if strcmp('png', ImagePath(end-2:end))
        fg.AlphaData = alphachannel;
    end

    pdfSettings.figureCount = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, pdfSettings, canvas);
end

function [strout] = xASL_qc_ParsePdfConfig_sub_WildcardReplace(strin, x)
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

function [pdfSettings] = xASL_qc_ParsePdfConfig_sub_PrintScan(input, x, currentFigure, pdfSettings)
    position = [xASL_str2num(input.position) xASL_str2num(input.size)];
    [canvas] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position(1:2), position(3:4), pdfSettings.canvas);

    if ~isfield(x.P, input.name)
        warning (['could not print', input.name, 'check if nifti exists in ExploreASL/Derivatives/Population']);
        return
    end

    ImIn = {x.P.(input.name)};
    header = input.name;
    
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
        [x.S.TraSlices] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(input.slice, 'TraSlice');
        [x.S.CorSlices] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(input.slice, 'CorSlice');
        [x.S.SagSlices] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(input.slice, 'SagSlice');       
    else
        x.S.TraSlices = [25, 50, 90];
        x.S.SagSlices = [25, 50, 90];
        x.S.CorSlices = [25, 50, 90];
    end

    % Create the image from the defined scans
    img = xASL_vis_CreateVisualFig(x, ImIn);   

    % Delete all previously defined graphics objects, and create new axes and image
    cla;
    ax = axes('Position', canvas, 'Visible', pdfSettings.axesVisible, 'Parent', currentFigure);
    fg = imshow(img);

    % Print the image header and update the figure count
    pdfSettings.figureCount = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, pdfSettings, canvas);
end

function [slice] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(struct, name)

    if isfield(struct, name) && ~isempty(struct.(name))
        slice = str2num(struct.(name));
    else 
        slice = [];
    end
end

function [newCanvas, line] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position, size, oldCanvas)
    newPosition = [oldCanvas(1) + position(1) * oldCanvas(3) oldCanvas(2) + position(2) * oldCanvas(4)];
    newSize = [oldCanvas(3:4) .* size];
    newCanvas = [newPosition, newSize];
    line = [newPosition(1)  newPosition(2) + newSize(2) newSize(1) 0];
end

function line = xASL_qc_ParsePdfConfig_sub_PrintText(input, currentFigure, line, pdfSettings)

    switch class(input)
        case {'string', 'char'} 
            String = input;
        case 'struct'
            String = input.text;
            if isfield(input, 'settings')
                pdfSettings = xASL_qc_ParsePdfConfig_sub_loadSettings(input.settings, pdfSettings);
            end
        otherwise
            class(input)
            error('xASL_qc_ParsePdfConfig_sub_PrintText couldnt find string of printable text')
    end

    ax = axes('Position', line , 'Visible', pdfSettings.axesVisible, 'Parent', currentFigure);
    text(0, 0, String, 'Parent', ax, 'FontSize', pdfSettings.fontSize, 'FontWeight', pdfSettings.fontWeight, 'Color', pdfSettings.color, 'FontName', pdfSettings.fontName, 'Interpreter', 'none', 'VerticalAlignment', 'top');
    line = xASL_qc_ParsePdfConfig_sub_NewLine(line, pdfSettings);
end

function line = xASL_qc_ParsePdfConfig_sub_NewLine(line, pdfSettings)
    line(2) = line(2) - (pdfSettings.lineSpacing) - (pdfSettings.fontSize * 0.001);

    if line(2) < 0 
        warning('No space left on page!');
    elseif line(2) < pdfSettings.canvas(2) 
        warning('Printing outside canvas, check block settings.');
    end
end

function [figureCount] = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, pdfSettings, position)

    if isempty(header)
        figureCount = pdfSettings.figureCount;
        return
    else
        figureCount = pdfSettings.figureCount + 1;
    end

    position(4) = 0;
    pdfSettings.HorizontalAlignment = 'center';
    text = ['Fig ' num2str(figureCount) ': ' header];
    xASL_qc_ParsePdfConfig_sub_PrintText(text, currentFigure, position, pdfSettings);
end

function [pdfSettings] = xASL_qc_ParsePdfConfig_sub_loadSettings(json, pdfSettings)

    fields = fieldnames(json);
    
    for iField = 1:length(fields)
        if strcmp(fields{iField}, 'fontSize') || strcmp(fields{iField}, 'lineSpacing')
            pdfSettings.(fields{iField}) = xASL_str2num(json.(fields{iField}));
        else
            pdfSettings.(fields{iField}) = json.(fields{iField});
        end 
    end
end


function [pdfSettings] = xASL_qc_ParsePdfConfig_sub_defaultSettings()
    % This function sets the default settings for the PDF report.
    pdfSettings.color = 'k';
    pdfSettings.HorizontalAlignment = 'left';
    pdfSettings.fontWeight = 'normal';
    pdfSettings.axesVisible = 'off';
    pdfSettings.fontName = 'default';
    pdfSettings.lineSpacing = 0.005;
    pdfSettings.figureCount = 0;
    pdfSettings.canvas = [0 0 1 1];

    if ispc
        pdfSettings.fontSize = 10;
    elseif isunix || ismac
        pdfSettings.fontSize = 8;
    end 

end