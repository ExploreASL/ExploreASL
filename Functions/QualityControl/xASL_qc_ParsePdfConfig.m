function [settingsPDF] = xASL_qc_ParsePdfConfig(layoutStructure, x, currentFigure, line, settingsPDF)
% xASL_qc_ParsePdfConfig function used by xASL_qc_CreatePDF to parse the configuration file loaded by xASL_adm_LoadPdfConfig.
%   This function will run recursively, and will print all text, images, scans and other content specified in the json file.
%   The json file should contain a structure with fields specified in the manual.
%   
% FORMAT: xASL_qc_ParsePdfConfig(json, x, currentFigure, line, settingsPDF)
%
% INPUT:
%   layoutStructure - json structure containing all information to be printed (REQUIRED)
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   currentFigure      - currentFigure handle to print to (OPTIONAL, defaults to matlab main figure)
%   line        - line to print to (OPTIONAL, defaults to [0 0.93 1 0] when it generates a new page)
%   settingsPDF    - settings used to print (OPTIONAL, will set default settings if not specified)
%
% OUTPUT: 
%   settingsPDF    - settings to print with
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

if nargin < 5 || isempty(settingsPDF)
    settingsPDF = xASL_qc_ParsePdfConfig_sub_defaultSettings();
end

fields = fieldnames(layoutStructure);

for iField = 1:length(fields)
    currentField = layoutStructure.(fields{iField});

    if strcmp(fields{iField}, 'content') || strcmp(fields{iField}, 'pages')
        for iContent = 1:length(currentField)
            if iscell(currentField(iContent))
                [settingsPDF, line] = xASL_qc_ParsePdfConfig_sub_parseContent(currentField{iContent}, x, currentFigure, line, settingsPDF);
            else
                [settingsPDF, line] = xASL_qc_ParsePdfConfig_sub_parseContent(currentField(iContent), x, currentFigure, line, settingsPDF);
            end
        end
    else
        [settingsPDF, line] = xASL_qc_ParsePdfConfig_sub_parseContent(currentField, x, currentFigure, line, settingsPDF);
    end 
end

end

% ====================================================================================================================================================
% ============================================================Content Parsing Functions===============================================================
% ====================================================================================================================================================

function [settingsPDF, line] = xASL_qc_ParsePdfConfig_sub_parseContent(currentField, x, currentFigure, line, settingsPDF)

    if ~isstruct(currentField)
%        fprintf('xASL_qc_ParsePdfConfig_sub_parseContent couldnt find struct of printable content \n');
        return;
    end

    if ~isfield(currentField, 'type')
        error([currentField, 'contains no field specifying the type of content, add field "type": "text" or "type": "image" to this json field for example']);
    end

    switch currentField.type
        case 'text' 
            line = xASL_qc_ParsePdfConfig_sub_PrintText(currentField, currentFigure, line, settingsPDF);
        case 'qc'
            line = xASL_qc_ParsePdfConfig_sub_PrintQC(currentField, x, currentFigure, line, settingsPDF);
        case 'image'
            settingsPDF = xASL_qc_ParsePdfConfig_sub_PrintImage(currentField, x, currentFigure, settingsPDF);  
        case 'scan'
            settingsPDF = xASL_qc_ParsePdfConfig_sub_PrintScan(currentField, x, currentFigure, settingsPDF);
        case 'page'
            xASL_qc_ParsePdfConfig_sub_printPage(currentField, x, settingsPDF);  
        case 'block'
            xASL_qc_ParsePdfConfig_sub_printBlock(currentField, x, currentFigure, settingsPDF);
        case 'settings'
            settingsPDF = xASL_qc_ParsePdfConfig_sub_loadSettings(currentField, settingsPDF);  
        case 'patients' 
            line = xASL_qc_ParsePdfConfig_sub_PrintPatient(x, currentFigure, line, settingsPDF);
    end  
end

% ====================================================================================================================================================

function  xASL_qc_ParsePdfConfig_sub_printPage(pageStruct, x, settingsPDF)

    %% Create the figure defaults
    figPrimary = figure('visible', 'off', 'Units', 'centimeters', 'Position', [0 0 21 29.7]);
    ax = axes('Position', [0 0 1 1], 'Visible', 'off', 'Parent', figPrimary);

    %% Print the Title
    pathLogo = fullfile(x.opts.MyPath, 'Design', 'ExploreASL_logoHeader.png');
    xASL_qc_ParsePdfConfig_sub_PrintImage(pathLogo, [], figPrimary, settingsPDF, [0 0.96 1 0.04]);

    %% Print the Footer
    xASL_qc_ParsePdfConfig_sub_PrintText('This report was automatically generated by ExploreASL', figPrimary, [0 0.02 1 0], settingsPDF);

    %% Parse pageStruct and create the page as defined in the json file
    xASL_qc_ParsePdfConfig(pageStruct, x, figPrimary, [0 0.93 1 0], settingsPDF);

    fileName = ['xASL_Report_', pageStruct.pageIdentifier];
    printPathFile = fullfile(x.dir.xASLDerivatives, x.SUBJECT, fileName); 

    print(figPrimary, printPathFile, '-dpdf', '-bestfit');

end

% ====================================================================================================================================================

function xASL_qc_ParsePdfConfig_sub_printBlock(blockStruct, x, pageFig, settingsPDF)
    position = xASL_str2num(blockStruct.position);
    size = xASL_str2num(blockStruct.size);
    [settingsPDF.canvas, line] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position, size, settingsPDF.canvas);
    xASL_qc_ParsePdfConfig(blockStruct, x, pageFig, line, settingsPDF);
end

% ====================================================================================================================================================

function [newCanvas, line] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position, size, oldCanvas)
    newPosition = [oldCanvas(1) + position(1) * oldCanvas(3) oldCanvas(2) + position(2) * oldCanvas(4)];
    newSize = [oldCanvas(3:4) .* size];
    newCanvas = [newPosition, newSize];
    line = [newPosition(1)  newPosition(2) + newSize(2) newSize(1) 0];
end

% ====================================================================================================================================================
% ==============================================================Image Based Functions=================================================================
% ====================================================================================================================================================

function [settingsPDF] = xASL_qc_ParsePdfConfig_sub_PrintImage(input, x, currentFigure, settingsPDF, position)
    header = '';
    switch nargin
    case 4
        imageStruct = input;
        position = [xASL_str2num(imageStruct.position) xASL_str2num(imageStruct.size)];

        if isfield(imageStruct, 'xPath') && isfield(x.P, imageStruct.name)
            ImagePath = x.P.(imageStruct.name);
        elseif isfield(imageStruct, 'absolutePath')
            ImagePath = imageStruct.absolutePath;
            ImagePath = xASL_qc_ParsePdfConfig_sub_WildcardReplace(ImagePath, x);
        elseif isfield(imageStruct, 'popPath')
            ImagePath = fullfile(x.dir.xASLDerivatives, 'Population', imageStruct.popPath);
            ImagePath = xASL_qc_ParsePdfConfig_sub_WildcardReplace(ImagePath, x);
        elseif isfield(imageStruct, 'subjPath')
            ImagePath = fullfile(x.dir.xASLDerivatives, x.SUBJECT, imageStruct.subjPath);
        else
            warning('xASL_qc_ParsePdfConfig_sub_PrintImage didnt have a defined path') ;
            ImagePath = xASL_qc_ParsePdfConfig_sub_WildcardReplace(ImagePath, x);
            return
        end

        if isfield(imageStruct, 'header')
            header = imageStruct.header;
        end
    case 5
        ImagePath = input;
    end

    [canvas] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position(1:2), position(3:4), settingsPDF.canvas);

    cla;
    ax = axes('Position', canvas, 'Visible', settingsPDF.axesVisible, 'Parent', currentFigure);
    [img, ~, alphachannel] = imread(ImagePath);
    fg = imshow(img);

    if strcmp('png', ImagePath(end-2:end))
        fg.AlphaData = alphachannel;
    end

    settingsPDF.figureCount = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, settingsPDF, canvas);
end

% ====================================================================================================================================================

function [settingsPDF] = xASL_qc_ParsePdfConfig_sub_PrintScan(scanStruct, x, currentFigure, settingsPDF)
    position = [xASL_str2num(scanStruct.position) xASL_str2num(scanStruct.size)];
    [canvas] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position(1:2), position(3:4), settingsPDF.canvas);

    if ~isfield(x.P, scanStruct.name)
        warning (['could not print', scanStruct.name, 'check if nifti exists in ExploreASL/Derivatives/Population']);
        return
    end

    ImIn = {x.P.(scanStruct.name)};
    header = scanStruct.name;
    
    if isfield(scanStruct, 'overlay')
        fields = fieldnames(scanStruct.overlay);
        for iField = 1:length(fields)
            if isfield(x.P, fields(iField))
                ImIn(end+1) = {x.P.(fields{iField})};
                header = [header ' + ' fields{iField}];
            end
        end
    end

    if isfield(scanStruct, 'slice')
        [x.S.TraSlices] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(scanStruct.slice, 'TraSlice');
        [x.S.CorSlices] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(scanStruct.slice, 'CorSlice');
        [x.S.SagSlices] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(scanStruct.slice, 'SagSlice');       
    else
        x.S.TraSlices = [25, 50, 90];
        x.S.SagSlices = [25, 50, 90];
        x.S.CorSlices = [25, 50, 90];
    end

    % Create the image from the defined scans
    img = xASL_vis_CreateVisualFig(x, ImIn);   

    % Delete all previously defined graphics objects, and create new axes and image
    cla;
    ax = axes('Position', canvas, 'Visible', settingsPDF.axesVisible, 'Parent', currentFigure);
    fg = imshow(img);

    % Print the image header and update the figure count
    settingsPDF.figureCount = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, settingsPDF, canvas);
end

% ====================================================================================================================================================

function [slice] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(struct, name)

    if isfield(struct, name) && ~isempty(struct.(name))
        slice = str2num(struct.(name));
    else 
        slice = [];
    end
end

% ====================================================================================================================================================
% ==============================================================Text Based Functions==================================================================
% ====================================================================================================================================================

function line = xASL_qc_ParsePdfConfig_sub_PrintText(input, currentFigure, line, settingsPDF)

    switch class(input)
        case {'string', 'char'} 
            String = input;
        case 'struct'
            textStruct = input;
            String = textStruct.text;
            if isfield(textStruct, 'settings')
                settingsPDF = xASL_qc_ParsePdfConfig_sub_loadSettings(textStruct.settings, settingsPDF);
            end
        otherwise
            class(input)
            error('xASL_qc_ParsePdfConfig_sub_PrintText couldnt find string of printable text')
    end

    ax = axes('Position', line , 'Visible', settingsPDF.axesVisible, 'Parent', currentFigure);
    text(0, 0, String, 'Parent', ax, 'FontSize', settingsPDF.fontSize, 'FontWeight', settingsPDF.fontWeight, 'Color', settingsPDF.color, 'FontName', settingsPDF.fontName, 'Interpreter', 'none', 'VerticalAlignment', 'top');
    line = xASL_qc_ParsePdfConfig_sub_NewLine(line, settingsPDF);
end

% ====================================================================================================================================================

function line = xASL_qc_ParsePdfConfig_sub_NewLine(line, settingsPDF)
    line(2) = line(2) - (settingsPDF.lineSpacing) - (settingsPDF.fontSize * 0.001);

    if line(2) < 0 
        warning('No space left on page!');
    elseif line(2) < settingsPDF.canvas(2) 
        warning('Printing outside canvas, check block settings.');
    end
end


% ====================================================================================================================================================

function [figureCount] = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, settingsPDF, position)

    if isempty(header)
        figureCount = settingsPDF.figureCount;
        return
    else
        figureCount = settingsPDF.figureCount + 1;
    end

    position(4) = 0;
    settingsPDF.HorizontalAlignment = 'center';
    text = ['Fig ' num2str(figureCount) ': ' header];
    xASL_qc_ParsePdfConfig_sub_PrintText(text, currentFigure, position, settingsPDF);
end

% ====================================================================================================================================================

function line = xASL_qc_ParsePdfConfig_sub_PrintPatient(x, currentFigure, line, settingsPDF)
    % Check if participant data exists.
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
    
    line = xASL_qc_ParsePdfConfig_sub_PrintText('Participant Information', currentFigure, line, settingsPDF );
    for iEntry = 1:size(PatientInfo, 2)
        if strcmp(PatientInfo{1, iEntry}, 'participant_id')
            line = xASL_qc_ParsePdfConfig_sub_PrintText(['Participant: ', PatientInfo{2, iEntry}], currentFigure, line, settingsPDF );
        elseif strcmp(PatientInfo{1, iEntry}, 'age')
            line = xASL_qc_ParsePdfConfig_sub_PrintText(['Participant age: ', num2str(PatientInfo{2, iEntry})], currentFigure, line, settingsPDF );
        elseif strcmp(PatientInfo{1, iEntry}, 'sex')
            line = xASL_qc_ParsePdfConfig_sub_PrintText(['Participant sex: ', PatientInfo{2, iEntry}], currentFigure, line, settingsPDF );
        end
    end

end

function [strout] = xASL_qc_ParsePdfConfig_sub_WildcardReplace(strin, x)
    strout = strin;
    substring = regexp(strin, '<\w*>', 'match');
    for substringIndex=1:length(substring)
        if ~isfield(x, substring{substringIndex}(2:end-1))
            warning(['Could not replace ', substring{substringIndex}, ' check if file exists in ExploreASL/Derivatives/Population']);
        else
            strout = strrep(strout, substring{substringIndex}, x.(substring{substringIndex}(2:end-1)));
        end
    end
end


% ====================================================================================================================================================

function line = xASL_qc_ParsePdfConfig_sub_PrintQC(qcStruct, x, currentFigure, line, settingsPDF)
    % Stop if module and field don't exists in output
    if ~isfield(x.Output, (qcStruct.module)) && ~isfield(x.Output.(qcStruct.module), qcStruct.parameter)
        return     
    end

    % Use field specific settings if they exists, otherwise default to monospace for QC values. 
    settingsPDF.fontName = 'monospace';
    if isfield(qcStruct, 'settings')
        settingsPDF = xASL_qc_ParsePdfConfig_sub_loadSettings(qcStruct.settings, settingsPDF);
    end

    TempValue = x.Output.(qcStruct.module).(qcStruct.parameter);
    
    if ~isfield(qcStruct, 'alias') 
        qcStruct.alias = qcStruct.parameter;
    end

    if ~isfield(qcStruct, 'unit') 
        qcStruct.unit = '';
    end

    if isfield(qcStruct, 'range') 
        [range] = strsplit(qcStruct.range, '-');
        if (TempValue < xASL_str2num(range{1})) || (TempValue > xASL_str2num(range{2}))
            settingsPDF.color = 'r';
        end
        qcStruct.range = ['(' qcStruct.range ')'];
    else
        qcStruct.range = '';
    end

    if isnumeric(TempValue)
        TempValue = xASL_num2str(TempValue);
    end

    if size(TempValue, 1) == 1
        TextString = sprintf([sprintf('%-20s', [qcStruct.alias, ':']), sprintf('%8s', TempValue), sprintf('%-12s', [qcStruct.unit, ' ' , qcStruct.range]), ' \n']);
    end

    line = xASL_qc_ParsePdfConfig_sub_PrintText(TextString, currentFigure, line, settingsPDF);

end

% ====================================================================================================================================================
% ==============================================================Settings Based Functions==============================================================
% ====================================================================================================================================================

function [settingsPDF] = xASL_qc_ParsePdfConfig_sub_loadSettings(json, settingsPDF)

    fields = fieldnames(json);
    
    for iField = 1:length(fields)
        if strcmp(fields{iField}, 'fontSize') || strcmp(fields{iField}, 'lineSpacing')
            settingsPDF.(fields{iField}) = xASL_str2num(json.(fields{iField}));
        else
            settingsPDF.(fields{iField}) = json.(fields{iField});
        end 
    end
end

% ===================================================================================================================================================

function [settingsPDF] = xASL_qc_ParsePdfConfig_sub_defaultSettings()
    % This function sets the default settings for the PDF report.
    settingsPDF.color = 'k';
    settingsPDF.HorizontalAlignment = 'left';
    settingsPDF.fontWeight = 'normal';
    settingsPDF.axesVisible = 'off';
    settingsPDF.fontName = 'default';
    settingsPDF.lineSpacing = 0.005;
    settingsPDF.figureCount = 0;
    settingsPDF.canvas = [0 0 1 1];

    if ispc
        settingsPDF.fontSize = 10;
    elseif isunix 
        settingsPDF.fontSize = 8;
    else 
        error('xASL_qc_ParsePdfConfig couldnt find OS')
    end 

end