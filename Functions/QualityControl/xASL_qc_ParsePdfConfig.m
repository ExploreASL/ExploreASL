function [settingsPDF] = xASL_qc_ParsePdfConfig(layoutStructure, x, currentFigure, line, settingsPDF)
% xASL_qc_ParsePdfConfig function used by xASL_qc_GenerateReport to parse the configuration file loaded by xASL_adm_LoadPdfConfig.
%   
%   
% FORMAT: xASL_qc_ParsePdfConfig(layoutStructure, x[, currentFigure, line, settingsPDF])
%
% INPUT:
%   layoutStructure     - json structure containing all information to be printed (REQUIRED)
%   x                   - structure containing fields with all information required to run this submodule (REQUIRED)
%   currentFigure       - currentFigure handle to print to (OPTIONAL, defaults to matlab main figure)
%   line                - line to print to (OPTIONAL, defaults to [0 0.93 1 0] when it generates a new page)
%   settingsPDF         - settings used to print (OPTIONAL, will set default settings if not specified)
%
% OUTPUT: 
%   settingsPDF    - settings to print with
%
% OUTPUTFILE:
%   xASL_Report_SubjectName.pdf - printed PDF rapport containing QC images & values
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  xASL_qc_ParsePdfConfig function used by xASL_qc_GenerateReport to parse the configuration file loaded by xASL_adm_LoadPdfConfig.
%               This function will run recursively, and will print all text, images, scans and other content specified in the json file.
%               The json file should contain a structure with fields specified in the manual.
% 
% EXAMPLE: xASL_qc_ParsePdfConfig(layoutStructure, x);
% __________________________________
% Copyright (C) 2015-2023 ExploreASL

if nargin < 2 
    error('xASL_qc_ParsePdfConfig requires at least 2 input arguments');
end

if nargin < 3 || isempty(currentFigure)
    currentFigure = figure('Visible','off');
end

if nargin < 4 || isempty(line)
    line = [0 0.93 1 1];
end

if nargin < 5 || isempty(settingsPDF)
    settingsPDF = xASL_qc_ParsePdfConfig_sub_defaultSettings(x);
end

fields = fieldnames(layoutStructure);

for iField = 1:length(fields)
    % xASL_TrackProgress(iField, length(fields));
    % This progress counting doesnt work, because the loop loops in itself
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
% This function parses the content of the json file, and calls the appropriate function to print the content.

    % It first checks if the currentField is a struct, and if it contains a field "type" which specifies the type of content to be printed.
    % If the currentField is not a struct, it will ignore the currentField and the next iteration will start.
    if ~isstruct(currentField)
        return;
    end

    if ~isfield(currentField, 'type') || ~isfield(currentField, 'category')
        error([currentField, 'contains no field specifying the category and type of content, add field "category": "content" or "type": "image2D" to this json field for example']);
    end

    % Depending on the type of content, it will call the appropriate function to print the content.
    % With exception of the type "settings", which is used to change variables like the font instead.
    switch currentField.category
        case 'content' 
            switch currentField.type
            case 'text' 
                line = xASL_qc_ParsePdfConfig_sub_PrintText(currentField, currentFigure, line, settingsPDF);
            case 'QCValues'
                line = xASL_qc_ParsePdfConfig_sub_PrintQC(currentField, x, currentFigure, line, settingsPDF);
            case 'image2D'
                settingsPDF = xASL_qc_ParsePdfConfig_sub_PrintImage(currentField, x, currentFigure, settingsPDF);  
            case 'image3D'
                settingsPDF = xASL_qc_ParsePdfConfig_sub_PrintScan(currentField, x, currentFigure, settingsPDF);
            case 'patients' 
                line = xASL_qc_ParsePdfConfig_sub_PrintPatient(x, currentFigure, line, settingsPDF);
        end  
        case 'metadata'
            switch currentField.type
            case 'page'
                xASL_qc_ParsePdfConfig_sub_printPage(currentField, x, settingsPDF);  
            case 'block'
                xASL_qc_ParsePdfConfig_sub_printBlock(currentField, x, currentFigure, settingsPDF);
            case 'textSettings'
                settingsPDF = xASL_qc_ParsePdfConfig_sub_loadSettings(currentField, settingsPDF);  
        end  
    end
end

% ====================================================================================================================================================

function  xASL_qc_ParsePdfConfig_sub_printPage(pageStruct, x, settingsPDF)
% This function prints pages using the layout defined in the json file.
% It first creates a new figure, and then iterates over and prints all content in the pageStruct
% The pageStruct should contain a field "content" which contains all content to be printed on the page.
% The pageStruct should also contain a field "pageIdentifier" which is used to name the printed PDF file.

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

    % Finally it prints the page to a PDF file using the pageIdentifier as filename in the subject directory.
    fileName = ['xASL_Report_', pageStruct.pageIdentifier];
    printPathFile = fullfile(x.dir.xASLDerivatives, x.SUBJECT, fileName); 
    print(figPrimary, printPathFile, '-dpdf', '-bestfit');

end

% ====================================================================================================================================================

function xASL_qc_ParsePdfConfig_sub_printBlock(blockStruct, x, pageFig, settingsPDF)
% This function prints a content blocks using the layout defined in the json file.
% Content blocks can be used when you have a predefined layout that you want to use multiple times.
% The blockStruct should contain a field "content" which contains all content to be printed in the block.
% Using the function xASL_qc_ParsePdfConfig_sub_createNewCanvas it will create a new subcanvas for the block to be printed in.

    position = xASL_str2num(blockStruct.position);
    size = xASL_str2num(blockStruct.size);
    [settingsPDF.canvas, line] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position, size, settingsPDF.canvas);
    xASL_qc_ParsePdfConfig(blockStruct, x, pageFig, line, settingsPDF);
end

% ====================================================================================================================================================

function [newCanvas, line] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position, size, oldCanvas)
% This function creates a new canvas based on the old canvas, and the position and size of the new canvas.
% The canvas is used to define the position and size where content is printed.
% This is for example used to define where a new image, scan or block is printed on the page.
% This function ensures that the new canvas is within the old canvas, and scaled appropriately.
% This is for example important when a predefined block contains images, which need to be scaled appropriately.
    
    % A canvas is defined as a 4 element vector [x y width height] where x and y are the lower left corner of the canvas.
    % This function also returns the variable line, which is the position of the next line to be printed.
    newPosition = [oldCanvas(1) + position(1) * oldCanvas(3) oldCanvas(2) + position(2) * oldCanvas(4)];
    newSize = [oldCanvas(3:4) .* size];
    newCanvas = [newPosition, newSize];
    line = [newPosition(1)  newPosition(2) + newSize(2) newSize(1) 0];
end

% ====================================================================================================================================================
% ==============================================================Image Based Functions=================================================================
% ====================================================================================================================================================

function [settingsPDF] = xASL_qc_ParsePdfConfig_sub_PrintImage(input, x, currentFigure, settingsPDF, position)
% This function prints images using the layout defined in the json file.

    % It first checks what the type is of the input.
    % If the image is defined as a struct, it will use the fields in the struct to define the image.
    % If the image is defined as a string, it will use the string as the path to the image.
    header = '';
    switch nargin
    case 4
        imageStruct = input;
        position = [xASL_str2num(imageStruct.position) xASL_str2num(imageStruct.size)];

        % It then checks what the user used to define the path to the image.
        % In the json file, the user can define the path to the image in 4 different ways:
        % If the path isnt defined in the json file, it will use throw a warning and return.
        if isfield(imageStruct, 'xPath') && isfield(x.P, imageStruct.name)
            ImagePath = x.P.(imageStruct.name);
        elseif isfield(imageStruct, 'absolutePath')
            ImagePath = imageStruct.absolutePath;
            ImagePath = xASL_qc_ParsePdfConfig_sub_WildcardReplace(ImagePath, x, settingsPDF);
        elseif isfield(imageStruct, 'popPath')
            ImagePath = fullfile(x.dir.xASLDerivatives, 'Population', imageStruct.popPath);
            ImagePath = xASL_qc_ParsePdfConfig_sub_WildcardReplace(ImagePath, x, settingsPDF);
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

    % First it calculates the size of the canvas for the image to be printed in.
    [canvas] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position(1:2), position(3:4), settingsPDF.canvas);

    % Finally it prints the image to the current figure, and updates the figure count.
    ax = axes('Position', canvas, 'Visible', settingsPDF.axesVisible, 'Parent', currentFigure);
    [img, ~, alphachannel] = imread(ImagePath);
    fg = imshow(img);
    %imagesc(img);

    % add alpha channel if png.
    if strcmp('png', ImagePath(end-2:end))
        fg.AlphaData = alphachannel;
    end

    settingsPDF.figureCount = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, settingsPDF, canvas);
end

% ====================================================================================================================================================

function [settingsPDF] = xASL_qc_ParsePdfConfig_sub_PrintScan(scanStruct, x, currentFigure, settingsPDF)
% This function prints scans using the layout defined in the json file.

    % It first checks if the requected scan exists in the x.P structure, and if it doesnt it will throw a warning and return.
    if ~isfield(x.P, scanStruct.name)
        warning (['could not print ', scanStruct.name, ', check if NIfTI exists in ExploreASL/Derivatives/Population']);
        return
    end

    % It then makes a list of all scans to be printed, and creates a header for the image.
    ImIn = {x.P.(scanStruct.name)};
    header = scanStruct.name;
    
    % Any overlay scans to be printed are added to the list of scans to be printed, and added to the header.
    if isfield(scanStruct, 'overlay')
        fields = fieldnames(scanStruct.overlay);
        for iField = 1:length(fields)
            if isfield(x.P, fields(iField))
                ImIn(end+1) = {x.P.(fields{iField})};
                header = [header ' + ' fields{iField}];
            end
        end
    end

    % The slices to be printed are defined in the json file, and are added to the x.S structure for printing.
    if isfield(scanStruct, 'slice')
        [x.S.TraSlices] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(scanStruct.slice, 'TraSlice');
        [x.S.CorSlices] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(scanStruct.slice, 'CorSlice');
        [x.S.SagSlices] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(scanStruct.slice, 'SagSlice');       
    else
        x.S.TraSlices = [25, 50, 90];
        x.S.SagSlices = [25, 50, 90];
        x.S.CorSlices = [25, 50, 90];
    end

    % Create the image from the scans defined in the json file.
    imageToPrint = xASL_vis_CreateVisualFig(x, ImIn, [], [], [], [], [], [], [], [], [], 0); % no verbosity   

    % First it calculates the size of the canvas for the image to be printed in.
    position = [xASL_str2num(scanStruct.position) xASL_str2num(scanStruct.size)];
    [canvas] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position(1:2), position(3:4), settingsPDF.canvas);

    % Finally it prints the image to the current figure, and updates the figure count.
    ax = axes('Position', canvas, 'Visible', settingsPDF.axesVisible, 'Parent', currentFigure);
    fg = imshow(imageToPrint);

    % Print the image header and update the figure count
    settingsPDF.figureCount = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, settingsPDF, canvas);
end

% ====================================================================================================================================================

function [slice] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(struct, name)
% This function changes the strings in the json file to numbers, and returns the slice to be printed.
% This specifically cannot use xASL_str2num, because the way that function returns NaNs is not compatible the image generation function.

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
% This function prints text using the layout defined in the json file.

    % It first checks what the type is of the input.
    % If the text is defined as a struct, it will use the fields in the struct to define the text.
    % If the text is defined as a string, it will use the string as the text to be printed.
    switch class(input)
        case {'string', 'char'} 
            String = input;
        case 'struct'
            textStruct = input;
            String = textStruct.text;
            % It then checks if the textStruct contains a field "settings" which is used to change variables like the font.
            if isfield(textStruct, 'textSettings')
                settingsPDF = xASL_qc_ParsePdfConfig_sub_loadSettings(textStruct.textSettings, settingsPDF);
            end
            % If the input cannot be parsed, it will throw an error.
        otherwise
            class(input)
            error('xASL_qc_ParsePdfConfig_sub_PrintText couldnt find string of printable text')
    end

    % It then prints the text to the current figure.
    ax = axes('Position', line , 'Visible', settingsPDF.axesVisible, 'Parent', currentFigure);
    String = strrep(String, '_', ' ');
    text(0, 0, String, 'Parent', ax, 'FontSize', settingsPDF.fontSize, 'FontWeight', settingsPDF.fontWeight, 'Color', settingsPDF.color, 'FontName', settingsPDF.fontName, 'Interpreter', 'none', 'VerticalAlignment', 'top');
    
    % And finally it updates the line position for the next line to be printed.
    line = xASL_qc_ParsePdfConfig_sub_NewLine(line, settingsPDF);
end

% ====================================================================================================================================================

function line = xASL_qc_ParsePdfConfig_sub_NewLine(line, settingsPDF)
% This function updates the line position for the next line to be printed.

    % The "newline" distance has so far been hardcoded based on experience, but this can be improved in the future.
    line(2) = line(2) - (settingsPDF.lineSpacing) - (settingsPDF.fontSize * 0.001);

    % It simply checks if the line position is lower than 0, and if it is it will throw a warning that the text is printed outside the canvas.
    if line(2) < 0 
        warning('No space left on page!');
    elseif line(2) < settingsPDF.canvas(2) 
        warning('Printing outside canvas, check block settings.');
    end
end


% ====================================================================================================================================================

function [figureCount] = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, settingsPDF, position)
% This function prints a header underneith the image to be printed.

    % If no header is specified, it will exit and not iterate the figure count.
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
% This function prints the patient information to the PDF report.
% The patient information is extracted from the participants.tsv file in the derivatives directory.
% The patient information is printed similar to how a newline is printed.

    % Check if participant data exists.
    ParticipantsTSV = xASL_adm_GetFileList(x.D.ROOT, 'participants.tsv');
    PatientInfo = {};


    % Check if participants.tsv exists
    if isempty(ParticipantsTSV)
        fprintf ('No participants.tsv found in the ExploreASL derivatives folder.\n');
        return
    end

    % Extract participant information from participants.tsv
    structParticipants = xASL_tsvRead(ParticipantsTSV{1});
    nParticipants = size(structParticipants, 1);
    for iPar = 2:nParticipants
        if structParticipants{iPar, 1} == x.SUBJECT
            PatientInfo = [structParticipants(1, :); structParticipants(iPar, :)];
        end
    end

    % Check if participant information exists in the participants.tsv
    if isempty(PatientInfo)
        fprintf ('No patient information found in participants.tsc in the ExploreASL derivatives folder.\n');
        return
    end
    
    % Print the patient information to the PDF report. (WIP make it oneline?)
    settingsTitle = settingsPDF;
    settingsTitle.fontWeight = 'bold';
    line = xASL_qc_ParsePdfConfig_sub_PrintText('Participant Information', currentFigure, line, settingsTitle );

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

function [strout] = xASL_qc_ParsePdfConfig_sub_WildcardReplace(strin, x, settingsPDF)
% This function replaces wildcards in the path to the image.
% Wildcards are defineds as <wildcard> in the json file, and are replaced with the corresponding field in the x structure.
% E.g., <SUBJECT> is replaced with the the value in the x.SUBJECT field

    strout = strin;
    substring = regexp(strin, '<\w*>', 'match');
    if  settingsPDF.BIDS_Translation
        substring = xASL_qc_ParsePdfConfig_sub_BIDS_Translation(substring);
    end

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
% This function prints QC values using the layout defined in the json file.
% QC Values are extracted from the x.Output structure, and printed to the PDF report in a single line.
% Certain settings can be applied to the QC values, for example a range can be specified.
% If the QC value is outside the range, it will be printed in red.

    % Stop if module and field don't exists in output
    if ~isfield(x.Output, (qcStruct.module)) || ~isfield(x.Output.(qcStruct.module), qcStruct.parameter)
        return     
    end

    % Use field specific settings if they exists, otherwise default to monospace for QC values. 
    settingsPDF.fontName = 'monospace';
    if isfield(qcStruct, 'textSettings')
        settingsPDF = xASL_qc_ParsePdfConfig_sub_loadSettings(qcStruct.textSettings, settingsPDF);
    end

    % Generate the string to be written in the pdf Report
    string = xASL_qc_ParsePdfConfig_sub_Generate_QC_String(qcStruct, x, settingsPDF);

    % Print the string to the PDF report.
    line = xASL_qc_ParsePdfConfig_sub_PrintText(string, currentFigure, line, settingsPDF);

end

% ====================================================================================================================================================
% ==============================================================Settings Based Functions==============================================================
% ====================================================================================================================================================

function [settingsPDF] = xASL_qc_ParsePdfConfig_sub_loadSettings(json, settingsPDF)
% This function replaces existing settings with new ones from the json file.
% If the json file contains a field "fontSize" or "lineSpacing" it will convert the string to a number.
% Otherwise it will simply replace the existing setting with the new one.

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

function [string] = xASL_qc_ParsePdfConfig_sub_BIDS_Translation(string)
    % This function replaces BIDS terminolgy into ExploreASL Legacy terminology.
    string = strrep(string, 'RUN', 'SESSION');
    string = strrep(string, 'SESSION', 'VISIT');       
end

function [string] = xASL_qc_ParsePdfConfig_sub_Generate_QC_String(qcStruct, x, settingsPDF)
    
    if ~isfield(x.Output, (qcStruct.module)) 
        fprintf (['QC Value ' x.Output.(qcStruct.module) ' not found, skipping. \n']);
        return
    end 

    if ~isfield(x.Output.(qcStruct.module),( qcStruct.parameter)) 
        fprintf (['QC Value ' x.Output.(qcStruct.module).(qcStruct.parameter) ' not found, skipping. \n']);
        return
    end

    TempValue = x.Output.(qcStruct.module).(qcStruct.parameter);
    
    if settingsPDF.QC_Value_alias
        qcStruct = xASL_qc_ParsePdfConfig_sub_QC_Translation(qcStruct, settingsPDF);
    end

    % Check if the QC has an alias, if so replaced the parameter with the alias in the printed text.
    if ~isfield(qcStruct, 'alias') || isempty(qcStruct.alias)
        qcStruct.alias = qcStruct.parameter;
    end

    % Check if the QC has a unit, if not set it to an empty string.
    if ~isfield(qcStruct, 'unit') 
        qcStruct.unit = '';
    end

    % Check if the QC has a range, if so check if the value is within the range, else print in red.
    if isfield(qcStruct, 'range') 
        if ~isempty(qcStruct.range)
            [range] = strsplit(qcStruct.range, '-');
            if (TempValue < xASL_str2num(range{1})) || (TempValue > xASL_str2num(range{2}))
                settingsPDF.color = 'r';
            end
            qcStruct.range = ['(' qcStruct.range ')'];
        end
    else
        qcStruct.range = '';
    end

    % Convert the value to a string, and print it to the PDF report.
    if isnumeric(TempValue)
        TempValue = xASL_num2str(TempValue);
    end

    % Combine name, value, unit and range into a single string for printing.
    if size(TempValue, 1) == 1
        string = sprintf([sprintf('%-20s', [qcStruct.alias, ':']), sprintf('%8s', TempValue), sprintf('%-12s', [qcStruct.unit, ' ' , qcStruct.range]), ' \n']);
    end

end

function [struct] = xASL_qc_ParsePdfConfig_sub_QC_Translation(struct, settingsPDF)
    % This function replaces QC values with long names with easier to read names with units.
    name = struct.parameter;
    index = find(strcmp(settingsPDF.QC_Translation(:,1), name));
    if ~isempty(index)
        struct.alias = char(settingsPDF.QC_Translation(index, 2));
        struct.unit  = char(settingsPDF.QC_Translation(index, 3));
        struct.range = char(settingsPDF.QC_Translation(index, 4));
    end
end

function [settingsPDF] = xASL_qc_ParsePdfConfig_sub_defaultSettings(x)
% This function sets the default settings for the PDF report.

    settingsPDF.color = 'k';
    settingsPDF.HorizontalAlignment = 'left';
    settingsPDF.fontWeight = 'normal';
    settingsPDF.axesVisible = 'off';
    settingsPDF.fontName = 'default';
    settingsPDF.lineSpacing = 0.005;
    settingsPDF.figureCount = 0;
    settingsPDF.canvas = [0 0 1 1];
    settingsPDF.BIDS_Translation = 0;
    settingsPDF.QC_Value_alias = 1; %0 = no translation, 1=translation

    settingsPDF.QC_Translation = xASL_tsvRead(fullfile(x.opts.MyPath, 'Functions', 'QualityControl', 'qc_translation.tsv'));


    if ispc
        settingsPDF.fontSize = 10;
    elseif isunix 
        settingsPDF.fontSize = 8;
    else 
        error('xASL_qc_ParsePdfConfig couldnt find OS')
    end 

end