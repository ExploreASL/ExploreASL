function [x,config] = xASL_qc_GenerateReport(x, subject)
% xASL_qc_GenerateReport Print output QC
%
% FORMAT: xASL_qc_GenerateReport(x)
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%
% OUTPUT: n/a
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
    for iFile=1:size(ExistingPrintFiles)
        xASL_delete(ExistingPrintFiles{iFile});
    end
end
% Delete existing xASL report files
ExistingPrintFiles = xASL_adm_GetFileList(PrintDir, '^xASL_subImage.+$');
if ~isempty(ExistingPrintFiles)
    for iFile=1:size(ExistingPrintFiles)
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
[~] = xASL_sub_parseJson(config, x, [], [], settings);

end


%% ================================================================================================================================================
%% ================================================================================================================================================
function  xASL_sub_printPage(json, x, settings)

    clear primaryFig ax  
    primaryFig = figure();
    set(primaryFig, 'visible', 'off');
    set(primaryFig, 'Units', 'centimeters' );
    set(primaryFig, 'Position', [0 0 21 29.7] );


    %% Print the Title
    logoPath = fullfile(x.opts.MyPath, 'Design/ExploreASL_logoHeader.png');
    xASL_sub_PrintImage(logoPath, [], primaryFig, [0 0.96 1 0.04]);
    line = xASL_sub_PrintText('ExploreASL QC summary report ', primaryFig, [0 0.96 1 0], settings);

    %% Print the Footer
    xASL_sub_PrintText('This report was automatically generated by ExploreASL', primaryFig,  [0 0.02 1 0], settings);

    %% Parse JSON and create the file
    [~] = xASL_sub_parseJson(json, x, primaryFig, line, settings);

    PrintFile = ['xASL_Report', json.number];
    PrintPath = fullfile(x.dir.xASLDerivatives, x.SUBJECT, PrintFile); 

    print(primaryFig, PrintPath, '-dpdf', '-bestfit');

    clear primaryFig ax
end


function xASL_sub_parseBlock(input, x, pageFig, settings)
    clear local_image
    position = [str2num(input.position) str2num(input.size)];
    local_image = figure();
    set(local_image, 'visible', 'off');
    set(local_image, 'Units', 'normalized' );
    set(local_image, 'Position', position );

    xASL_sub_parseJson(input, x, local_image, [0 0.96 1 0.04], settings);

    PrintFile = 'xASL_subImage';
    PrintPath = fullfile(x.dir.xASLDerivatives, x.SUBJECT, PrintFile); 
    
    print(local_image, PrintPath, '-djpeg')
    clear local_image
    
    figure(pageFig);
    set(pageFig, 'visible', 'off');
    xASL_sub_PrintImage([PrintPath '.jpg'], [], pageFig, position);    
end


function [settings] = xASL_sub_parseJson(json, x, figure, line, settings)
    fields = fieldnames(json);
    for iField=1:length(fields)
        currentField=json.(fields{iField});

        if ~isstruct(currentField)
            continue;
        end

        if ~isfield(currentField, 'type')
            error([currentField, 'containts no field type']);
        end

        switch currentField.type
            case 'text' 
                line = xASL_sub_PrintText(currentField, figure, line, settings);
            case 'qc'
                line = xASL_sub_PrintQC(currentField, fields{iField}, x, figure, line, settings);
            case 'image'
                xASL_sub_PrintImage(currentField, x, figure);  
            case 'scan'
                xASL_sub_PrintScan(currentField, x, figure);
            case 'page'
                xASL_sub_printPage(currentField, x, settings);  
            case 'block'
                xASL_sub_parseBlock(currentField, x, figure, settings);
            case 'settings'
                settings = xASL_sub_loadSettings(currentField, settings);  
            case 'patients' 
                PatientInfo = xASL_sub_readParticipantsTSV(x);
                line = xASL_sub_PrintPatient(PatientInfo, figure, line, settings);
        end  

    end
end

function PatientInfo = xASL_sub_readParticipantsTSV(x)
    %Check if participant data exists.
    ExistsParticipants = xASL_adm_GetFileList(x.D.ROOT, 'participants.tsv');
    PatientInfo = {};

    if isempty(ExistsParticipants)
        return
    end

    sParticipants = xASL_tsvRead(ExistsParticipants{1});
    nParticipants = size(sParticipants, 1);
    for iPar=2:nParticipants
        if sParticipants{iPar, 1} == x.SUBJECT
            %Extract Participants table, and participants info
            PatientInfo = [sParticipants(1,:); sParticipants(iPar, :)];
        end
    end
    
end

function line = xASL_sub_PrintPatient(PatientInfo, figure, line, settings)
    if isempty(PatientInfo)
        fprintf ('no patient info found to be printed \n');
        return
    end
    line = xASL_sub_PrintText('Participant Information', figure, line, settings );
    for iEntry=1:size(PatientInfo,2)
        if strcmp(PatientInfo{1, iEntry},'participant_id')
            line = xASL_sub_PrintText(['Participant: ', PatientInfo{2, iEntry}], figure, line, settings );
        elseif strcmp(PatientInfo{1, iEntry},'age')
            line = xASL_sub_PrintText(['Participant age: ', num2str(PatientInfo{2, iEntry})], figure, line, settings );
        elseif strcmp(PatientInfo{1, iEntry},'sex')
            line = xASL_sub_PrintText(['Participant sex: ', PatientInfo{2, iEntry}], figure, line, settings );
        end
    end

end

function line = xASL_sub_PrintQC(json, name, x, figure, line, settings)
    % Printing header/title
    json.name = name;

    % Stop if module does not exists in output
    if ~isfield(x.Output, (json.module))
        return     
    end

    % Stop if field does not exists in output
    if ~isfield(x.Output.(json.module), json.name)
        return     
    end

    % Use field specific settings if they exists
    settings.fontName = 'monospace';
    if isfield(json, 'settings')
        settings = xASL_sub_loadSettings(json.settings, settings);
    end

    % Print only if visibility exists and is true. 
    if isfield(json, "visible")  && json.visible 
        
        TempValue = x.Output.(json.module).(json.name);
        
        if ~isfield(json, 'alias') 
            json.alias = name;
        end

        if ~isfield(json, 'unit') 
            json.unit = '';
        end

        if isfield(json, 'range') 
            [range] = strsplit(json.range, '-');
            if (TempValue < str2double(range{1})) || (TempValue > str2double(range{2}))
                settings.color = 'r';
            end
        else
            json.range = '';
        end

        if isnumeric(TempValue)
            TempValue = num2str(TempValue);
        end

        if size(TempValue, 1) == 1
            TextString = sprintf([sprintf('%-20s',[json.alias,':']), sprintf('%8s',TempValue),sprintf('%-12s', [json.unit, ' ' , json.range]), ' \n']);
        end

        line = xASL_sub_PrintText(TextString, figure, line, settings);
    end
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
    line = xASL_sub_NewLine(line, settings);
    ax=axes('Position', line ,'Visible', settings.axesVisible, 'Parent', figure);
    text(0,0, String, 'Parent', ax, 'FontSize', settings.fontSize,'FontWeight', settings.fontWeight, 'Color', settings.color, 'FontName', settings.fontName, 'Interpreter', 'none', 'VerticalAlignment','top');
end

function xASL_sub_PrintImage(input, x, figure, position)
    switch nargin
    case 3
        position = [str2num(input.position) str2num(input.size)];
        ImagePath = input.filePath;
    case 4
        ImagePath = input;
    end
    ax=axes('Position', position, 'Visible', 'off', 'Parent', figure);
    [img, ~, alphachannel] = imread(ImagePath);
    fg= imshow(img);
    if strcmp('png', ImagePath(end-2:end))
        fg.AlphaData=alphachannel;
    end
    clear fg
end

function xASL_sub_PrintScan(input, x, figure)
    position = [str2num(input.position) str2num(input.size)];
    if ~isfield(x.P, input.name)
        warning (['could not print', input.name, 'check if nifti exists in ExploreASL/Derivatives/Population']);
        return
    end
    ImIn = {x.P.(input.name)};

    ax=axes('Position', position, 'Visible', 'off', 'Parent', figure);
    if isfield(input, 'overlay')
        fields = fieldnames(input.overlay);
        for iField=1:length(fields)
            if isfield(x.P, fields(iField))
                ImIn(end+1) = {x.P.(fields{iField})};
            end
        end
    end
    if isfield(input, 'slice')
        x.S.TraSlices = [str2num(input.slice.TraSlice)];
        x.S.CorSlices = [str2num(input.slice.CorSlice)];
        x.S.SagSlices = [str2num(input.slice.SagSlice)];
    else
        x.S.TraSlices = [25, 50, 90];
        x.S.SagSlices = [25, 50, 90];
        x.S.CorSlices = [25, 50, 90];
    end

    img = xASL_vis_CreateVisualFig(x, ImIn, [], [], [], [], [], [], [], [], [], [], []);
    fg= imshow(img);
    clear fg
end

function line = xASL_sub_NewLine(line, settings);
    line(2) = line(2) - (settings.fontSize*0.0022);
    if line(2) < 0 
        warning('No space left on page!');
    end
end

function [settings] = xASL_sub_loadSettings(json, settings)
    if isfield(json, 'color')
        settings.color = json.color;
    end  
    if isfield(json, 'HorizontalAlignment')
        settings.HorizontalAlignment = json.HorizontalAlignment;
    end  
    if isfield(json, 'fontWeight')
        settings.fontWeight = json.fontWeight;
    end  
    if isfield(json, 'fontSize')
        settings.fontSize = str2double(json.fontSize);
    end   
    if isfield(json, 'axesVisible')
        settings.axesVisible = json.axesVisible;
    end   

end

function [settings] = xASL_sub_defaultSettings()
    settings.color = 'k';
    settings.HorizontalAlignment = 'left';
    settings.fontWeight = 'normal';
    settings.axesVisible = 'off';
    settings.fontName = 'default';
    if ispc
        settings.fontSize = 10;
    elseif isunix || ismac
        settings.fontSize = 8;
    end 
end
