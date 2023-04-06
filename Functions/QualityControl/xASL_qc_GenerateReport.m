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

x = xASL_adm_LoadX(x, PathX, false); % Assume memory x is newer than x.mat

% Make sure that the directory exists
PrintDir = fullfile(x.dir.xASLDerivatives, subject);
xASL_adm_CreateDir(PrintDir);

% Delete existing xASL report files
ExistingPrintFiles = xASL_adm_GetFileList(pwd, '^xASL_Report_.+$');
if ~isempty(ExistingPrintFiles)
    xASL_delete(ExistingPrintFiles{1});
end

% Set defaults
settings = xASL_sub_defaultSettings();

%% Print the title
fprintf('Printing ExploreASL PDF report:   \n');

clear Str spm_fig ax OutputFields StrucFields
spm_fig = spm_figure('Create','Graphics','visible','off');
set(spm_fig,'windowstyle','normal'); 
spm_figure('Clear','Graphics');

%% Print the Title
logoPath = fullfile(x.opts.MyPath, 'Design/ExploreASL_logoHeader.png');
xASL_sub_PrintImage(logoPath, spm_fig, [0 0.96 1 0.04]);
line = xASL_sub_PrintText('ExploreASL QC summary report ', spm_fig, [0 0.96 1 0], settings);

%% Print the Footer
xASL_sub_PrintText('This report was automatically generated with ExploreASL', spm_fig,  [0 0.02 1 0], settings);

%% Load Pdf configuration
config = xASL_adm_LoadPDFConfig(x);

% Load participants data
PatientInfo = xASL_sub_readParticipantsTSV(x);
line = xASL_sub_PrintPatient(PatientInfo, spm_fig, line, settings);

%% Parse JSON and create the file
[line, settings] = xASL_sub_parseJson(config, x, spm_fig, line, settings);

PrintFile = 'xASL_Report';
PrintPath = fullfile(PrintDir, PrintFile); 

% Spm print adds _001 to the name of the file, I havent figured out how to get rid of this yet(ugly);
xASL_delete(strcat(PrintPath, '_001.pdf')); 
spm_figure('Print',spm_fig, PrintPath, 'pdf');

clear ax spm_fig

end


%% ================================================================================================================================================
%% ================================================================================================================================================
function [line, settings] = xASL_sub_parseJson(json, x, figure, line, settings)
    fields = fieldnames(json);
    for iField=1:length(fields)
        currentField=json.(fields{iField});

        if ~isfield(currentField, 'type')
            error([currentField, 'containts no field type']);
        end

        switch currentField.type
            case 'text' 
                line = xASL_sub_PrintText(currentField, figure, line, settings);
            case 'qc'
                line = xASL_sub_PrintQC(currentField, fields{iField}, x, figure, line, settings);
            case 'image'
                xASL_sub_PrintImage(currentField, figure);  
            case 'settings'
                settings = xASL_sub_loadSettings(currentField, settings);  
        end  

    end
end

function PatientInfo = xASL_sub_readParticipantsTSV(x)
    %Check if participant data exists.
    ExistsParticipants = xASL_adm_GetFileList(x.D.ROOT, 'participants.tsv');

    if isempty(ExistsParticipants)
        PatientInfo = {};
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
        fprintf ("no patient info found to be printed")
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
            if (TempValue < str2double(range{1})) | (TempValue > str2double(range{2}))
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
    line = xASL_sub_NewLine(line);
    ax=axes('Position', line ,'Visible', settings.axesVisible, 'Parent', figure);
    text(0,0, String, 'Parent', ax, 'FontSize', settings.fontSize,'FontWeight', settings.fontWeight, 'Color', settings.color, 'FontName', settings.fontName, 'Interpreter', 'none', 'VerticalAlignment','top');
end

function xASL_sub_PrintImage(input, figure, position)
    switch nargin
    case 2
        position = [str2num(input.position) str2num(input.size)];
        ImagePath = input.filePath;
    case 3
        ImagePath = input;
    end
    ax=axes('Position', position, 'Visible', 'off', 'Parent', figure);
    [img, ~, alphachannel] = imread(ImagePath);
    fg= imshow(img);
    fg.AlphaData=alphachannel;
    clear fg
end

function line = xASL_sub_NewLine(line);
    line(2) = line(2) - 0.016;
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
        settings.fontSize = str2num(json.fontSize);
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
