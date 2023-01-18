function xASL_qc_GenerateReport(x)
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
% Copyright (C) 2015-2022 ExploreASL

% Determine x.mat file
PathX = fullfile(x.dir.SUBJECTDIR,'x.mat');

% Check if x.mat file exists already
if ~exist(PathX, 'file')
    warning([PathX ' didnt exist, skipping xASL_qc_CreateOutputPDF']);
    return;
end

x = xASL_adm_LoadX(x, PathX, false); % Assume memory x is newer than x.mat

% Determine subject and session
SuSeID = [x.SUBJECT '_' x.SESSIONS{1}];
iSubjSess = (1)*x.dataset.nSessions+1;

% Make sure that the directory exists
PrintDir = fullfile(x.D.ROOT, x.SUBJECT);
xASL_adm_CreateDir(PrintDir);

% Delete existing xASL report files
ExistingPrintFiles = xASL_adm_GetFileList(pwd, '^xASL_Report_.+$');
if ~isempty(ExistingPrintFiles)
    xASL_delete(ExistingPrintFiles{1});
end

% Determine font size
if ispc
    fontsize = 10;
elseif isunix || ismac
    fontsize = 8;
end

% Load participants data
PatientInfo = xASL_adm_readParticipantsTSV(x);

% Debugg print participants
if ~isempty(PatientInfo)
    PatientInfo
end

%% Print the title
fprintf('Printing ExploreASL PDF report:   \n');

clear Str spm_fig ax OutputFields StrucFields
spm_fig = spm_figure('Create','Graphics','visible','off');
set(spm_fig,'windowstyle','normal'); 
spm_figure('Clear','Graphics');

%% Print the Title
PrintImage('Design/ExploreASL_logoHeader.png',spm_fig,[0 0.96 1 0.04]);
NewLine = PrintBold(['ExploreASL QC summary report ', SuSeID], spm_fig, [0 0.96 1 0], fontsize+2);

%% Print the Footer
PrintBold('This report was automatically generated with ExploreASL', spm_fig,  [0 0.02 1 0], fontsize+2);

%% Check if config file exists, else make one
%% Load Standard Configuration
cPath = strcat(x.dir.xASLDerivatives, '/ConfigReportPDF.json');
if isfile(cPath)
    config = xASL_adm_LoadPDFConfig(cPath);
else
    fprintf('configReportPDF.json does not exists, using default config \n');
    sPath = strcat(x.opts.MyPath, '/Functions/QualityControl/templateConfigReportPDF.json');
    config = xASL_adm_LoadPDFConfig(sPath);
end

%% -----------------------------------------------------------------------------------------------
%% Collect field name & field values to print
OutputFields = fieldnames(x.Output);
nOutputFields = length(OutputFields);

if nOutputFields>5
    fprintf('%s\n',['Printing QC for ' num2str(nOutputFields) ' modalities, not sure if this fits']);
end 

for iField=1:nOutputFields
    if length(x.Output.(OutputFields{iField}))<iSubjSess
        Ind = 1; 
    else
        Ind = iSubjSess;
    end
    StrucFields{iField} = sort(fieldnames(x.Output.(OutputFields{iField})(Ind)));
    for iV=1:length(StrucFields{iField})     
            Paragraphs{iField}(iV).name = StrucFields{iField}{iV};
            Paragraphs{iField}(iV).value = x.Output.(OutputFields{iField})(Ind).(StrucFields{iField}{iV});
    end
end 

for iField=1:length(Paragraphs) 
    NewLine = PrintParagraph(x, config, spm_fig, Paragraphs{iField}, OutputFields{iField}, fontsize, NewLine);
end

PrintFile = 'xASL_Report';
PrintPath = fullfile(PrintDir, PrintFile); 

% Spm print adds _001 to the name of the file, I havent figured out how to get rid of this yet(ugly);
xASL_delete(strcat(PrintPath, '_001.pdf')); 
spm_figure('Print',spm_fig, PrintPath, 'pdf');

clear ax spm_fig

end


%% ========================================================================
%% ========================================================================
function sResult = xASL_qc_CheckConfigStruct(ConfigStruct, subModule, QCParm, QCOption, AlternativeString)
    % Checks if value exists in config struct, then compares it 
    if nargin<5 || isempty(AlternativeString)
        AlternativeString= '';
    end

    if isfield(ConfigStruct.(subModule).(QCParm), (QCOption)) 
        sResult = ConfigStruct.(subModule).(QCParm).(QCOption);
    else
        sResult = AlternativeString;
    end
end

function PatientInfo = xASL_adm_readParticipantsTSV(x)
    %Check if participant data exists.
    ExistsParticipants = xASL_adm_GetFileList(x.D.ROOT, 'participants.tsv');
    if ~isempty(ExistsParticipants)
        sParticipants = xASL_tsvRead(ExistsParticipants{1});
        nParticipants = size(sParticipants, 1);
        for iPar=2:nParticipants
            if sParticipants{iPar, 1} == x.SUBJECT
                %Extract Participants table, and participants info
                PatientInfo = [sParticipants(1,:); sParticipants(iPar, :)];
            end
        end
    else
        PatientInfo = {};
    end
end

function graph = xASL_Average_Intensity(x)
    if ~xASL_exist(x.D.ROOT); return; end
    
    graph = [];
end

function SlabPlaced = xASL_vis_PlotSlab(x)
    if ~xASL_exist(x.D.ROOT); return; end
    Dlist = xASL_adm_GetFileList(x.D.ROOT,'^.*$','List',[0 Inf], true);
    Flist = xASL_adm_GetFileList(x.D.ROOT,'^.nii*$','List',[0 Inf], false);
    bAnat = false;
    NiftiList = {};

    %Load T1W image
    if contains(fieldnames(Flist), 'T1*')
        bAnat = true;
        NiftiList{1,1} = fullfile(x.D.ROOT, 'T1.nii.gz');
    end

    if contains(fieldnames(Flist), 'T2*')
        bAnat = true;
        NiftiList{end,1} = fullfile(x.D.ROOT, 'T2.nii.gz');
    end

    if contains(fieldnames(Flist), 'FLAIR*')
        bAnat = true;
        NiftiList{end,1} = fullfile(x.D.ROOT, 'Flair.nii.gz');
    end

    if ~(bAnat)
        warning('No anatomical scans found, cannot plot orientation.');
    end

    %Load ASL image
    if contains(fieldnames(Dlist), 'ASL')
        
    end

    xASL_qc_PrintOrientation(NiftiList, x.D.ROOT, '');
    SlabPlaced = [];
end

function line = PrintParagraph(x, config, figure, Contents, ModName, fontsize, line)
    % Printing header/title
    line = PrintBold(ModName, figure, line, fontsize);    
    for ii=1:size(Contents,2) % iterating over text lines

        TempName = Contents(ii).name;
    
        % Prepare text to print
        if ~isstruct(Contents(ii).value) % no extra layer
            TempValue  = Contents(ii).value;
            if isnumeric(TempValue) || islogical(TempValue)
                TempValue = num2str(TempValue, 4); 
            end

        else % add extra layer
            TempFields = fieldnames(x.Output.(ModName).(Contents(ii).name));
            TempValue = '';
            for iK=1:length(TempFields)
                TempName = [TempName  '/' TempFields{iK}];
                TV2 = x.Output.(ModName).(Contents(ii).name).(TempFields{iK});
                if  isnumeric(TV2)
                    TV2 = num2str(TV2); 
                end
                TempValue = [TempValue ' ' TV2];
            end
        end

        if sum(strcmp(fieldnames(config.(ModName)),TempName))
            % Print only if visibility exists and is true. 
            if isfield(config.(ModName).(TempName), "Visible")  && config.(ModName).(TempName).Visible 

                % Define the variables from Output
                TempRange = xASL_qc_CheckConfigStruct(config, ModName, TempName, "Range");
                TempUnit  = xASL_qc_CheckConfigStruct(config, ModName, TempName, "Unit");
                TempAlias = xASL_qc_CheckConfigStruct(config, ModName, TempName, "Alias", TempName);

                if size(TempValue, 1) == 1
                    TextString = sprintf([sprintf('%-20s',[TempAlias,':']), sprintf('%8s',TempValue),sprintf('%-12s', [TempUnit, ' ' , TempRange]), ' \n']);
                end
                line = PrintText(TextString, figure, line, fontsize);
            end
        end
    end
end



function line = PrintBold(String, figure, line, fontsize)
    line = NewLine(line);
    ax=axes('Position', line ,'Visible','off','Parent',figure);
    text(0,1, String, 'FontSize', fontsize, 'FontWeight', 'bold', 'Interpreter', 'none', 'Parent', ax, 'VerticalAlignment','top');
end

function line = PrintText(String, figure, line, fontsize)
    line = NewLine(line);
    ax=axes('Position', line ,'Visible','off','Parent',figure);
    text(0,1, String, 'FontSize', fontsize, 'Interpreter', 'none', 'Parent', ax, 'VerticalAlignment','top');
end

function line = NewLine(line);
    line(2) = line(2) - 0.016;
    if line(2) < 0 
        warning('No space left on page!');
    end
end

function PrintImage(ImagePath, figure, Canvas)
    ax=axes('Position', Canvas, 'Visible', 'off', 'Parent', figure);
    [img, ~, alphachannel] = imread(ImagePath);
    fg= imshow(img);
    fg.AlphaData=alphachannel;
    clear fg
end