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

%% Print the title
fprintf('Printing ExploreASL PDF report:   \n');

clear Str spm_fig ax OutputFields StrucFields
spm_fig = spm_figure('Create','Graphics','visible','off');
set(spm_fig,'windowstyle','normal'); 
spm_figure('Clear','Graphics');

%% Print the Title
ax=axes('Position',[0 0.96 1 0.04], 'Visible', 'off', 'Parent', spm_fig);
[img, ~, alphachannel] = imread('Design/ExploreASL_logoHeader.png');
fg= imshow(img);
fg.AlphaData=alphachannel;
clear fg
ax=axes('Position',[0.01 0.91 1 0.05],'Visible','off','Parent',spm_fig);
text(0,1, ['ExploreASL QC summary report ', SuSeID], 'FontSize',fontsize+2,'FontWeight','bold','Interpreter','none','Parent',ax, 'VerticalAlignment','top');
ax=axes('Position',[0 0 1 0],'Visible','off','Parent',spm_fig);
text(0,1, 'This report was automatically generated with ExploreASL', 'FontSize',fontsize+2,'FontWeight','bold','Interpreter','none','Parent',ax, 'VerticalAlignment','top');
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
            Str{iField}(iV).name = StrucFields{iField}{iV};
        Str{iField}(iV).value = x.Output.(OutputFields{iField})(Ind).(StrucFields{iField}{iV});
    end
end 

ImFieldNames = fieldnames(x.Output_im);
nNext = 1;
for iField=1:nOutputFields
    for iImage=1:length(ImFieldNames)
        if  strcmp(OutputFields{iField},ImFieldNames{iImage})
            ImFieldsOrder(nNext) = iImage;
            nNext = nNext+1;
        end
    end
end

for iField=1:length(Str) % iterating over categories
    % Printing header/title
    ax=axes('Position',[(1-iField/length(Str)) 0.89 (1/length(Str)) 0.05],'Visible','off','Parent',spm_fig);
    text(0,1,[OutputFields{iField}], 'FontSize',fontsize+1,'FontWeight','bold','Interpreter','none','Parent' , ax ,'VerticalAlignment','top'); 

    sPrintArray = [];
    for ii=1:size(Str{iField},2) % iterating over text lines

        TempName = Str{iField}(ii).name;
    
        % Prepare text to print
        if ~isstruct(Str{iField}(ii).value) % no extra layer
            TempValue  = Str{iField}(ii).value;
            if isnumeric(TempValue) || islogical(TempValue)
                TempValue = num2str(TempValue, 4); 
            end

        else % add extra layer
            TempFields = fieldnames(x.Output.(OutputFields{iField}).(Str{iField}(ii).name));
            TempValue = '';
            for iK=1:length(TempFields)
                TempName = [TempName  '/' TempFields{iK}];
                TV2 = x.Output.(OutputFields{iField}).(Str{iField}(ii).name).(TempFields{iK});
                if  isnumeric(TV2)
                    TV2 = num2str(TV2); 
                end
                TempValue = [TempValue ' ' TV2];
            end
        end

        f = '';
        if sum(strcmp(fieldnames(config.(OutputFields{iField})),TempName))
            % Print only if visibility exists and is true. 
            if isfield(config.(OutputFields{iField}).(TempName), "Visible")  && config.(OutputFields{iField}).(TempName).Visible 

                % Define the variables from Output
                TempRange = xASL_qc_CheckConfigStruct(config, OutputFields{iField}, TempName, "Range");
                TempUnit  = xASL_qc_CheckConfigStruct(config, OutputFields{iField}, TempName, "Unit");
                TempAlias = xASL_qc_CheckConfigStruct(config, OutputFields{iField}, TempName, "Alias", TempName);

                if size(TempValue, 1) == 1
                    f = sprintf([sprintf('%-20s',[TempAlias,':']), sprintf('%8s',TempValue),sprintf('%-12s', [TempUnit, ' ' , TempRange]), ' \n']);
                end
            end
        end
    sPrintArray = [sPrintArray, f];
    end
    text(0, 0.70, sPrintArray , 'FontName', 'FixedWidth','FontSize',fontsize,'Interpreter','none','Parent',ax, 'VerticalAlignment','top');
end

%% -----------------------------------------------------------------------------------------------
%% 2) Print bottom half with image fields

% Print images
% ====== WIP to be replaced from manual image generation =========================================
for iField=1:length(ImFieldsOrder)
    imsCurrField   = x.Output_im.(ImFieldNames{ImFieldsOrder(iField)});
    for iImage=1:min(length(imsCurrField),32) % iterate over images in this field, max 16 images allowed
        CurrentIm  = double(imsCurrField{iImage});
        % First treat the image
        if  size(CurrentIm,3)==1 % convert to color
            CurrentIm = repmat(CurrentIm,[1 1 3]);     
        end

        CurrentIm = CurrentIm./max(CurrentIm(:)); % rescale

        % Then print the image at respective position
        iRows = ceil(sqrt(length(imsCurrField)));
        Spacing = 1/length(Str)/iRows;
        % set the ax handle, scales of PDF size and structure amount
        ax=axes('Position',[((1-iField/length(Str))+(mod(iImage-1,iRows)*Spacing)) (0.03 + Spacing*floor((iImage-1)/iRows)) Spacing Spacing],'HandleVisibility','off','Parent',spm_fig); % x(1) position, y(1) position, x-size, y-size
        image([0.01 0.50],[0.01 0.50],CurrentIm,'HandleVisibility','off','Parent',ax); % print the image
        set(ax,'handlevisibility','off','visible','off');
    end
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