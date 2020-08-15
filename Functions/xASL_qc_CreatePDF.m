function xASL_qc_CreatePDF(x, DoSubject)
%xASL_qc_CreatePDF Print output QC
%
% FORMAT: xASL_qc_CreatePDF(x[, DoSubject])
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   DoSubject   - index of subject to process (REQUIRED)
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
%              Further code explanation:
%              Below, using the Matlab & SPM Figure tools we create an image, which is
%              then printed to a PDF file
%              fg = the main Figure handle
%              ax = "axes" handles, these are objects containing either 1) text or 2)
%              images, with fg as "parent" (1) & (2) images have ax as "parent"
%              Positions are calculated in such a way that 4 categories can be printed,
%              which will be the first 4 fields found in x.Output
%              then allowing 8 single slice images, and 15 text lines (name & value
%              columns)
% 
% EXAMPLE: xASL_qc_CreatePDF(x);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


if ~usejava('jvm') % only if JVM loaded
    fprintf('Warning: skipping PDF report, JVM missing\n');
    return;
end

fprintf('Printing ExploreASL PDF report:   ');


%% -----------------------------------------------------------------------------------------------
%  Admin
if length(DoSubject)~=1
    warning('Too many indices provided for PDF creation! Using first subject only');
end

PathX = fullfile(x.SUBJECTDIR,'x.mat');

if ~exist(PathX, 'file')
    warning([PathX ' didnt exist, skipping xASL_qc_CreateOutputPDF']);
    return;
end
x = xASL_adm_LoadX(x, PathX, false); % assume memory x is newer than x.mat

SuSeID = [x.SUBJECTS{DoSubject(1)} '_' x.SESSIONS{1}];
iSubjSess = (DoSubject(1)-1)*x.nSessions+1;

PrintDir = fullfile(x.D.ROOT, x.SUBJECTS{DoSubject(1)});
xASL_adm_CreateDir(PrintDir);
PrintFile = ['xASL_Report_' x.SUBJECTS{DoSubject(1)} '.pdf'];
xASL_delete(PrintFile);
PrintPath = fullfile(PrintDir, PrintFile);        

try
    clear Str fg ax OutputFields StrucFields

    %% Create the Figure
    fg = spm_figure('FindWin','Graphics'); 
    set(0,'CurrentFigure',fg);

    fg = spm_figure('Create','Graphics','visible','off');

    set(fg,'windowstyle','normal'); 
    spm_figure('Clear','Graphics'); 
    % Determine font size
    if ispc
        fontsize = 8;
    elseif isunix || ismac
        fontsize = 6.5;
    end

    %% Print the title
    ax=axes('Position',[0.01 0.75 1 0.24],'Visible','off','Parent',fg);
    text(0,1,  ['xASL report: ' x.name ', ' SuSeID],...
      'FontSize',fontsize+1,'FontWeight','bold','Interpreter','none','Parent',ax);

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
            % Ensure prefix=ScanType
            if ~strcmp(OutputFields{iField},'Structural') && ~strcmp(OutputFields{iField},'SoftwareVersion')
                CheckIs = regexp(StrucFields{iField}{iV},[OutputFields{iField} '_']);
                if ~isempty(CheckIs) && CheckIs==1
                    Str{iField}(iV).name = StrucFields{iField}{iV};
                else
                    Str{iField}(iV).name = StrucFields{iField}{iV}; % avoid redundancy
%                     Str{iField}(iV).name = [OutputFields{iField} '_' StrucFields{iField}{iV}];
                end
            else
                Str{iField}(iV).name = StrucFields{iField}{iV};
            end
            Str{iField}(iV).value = x.Output.(OutputFields{iField})(Ind).(StrucFields{iField}{iV});
        end
    end 


    %% -----------------------------------------------------------------------------------------------
    %% 1) Print first column with text fields
    % Get maximal length of Str's to print
    for iField=1:length(Str)
        lengthStr(iField) = length(Str{iField}); 
    end

    % Text positions
    iY  = 0.11; % same with image positions below
    tI  = 0.055; % text line shift
    tX  = 0; % text x positions
    tX2 = 0.22; % text x positions

    %% Determine vertical starting positions for each output category
    InitialVertPos = 0.97;
    for iL=1:nOutputFields % lines per output category
        nLines(iL) = length(StrucFields{iL});
    end
    for ii=1:nOutputFields
        InitialVertPos(ii+1) = InitialVertPos(ii) - 1*nLines(ii)*tI-0.1;
    end

    for iField=1:length(Str) % iterating over categories
        VertPos = InitialVertPos(iField)-tI;
        % Printing header/title
        text(0,VertPos, OutputFields{iField} ,'FontWeight','bold','FontSize',fontsize, 'Interpreter','none','Parent',ax);
        VertPos = VertPos - tI;
        for ii=1:size(Str{iField},2) % iterating over text lines

            TempName = Str{iField}(ii).name;

            % Prepare text to print
            if ~isstruct(Str{iField}(ii).value) % no extra layer
                TempValue  = Str{iField}(ii).value;
                if  isnumeric(TempValue); TempValue = num2str(TempValue,4); end

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

            % Print the text (Name & Value columns)
            htext(1,ii,1) = text( tX,VertPos, TempName,'FontSize',fontsize, 'Interpreter','none','Parent',ax);
            htext(1,ii,2) = text(tX2,VertPos, TempValue,'FontSize',fontsize, 'Interpreter','none','Parent',ax);
            VertPos = VertPos - tI; % Vertical Position
        end
    end



    %% -----------------------------------------------------------------------------------------------
    %% 2) Print second column with image fields
    %  First synchronize the text fields & image fields
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

    % The images are square (122x122 for 1.5 mm MNI), but
    % the PDF output has a ratio 0.647, 0.16 0.11 works

    % Create image positions
    SX = 0.16; SY = 0.11; % dimensions single image slice
    iX = 0.16; % horizontal difference between images
    bY = 1-1*iY-0.004; % vertical position
    bX = 0.36; % horizontal start position

    for iI=1:9 % Define image positions in PDF
        Pos((iI-1)*8+[1:4],1:4) = [[bX:iX:bX+3*iX]' repmat(bY,4,1) repmat(SX,4,1) repmat(SY,4,1)]; 
        bY = bY-iY;
        Pos((iI-1)*8+[5:8],1:4) = [[bX:iX:bX+3*iX]' repmat(bY,4,1) repmat(SX,4,1) repmat(SY,4,1)];
        bY = bY-iY-0.01;
    end


    % Print images
    for iField=1:length(ImFieldsOrder)
        imsCurrField   = x.Output_im.(ImFieldNames{ImFieldsOrder(iField)});
        for iImage=1:min(length(imsCurrField),32) % iterate over images in this field, max 16 images allowed
            CurrentIm  = double(imsCurrField{iImage});
            % First treat the image
            if  size(CurrentIm,3)==1 % convert to color
                CurrentIm = repmat(CurrentIm,[1 1 3]);     
            end
            CurrentIm = CurrentIm+eps;
            CurrentIm = CurrentIm./max(CurrentIm(:)); % rescale

            % Then print the image at respective position
            iI=(iField-1)*16+iImage; % determine the position
            % set the ax handle
            ax=axes('Position',Pos(iI,1:4),'HandleVisibility','off','Parent',fg); % x(1) position, y(1) position, x-size, y-size
            % NB: scales with PDF size
            % NB: sizes have no effect when too close to border
            % (which is the other side with y)
            image([0.01 0.50],[0.01 0.50],CurrentIm,'HandleVisibility','off','Parent',ax); % print the image
            set(ax,'handlevisibility','off','visible','off');
        end
    end

    % Create PDF file
    xASL_delete(PrintPath);
    print(fg, '-dpdf', '-r600', PrintPath);
catch ME
    fprintf('%s\n', ['Creation of ' PrintFile ' failed:']);
    warning(ME.message);
end

fprintf('\n');
close all; % close the Matlab figure


end