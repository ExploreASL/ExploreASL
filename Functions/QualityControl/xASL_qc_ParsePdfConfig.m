function [settingsPDF] = xASL_qc_ParsePdfConfig(layoutStructure, x, currentFigure, line, settingsPDF)
% xASL_qc_ParsePdfConfig function used by xASL_qc_GenerateReport to parse the configuration file loaded by xASL_qc_LoadPdfConfig.
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
%   settingsPDF         - structure containing all settings used how to print information, used in recursive calls.
%
% OUTPUTFILE:
%   xASL_Report_SubjectName.pdf - printed PDF rapport containing QC images & values
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  xASL_qc_ParsePdfConfig function used by xASL_qc_GenerateReport to parse the configuration JSON file loaded by xASL_qc_LoadPdfConfig,
%               which it will use to generate a PDF report.
%
%               This function will loop recursively over all elements/fields specified to in the JSON file,
%               as recursive calls to this same function, as each subfunction manages what to do with the individual nested json elements (fields).
%               These elements include text, images, scans and other content, and add them iteratively to a figure, that will be printed to a PDF file.
%
%               The configuration json file should contain a structure with these fields (elements) as specified in the manual.
%
%               This function is separated in the following subfunctions dealing with general content, specific images or text, and general settings:
%
% ============================================================Content Parsing Functions===============================================================
% xASL_qc_ParsePdfConfig_sub_parseContent       call appropriate function to print the content
%                                               content can be text, QCvalues, image2D, image3D, patients
%                                               or module, page, block, textSettings
% xASL_qc_ParsePdfConfig_sub_printPage          initialize a page, loop its content, and print it to a file
% xASL_qc_ParsePdfConfig_sub_printBlock         print a block, i.e. a content layout that is used multiple times
% xASL_qc_ParsePdfConfig_sub_createNewCanvas    define the positioning and scaling of a block 
%
% ==============================================================Image Based Functions=================================================================
% xASL_qc_ParsePdfConfig_sub_PrintImage             Add an image file to the PDF
% xASL_qc_ParsePdfConfig_sub_printQCImages          Add an ExploreASL QC image to the PDF
% xASL_qc_ParsePdfConfig_sub_PrintScan              Add a NIfTI file as images to the PDF
% xASL_qc_ParsePdfConfig_sub_getSliceFromStruct     Manage slices from the NIfTI file
% xASL_qc_ParsePdfConfig_sub_PrintHeader            Add a header underneath an image that is added
% xASL_qc_ParsePdfConfig_sub_WildcardReplace        Manage wildcards (e.g. <SUBJECT>) for adding an image
%
% ==============================================================Text Based Functions==================================================================
% xASL_qc_ParsePdfConfig_sub_PrintText              Add text to the PDF
% xASL_qc_ParsePdfConfig_sub_PrintPatient           Add patient information from participants.tsv
% xASL_qc_ParsePdfConfig_sub_PrintQC                Add key-values from x.Output
% xASL_qc_ParsePdfConfig_sub_Generate_QC_String     Manage the key, value, unit, range from x.Output
% xASL_qc_ParsePdfConfig_sub_PaddedString           Manage the text position by padding empty strings
%
% ==============================================================Settings Based Functions==============================================================
% xASL_qc_ParsePdfConfig_sub_loadSettings       Load configuration JSON settings as string or number to settingsPDF
% xASL_qc_ParsePdfConfig_sub_defaultSettings    Set default settings for the PDF report creation
%
% EXAMPLE: xASL_qc_ParsePdfConfig(layoutStructure, x);
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


%% ====================================================================================================================================================
% Admin
if nargin < 2 
    error('At least 2 input arguments required');
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

    currentField = layoutStructure.(fields{iField});

    if strcmp(fields{iField}, 'content') || strcmp(fields{iField}, 'pages') || strcmp(fields{iField}, 'modules')
        % For content, pages, or modules, we need to loop
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
%
% xASL_qc_ParsePdfConfig_sub_parseContent       call appropriate function to print the content
%                                               content can be text, QCvalues, image2D, image3D, patients
%                                               or module, page, block, textSettings
% xASL_qc_ParsePdfConfig_sub_printPage          initialize a page, loop its content, and print it to a file
% xASL_qc_ParsePdfConfig_sub_printBlock         print a block, i.e. a content layout that is used multiple times
% xASL_qc_ParsePdfConfig_sub_createNewCanvas    define the positioning and scaling of a block


%% ====================================================================================================================================================
function [settingsPDF, line] = xASL_qc_ParsePdfConfig_sub_parseContent(currentField, x, currentFigure, line, settingsPDF)
% This function parses the content of the json file, and calls the appropriate function to print the content.

    % It first checks if the currentField is a struct, and if it contains a field "type" which specifies the type of content to be printed.
    % If the currentField is not a struct, it will ignore the currentField and the next iteration will start.
    if ~isstruct(currentField)
        return;
    end

    if ~isfield(currentField, 'type') || ~isfield(currentField, 'category')
        error([currentField ' contains no field specifying the category and type of content, add field "category": "content" or "type": "image2D" to this json field for example']);
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
                case 'module'
                    % print multiple pages
                    for iPage = 1:length(currentField.content)
                        xASL_qc_ParsePdfConfig_sub_printPage(currentField.content(iPage), x, settingsPDF);
                    end

                case 'page'
                    xASL_qc_ParsePdfConfig_sub_printPage(currentField, x, settingsPDF);  
                case 'block'
                    xASL_qc_ParsePdfConfig_sub_printBlock(currentField, x, currentFigure, settingsPDF);
                case 'textSettings'
                    settingsPDF = xASL_qc_ParsePdfConfig_sub_loadSettings(currentField, settingsPDF);
            end  
    end
end


%% ====================================================================================================================================================
function  xASL_qc_ParsePdfConfig_sub_printPage(pageStruct, x, settingsPDF)
% This function prints pages using the layout defined in the json file.
% It first creates a new figure, and then iterates over and prints all content in the pageStruct
% The pageStruct should contain a field "content" which contains all content to be printed on the page.
% The pageStruct should also contain a field "identifier" which is used to name the printed PDF file.

    %% Create the PDF figure defaults
    figPrimary = figure('visible', 'off', 'Units', 'centimeters', 'Position', [0 0 21 29.7]);
    ax = axes('Position', [0 0 1 1], 'Visible', 'off', 'Parent', figPrimary);

    %% Print the Title
    pathLogo = fullfile(x.opts.MyPath, 'Design', 'ExploreASL_logoHeader.png');
    xASL_qc_ParsePdfConfig_sub_PrintImage(pathLogo, [], figPrimary, settingsPDF, [0 0.96 1 0.04]); % [0 0.96 1 0.04]

    %% Print the date and time
    stringDateTime = char(datetime('now','TimeZone','local','Format','d MMMM y, HH:mm:ss'));
    xASL_qc_ParsePdfConfig_sub_PrintText(stringDateTime, figPrimary, [0 0.99 1 0], settingsPDF); % [0 0.96 1 0.08]

    %% Print the Footer
    % xASL_qc_ParsePdfConfig_sub_PrintText('This report was automatically generated by ExploreASL', figPrimary, [0 0.02 1 0], settingsPDF);

    %% Parse pageStruct and create the page as defined in the configuration json file
    xASL_qc_ParsePdfConfig(pageStruct, x, figPrimary, [0 0.93 1 0], settingsPDF);

    % Finally it prints the page to a PDF file using the identifier as filename in the subject directory.
    fileName = ['xASL_Report_' x.SUBJECT '_' pageStruct.identifier];
    printPathFile = fullfile(x.dir.xASLDerivatives, x.SUBJECT, fileName);
    xASL_delete(printPathFile);
    fprintf('%s\n', ['Printing ' fileName '.pdf']);
    print(figPrimary, printPathFile, '-dpdf', '-fillpage'); % -bestfit

end


%% ====================================================================================================================================================
function xASL_qc_ParsePdfConfig_sub_printBlock(currentField, x, currentFigure, settingsPDF)
% This function prints a content blocks using the layout defined in the json file.
% Content blocks can be used when you have a predefined layout that you want to use multiple times.
% The blockStruct should contain a field "content" which contains all content to be printed in the block.
% Using the function xASL_qc_ParsePdfConfig_sub_createNewCanvas it will create a new subcanvas for the block to be printed in.

    position = xASL_str2num(currentField.position);
    size = xASL_str2num(currentField.size);
    [settingsPDF.canvas, line] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position, size, settingsPDF.canvas);
    xASL_qc_ParsePdfConfig(currentField, x, currentFigure, line, settingsPDF);
end


%% ====================================================================================================================================================
function [newCanvas, line] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position, size, oldCanvas)
% This function creates a new canvas based on the old canvas, and the position and size of the new canvas.
% The canvas is used to define the position and size where content is printed.
% This is for example used to define where a new image, scan or block is printed on the page.
% This function ensures that the new canvas is within the old canvas, and scaled appropriately.
% This is for example important when a predefined block contains images, which need to be scaled appropriately.
    
    % A canvas is defined as a 4 element vector [x y width height] where x and y are the lower left corner of the canvas.
    % This function also returns the variable line, which is the position of the next line to be printed.
    newPosition = [oldCanvas(1) + position(1) * oldCanvas(3) oldCanvas(2) + position(2) * oldCanvas(4)];
    newSize = oldCanvas(3:4) .* size;
    newCanvas = [newPosition, newSize];
    line = [newPosition(1)  newPosition(2) + newSize(2) newSize(1) 0];
end







% ====================================================================================================================================================
% ==============================================================Image Based Functions=================================================================
% ====================================================================================================================================================
%
% xASL_qc_ParsePdfConfig_sub_PrintImage             Add an image file to the PDF
% xASL_qc_ParsePdfConfig_sub_printQCImages          Add an ExploreASL QC image to the PDF
% xASL_qc_ParsePdfConfig_sub_PrintScan              Add a NIfTI file as images to the PDF
% xASL_qc_ParsePdfConfig_sub_getSliceFromStruct     Manage slices from the NIfTI file
% xASL_qc_ParsePdfConfig_sub_PrintHeader            Add a header underneath an image that is added
% xASL_qc_ParsePdfConfig_sub_WildcardReplace        Manage wildcards (e.g. <SUBJECT>) for adding an image

%% ====================================================================================================================================================
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
                warning('input needs to have either xPath, absolutePath, popPath, or subjPath') ;
                % ImagePath = xASL_qc_ParsePdfConfig_sub_WildcardReplace(ImagePath, x);
                return;
            end
    
            if isfield(imageStruct, 'header')
                header = imageStruct.header;
            end
        case 5
            ImagePath = input;
    end

    % First it calculates the size of the canvas for the image to be printed in.
    canvas = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position(1:2), position(3:4), settingsPDF.canvas);

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


%% ====================================================================================================================================================
function [settingsPDF] = xASL_qc_ParsePdfConfig_sub_printQCImages(qcStruct, x, currentFigure, settingsPDF)
    % This function prints images using the layout defined in the json file.

    if ~isfield(qcStruct, 'module') ||  ~isfield(x.Output, (qcStruct.module)) % If the module is missing
        warning(['Derivatives for PDF missing for module: ' qcStruct.module]);
        return
    elseif ~isfield(qcStruct, 'session') || qcStruct.session == "" % If we don't have sessions
        allImages = x.Output_im.(qcStruct.module);
    elseif isfield(x.Output_im.(qcStruct.module), qcStruct.session) % If we have sessions, use the session
        allImages = x.Output_im.(qcStruct.module).(qcStruct.session);
    else
        warning('No derivatives found for PDF');
        return
    end

    if ~isfield(qcStruct, 'position') ||~isfield(qcStruct, 'size')
        warning('Position and size parameters missing for PDF printing');
        return
    end
    
    % First calculate the size of the canvas for the image to be printed in.
    position = [xASL_str2num(qcStruct.position) xASL_str2num(qcStruct.size)];
    canvas = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position(1:2), position(3:4), settingsPDF.canvas);
    
    % Dimensions, at most 24 images will be printed
    imageFields = fields(allImages);
    nImages = length(imageFields);
    nImages  = min(nImages, 24); % Allow maximum of 24 images to be printed
    imPerRow = ceil(sqrt(nImages));
    imSize   = 1/imPerRow;

    % Final image position, to be used by headers below
    canvasHeader = canvas;
    xPos = mod(nImages - 1, imPerRow) * imSize;
    yPos = (imPerRow - ceil(nImages / imPerRow)) * imSize;
    position = xASL_qc_ParsePdfConfig_sub_createNewCanvas([xPos, yPos], [imSize, imSize], canvas);
    canvasHeader(2) = position(2); % the y-position is needed to print headers right below the QC images

    % Print images
    for iImage = 1:nImages
        imageName = imageFields{iImage};
        CurrentIm = allImages.(imageName);
        % This should have only a single image per imageName, per xASL_vis_AddIM2QC
        if isempty(CurrentIm)
            warning(['Something went wrong in xASL_vis_AddIM2QC, image missing: ' imageFields{iImage}]);
        elseif iscell(CurrentIm)
            warning(['Something went wrong in xASL_vis_AddIM2QC, too many images: ' imageFields{iImage}]);
        end
        CurrentIm = double(allImages.(imageFields{iImage}));

        % Convert grayscale images to color
        if  size(CurrentIm,3) == 1 
            CurrentIm = repmat(CurrentIm,[1 1 3]);     
        end

        % Rescale image to fit between 0-1
        CurrentIm = CurrentIm ./ max(CurrentIm(:)); %

        % Determine the position and size of the image
        xPos = mod(iImage - 1, imPerRow) * imSize;
        yPos = (imPerRow - ceil(iImage / imPerRow)) * imSize;
        position = xASL_qc_ParsePdfConfig_sub_createNewCanvas([xPos, yPos], [imSize, imSize], canvas);

        % Finally print the image to the current figure, and updates the figure count.
        ax = axes('Position', position , 'Visible', settingsPDF.axesVisible, 'Parent', currentFigure);
        fg = imshow(CurrentIm);

        % Print the image header and update the figure count
        header = imageName;

        % Add the NIfTI file-name
        % here, the later the filename to check is placed in the vector, 
        % the higher its priority to print its NIfTI file.
        % E.g., For 'rT1_with_rc2T1', it should print rc2T1.nii and not rT1.nii
        fileName = 'something.nii'; % default for NIfTI files that haven't been added here

        namesToCheck = {'qCBF' 'M0' 'Tex' 'ITT' 'ATT' 'SD' 'SNR' 'mean_control' 'noSmooth_M0' 'rFLAIR'...
            'rT1' 'rT1_ORI' 'rFLAIR_ORI' 'rc1T1' 'rc2T1' 'rc3T1' 'PV_pGM' 'PV_pWM' 'PV_pCSF' 'CentralWM_QC' 'rWMH_SEGM'};
        for iName=1:length(namesToCheck)
            if contains(imageName, namesToCheck{iName})
                fileName = [namesToCheck{iName} '.nii'];
            end
        end

        % Make the header a bit more human readable
        header = strrep(header, 'Tra', 'Transversal');
        header = strrep(header, 'Cor', 'Coronal');
        header = strrep(header, 'Sag', 'Sagittal');

        header = strrep(header, '_ORI', '_before lesion filling');
        header = strrep(header, '_rWMH_SEGM', '_WMH segmentation');
        header = strrep(header, '_Reg_', ' ');
        header = strrep(header, '_with_', ' overlaid with ');

        header = strrep(header, '_', ' '); % remove underscores

        header = strrep(header, 'noSmooth M0', 'M0 without smoothing');

        header = strrep(header, 'rc1T1', 'GM segmentation');
        header = strrep(header, 'rc2T1', 'WM segmentation');
        header = strrep(header, 'rc3T1', 'CSF segmentation');

        header = strrep(header, 'T1', 'T1w');
        header = strrep(header, 'rT1', 'T1');
        header = strrep(header, 'rFLAIR', 'FLAIR');
        header = strrep(header, 'qCBF', 'CBF');

        header = strrep(header, 'CentralWM QC', 'central WM QC region');
                

        % Add the filename as well
        header = [header ' (' fileName ')'];
        header = strrep(header, 'M0 (', 'M0 after smoothing (');

        settingsPDF.figureCount = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, settingsPDF, canvasHeader);

    end

end


%% ====================================================================================================================================================
function [settingsPDF] = xASL_qc_ParsePdfConfig_sub_PrintScan(scanStruct, x, currentFigure, settingsPDF)
% This function prints scans using the layout defined in the json file.

    % It first checks if the requected scan exists in the x.P structure, and if it doesnt it will throw a warning and return.
    if ~isfield(x.P, scanStruct.name)
        warning (['could not print ' scanStruct.name ', check if NIfTI exists in ExploreASL/Derivatives/Population']);
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
    x.S.ConcatSliceDims = 0;
    % Create the image from the scans defined in the json file.
    imageToPrint = xASL_vis_CreateVisualFig(x, ImIn, [], [], [], [], [], [], [], [], [], 0); % no verbosity   

    % Calculate the size of the canvas for the image to be printed in
    position = [xASL_str2num(scanStruct.position) xASL_str2num(scanStruct.size)];
    [canvas] = xASL_qc_ParsePdfConfig_sub_createNewCanvas(position(1:2), position(3:4), settingsPDF.canvas);

    % Print the image to the current figure, and update the figure count.
    ax = axes('Position', canvas, 'Visible', settingsPDF.axesVisible, 'Parent', currentFigure);
    fg = imshow(imageToPrint);

    % Print the image header and update the figure count
    settingsPDF.figureCount = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, settingsPDF, canvas);
end


%% ====================================================================================================================================================
function [slice] = xASL_qc_ParsePdfConfig_sub_getSliceFromStruct(struct, name)
% This function changes the strings in the json file to numbers, and returns the slice to be printed.
% This specifically cannot use xASL_str2num, because the way that function returns NaNs is not compatible the image generation function.

    if isfield(struct, name) && ~isempty(struct.(name))
        slice = str2num(struct.(name));
    else 
        slice = [];
    end
end


%% ====================================================================================================================================================
function [figureCount] = xASL_qc_ParsePdfConfig_sub_PrintHeader(header, currentFigure, settingsPDF, position)
% This function prints a "header" underneath the image to be printed (header == legend for the contents of the image)

    % If no header is specified, it will exit and not iterate the figure count.
    if isempty(header) || ~settingsPDF.imageHeaders
        figureCount = settingsPDF.figureCount;
        return
    else
        figureCount = settingsPDF.figureCount + 1;
    end

    % position is the canvas obtained from the QC figures to be printed
    position(2) = position(2) - figureCount*0.01; % this is the y position in the PDF, the legend should be printed right below the image
    position(4) = 0;
    settingsPDF.HorizontalAlignment = 'center';
    text = ['Figure ' num2str(figureCount) ': ' header];
    xASL_qc_ParsePdfConfig_sub_PrintText(text, currentFigure, position, settingsPDF);
end


%% ====================================================================================================================================================
function [strout] = xASL_qc_ParsePdfConfig_sub_WildcardReplace(strin, x, settingsPDF)
% This function replaces wildcards in the path to the image.
% Wildcards are defineds as <wildcard> in the json file, and are replaced with the corresponding field in the x structure.
% E.g., <SUBJECT> is replaced with the the value in the x.SUBJECT field

    strout = strin;
    substring = regexp(strin, '<\w*>', 'match');
    if  settingsPDF.BIDS_Translation
        % Replace BIDS nomenclature to ExploreASL Legacy nomenclature
        substring = strrep(substring, 'RUN', 'SESSION');
        substring = strrep(substring, 'SESSION', 'VISIT');
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
% ==============================================================Text Based Functions==================================================================
% ====================================================================================================================================================
%
% xASL_qc_ParsePdfConfig_sub_PrintText              Add text to the PDF
% xASL_qc_ParsePdfConfig_sub_PrintPatient           Add patient information from participants.tsv
% xASL_qc_ParsePdfConfig_sub_PrintQC                Add key-values from x.Output
% xASL_qc_ParsePdfConfig_sub_Generate_QC_String     Manage the key, value, unit, range from x.Output
% xASL_qc_ParsePdfConfig_sub_PaddedString           Manage the text position by padding empty strings



% ====================================================================================================================================================
function line = xASL_qc_ParsePdfConfig_sub_PrintText(input, currentFigure, line, settingsPDF)
% This function prints text using the layout defined in the json file.

    % First check the input type:
    % If the text is defined as a string, it will use the string as the text to be printed.
    % If the text is defined as a struct, it will use the fields in the struct to define the text.    
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
            % If the input cannot be parsed, it will issue an error.
        otherwise
            class(input)
            error('Could not find a string of printable text');
    end

    % It then prints the text to the current figure.
    ax = axes('Position', line , 'Visible', settingsPDF.axesVisible, 'Parent', currentFigure);
    % String = strrep(String, '_', ' '); % now we do this specifically separated for keys and values
    text(0, 0, String, 'Parent', ax, 'FontSize', settingsPDF.fontSize, 'FontWeight', settingsPDF.fontWeight, 'Color', settingsPDF.color, 'FontName', settingsPDF.fontName, 'Interpreter', 'none', 'VerticalAlignment', 'top');
    
    %% Update the line position for the next line to be printed

    % The "newline" distance has so far been hardcoded based on experience, but this can be improved in the future.
    line(2) = line(2) - (settingsPDF.lineSpacing) - (settingsPDF.fontSize * 0.001);

    % It simply checks if the line position is lower than 0, and if it is it will throw a warning that the text is printed outside the canvas.
    if line(2) < 0 
        warning('No space left on page!');
    elseif line(2) < settingsPDF.canvas(2) 
        warning('Printing outside canvas, check block settings.');
    end


end


%% ====================================================================================================================================================
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
        fprintf ('No patient information found in participants.tsv in the ExploreASL derivatives folder.\n');
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


%% ====================================================================================================================================================
function line = xASL_qc_ParsePdfConfig_sub_PrintQC(qcStruct, x, currentFigure, line, settingsPDF)
% This function prints QC values using the layout defined in the json file.
% QC Values are extracted from the x.Output structure, and printed to the PDF report in a single line.
% Certain settings can be applied to the QC values, for example a range can be specified.
% If the QC value is outside the range, it will be printed in red.
    
    if ~isfield(qcStruct, 'parameter') || ~isfield(qcStruct, 'module')
        fprintf('QC content field not properly defined in JSON, skipping printing of QC parameter. \n');
        return
    elseif strcmp(qcStruct.parameter, 'qc_images')
        xASL_qc_ParsePdfConfig_sub_printQCImages(qcStruct, x, currentFigure, settingsPDF);
        return
    end

    % Stop if module and field don't exists in output
    if ~isfield(qcStruct, 'session') || qcStruct.session == ""
        if ~isfield(x.Output, (qcStruct.module)) || ~isfield(x.Output.(qcStruct.module), qcStruct.parameter) 
            return
        end
    elseif ~isfield(x.Output, (qcStruct.module)) || ~isfield(x.Output.(qcStruct.module), qcStruct.session)  || ~isfield(x.Output.(qcStruct.module).(qcStruct.session), qcStruct.parameter) 
        return
    end

    % Use field specific settings if they exists, otherwise default to monospace for QC values. 
    settingsPDF.fontName = 'monospace';
    if isfield(qcStruct, 'textSettings')
        settingsPDF = xASL_qc_ParsePdfConfig_sub_loadSettings(qcStruct.textSettings, settingsPDF);
    end

    % Generate the string to be written in the pdf Report
    % while managing the key, value, unit, range from x.Output
    [string, settingsPDF] = xASL_qc_ParsePdfConfig_sub_Generate_QC_String(qcStruct, x, settingsPDF);

    % Print the string to the PDF report
    parametersToSkip = {'ID'};
    
    if sum(strcmp(parametersToSkip, qcStruct.parameter))>0
        % skip this parameter
    else
        line = xASL_qc_ParsePdfConfig_sub_PrintText(string, currentFigure, line, settingsPDF);
    end

end


%% ===================================================================================================================================================
function [string, settingsPDF] = xASL_qc_ParsePdfConfig_sub_Generate_QC_String(qcStruct, x, settingsPDF)
% This function creates QC lines from the key names and values
% It performs the following tasks:
% 1. Translate from the provided translation tsv
% 2. Manage key/alias
% 3. Manage unit
% 4. Manage range
% 5. Manage values
% 6. Equalize string lengths
% 7. Combine name, value, unit and range into a single string for printing


    if ~isfield(qcStruct, 'module')
        fprintf('%s\n', 'No module specified to be printed in the PDF printing, skipping');
        return;
    elseif ~isfield(x.Output, (qcStruct.module))
        fprintf('%s\n', ['Module ' qcStruct.module ' was not processed, not printing to the PDF']);
        return;
    % if we parse a module without sessions, e.g., Structural module
    elseif ~isfield(qcStruct, 'session') || qcStruct.session == "" 

        if ~isfield(qcStruct, 'parameter')
            fprintf('%s\n', 'No value specified to be printed in the PDF printing, skipping');
            return;
        elseif ~isfield(x.Output.(qcStruct.module), qcStruct.parameter) 
            fprintf('%s\n', ['Parameter ' qcStruct.parameter ' was not processed, skipping']);
            return;
        else
            TempValue = x.Output.(qcStruct.module).(qcStruct.parameter);
        end
    % if we parse a module with sessions, e.g., ASL module
    elseif isfield(x.Output.(qcStruct.module), qcStruct.session)
        if ~isfield(qcStruct, 'parameter')
            fprintf('%s\n', 'No parameters found for printing to PDF, skipping');
            return;
        elseif ~isfield(x.Output.(qcStruct.module).(qcStruct.session), qcStruct.parameter) 
            fprintf('%s\n', ['Parameter ' qcStruct.parameter ' was not processed, skipping printing to PDF']);
            return;
        else
            TempValue = x.Output.(qcStruct.module).(qcStruct.session).(qcStruct.parameter);
        end
    end


    %% 1. Translate from the provided translation tsv if that's enabled
    % replacing QC keys with long names with easier to read names, units, and range
    if settingsPDF.QC_TSV_Translations

        indexTranslation = find(strcmp(settingsPDF.QC_Translation(:,1), qcStruct.parameter));
        if ~isempty(indexTranslation)
            qcStruct.alias = char(settingsPDF.QC_Translation(indexTranslation, 2));
            qcStruct.unit  = char(settingsPDF.QC_Translation(indexTranslation, 4));
            qcStruct.range = char(settingsPDF.QC_Translation(indexTranslation, 5));
        end

    end


    %% 2. Manage key/alias
    % Check if the QC has an alias, if so replaced the parameter with the alias in the printed text.
    if ~isfield(qcStruct, 'alias') || isempty(qcStruct.alias)
        qcStruct.alias = qcStruct.parameter;
    end

    % Remove '_' only from the key, not from the value
    % (previously, this was removed from all strings)
    qcStruct.alias = strrep(qcStruct.alias, '_', ' ');


    % Ratios -> percentages
    % volumetric ratios or contrast-to-noise (e.g. GM-WM)/SD(WM) ) are usually better represented as percentages
    bRatio2Percentage = contains(qcStruct.alias, 'volume (ratio)') || contains(qcStruct.alias, 'CNR (ratio)');
    if bRatio2Percentage
        qcStruct.alias = strrep(qcStruct.alias, 'ratio', '%');
    end


    %% 3. Manage unit
    % (currently, we are not using this, but printing the unit inside the key/alias as that is a more enonomical use of the PDF space)
    % Check if the QC has a unit, if not set it to an empty string.
    if ~isfield(qcStruct, 'unit') 
        qcStruct.unit = '';
    end


    %% 4. Manage range
    % Check if the QC has a range, if so check if the value is within the range, else print in red.
    if isfield(qcStruct, 'range') && isnumeric(TempValue)
        if ~isempty(qcStruct.range)
            range = xASL_str2num(strsplit(qcStruct.range, '-'));

            if ~isequal(size(range),[1 2]) % check first that range has the correct size, should have 2 values
                nRangeValues = numel(range);
                warning(['range parameter should have 2 values in configReportPDF.json but had ' xASL_num2str(nRangeValues) ' values']);
            elseif TempValue < range(1) || TempValue > range(2) 
                % check if the value is too low or too high (i.e. outside the allowed range)
                settingsPDF.color = 'r'; % then color the value red
            end
            
            % Ratio -> percentage
            if bRatio2Percentage
                range = range.*100;
            end

            qcStruct.range = ['(' num2str(range(1)) '—' num2str(range(2)) ')'];
        end
    else
        qcStruct.range = '';
    end


    %% 5. Manage values

    % Convert the value to a string.
    if ~strcmp(qcStruct.alias, 'Version ExploreASL Git commit')
        % ensure a Git commit is not managed as a number
        TempValue = xASL_num2str(TempValue, settingsPDF.numberFormat);
    end

    % Manage cells
    if iscell(TempValue)
        TempValue = cell2mat(TempValue);
    end

    % Manage specific value cases:
    % 1. Manage subject ID:
    if strcmpi(qcStruct.alias, 'subject ID')
        TempValue = strrep(TempValue, 'sub-', '');
        TempValue = regexprep(TempValue, '_ASL_\d*', '');
    end

    numericalValue = xASL_str2num(TempValue); % outputs NaN for strings
    if strcmp(qcStruct.alias, 'Version ExploreASL Git commit')
        % print the abbreviated commit number
        % & ensure this is not managed as a number
        TempValue = TempValue(1:7);

    elseif ~isnan(numericalValue(1)) % when we have a numeric value

        % 2. Multiple identical values will only be printed once
        uniqueValues = unique(numericalValue);
        if numel(numericalValue)>1 && isscalar(uniqueValues)
            numericalValue = numericalValue(1);
        end

        % 3. Ratio -> percentages
        if bRatio2Percentage
            numericalValue = numericalValue.*100;
        end

        % 4. Print numerical values only with 4 significant digits
        for iNum=1:length(numericalValue)
            TempValue2{iNum} = xASL_adm_formatWithRoundedMagnitude(numericalValue(iNum), 4);
        end

        TempValue = strjoin(TempValue2, ' ');
    end



    %% 6. Equalize string lengths
    qcStruct.alias  = xASL_qc_ParsePdfConfig_sub_PaddedString( qcStruct.alias, 35); % 25 % this is the key
    try
        TempValue       = xASL_qc_ParsePdfConfig_sub_PaddedString( TempValue, 18, 'right'); % this is the value
    catch

        disp('piet');

    end

    UnitRange       = xASL_qc_ParsePdfConfig_sub_PaddedString( [qcStruct.unit, qcStruct.range], 20);


    %% 7. Combine name, value, unit and range into a single string for printing
    if size(TempValue, 1) == 1
        % string = sprintf([qcStruct.alias ':' TempValue ' ' UnitRange ' \n']);
        try
            string = [qcStruct.alias ':' TempValue ' ' UnitRange];
        catch
            string = 'qcStruct.alias: complicated value';
            fprintf('%s\n', ['Cannot print this parameter to PDF: ' qcStruct.alias]);
        end
    end


end


%% ===================================================================================================================================================
function resultText = xASL_qc_ParsePdfConfig_sub_PaddedString(textToPrint, textWidth, align, SymbolToFill)
% This function manages the text position by padding empty strings
    if nargin < 3 || isempty(align)
        align = 'left';
    end

    if nargin < 4 || isempty(SymbolToFill)
        SymbolToFill = ' ';
    end

    %% Create default string
    resultText = repmat(SymbolToFill, 1, textWidth);
    
    % Get string sizes
    [xSize, ySize] = size(textToPrint);

    % First check if textToPrint is really text
    if isstruct(textToPrint)
        nFields = length(fields(textToPrint));
        error(['textToPrint is a struct with ' num2str(nFields) ' fields']);
    end

    % If a value contains Two (or more) entries, add a space and concatenate them horizontally
    if xSize > 1
        textToPrint(1, ySize + 1) = ' ';  % Add space to each element
        concatenated = '';
        
        for element=1 : xSize
            concatenated = append(concatenated, textToPrint(element, :));
        end
        
        textToPrint = concatenated(1 : end-1);
        [~ , ySize] = size(textToPrint);
    end

    % If the string size is smaller than the allotted size, align the text left or right, otherwise replace final characters with elipses
    if ySize < textWidth
        if strcmp(align,'left')
            resultText(1:ySize) = textToPrint;
        elseif strcmp(align,'right')
            resultText(end - ySize+1:textWidth) = textToPrint;
        end
    else
        resultText(1:textWidth-2) = textToPrint(1:textWidth-2);
        resultText(textWidth-2:textWidth) = '...';
    end
   
end



% ====================================================================================================================================================
% ==============================================================Settings Based Functions==============================================================
% ====================================================================================================================================================
%
% xASL_qc_ParsePdfConfig_sub_loadSettings       Load configuration JSON settings as string or number to settingsPDF
% xASL_qc_ParsePdfConfig_sub_defaultSettings    Set default settings for the PDF report creation


%% ===================================================================================================================================================
function [settingsPDF] = xASL_qc_ParsePdfConfig_sub_loadSettings(json, settingsPDF)
% This function replaces existing settings with new ones from the json file.
% If the json file contains a field "fontSize" or "lineSpacing" it will convert the string to a number.
% Otherwise it will simply replace the existing setting with the new one.

    fields = fieldnames(json);
    
    for iField = 1:length(fields)
        if strcmp(fields{iField}, 'fontSize') || strcmp(fields{iField}, 'lineSpacing') || strcmp(fields{iField}, 'imageHeaders')
            settingsPDF.(fields{iField}) = xASL_str2num(json.(fields{iField}));
        else
            settingsPDF.(fields{iField}) = json.(fields{iField});
        end 
    end
end


%% ===================================================================================================================================================
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
    settingsPDF.QC_TSV_Translations = 1; 
    settingsPDF.imageHeaders = 1;
    settingsPDF.numberFormat = '%.2f';

    settingsPDF.QC_Translation = xASL_tsvRead(fullfile(x.opts.MyPath, 'Functions', 'QualityControl', 'qc_glossary.tsv'));


    if ispc
        settingsPDF.fontSize = 9; % 10
    elseif isunix 
        settingsPDF.fontSize = 6.5; % 8
    else 
        error('Could not find OS');
    end 

end