function reportFile = xASL_report_Main(x)
%xASL_report_Main Create the population based multi-page xASL report PDF
%
% FORMAT:       reportFile = xASL_report_Main(x)
% 
% INPUT:        x          - ExploreASL x struct (STRUCT, REQUIRED)
%
% OUTPUT:       reportFile - Path to exported PDF
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Create the population based multi-page xASL report PDF.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      reportFile = xASL_report_Main(x);
%
% __________________________________
% Copyright 2015-2021 ExploreASL

    % Get template and output path
    reportTemplate = fullfile(x.opts.MyPath, 'Functions', 'xASL_report_Template.m');
    populationDir = fullfile(x.dir.xASLDerivatives, 'Population');
    htmlDir = fullfile(populationDir, 'html');
    reportPDFhtml = fullfile(htmlDir, 'xASL_report.pdf');
    reportOutput = fullfile(populationDir, 'xASL_report.m');
    reportPDFpop = fullfile(populationDir, 'xASL_report.pdf');
    
    % Change to population directory
    cd(populationDir);
    
    % Read template and replace absolute path placeholders
    textArray = xASL_report_ReplacePlaceholders(reportTemplate, reportOutput, populationDir);
    
    % Copy the logo
    xASL_Copy(fullfile(x.opts.MyPath,'Design','ExploreASL_logoHeader.png'),fullfile(htmlDir,'logo.png'));
    
    % Add CBFs
    textArray = xASL_report_AddImageRow(populationDir,htmlDir,reportOutput,'qCBF_sub-',textArray,'CBF Images');
    
    % Add M0s
    textArray = xASL_report_AddImageRow(populationDir,htmlDir,reportOutput,'M0_sub-',textArray,'M0 Images');
    
    % Add T1ws
    textArray = xASL_report_AddImageRow(populationDir,htmlDir,reportOutput,'rT1_sub-',textArray,'T1w Images');
    
    % Add FLAIRs
    textArray = xASL_report_AddImageRow(populationDir,htmlDir,reportOutput,'rFLAIR_sub-',textArray,'FLAIR Images');
    
    % Add Slice Gradients
    textArray = xASL_report_AddImageRow(populationDir,htmlDir,reportOutput,'SliceGradient_sub-',textArray,'Slice Gradients');
    
    % Add PV GM Images
    textArray = xASL_report_AddImageRow(populationDir,htmlDir,reportOutput,'PV_pGM_sub-',textArray,'PV GM Images');
    
    % Add PV WM Images
    textArray = xASL_report_AddImageRow(populationDir,htmlDir,reportOutput,'PV_pWM_sub-',textArray,'PV WM Images');
    
    % Add SD Images
    textArray = xASL_report_AddImageRow(populationDir,htmlDir,reportOutput,'SD_sub-',textArray,'SD Images');
    
    % Add SNR Images
    textArray = xASL_report_AddImageRow(populationDir,htmlDir,reportOutput,'SNR_sub-',textArray,'SNR Images');
    
    
    % Create the report
    reportFile = publish(reportOutput, 'pdf');
    xASL_delete(reportOutput);
    
    % Move report
    xASL_Move(reportPDFhtml,reportPDFpop,1);
    
    % Remove html directory
    xASL_delete(htmlDir,1);
    
    % Change back to xASL directory
    cd(x.opts.MyPath);


end


%% Replace placeholders
function textArray = xASL_report_ReplacePlaceholders(reportTemplate, reportOutput, populationDir)
    
    % Read text file
    textArray = xASL_io_ReadTextFileLineByLine(reportTemplate);
    
    % Get placeholder
    placeholder = 'absolute_path_to_html_dir/';
    
    % Get absolute path
    absolutePath = [fullfile(populationDir,'html'),filesep];
    
    % Iterate over lines
    for iLine = 1:numel(textArray)
        % Get current line
        thisLine = textArray{iLine};
        if ~isempty(regexpi(thisLine,placeholder))
            textArray{iLine} = strrep(thisLine,placeholder,absolutePath);
        end
    end
    
    % Write the text to a PDF line by line
    xASL_report_WriteTextLineByLine(reportOutput,textArray);

end


%% Write the text to a PDF line by line
function xASL_report_WriteTextLineByLine(outputPath,textArray)

    % Write new array to output file
    fid = fopen(outputPath, 'w');
    for iLine = 1:numel(textArray)
        fprintf(fid, '%s\n', textArray{iLine});
    end
    fclose(fid);

end


%% Add a text section
function textArray = xASL_report_AddTextSection(textArray,newLines,outputPath)

    textArray = vertcat(textArray,newLines);
    xASL_report_WriteTextLineByLine(outputPath,textArray)

end


%% Create a horizontal figure and export the path
function figurePath = xASL_report_CreateHorizontalFigure(filepathsImages,htmlDir,name,filePrefix)

    % Number of images
    numImages = 4;

    % Create the figure
    fig = figure('visible', 'off', 'Position', [10 10 4*250 250]);

    %% Iterate over file paths
    
    % Currently we only print up to four images
    for iFile=1:numImages
        
        % Get file name
        if iFile<=numel(filepathsImages)
            thisFilePath = filepathsImages{iFile};
            [~,thisFile] = xASL_fileparts(thisFilePath);
            imageSlice = xASL_report_ExplortMiddleSlice(thisFilePath);
        else
            thisFile = '';
            imageSlice = zeros(10,10);
        end
        % Limits of the image
        limitMax = max(imageSlice(:));
        limitMin = min(imageSlice(:));
        if limitMin<limitMax
            limitsImage = [limitMin limitMax];
        else
            limitsImage = [0 100];
        end
        % Custom colormaps
        if ~isempty(regexpi(filePrefix,'cbf'))
            limitsImage = [0 150];
        elseif ~isempty(regexpi(filePrefix,'snr'))
            limitsImage = [0 20];
        elseif ~isempty(regexpi(filePrefix,'pv_'))
            limitsImage = [0 1];
        end
        % Create the subplot
        ax(iFile) = subplot(1,numImages,iFile);
        imagesc(rot90(imageSlice),limitsImage);
        set(gca,'XColor', 'none','YColor','none')
        if iFile<=numel(filepathsImages)
            if ~isempty(regexpi(filePrefix,'cbf'))
                colormap(ax(iFile),'hot');
            elseif ~isempty(regexpi(filePrefix,'slicegradient'))
                colormap(ax(iFile),'jet');
            else
                colormap(ax(iFile),'gray');
            end
        else
            colormap(ax(iFile),'white');
        end
        titleName = strrep(thisFile,filePrefix,'');
        title(titleName);
        % Print colorbar for last image
        if iFile==numImages
            width = 0.02;
            offset = 0.92;
            h1 = get(subplot(1,numImages,1),'Position');
            cb = colorbar(ax(1), 'Position', [offset  h1(2)  width  h1(4)]);
        end
        % Custom colorbars
        if ~isempty(regexpi(filePrefix,'cbf'))
            cb.Limits = [0 150];
        elseif ~isempty(regexpi(filePrefix,'snr'))
            cb.Limits = [0 20];
        elseif ~isempty(regexpi(filePrefix,'pv_'))
            cb.Limits = [0 1];
        end
        
    end
    
    % Save and close the figure
    figurePath = [fullfile(htmlDir,name) '.png'];
    print(fig,figurePath,'-dpng');
    close all
    

end


%% Extract middle slice (untouched)
function middleSlice = xASL_report_ExplortMiddleSlice(niftiPath)

    % Open image untouched
    niftiPathCopy = strrep(niftiPath,'.nii','_copy.nii');
    xASL_Copy(niftiPath,niftiPathCopy);
    imageNifti = xASL_io_Nifti2Im(niftiPathCopy);
    niftiPathCopy = strrep(niftiPathCopy,'.gz','');
    xASL_delete(niftiPathCopy);
    
    % Get middle z-slice
    numZ = int32(size(imageNifti,3)/2);
    middleSlice = imageNifti(:,:,numZ);

end


%% Add an image row
function textArray = xASL_report_AddImageRow(populationDir,htmlDir,reportOutput,filePrefix,textArray,imageTitle)

    filepathsImages = xASL_adm_GetFileList(populationDir,['^' filePrefix '.+'],true);
    if ~isempty(filepathsImages)
        figurePath = xASL_report_CreateHorizontalFigure(filepathsImages,htmlDir,filePrefix,filePrefix);
        newText = {['%% ' imageTitle], '% ', ['% <<' figurePath '>>'], '% '}';
        textArray = xASL_report_AddTextSection(textArray,newText,reportOutput);
    end

end


