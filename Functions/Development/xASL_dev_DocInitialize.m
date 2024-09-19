function xASL_dev_DocInitialize(baseOutputFolder)
%xASL_dev_DocInitialize This function prepares the documentation for deployment
% functions.
%
% FORMAT:       xASL_dev_DocInitialize
% 
% INPUT:        baseOutputFolder     - Folder where the generated markdown files are stored (OPTIONAL, DEFAULT = ExploreASL/Documentation repository)
%
% OUTPUT:       n/a
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function generates all markdown files, which are necessary for the mkdocs documentation - i.e. it initializes the MD files needed
%               for deploying the documentation. It copies the manually edited MD files and crawls through function headers to ready everything for the
%               deployment
%
%               1. Administration
%               2. Copy MD-files from ExploreASL source-code
%               3. Copy the manually edited MD-files 
%               4. Crawl through function headers to generate reference manual
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_dev_DocInitialize(fullfile(x.opts.MyPath,'Development','Documentation_GitHub'));
% __________________________________
% Copyright (c) 2015-2022 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________



    %% 1. Administration

    % Initialize ExploreASL
    x = ExploreASL_Initialize;
    
    % Reminder on the structures of ASL source code headers
    xASL_adm_BreakString('REMINDER','=',[],1);
    fprintf('A correct ExploreASL header should include the following tags in the correct order:\n');
    fprintf('"FORMAT:", "INPUT:", "OUTPUT:", "DESCRIPTION:" and "EXAMPLE:"\n');
    xASL_adm_BreakString('','=',[],1);
    
    % Try to find the ExploreASL/Documentation repository on the same folder level
    [adminPath,~,~] = fileparts(mfilename('fullpath'));
    [functionsPath,~,~] = fileparts(adminPath);
    potentialDocumentationDirectory = strrep(functionsPath,fullfile('ExploreASL','Functions'),'Documentation');
    if ~isempty(xASL_adm_GetFileList(potentialDocumentationDirectory))
        fprintf('Documentation repository found...\n');
        baseOutputFolder = potentialDocumentationDirectory;
    end
    
    % Ask user for Documentation repository
    if nargin<1 && isempty(baseOutputFolder)
        if ~usejava('desktop') || ~usejava('jvm') || ~feature('ShowFigureWindows')
            baseOutputFolder = input('Insert Documentation directory: ');
        else
            baseOutputFolder = uigetdir('Select Documentation repository...');
        end
    end
    
    % Set up output folder
    outputFolder = fullfile(baseOutputFolder,'docs');
    templatesDir = fullfile(baseOutputFolder,'source');
    
	%% 2. Copy MD-files from ExploreASL source-code
	
    % Copy CHANGES.md
    xASL_Copy(fullfile(x.opts.MyPath,'CHANGES.md'),fullfile(outputFolder,'Changes.md'),1);
    
    % Copy and modify the index README
    xASL_Copy(fullfile(x.opts.MyPath,'README.md'),fullfile(outputFolder,'index.md'),1);

	% Update figure links in these files
	
	% Logo
    swapTextInFile(fullfile(outputFolder,'index.md'),...
                  '(Design/ExploreASL_logoHeader.png)',...
                  '(./img/ExploreASL_logoHeader.png)');
    % Workflow
    swapTextInFile(fullfile(outputFolder,'index.md'),...
                  '(Design/WorkflowUpdate.png "Workflow of ExploreASL")',...
                  '(./img/ExploreASL_Workflow.png "Workflow ExploreASL")');
	
	%% 3. Copy the manually edited MD-files
	% from the Documentation/source directory where they can be edited to the Documentation/docs directory, 
	% where they are ready for the deployment.
	
    % Copy the ProcessingParameters
    xASL_Copy(fullfile(templatesDir,'ProcessingParameters.md'),fullfile(outputFolder,'ProcessingParameters.md'),1);
                  
    % Copy the REQUIREMENTS file
    xASL_Copy(fullfile(templatesDir,'Requirements.md'),fullfile(outputFolder,'Requirements.md'),1);
    
    % Copy the ABOUT file
    xASL_Copy(fullfile(templatesDir,'About.md'),fullfile(outputFolder,'About.md'),1);
    
    % Copy the TUTORIALS files
    xASL_Copy(fullfile(templatesDir,'FAQ.md'),fullfile(outputFolder,'FAQ.md'),1);
    xASL_Copy(fullfile(templatesDir,'Tutorials-Install.md'),fullfile(outputFolder,'Tutorials-Install.md'),1);
	xASL_Copy(fullfile(templatesDir,'Tutorials-Execution.md'),fullfile(outputFolder,'Tutorials-Execution.md'),1);
	xASL_Copy(fullfile(templatesDir,'Tutorials-Import.md'),fullfile(outputFolder,'Tutorials-Import.md'),1);
	xASL_Copy(fullfile(templatesDir,'Tutorials-Processing.md'),fullfile(outputFolder,'Tutorials-Processing.md'),1);
	xASL_Copy(fullfile(templatesDir,'Tutorials-QC.md'),fullfile(outputFolder,'Tutorials-QC.md'),1);
	xASL_Copy(fullfile(templatesDir,'Tutorials-Developer.md'),fullfile(outputFolder,'Tutorials-Developer.md'),1);	
	
	%% 4. Crawl through function headers to generate reference manual
    % Create the functions markdown file
    xASL_dev_DocCrawler(fullfile(x.opts.MyPath,'Functions'), fullfile(outputFolder,'Functions.md'),'Functions');
    
    % Create the functions markdown file
    xASL_dev_DocCrawler(fullfile(x.opts.MyPath,'External','SPMmodified','xASL'), fullfile(outputFolder,'SPMxASL.md'),'SPMxASL');

    % Convert and copy lincense file
    convertLicenseToMarkdown(fullfile(x.opts.MyPath,'LICENSE-EXPLOREASL'),fullfile(outputFolder,'License.md'));
    
    % Use documentation crawler for modules
    xASL_dev_DocCrawler({...
                        fullfile(x.opts.MyPath,'Modules','xASL_module_Import.m'),...
                        fullfile(x.opts.MyPath,'Modules','xASL_module_Structural.m'),...
                        fullfile(x.opts.MyPath,'Modules','xASL_module_ASL.m'),...
                        fullfile(x.opts.MyPath,'Modules','xASL_module_Population.m')...
                        },...
                        fullfile(outputFolder,'Modules.md'),'Modules');
    
    % Use documentation crawler for submodules
    xASL_dev_DocCrawler(fullfile(x.opts.MyPath,'Modules','Module_Import'), fullfile(outputFolder,'Import_Module.md'),'ImportModule');
    xASL_dev_DocCrawler(fullfile(x.opts.MyPath,'Modules','Module_Processing','Module_Structural'), fullfile(outputFolder,'Structural_Module.md'),'StructuralModule');
    xASL_dev_DocCrawler(fullfile(x.opts.MyPath,'Modules','Module_Processing','Module_ASL'), fullfile(outputFolder,'ASL_Module.md'),'ASLModule');
    xASL_dev_DocCrawler(fullfile(x.opts.MyPath,'Modules','Module_Processing','Module_Population'), fullfile(outputFolder,'Population_Module.md'),'PopulationModule');
    

end

%% Function to swap text segments of documentation files
function swapTextInFile(filePath,text2swap,newText,onlyFirst)
    changeText = true;
    % Check nargin
    if nargin < 4
        onlyFirst = false;
    end
    % Open file
    file_id=fopen(filePath,'r');
    text_cell=cell(1);
    while 1
        text_line_read=fgetl(file_id);
        if text_line_read == -1
            break
        else
            text_cell(end,1)=cellstr(text_line_read);
            text_cell(end+1,1)=cell(1);
        end
    end
    fclose(file_id);
    
    %% Change text
    for curLine=1:numel(text_cell)
        if ~isempty(text_cell{curLine,1})
            if ~isempty(contains(text_cell{curLine,1},text2swap))
                curLineText = char(text_cell{curLine,1});
                if changeText && contains(curLineText,text2swap)
                    text_cell{curLine,1}= strrep(curLineText,text2swap,newText);
                    if onlyFirst
                        changeText = false;
                    end
                end
            end
        end
    end
    
    %% Change file
    file_id=fopen(filePath,'w');
    for i=1:length(text_cell)    
        fprintf(file_id,'%s\n', text_cell{i});
    end
    fclose(file_id);

end

%% Convert license file to markdown file
function convertLicenseToMarkdown(filePath,newPath)
    %% Open file
    file_id=fopen(filePath,'r');
    text_cell=cell(1);
    while 1
        text_line_read=fgetl(file_id);
        if text_line_read == -1
            break
        else
            text_cell(end,1)=cellstr(text_line_read);
            text_cell(end+1,1)=cell(1);
        end
    end
    fclose(file_id);

    %% Markdown styling
    for curLine=1:numel(text_cell)
        if ~isempty(text_cell{curLine,1})
            % Convert "ExploreASL" to bold text
            if ~isempty(contains(text_cell{curLine,1},'ExploreASL'))
                text_cell{curLine,1}= strrep(char(text_cell{curLine,1}),'ExploreASL','**ExploreASL**');
            end
            % Convert "the Software" to bold text
            if ~isempty(contains(text_cell{curLine,1},'"the Software'))
                text_cell{curLine,1}= strrep(char(text_cell{curLine,1}),'the Software','**the Software**');
            end
            % Convert "the University" to bold text
            if ~isempty(contains(text_cell{curLine,1},'the University'))
                text_cell{curLine,1}= strrep(char(text_cell{curLine,1}),'the University','**the University**');
            end
        end
    end

    %% Save new file
    file_id=fopen(newPath,'w');
    for i=1:length(text_cell)    
        fprintf(file_id,'%s\n', text_cell{i});
    end
    fclose(file_id);

end
