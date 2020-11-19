function xASL_adm_DocInitialize(outputFolder,updateExploreASL)
%xASL_adm_DocInitialize Script to call the separate documentation crawler
% functions.
%
% FORMAT:       xASL_adm_DocInitialize
% 
% INPUT:        outputFolder     - Folder where the generated markdown files are stored (OPTIONAL, DEFAULT = fullfile(x.MyPath,'Development','Documentation_GitHub'))
%               updateExploreASL - Update README files in ExploreASL structure (OPTIONAL, DEFAULT = true)
%
% OUTPUT:       n/a
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function generates all markdown files, which are
%               necessary for the mkdocs documentation.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_adm_DocInitialize(fullfile(x.MyPath,'Development','Documentation_GitHub'),true);
% __________________________________
% Copyright 2015-2020 ExploreASL


    %% Workflow

    % Initialize ExploreASL
    x = ExploreASL_Initialize([],0);
    
    % Define output folder
    if nargin < 1
        outputFolder = fullfile(x.MyPath,'Development','Documentation_GitHub');
    end

    if nargin < 2
        updateExploreASL = true;
    end
    
    % Copy and modify the index README
    copyfile(fullfile(x.MyPath,'README.md'),fullfile(outputFolder,'index.md'));
    % Logo
    swapTextInFile(fullfile(outputFolder,'index.md'),...
                  '(https://github.com/ExploreASL/ExploreASL/blob/develop/ExploreASL_logoSmall.png)',...
                  '(./img/title.png "ExploreASL")');
    % Workflow
    swapTextInFile(fullfile(outputFolder,'index.md'),...
                  '(https://www.researchgate.net/profile/Andrew_Robertson7/publication/337328693/figure/fig1/AS:826578854481921@1574083164220/Schematic-diagram-of-ExploreASL-processing-steps-Steps-marked-with-a-are-optional.ppm "Workflow of ExploreASL")',...
                  '(./img/ExploreASL_Workflow.jpg "Workflow ExploreASL")');
              
    % Copy the REQUIREMENTS file
    copyfile(fullfile(x.MyPath,'REQUIREMENTS.md'),fullfile(outputFolder,'Requirements.md'));
    
    % Copy the ABOUT file
    copyfile(fullfile(x.MyPath,'ABOUT.md'),fullfile(outputFolder,'About.md'));
    
    % Create the functions markdown file
    xASL_adm_DocCrawler(fullfile(x.MyPath,'Functions'), fullfile(outputFolder,'Functions.md'),'Functions');

    % Convert and copy lincense file
    convertLicenseToMarkdown(fullfile(x.MyPath,'LICENSE-EXPLOREASL'),fullfile(outputFolder,'License.md'));
    
    % Use documentation crawler for modules
    xASL_adm_DocCrawler({fullfile(x.MyPath,'ExploreASL_Import.m'),...
                        fullfile(x.MyPath,'Modules','xASL_module_Structural.m'),...
                        fullfile(x.MyPath,'Modules','xASL_module_ASL.m'),...
                        fullfile(x.MyPath,'Modules','xASL_module_Population.m')},...
                        fullfile(outputFolder,'Modules.md'),'Modules');
    
    % Use documentation crawler for submodules
    xASL_adm_DocCrawler(fullfile(x.MyPath,'ExploreASL_ImportConfig.m'), fullfile(outputFolder,'Import_Module.md'),'ImportModule');
    xASL_adm_DocCrawler(fullfile(x.MyPath,'Modules','SubModule_Structural'), fullfile(outputFolder,'Structural_Module.md'),'StructuralModule');
    xASL_adm_DocCrawler(fullfile(x.MyPath,'Modules','SubModule_ASL'), fullfile(outputFolder,'ASL_Module.md'),'ASLModule');
    xASL_adm_DocCrawler(fullfile(x.MyPath,'Modules','SubModule_Population'), fullfile(outputFolder,'Population_Module.md'),'PopulationModule');
    
    
    %% Update documentation read me files in the ExploreASL folder structure
    if updateExploreASL
        copyfile(fullfile(outputFolder,'Functions.md'),fullfile(x.MyPath,'Functions','README.md'));                                     % Functions
        copyfile(fullfile(outputFolder,'Modules.md'),fullfile(x.MyPath,'Modules','README.md'));                                         % Modules
        copyfile(fullfile(outputFolder,'Structural_Module.md'),fullfile(x.MyPath,'Modules','SubModule_Structural','README.md'));        % SubModules (Structural)
        copyfile(fullfile(outputFolder,'ASL_Module.md'),fullfile(x.MyPath,'Modules','SubModule_ASL','README.md'));                  % SubModules (ASL)
        copyfile(fullfile(outputFolder,'Population_Module.md'),fullfile(x.MyPath,'Modules','SubModule_Population','README.md'));    % SubModules (Population)
    end



end

%% Function to swap text segments of documentation files
function swapTextInFile(filePath,text2swap,newText)
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
    
    %% Change text
    for curLine=1:numel(text_cell)
        if ~isempty(text_cell{curLine,1})
            if ~isempty(contains(text_cell{curLine,1},text2swap))
                text_cell{curLine,1}= strrep(char(text_cell{curLine,1}),text2swap,newText);
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




