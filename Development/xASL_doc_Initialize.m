function xASL_doc_Initialize
%xASL_doc_Initialize Script to call the separate documentation crawler
% functions.
%
% FORMAT:       xASL_docu_Initialize
% 
% INPUT:        None
%
% OUTPUT:       None
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function generates all markdown files, which are
%               necessary for the mkdocs documentation.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_doc_Initialize
% __________________________________
% Copyright 2015-2020 ExploreASL



    %% Workflow

    % Initialize ExploreASL
    x = ExploreASL_Initialize([],0);

    % Copy and modify the index README
    copyfile(fullfile(x.MyPath,'README.md'),fullfile(x.MyPath,'Development','Documentation_GitHub','index.md'));
    % Logo
    swapTextInFile(fullfile(x.MyPath,'Development','Documentation_GitHub','index.md'),...
                  '(https://github.com/ExploreASL/ExploreASL/blob/develop/ExploreASL_logoSmall.png)',...
                  '(./img/title.png "ExploreASL")');
    % Workflow
    swapTextInFile(fullfile(x.MyPath,'Development','Documentation_GitHub','index.md'),...
                  '(https://www.researchgate.net/profile/Andrew_Robertson7/publication/337328693/figure/fig1/AS:826578854481921@1574083164220/Schematic-diagram-of-ExploreASL-processing-steps-Steps-marked-with-a-are-optional.ppm "Workflow of ExploreASL")',...
                  '(./img/ExploreASL_Workflow.jpg "Workflow ExploreASL")');
    
    % Create the functions markdown file
    xASL_doc_Crawler(fullfile(x.MyPath,'Functions'), fullfile(x.MyPath,'Development','Documentation_GitHub','Functions.md'),'Functions');

    % Convert and copy lincense file
    convertLicenseToMarkdown(fullfile(x.MyPath,'LICENSE-EXPLOREASL'),fullfile(x.MyPath,'Development','Documentation_GitHub','License.md'));
    
    % Use documentation crawler for modules
    xASL_doc_Crawler({fullfile(x.MyPath,'ExploreASL_Import.m'),...
                      fullfile(x.MyPath,'Modules','xASL_module_Structural.m'),...
                      fullfile(x.MyPath,'Modules','xASL_module_ASL.m'),...
                      fullfile(x.MyPath,'Modules','xASL_module_Population.m')},...
                      fullfile(x.MyPath,'Development','Documentation_GitHub','Modules.md'),'Modules');
    
    % Use documentation crawler for submodules
    xASL_doc_Crawler(fullfile(x.MyPath,'Modules','SubModule_Structural'), fullfile(x.MyPath,'Development','Documentation_GitHub','Structural_Module.md'),'StructuralModule');
    xASL_doc_Crawler(fullfile(x.MyPath,'Modules','SubModule_ASL'), fullfile(x.MyPath,'Development','Documentation_GitHub','ASL_Module.md'),'ASLModule');
    xASL_doc_Crawler(fullfile(x.MyPath,'Modules','SubModule_Population'), fullfile(x.MyPath,'Development','Documentation_GitHub','Population_Module.md'),'PopulationModule');




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




