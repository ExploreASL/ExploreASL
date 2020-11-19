function xASL_adm_DocCrawler(inputPath,mdoutput,content)
%xASL_doc_Crawler Script to get information from the file headers and
% convert the information into a markdown file.
%
% FORMAT:       xASL_adm_DocCrawler(inputPath)
% 
% INPUT:        inputPath     - input folder (or file)
%               mdoutput      - result file
%               content       - "Functions", "StructuralModule", ...
%
% OUTPUT:       n/a
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function checks each individual file header and
%               extracts the information. The results is saved as a
%               markdown file.
%
%               If you want to use star symbols (*testFile.m e.g.) we
%               recommend not to use them in the same line with bold text
%               (which is written like this: {{bold text}}).
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_adm_DocCrawler('M:\...\Functions', 'M:\...\Output.md','Functions');
% __________________________________
% Copyright 2015-2020 ExploreASL

    %% Reminder
    fprintf('============================================= REMINDER =============================================\n');
    fprintf('A correct ExploreASL header should include the following tags in the correct order:\n');
    fprintf('"FORMAT:", "INPUT:", "OUTPUT:", "DESCRIPTION:" and "EXAMPLE:"\n');

    %% Defaults
    SeparatorLine = repmat('-',1,149);
    isFileList = false;
    
    % Sections for "Functions" folder
    SECTION = {'adm', 'bids', 'fsl', 'im', 'init', 'io', 'qc', 'quant', 'spm', 'stat', 'vis'}';
    SECTION_NAMES = {'Administration', 'BIDS', 'FSL', 'Imaging', 'Initialization', 'Input and Output', 'QC', 'Quantization', 'SPM', 'Statistics', 'Visualization'}';

    %% Input Check
    if nargin < 1
        error('Input path not defined...')
    end
    
    if nargin < 2
        error('Output file not defined...')
    end
    
    if nargin < 3
        error('Content not defined...')
    end

    %% Check directory
    if ~iscell(inputPath)
        listing = dir(inputPath);
    else
        listing = [];
    end

    %% Iterate over files

    % Check if input is folder or file
    if ~isempty(listing)
        % Remove folders from list
        folderList = [listing.isdir]';
        listing(folderList) = [];
    else
        % Input path is a list of files
        isFileList = true;
        for individualPath=1:numel(inputPath)
            listing(individualPath).name = inputPath{individualPath};
        end
    end

    % Current section and iterator
    cS = 0;
    it = 1;
    
    % Predefine cell array
    TEXT = cell(1,1);

    % Get header information from each file
    for i = 1:numel(listing)
        % Get filename
        fileName = listing(i).name;
        if isFileList
            [~, fileName, ~] = fileparts(fileName);
        end

        % Get information from header if it's not an .md or .mat file
        [extractHeaderInfo,formatText,descriptionText] = analyzeHeader(fileName);

        % Extract and save information
        if extractHeaderInfo

            % Write data to markdown

            % Length of format text
            lF = size(formatText);
            lF = lF(1);

            % Length of description text
            lD = size(descriptionText);
            lD = lD(1);

            % Start of function description
            if cS==0
                if strcmp(content,'Functions')
                    TEXT{it,1} = '# Functions';  it = it+1;
                    TEXT{it,1} = ' ';  it = it+1;
                    TEXT{it,1} = '----';  it = it+1;
                    TEXT{it,1} = '## General Functions';  it = it+1;
                    TEXT{it,1} = ' ';  it = it+1;
                elseif strcmp(content,'Modules')
                    TEXT{it,1} = '# Modules';  it = it+1;
                    TEXT{it,1} = ' ';  it = it+1;
                elseif strcmp(content,'ImportModule')
                    TEXT{it,1} = '# Submodules of the Import Module';  it = it+1;
                    TEXT{it,1} = ' ';  it = it+1;
                elseif strcmp(content,'StructuralModule')
                    TEXT{it,1} = '# Submodules of the Structural Module';  it = it+1;
                    TEXT{it,1} = ' ';  it = it+1;
                elseif strcmp(content,'ASLModule')
                    TEXT{it,1} = '# Submodules of the ASL Module';  it = it+1;
                    TEXT{it,1} = ' ';  it = it+1;
                elseif strcmp(content,'PopulationModule')
                    TEXT{it,1} = '# Submodules of the Population Module';  it = it+1;
                    TEXT{it,1} = ' ';  it = it+1;
                end
                cS = cS+1; 
            end
            
            if strcmp(fileName,'ExploreASL_Import')
                TEXT{it,1} = '----';  it = it+1;
                TEXT{it,1} = '## 1. Import Module';  it = it+1;
                TEXT{it,1} = ' ';  it = it+1;
            elseif strcmp(fileName,'xASL_module_Structural')
                TEXT{it,1} = '----';  it = it+1;
                TEXT{it,1} = '## 2. Structural Module';  it = it+1;
                TEXT{it,1} = ' ';  it = it+1;
            elseif strcmp(fileName,'xASL_module_ASL')
                TEXT{it,1} = '----';  it = it+1;
                TEXT{it,1} = '## 3. ASL Module';  it = it+1;
                TEXT{it,1} = ' ';  it = it+1;
            elseif strcmp(fileName,'xASL_module_Population')
                TEXT{it,1} = '----';  it = it+1;
                TEXT{it,1} = '## 4. Population Module';  it = it+1;
                TEXT{it,1} = ' ';  it = it+1;
            end
            
            % Get the current section
            if strcmp(content,"Functions")
                if cS <= length(SECTION)
                    if contains(fileName,['xASL_', char(SECTION{cS,1})])
                        TEXT{it,1} = ['## ', char(SECTION_NAMES{cS,1})];  it = it+1;
                        TEXT{it,1} = '';  it = it+1;
                        cS = cS+1;
                    end
                end
            end

            % Filenames
            TEXT{it,1} = '----'; it = it+1;
            escapedName = strrep(fileName,'_','\_');
            TEXT{it,1} = char(['### ' escapedName]); it = it+1;
            TEXT{it,1} = ''; it = it+1;

            % Create format description
            TEXT{it,1} = '#### Format'; it = it+1; 
            TEXT{it,1} = ''; it = it+1;
            TEXT{it,1} = '```matlab'; it = it+1;
            for iter=1:lF
                TEXT{it,1} = formatText(iter,:); it = it+1;
            end
            % Remove empty char array after FORMAT description
            if strcmp(strtrim(TEXT{it-1,1}),'')
                it=it-1;
            end
            TEXT{it,1} = '```'; it = it+1;
            TEXT{it,1} = ''; it = it+1;

            % Create description
            TEXT{it,1} = '#### Description';
            it = it+1;
            for iter=1:lD
                TEXT{it,1} = descriptionText(iter,:);
                % WITH THIS SYNTAX WE CANT USE STAR SYMBOLS AND BOLD TEXT IN ONE LINE
                % OUR NEW BOLD TEXT SYMBOLS ARE {{ AND }} (curly brackets)
                if contains(TEXT{it,1},' {') || contains(TEXT{it,1},'} ')
                    TEXT{it,:} = strrep(TEXT{it,:},'{{','**');
                    TEXT{it,:} = strrep(TEXT{it,:},'}}','**');
                elseif contains(TEXT{it,1},'*')
                    % Escape stars in the description
                    TEXT{it,:} = strrep(TEXT{it,:},'*','\*'); 
                end

                if contains(TEXT{it,1},'_')
                    % Escape underscores in the description
                    TEXT{it,:} = strrep(TEXT{it,:},'_','\_'); 
                end
                it = it+1;
            end

            % Empty lines
            TEXT{it,1} = ''; it = it+1;
        else
            fprintf('File %s does not fulfill the documentation requirements...\n', fileName);
        end

    end

    % Remove spaces and separator lines
    for i=1:numel(TEXT)
        % Trim char arrays
        TEXT{i,:} = strtrim(TEXT{i,:});

        % Remove line separators
        if strcmp(TEXT(i,:),SeparatorLine)
            TEXT{i,:} = strrep(TEXT{i,:},SeparatorLine,'');
        end
    end

    % Print information to markdown file
    fileID = fopen(mdoutput,'w');
    for i=1:numel(TEXT)
        fprintf(fileID,'%s\n',char(TEXT{i,:}));
    end
    fclose(fileID);

    % Final output
    fprintf('Markdown file generated...\n');

end

% Helper function
function [extractHeaderInfo,formatText,descriptionText] = analyzeHeader(fileName)

    % Fallback
    extractHeaderInfo = false;
    formatText = NaN;
    descriptionText = NaN;
    
    % Skip .md and .mat files
    if ~(contains(fileName,'.md'))
        if ~(contains(fileName,'.mat'))
            % Extract header
            try
                % Get header from help function
                header = help(fileName);

                % Convert header to cell array
                headerList = strsplit(header,'\n');
                headerList = headerList';

                % Get line with FORMAT description
                [stringsFORMAT, ~, ~]= regexpi(headerList, 'FORMAT:', 'match');

                % Get line with FORMAT description
                [stringsINPUT, ~, ~]= regexpi(headerList, 'INPUT:', 'match');

                % Get line with DESCRIPTION
                [stringsDESCRIPTION, ~, ~]= regexpi(headerList, 'DESCRIPTION:', 'match');

                % Get line with EXAMPLE
                [stringsEXAMPLE, ~, ~]= regexpi(headerList, 'EXAMPLE:', 'match');

                % Find indices
                indexFORMAT = find(~cellfun(@isempty,stringsFORMAT));
                indexINPUT = find(~cellfun(@isempty,stringsINPUT));
                indexDESCRIPTION = find(~cellfun(@isempty,stringsDESCRIPTION));
                indexEXAMPLE = find(~cellfun(@isempty,stringsEXAMPLE));

                % Check if all descriptions are there
                if isempty(indexFORMAT) || isempty(indexINPUT) || isempty(indexDESCRIPTION) || isempty(indexEXAMPLE)
                    extractHeaderInfo = false;
                else
                    % FORMAT text
                    formatText = headerList(indexFORMAT:indexINPUT-1);
                    formatText = strrep(formatText,'FORMAT:','');
                    formatText = char(strtrim(formatText));

                    % DESCRIPTION text
                    descriptionText = headerList(indexDESCRIPTION:indexEXAMPLE-1);
                    descriptionText = strrep(descriptionText,'DESCRIPTION:','');
                    descriptionText = char(strtrim(descriptionText));

                    % Give the okay to extract the information of the current file
                    extractHeaderInfo = true;
                end

            catch
                extractHeaderInfo = false;
            end 
        end
    end

end



