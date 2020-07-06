function xASL_docu_Crawler(folder,mdoutput)
%xASL_docu_Crawler Script to get information from the file headers and
% convert the information into a markdown file.
%
% FORMAT:       xASL_docu_Crawler(folder)
% 
% INPUT:        folder - input folder
%               mdoutput - result file
%
% OUTPUT:       None
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function checks each individual file header and
%               extracts the information. The results is saved as a
%               markdown file.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLES:     xASL_docu_Crawler('M:\...\Functions', 'M:\...\Output.md')
% __________________________________
% Copyright 2015-2020 ExploreASL

% Improve command window output
BreakString = [repmat('=',1,100),'\n'];
SeparatorLine = repmat('-',1,149);

%% Update GIT
fprintf(BreakString)
% try
%     fprintf('GIT: ');
%     system('git fetch','-echo');
%     system('git pull','-echo');
% catch
%     warning('It seems that your working directory is not the ExploreASL directory...');
% end

%% Input Check
if nargin < 1
    error('Input folder not defined...')
end

%% Check directory
listing = dir(folder);

%% Iterate over files

% Remove folders
folderList = [listing.isdir]';
listing(folderList) = [];

% Iterator
it = 1;

% Get header information from each file
for i = 1:numel(listing)
    
    fileName = listing(i).name;
    filePath = fullfile(folder, fileName);
    
    % Extract header
    try
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
            continue
        end
        
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
        
    catch
        extractHeaderInfo = false;
    end
    
    % Extract and save information
    if extractHeaderInfo
        
        % Write data to markdown
        
        % Length of format text
        lF = size(formatText);
        lF = lF(1);
        
        % Length of description text
        lD = size(descriptionText);
        lD = lD(1);
        
        % Filenames
        TEXT{it,1} = '----'; it = it+1;
        escapedName = strrep(fileName,'_','\_');
        TEXT{it,1} = char(['### ' escapedName]); it = it+1;
        TEXT{it,1} = ''; it = it+1;
        
        % Create format description
        TEXT{it,1} = '#### Format'; it = it+1; 
        TEXT{it,1} = ''; it = it+1;
        TEXT{it,1} = '```matlab'; it = it+1;
        for i=1:lF
            TEXT{it,1} = formatText(i,:); it = it+1;
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
        for i=1:lD
            TEXT{it,1} = descriptionText(i,:);
            if contains(TEXT{it,1},'_')
                TEXT{it,:} = strrep(TEXT{it,:},'_','\_'); % Escape underscores in the description
                TEXT{it,:} = strrep(TEXT{it,:},'*','\*'); % Escape stars in the description
            end
            it = it+1;
        end
        
        % Empty lines
        TEXT{it,1} = ''; it = it+1;
        
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
fprintf(BreakString);









