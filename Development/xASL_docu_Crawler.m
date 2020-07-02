function xASL_docu_Crawler(folder)
%xASL_docu_Crawler Script to get information from the file headers and
% convert the information into a markdown file.
%
% FORMAT:       xASL_docu_Crawler(folder)
% 
% INPUT:        folder
%
% OUTPUT:       None
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function checks each individual file header and
%               extracts the information. The results is saved as a
%               markdown file.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLES:     xASL_docu_Crawler('C:\...\ExploreASL\Functions')
% __________________________________
% Copyright 2015-2020 ExploreASL

% Improve command window output
BreakString = [repmat('=',1,100),'\n'];

%% Update GIT
fprintf(BreakString)
try
    fprintf('GIT: ');
    system('git fetch','-echo');
    system('git pull','-echo');
catch
    warning('It seems that your working directory is not the ExploreASL directory...');
end

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
        
        
        % RIGHT NOW IT STILL OVERWRITES THE TEXT EACH TIME ETC. BUT IT
        % SEEMS PRETTY EASY SO FAR...
        
        fileID = fopen('C:\...\test.md','w');
        
        fprintf(fileID,'----\nFORMAT:\n');
        
        for i=1:lF
            fprintf(fileID,'%s\n',formatText(i,:));
        end
        
        fprintf(fileID,'----\nDESCRIPTION:\n');
        
        for i=1:lD
            fprintf(fileID,'%s\n',descriptionText(i,:));
        end
        
        
        fclose(fileID);
        
    end
    
end



