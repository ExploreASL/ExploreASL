% Copyright 2015-2024 ExploreASL (Works In Progress code)
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

%xASL_adm_ChangeHeaders Change the header of each ExploreASL Matlab .m file
% e.g., to add a license statement
% Initialize ExploreASL;
x = ExploreASL;
% Define the directory containing your .m files
x.opts.MyPath
% Get the length of MyPath to print shorter paths
myPathLength = length(x.opts.MyPath);
% Get a list of all .m files in the directory
fileList = xASL_adm_GetFileList(x.opts.MyPath, '.*.m$', 'FPListRec');
% Define the text to add after the first '% Copyright' line
headerText = {
    '% Licensed under Apache 2.0, see permissions and limitations at'
    '% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE'
    '% you may only use this file in compliance with the License.'
    '% __________________________________';
    '';
};

count = 0; % count number of files without copyright line

% Loop through each file
for iFile = 1:length(fileList)
    % Skip 
    [fPath] = fileparts(fileList{iFile});
    if ~contains(fPath, 'External')
        % Only process files that are not within the external folder
        % i.e. subjective to the ExploreASL license
        % Read the existing content of the file into a cell array, one line per cell
        fileContent = fileread(fileList{iFile});
        fileLines = strsplit(fileContent, '\n', 'CollapseDelimiters', false)';
        
        % Find lines that starts with '% Copyright'
        copyrightIndex = find(~cellfun(@isempty, regexp(fileLines, '^\s*%\s*Copyright')));

        % First we fix the files that don't have a copyright statement
        if length(copyrightIndex)<1
            fprintf('%s\n', ['No copyright line found: ' fileList{iFile}(myPathLength+1:end)]);
            count = count+1;
            % Add the copyright statement at the beginning of the file
            copyrightStatement = '% Copyright 2015-2024 ExploreASL (Works In Progress code)';
            fileLines = [{copyrightStatement}; fileLines];
        end

        % Now we get the lines with copyright statement again, including the fixed files
        copyrightIndex = find(~cellfun(@isempty, regexp(fileLines, '^\s*%\s*Copyright')));        

        if length(copyrightIndex)>1
            fprintf('%s\n', [xASL_num2str(length(copyrightIndex)) ' copyright lines found: ' fileList{iFile}(myPathLength+1:end)]);
        else
            % Insert headerText after the first '% Copyright' line
            fileLines = [fileLines(1:copyrightIndex); headerText; fileLines(copyrightIndex+1:end)];
        end
       
        % Rejoin the file lines into a single string
        newFileContent = strjoin(fileLines, '\n');
        % Open the file for writing and overwrite with new content
        fid = fopen(fileList{iFile}, 'w');
        fprintf(fid, '%s', newFileContent);
        fclose(fid);
     end
end