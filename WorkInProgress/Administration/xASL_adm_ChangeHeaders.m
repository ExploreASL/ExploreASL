%xASL_adm_ChangeHeaders Change the header of each ExploreASL Matlab .m file
% e.g., to add a license statement

% Initialize ExploreASL;
x = ExploreASL;

% Define the directory containing your .m files
x.opts.MyPath

% Get the length of MyPath to print shorter paths
myPathLength = length(x.opts.MyPath);

folderPath = 'C:\path\to\your\folder';

% Get a list of all .m files in the directory
fileList = xASL_adm_GetFileList(x.opts.MyPath, '.*.m$', 'FPListRec');

% Define the text to add after the first '% Copyright' line
headerText = {
    '% Licensed under Apache 2.0, see permissions and limitations at'
    '% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE-EXPLOREASL'
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
        fileLines = strsplit(fileContent, '\n')';
        
        % Find lines that starts with '% Copyright'
        copyrightIndex = find(startsWith(strtrim(fileLines), '% Copyright'));
        if length(copyrightIndex)<1
            fprintf('%s\n', ['No copyright line found: ' fileList{iFile}(myPathLength+1:end)]);
            count = count+1;
        elseif length(copyrightIndex)>1
            fprintf('%s\n', [xASL_num2str(length(copyrightIndex)) ' copyright lines found: ' fileList{iFile}(myPathLength+1:end)]);
        else
    
            % Insert headerText after the first '% Copyright' line
            fileLines = [fileLines(1:copyrightIndex); headerText; fileLines(copyrightIndex+1:end)];
    
            % If no '% Copyright' line is found, you might decide to add the header at the top
            % Uncomment the following line if you want to add it at the top when not found:
            % fileLines = [headerText; fileLines];
            
            % Rejoin the file lines into a single string
            newFileContent = strjoin(fileLines, '\n');
            
            % Open the file for writing and overwrite with new content
            fid = fopen(fileList{iFile}, 'w');
            fprintf(fid, '%s', newFileContent);
            fclose(fid);
        end
    end
end