function [identical,results,reportTable] = xASL_bids_CompareStructures(pathDatasetA,pathDatasetB,bPrintReport,threshRmseNii,detailedOutput,printWarnings,ignoreLogs)
%xASL_bids_CompareStructures Function that compares two BIDS folders with several subfolders and studies and prints the differences.
%
% FORMAT: [identical,results,reportTable] = xASL_bids_CompareStructures(pathDatasetA,pathDatasetB,[bPrintReport,threshRmseNii,detailedOutput,printWarnings,ignoreLogs]);
%
% INPUT:
%        pathDatasetA       - path to first BIDS structure [char array] (REQUIRED)
%        pathDatasetB       - path to second BIDS structure [char array] (REQUIRED)
%        bPrintReport       - true or false to print console report (OPTIONAL, DEFAULT = true)
%        threshRmseNii      - normalized RMSE threshold for comparing NIFTI content (OPTIONAL, DEFAULT = 1e-5)
%        detailedOutput     - additional text ouput (also print that there are no missing files etc.) (OPTIONAL, DEFAULT = false)
%        printWarnings      - print differences as warnings (OPTIONAL, DEFAULT = true)
%        ignoreLogs         - ignore log files (OPTIONAL, DEFAULT = false)
%
% OUTPUT:
%        identical          - Returns 1 if both folder structures are identical and 0 if not
%        results            - structure containing (possible) differences of both folder structures
%        reportTable        - Report table (this is supposed to be the successor of "results", with a simple overview) (TABLE)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Function that compares two BIDS folders with several subfolders and studies and prints the differences.
%                   We recommend to set bPrintReport to true, because you otherwise can't see significant file content differences.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          pathDatasetA = '...\bids-examples\eeg_rest_fmri';
%                   pathDatasetB = '...\bids-examples\eeg_rest_fmri_exact_copy'
%                   [identical,results] = xASL_bids_CompareStructures(pathDatasetA,pathDatasetB,true,0.01);
%
% REFERENCES:       ...
% __________________________________
% Copyright (c) 2015-2021 ExploreASL


    %% Input Check

    % Check if both root folders are valid char arrays or strings
    if ~ischar(pathDatasetA)
        error('The path of structure A is not a char array');
    end
    if ~ischar(pathDatasetB)
        error('The path of structure B is not a char array');
    end
    if strcmp(pathDatasetA,pathDatasetB)
        warning('The path of dataset A is equal to the path of dataset B...');
    end

    % Check if both root folders exists
    if ~xASL_exist(pathDatasetA, 'dir')
        error('The root folder of structure A does not exist: %s',pathDatasetA);
    end
    if ~xASL_exist(pathDatasetB, 'dir')
        error('The root folder of structure B does not exist: %s',pathDatasetB);
    end
	
    % Default value for bPrintReport
	if nargin < 3 || isempty(bPrintReport)
		bPrintReport = true;
	end
    
    % Default value for RMSE threshold
    if nargin < 4 || isempty(threshRmseNii)
       threshRmseNii = 1e-5;
    end
    
    % Detailed output
    if nargin < 5 || isempty(detailedOutput)
       detailedOutput = false;
    end
    
    % Print differences as warnings
    if nargin < 6 || isempty(printWarnings)
       printWarnings = true;
    end
    
    % Ignore log files
    if nargin < 7 || isempty(ignoreLogs)
        ignoreLogs = false;
    end


    %% Defaults

    % Set identical to true (will be set to false as soon as a difference is found)
    identical = true;

    % Initialize results structure
    results = struct;
    
    % Report table
    % - dataset: datasetA or datasetB name
    % - name: file or folder name
    % - message: information about missing file or folder or difference in content
    reportTable = array2table(zeros(0,3), 'VariableNames',{'dataset','name','message'});

    %% Initialization
    [results,pathDatasetA,pathDatasetB,datasetA,datasetB,fileListA,fileListB,~,~] = ...
        xASL_bids_CompareStructures_Init(results,pathDatasetA,pathDatasetB,ignoreLogs);

    %% Checks
    
    % Identical check
    if ~isempty(results.(datasetA).missingFolders) || ~isempty(results.(datasetB).missingFolders)
        identical = false;
    end
    
    % Identical check
    if ~isempty(results.(datasetA).missingFiles) || ~isempty(results.(datasetB).missingFiles)
        identical = false;
    end
    
    % Full Report
    if bPrintReport
        xASL_bids_CompareStructures_PrintFullReport(results,datasetA,datasetB,detailedOutput);
    end
    
    % Compare file content
    [identical,results.differences] = xASL_bids_CompareStructuresCheckContent(...
        fileListA,fileListB,...
        pathDatasetA,pathDatasetB,...
        identical,bPrintReport,detailedOutput,threshRmseNii);
    
    % Add all entries (missing files, folders and differences in content) to the table
    reportTable = xASL_bids_CompareStructures_AddEntriesToTable(reportTable,results,datasetA,datasetB);
    
    % Print differences as warnings
    if printWarnings
        printMissingAsWarnings(results);
        printDifferencesAsWarnings(results.differences);
    end

end


%% Add all entries (missing files, folders and differences in content) to the table
function reportTable = xASL_bids_CompareStructures_AddEntriesToTable(reportTable,results,datasetA,datasetB)

    % Go through missing folders
    if size(results.(datasetA).missingFolders,1)>0 && ...
       size(results.(datasetA).missingFolders,2)>0 && ...
            ~isempty(results.(datasetA).missingFolders{1})
        for iElement=1:size(results.(datasetA).missingFolders,1)
            dataset = cellstr(datasetA);
            name = cellstr('Missing folder');
            message = cellstr(results.(datasetA).missingFolders{iElement,1});
            reportTable = xASL_bids_CompareStructures_AddTableRow(reportTable,dataset,name,message);
        end
    end
    if size(results.(datasetB).missingFolders,1)>0 && ...
       size(results.(datasetB).missingFolders,2)>0 && ...
            ~isempty(results.(datasetB).missingFolders{1})
        for iElement=1:size(results.(datasetB).missingFolders,1)
            dataset = cellstr(datasetB);
            name = cellstr('Missing folder');
            message = cellstr(results.(datasetB).missingFolders{iElement,1});
            reportTable = xASL_bids_CompareStructures_AddTableRow(reportTable,dataset,name,message);
        end
    end

    % Go through missing files
    if size(results.(datasetA).missingFiles,1)>0 && ...
       size(results.(datasetA).missingFiles,2)>0 && ...
            ~isempty(results.(datasetA).missingFiles{1})
        for iElement=1:size(results.(datasetA).missingFiles,1)
            dataset = cellstr(datasetA);
            name = cellstr('Missing file');
            message = cellstr(results.(datasetA).missingFiles{iElement,1});
            reportTable = xASL_bids_CompareStructures_AddTableRow(reportTable,dataset,name,message);
        end
    end
    if size(results.(datasetB).missingFiles,1)>0 && ...
       size(results.(datasetB).missingFiles,2)>0 && ...
            ~isempty(results.(datasetB).missingFiles{1})
        for iElement=1:size(results.(datasetB).missingFiles,1)
            dataset = cellstr(datasetB);
            name = cellstr('Missing file');
            message = cellstr(results.(datasetB).missingFiles{iElement,1});
            reportTable = xASL_bids_CompareStructures_AddTableRow(reportTable,dataset,name,message);
        end
    end

    % Go through differences
    if ~isempty(results.differences{1})
        for iElement=1:size(results.differences,1)
            dataset = cellstr('Both');
            name = cellstr('Different file content');
            message = cellstr(results.differences{iElement,1});
            reportTable = xASL_bids_CompareStructures_AddTableRow(reportTable,dataset,name,message);
        end
    end


end


%% Add row to table
function reportTable = xASL_bids_CompareStructures_AddTableRow(reportTable,dataset,name,message)
    
    % Create row
    thisStruct.dataset = dataset;
    thisStruct.name = name;
    thisStruct.message = message;
    thisRow = struct2table(thisStruct);
    % Add row
    reportTable = [reportTable;thisRow];

end


%% Compare Structures Initialization
function [results,pathDatasetA,pathDatasetB,datasetA,datasetB,fileListA,fileListB,folderListA,folderListB] = ...
    xASL_bids_CompareStructures_Init(results,pathDatasetA,pathDatasetB,ignoreLogs)

    % Remove last character if it is a slash
    if strcmp(pathDatasetA(end),'/') || strcmp(pathDatasetA(end),'\')
        pathDatasetA = pathDatasetA(1:end-1);
    end
    if strcmp(pathDatasetB(end),'/') || strcmp(pathDatasetB(end),'\')
        pathDatasetB = pathDatasetB(1:end-1);
    end

    % Convert to valid paths
    pathDatasetA = fullfile(pathDatasetA);
    pathDatasetB = fullfile(pathDatasetB);

    % Get dataset names
    [~, datasetA] = fileparts(pathDatasetA);
    [~, datasetB] = fileparts(pathDatasetB);

    % Make sure you have valid identifiers for the field names
    datasetA = matlab.lang.makeValidName(datasetA,'ReplacementStyle','delete');
    datasetB = matlab.lang.makeValidName(datasetB,'ReplacementStyle','delete');

    % Make sure that datasetA and datasetB are not exactly the same
    if strcmp(datasetA, datasetB)
        sameName = datasetA;
        datasetA = [sameName '_A'];
        datasetB = [sameName '_B'];
    end

    results.(datasetA) = struct;
    results.(datasetB) = struct;

    % Get files and folders of datasets A and B
    [~, dateVersion] = version;
    if xASL_str2num(dateVersion(end-4:end))>2016 % dir does not work recursively in older versions
        filesA = dir(fullfile(pathDatasetA, '**','*.*'));
        filesB = dir(fullfile(pathDatasetB, '**','*.*'));
    else
        [filesA, filesB] = xASL_bids_CompareStructures_GetFileListsUnix(pathDatasetA,pathDatasetB);
    end

    % Remove root path
    filesA = modifyFileList(filesA,pathDatasetA,ignoreLogs);
    filesB = modifyFileList(filesB,pathDatasetB,ignoreLogs);

    % Get lists
    fileListA = getListWithout('folders', filesA);
    fileListB = getListWithout('folders', filesB);
    folderListA = getListWithout('files', filesA);
    folderListB = getListWithout('files', filesB);

    % Get folder lists
    folderListA = unique({folderListA.folder}');
    folderListB = unique({folderListB.folder}');

    % Get real file lists
    fileListA = unique({fileListA.name}');
    fileListB = unique({fileListB.name}');

    % Missing Folders
    results.(datasetA).missingFolders = setdiff(folderListB,folderListA);
    results.(datasetB).missingFolders = setdiff(folderListA,folderListB);

    % Missing Files
    results.(datasetA).missingFiles = setdiff(fileListB,fileListA);
    results.(datasetB).missingFiles = setdiff(fileListA,fileListB);

end


%% Print full report
function xASL_bids_CompareStructures_PrintFullReport(results,datasetA,datasetB,detailedOutput)

    if detailedOutput
        fprintf(strcat(repmat('=',100,1)','\n'));
    end
    if detailedOutput || ~isempty(results.(datasetA).missingFolders) || ~isempty(results.(datasetA).missingFiles)
        fprintf('Dataset:   %s\n',datasetA)
    end
    printList(results.(datasetA).missingFolders)
    printList(results.(datasetA).missingFiles)

    if detailedOutput
        if isempty(results.(datasetA).missingFolders) && isempty(results.(datasetA).missingFiles)
            fprintf('           %s\n','No missing files');
        end
    end

    if detailedOutput
        fprintf(strcat(repmat('=',100,1)','\n'));
    end
    if detailedOutput || ~isempty(results.(datasetB).missingFolders) || ~isempty(results.(datasetB).missingFiles)
        fprintf('Dataset:   %s\n',datasetB)
    end
    printList(results.(datasetB).missingFolders)
    printList(results.(datasetB).missingFiles)

    if detailedOutput
        if isempty(results.(datasetB).missingFolders) && isempty(results.(datasetB).missingFiles)
            fprintf('           %s\n','No missing files');
        end
    end

    % End of report
    if detailedOutput
        fprintf(strcat(repmat('=',100,1)','\n'));
    end

end



%% Print differences as warnings
function printDifferencesAsWarnings(differences)

    % Iterate over differences
    if ~isempty(differences{1,1})
        for iT = 1:size(differences,1)
            warning([differences{iT,1} '  ']);
            fprintf('\n');
        end
    end
    
end


%% Print missing files and folders as warnings
function printMissingAsWarnings(results)

    fieldNames = fieldnames(results);
    for iT = 1:(length(fieldNames)-1)
        if ~strcmp(fieldNames(iT),'differences')
            % Folders
            if ~isempty(results.(fieldNames{iT}).missingFolders)
                for thisWarning = 1:size(results.(fieldNames{iT}).missingFolders,1)
                    warning('Missing folder %s  ', results.(fieldNames{iT}).missingFolders{thisWarning,1});
                    fprintf('\n');
                end                
            end
            % Files
            if ~isempty(results.(fieldNames{iT}).missingFiles)
                for thisWarning = 1:size(results.(fieldNames{iT}).missingFiles,1)
                    warning('Missing file %s  ', results.(fieldNames{iT}).missingFiles{thisWarning,1});
                    fprintf('\n');
                end  
            end            
        end
    end
    
end


%% Get list without files/folders
function returnList = getListWithout(thisType,List)
    
    element = 1;
    for iT = 1:size(List,1)
        if strcmp(thisType,'folders')
            if ~List(iT,1).isdir
                returnList(element,1) = List(iT,1);
                element = element+1;
            end
        end
        if strcmp(thisType,'files')
            if List(iT,1).isdir
                returnList(element,1) = List(iT,1);
                element = element+1;
            end
        end
    end
    
    % Make sure there is a list to return
    if ~exist('returnList','var')
        returnList.name = '';
        returnList.folder = '';
        returnList.date = '';
        returnList.bytes = '';
        returnList.isdir = 0;
        returnList.datenum = 0;
    end

end


%% Modify file lists
function fileList = modifyFileList(fileList,root,ignoreLogs)
    
    % Iterate over file list: change folder names
    for iFile=1:numel(fileList)
        fileList(iFile).folder = strrep(fileList(iFile).folder,root,'');
    end
    
    % Remove log files from list
    if ignoreLogs
        for iFile=1:numel(fileList)
            [~, ~, extension] = xASL_fileparts(fileList(iFile).name);
            if strcmp(extension,'.log')
                % The easiest way to remove them from the list was to make them an empty directory basically
                fileList(iFile).name = '.';
                fileList(iFile).folder = '';
                fileList(iFile).isdir = 1;
            end
        end
    end
    
    % Iterate over file list: change file names
    for iFile=1:numel(fileList)
        % Check that the current element is not a folder
        if ~strcmp(fileList(iFile).name,'.') && ~strcmp(fileList(iFile).name,'..')
            fileList(iFile).name = fullfile(fileList(iFile).folder,fileList(iFile).name);
        end
    end
    
    
    
end


%% Print list functions
function printList(currentList)
    % Iterate over list
    if ~isempty(currentList)
        for iFile=1:length(currentList)
            fprintf('Missing:   %s\n',currentList{iFile})
        end
    end
end





%% xASL_bids_CompareStructures_GetFileListsUnix
function [filesA, filesB] = xASL_bids_CompareStructures_GetFileListsUnix(pathDatasetA,pathDatasetB)

    fprintf('This method is not able to find empty directories right now...\n');

    % Get file lists
    onlyFilesA = xASL_adm_GetFileList(pathDatasetA,'^.+$','FPListRec');
    onlyFilesB = xASL_adm_GetFileList(pathDatasetB,'^.+$','FPListRec');
    
    % Create structure
    filesA = struct;
    
    % Iterate over file lists
    for iFile = 1:numel(onlyFilesA)
        [thisFolder, thisFile, thisExtension] = xASL_fileparts(onlyFilesA(iFile,1));
        filesA(iFile).name = [thisFile thisExtension];
        filesA(iFile).isdir = 0;
        filesA(iFile).folder = thisFolder;
    end
    filesB = struct;
    for iFile = 1:numel(onlyFilesB)
        [thisFolder, thisFile, thisExtension] = xASL_fileparts(onlyFilesB(iFile,1));
        filesB(iFile).name = [thisFile thisExtension];
        filesB(iFile).isdir = 0;
        filesB(iFile).folder = thisFolder;

    end
    
    % Add folder list
    addId = size(filesA,2)+1;
    for iFile = 1:size(filesA,2)
        relativeFolder = strrep(filesA(iFile).folder,pathDatasetA,'');
        if ~isempty(relativeFolder)
            filesA(addId).name = relativeFolder;
            filesA(addId).isdir = 1;
            filesA(addId).folder = relativeFolder;
            addId = addId+1;
        end
    end
    addId = size(filesB,2)+1;
    for iFile = 1:size(filesB,2)
        relativeFolder = strrep(filesB(iFile).folder,pathDatasetB,'');
        if ~isempty(relativeFolder)
            filesB(addId).name = relativeFolder;
            filesB(addId).isdir = 1;
            filesB(addId).folder = relativeFolder;
            addId = addId+1;
        end
    end
    
    % Transpose
    filesA = filesA';
    filesB = filesB';

end


