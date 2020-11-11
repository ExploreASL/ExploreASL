function [identical,results] = xASL_bids_CompareStructures(pathDatasetA,pathDatasetB,printReport,thresh)
%xASL_bids_CompareStructures Function that compares two BIDS folders with several subfolders and studies and prints the differences.
%
% FORMAT: [identical,results] = xASL_bids_CompareStructures(pathDatasetA,pathDatasetB,[printReport,thresh]);
%
% INPUT:
%        pathDatasetA       - path to first BIDS structure (REQUIRED)
%        pathDatasetB       - path to second BIDS structure (REQUIRED)
%        printReport        - true or false to print console report (OPTIONAL, DEFAULT = true)
%        thresh             - RMSE threshold for comparing NIFTI content (OPTIONAL, DEFAULT = 0.1)
%
% OUTPUT:
%        identical          - Returns 1 if both folder structures are identical and 0 if not
%        results            - structure containing (possible) differences of both folder structures
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Function that compares two BIDS folders with several subfolders and studies and prints the differences.
%                   We recommend to set printReport to true, because you otherwise can't see significant file content differences.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          pathDatasetA = '...\bids-examples\eeg_rest_fmri';
%                   pathDatasetB = '...\bids-examples\eeg_rest_fmri_exact_copy'
%                   [identical,results] = xASL_bids_CompareStructures(pathDatasetA,pathDatasetB,true);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2020 ExploreASL


    %% Input Check

    % Check if both root folders are valid char arrays or strings
    if ~(ischar(pathDatasetA) || isstring(pathDatasetA))
        error('The path of structure A is neither a char array not a string...');
    end
    if ~(ischar(pathDatasetB) || isstring(pathDatasetB))
        error('The path of structure A is neither a char array not a string...');
    end

    % Check if both root folders exists
    if ~(xASL_exist(pathDatasetA)==7)
        error('The root folder of structure A does not exist...');
    end
    if ~(xASL_exist(pathDatasetB)==7)
        error('The root folder of structure B does not exist...');
    end
	
    % Default value for printReport
	if nargin < 3 || isempty(printReport)
		printReport = true;
    end
    
    % Default value for RMSE threshold
    if nargin < 4
       thresh = 0.1; 
    end


    %% Defaults

    % Set identical to true (will be set to false as soon as a difference is found)
    identical = true;

    % Initialize results structure
    results = struct;

    %% Initialization
    
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
    [~,datasetA,~] = fileparts(pathDatasetA);
    [~,datasetB,~] = fileparts(pathDatasetB);

    % Make sure you have valid identifiers for the field names
    datasetA = matlab.lang.makeValidName(datasetA,'ReplacementStyle','delete');
    datasetB = matlab.lang.makeValidName(datasetB,'ReplacementStyle','delete');
    results.(datasetA) = struct;
    results.(datasetB) = struct;
    
    % Get files and folders of datasets A and B
    filesA = dir(fullfile(pathDatasetA, '**','*.*'));
    filesB = dir(fullfile(pathDatasetB, '**','*.*'));
    
    % Remove root path
    filesA = modifyFileList(filesA,pathDatasetA);
    filesB = modifyFileList(filesB,pathDatasetB);
    
    % Get folder lists
    folderListA = unique(string({filesA.folder}'));
    folderListB = unique(string({filesB.folder}'));
    
    % Get real file lists
    fileListA = unique(string({filesA.name}'));
    fileListB = unique(string({filesB.name}'));
    
    % Missing Folders
    results.(datasetA).missingFolders = setdiff(folderListB,folderListA);
    results.(datasetB).missingFolders = setdiff(folderListA,folderListB);
    
    % Missing Files
    results.(datasetA).missingFiles = setdiff(fileListB,fileListA);
    results.(datasetB).missingFiles = setdiff(fileListA,fileListB);
    
    % Identical check
    if ~isempty(results.(datasetA).missingFolders) || ~isempty(results.(datasetB).missingFolders)
        identical = false;
    end
    
    % Identical check
    if ~isempty(results.(datasetA).missingFiles) || ~isempty(results.(datasetB).missingFiles)
        identical = false;
    end
    
    % Full Report
    if printReport
        fprintf(strcat(repmat('=',100,1)','\n'));
        fprintf('Dataset:\t\t%s\n',datasetA)
        printList(results.(datasetA).missingFolders)
        printList(results.(datasetA).missingFiles)

        fprintf(strcat(repmat('=',100,1)','\n'));
        fprintf('Dataset:\t\t%s\n',datasetB)
        printList(results.(datasetB).missingFolders)
        printList(results.(datasetB).missingFiles)

        % End of report
        fprintf(strcat(repmat('=',100,1)','\n'));
    end
    
    % Compare file content
    identical = checkFileContents(fileListA,fileListB,pathDatasetA,pathDatasetB,identical,printReport,thresh);

    
    

end

%% Modify file lists
function fileList = modifyFileList(fileList,root)
    % Iterate over file list: change folder names
    for it=1:numel(fileList)
        fileList(it).folder = strrep(fileList(it).folder,root,'');
    end
    % Iterate over file list: change file names
    for it=1:numel(fileList)
        % Check that the current element is not a folder
        if ~strcmp(fileList(it).name,'.') && ~strcmp(fileList(it).name,'..')
            fileList(it).name = fullfile(fileList(it).folder,fileList(it).name);
        end
    end
end

% Print list functions
function printList(currentList)
    % Iterate over list
    if ~isempty(currentList)
        for it=1:length(currentList)
            fprintf('Missing:\t\t%s\n',currentList(it))
        end
    end
end

%% Check file contents
function identical = checkFileContents(filesDatasetA,filesDatasetB,pathDatasetA,pathDatasetB,identical,printReportFunction,thresh)
    
    % All files
    allFiles = unique([filesDatasetA',filesDatasetB']');
    
    % Iterate over list
    for it=1:length(allFiles)
        % Assign root directory of dataset A
        currentFileA = fullfile(pathDatasetA,allFiles(it));
        currentFileB = fullfile(pathDatasetB,allFiles(it));
        % Get extension
        [~,~,extension] = fileparts(allFiles(it));
        % Check extension
        if strcmp(extension,'.json')
            if printReportFunction
                fprintf('Checking:\t\t%s\n',allFiles(it))
            end
            % Compare JSON files on field basis
            if (isfile(currentFileA) && isfile(currentFileB)) % xASL_exist somehow didn't work here (again)
                % Import JSON files
                jsonA = xASL_import_json(char(currentFileA));
                jsonB = xASL_import_json(char(currentFileB));
                % Get JSON field names
                fieldNamesA = fieldnames(jsonA);
                fieldNamesB = fieldnames(jsonB);
                % Check if both JSON files have the same fields
                differenceA = setdiff(fieldNamesA,fieldNamesB);
                differenceB = setdiff(fieldNamesB,fieldNamesA);
                if isempty(differenceA) && isempty(differenceB)
                    % Both JSON files have the same fields, so we can compare the content
                    compareFieldLists(jsonA,jsonB,fieldNamesA,fieldNamesB,printReportFunction);
                else
                    % The JSON files have different fields
                    sharedFieldsAB = intersect(fieldNamesB,fieldNamesA);
                    allFields = union(fieldNamesA,fieldNamesB);
                    missingFields = setdiff(allFields,sharedFieldsAB);
                    % Print out and compare shared fields 
                    for sharedField=1:numel(sharedFieldsAB)
                        if printReportFunction
                            % fprintf('\t\t\t\tShared field %s ...\n',string(sharedFieldsAB{sharedField}));
                        end
                        % Now we can compare these fields like in the part above
                        compareFieldLists(jsonA,jsonB,sharedFieldsAB,sharedFieldsAB,printReportFunction);
                    end
                    % Print out missing fields 
                    for missingField=1:numel(missingFields)
                        if printReportFunction
                            fprintf('\t\t\t\tMissing field %s ...\n',string(missingFields{missingField}));
                        end
                    end
                end
            end
        elseif strcmp(extension,'.tsv') || strcmp(extension,'.txt') || strcmp(extension,'.csv')
            % Read files if they exist
            if (isfile(currentFileA) && isfile(currentFileB)) % xASL_exist somehow didn't work here (again)
                % Compare text files content directly
                currentFileTextA = fileread(currentFileA);
                currentFileTextB = fileread(currentFileB);
                if printReportFunction
                    fprintf('Checking:\t\t%s\n',allFiles(it))
                    if ~strcmp(currentFileTextA,currentFileTextB)
                        fprintf('\t\t\t\tDifferent file content!\n')
                        identical = false;
                    end
                end
            end
        elseif strcmp(extension,'.nii') || contains(allFiles(it),'.nii.gz') % Warning: This unzips your NIFTIs right now!
            % Read files if they exist
            if (isfile(currentFileA) && isfile(currentFileB)) % xASL_exist somehow didn't work here (again)
                % Check file size (there were some 0KB images in the bids-examples)
                tmpFileA = dir(currentFileA);
                tmpFileB = dir(currentFileB);
                sizeA = tmpFileA.bytes;
                sizeB = tmpFileB.bytes;
                if printReportFunction
                    fprintf('Checking:\t\t%s\n',allFiles(it));
                end
                % Check file size (images with a file size lower than 1 byte are corrupt)
                if (sizeA>1 && sizeB>1)
                    if contains(allFiles(it),'.nii.gz')
                        % Unzip first
                        tmpPathImageA = xASL_adm_UnzipNifti(char(currentFileA));
                        tmpPathImageB = xASL_adm_UnzipNifti(char(currentFileB));
                        % Read in image
                        imageA = xASL_io_Nifti2Im(tmpPathImageA);
                        imageB = xASL_io_Nifti2Im(tmpPathImageB);
                        % GZIP Niftis afterwards
                        pathImageA = xASL_adm_GzipNifti(tmpPathImageA);
                        pathImageB = xASL_adm_GzipNifti(tmpPathImageB);
                    else
                        imageA = xASL_io_Nifti2Im(char(currentFileA));
                        imageB = xASL_io_Nifti2Im(char(currentFileB));
                    end
                    % Report function which prints to the console
                    if printReportFunction
                        % differenceAB = sum(imageA-imageB,'all');
                        RMSE = sqrt(mean(imageA(:).^2 - imageB(:).^2)) ;
                        if (RMSE>thresh)
                            fprintf('RMSE of NIFTIs above threshold!\n')
                            identical = false;
                        end
                    end
                else
                    if printReportFunction
                        fprintf('\t\t\t\tFile is too small to be a real image!\n');
                    end
                end
            end
        end
    end
end

%% Compare field lists
function compareFieldLists(jsonStructA,jsonStructB,fieldListA,fieldListB,printReportNow)

    % Iterate over fields
    for curField=1:numel(fieldListA)
        curFieldNameA = string(fieldListA(curField));
        curFieldNameB = string(fieldListB(curField));
        fieldContentA = jsonStructA.(fieldListA{curField});
        fieldContentB = jsonStructB.(fieldListB{curField});
        if isnumeric(fieldContentA) && isnumeric(fieldContentB)
            % Compare numbers
            if ~(fieldContentA==fieldContentB)
                if printReportNow
                    fprintf('\t\t\t\tMismatch of JSON field %s ...\n', curFieldNameA);
                end
            end
        elseif (ischar(fieldContentA) || isstring(fieldContentA)) && (ischar(fieldContentB) || isstring(fieldContentB))
            % Compare char arrays and strings
            if ~(strcmp(fieldContentA,fieldContentB))
                if printReportNow
                    fprintf('\t\t\t\tMismatch of JSON field %s ...\n', curFieldNameA);
                end
            end
        elseif iscell(fieldContentA) && iscell(fieldContentB)
            % Compare cell arrays
            if ~(isempty(setdiff(fieldContentA,fieldContentB)) && isempty(setdiff(fieldContentB,fieldContentA)))
                if printReportNow
                    fprintf('\t\t\t\tMismatch of JSON field %s ...\n', curFieldNameA);
                end
            end
        elseif isstruct(fieldContentA) && isstruct(fieldContentB)
            % Compare cell arrays
            if ~(isequal(fieldContentA,fieldContentB))
                if printReportNow
                    fprintf('\t\t\t\tMismatch of JSON field %s ...\n', curFieldNameA);
                end
            end
        else
            % Neither number nor text
            if printReportNow
                fprintf('\t\t\t\tIt seems we are comparing neither numbers nor text ...\n');
            end
        end
    end
    
end



