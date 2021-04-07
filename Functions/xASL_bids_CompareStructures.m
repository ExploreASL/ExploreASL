function [identical,results] = xASL_bids_CompareStructures(pathDatasetA,pathDatasetB,bPrintReport,threshRmseNii,detailedOutput,printWarnings)
%xASL_bids_CompareStructures Function that compares two BIDS folders with several subfolders and studies and prints the differences.
%
% FORMAT: [identical,results] = xASL_bids_CompareStructures(pathDatasetA,pathDatasetB,[bPrintReport,threshRmseNii]);
%
% INPUT:
%        pathDatasetA       - path to first BIDS structure [char array] (REQUIRED)
%        pathDatasetB       - path to second BIDS structure [char array] (REQUIRED)
%        bPrintReport       - true or false to print console report (OPTIONAL, DEFAULT = true)
%        threshRmseNii      - normalized RMSE threshold for comparing NIFTI content (OPTIONAL, DEFAULT = 1e-5)
%        detailedOutput     - additional text ouput (also print that there are no missing files etc.) (OPTIONAL, DEFAULT = false)
%        printWarnings      - print differences as warnings (OPTIONAL, DEFAULT = true)
%
% OUTPUT:
%        identical          - Returns 1 if both folder structures are identical and 0 if not
%        results            - structure containing (possible) differences of both folder structures
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
% Copyright @ 2015-2021 ExploreASL


    %% Input Check

    % Check if both root folders are valid char arrays or strings
    if ~ischar(pathDatasetA)
        error('The path of structure A is not a char array');
    end
    if ~ischar(pathDatasetB)
        error('The path of structure B is not a char array');
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
    [~, datasetA] = fileparts(pathDatasetA);
    [~, datasetB] = fileparts(pathDatasetB);

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
        if detailedOutput
            fprintf(strcat(repmat('=',100,1)','\n'));
        end
        if detailedOutput || ~isempty(results.(datasetA).missingFolders) || ~isempty(results.(datasetA).missingFiles)
            fprintf('Dataset:\t\t%s\n',datasetA)
        end
        printList(results.(datasetA).missingFolders)
        printList(results.(datasetA).missingFiles)
        
        if detailedOutput
            if isempty(results.(datasetA).missingFolders) && isempty(results.(datasetA).missingFiles)
                fprintf('\t\t\t\t%s\n','No missing files');
            end
        end
        
        if detailedOutput
            fprintf(strcat(repmat('=',100,1)','\n'));
        end
        if detailedOutput || ~isempty(results.(datasetB).missingFolders) || ~isempty(results.(datasetB).missingFiles)
            fprintf('Dataset:\t\t%s\n',datasetB)
        end
        printList(results.(datasetB).missingFolders)
        printList(results.(datasetB).missingFiles)
        
        if detailedOutput
            if isempty(results.(datasetB).missingFolders) && isempty(results.(datasetB).missingFiles)
                fprintf('\t\t\t\t%s\n','No missing files');
            end
        end
        
        % End of report
        if detailedOutput
            fprintf(strcat(repmat('=',100,1)','\n'));
        end
    end
    
    
    % Compare file content
    [identical,results.differences] = checkFileContents(fileListA,fileListB,pathDatasetA,pathDatasetB,identical,bPrintReport,threshRmseNii);
    
    % Print differences as warnings
    if printWarnings
        printMissingAsWarnings(results);
        printDifferencesAsWarnings(results.differences);
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


end

%% Modify file lists
function fileList = modifyFileList(fileList,root)
    % Iterate over file list: change folder names
    for iFile=1:numel(fileList)
        fileList(iFile).folder = strrep(fileList(iFile).folder,root,'');
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
            fprintf('Missing:\t\t%s\n',currentList{iFile})
        end
    end
end

%% Check file contents
function [identical,differences] = checkFileContents(filesDatasetA,filesDatasetB,pathDatasetA,pathDatasetB,identical,bPrintReport,threshRmseNii)
    
    % All files
    allFiles = unique([filesDatasetA',filesDatasetB']');
    
    % Differences
    differences = cell(1,1);
    
    % Difference number
    dn = 1;
    
    % Iterate over list
    for iFile=1:length(allFiles)
        % Assign root directory of dataset A
        currentFileA = fullfile(pathDatasetA,allFiles{iFile});
        currentFileB = fullfile(pathDatasetB,allFiles{iFile});
        % Get extension
        [~,~,extension] = fileparts(allFiles{iFile});
        % Check extension
        if strcmp(extension,'.json')
			jsonErrorReport='';

			% Compare JSON files on field basis
			if (isfile(currentFileA) && isfile(currentFileB)) % xASL_exist somehow didn't work here (again)
				% Import JSON files
				jsonA = spm_jsonread(char(currentFileA));
				jsonB = spm_jsonread(char(currentFileB));
				
				% Get JSON field names
				fieldNamesA = fieldnames(jsonA);
				fieldNamesB = fieldnames(jsonB);
				
				% Check which fields are shared and which different
				sharedFieldsAB = intersect(fieldNamesB,fieldNamesA);

				% Fields that are in B, but missing in A
				missingFields = setdiff(fieldNamesB,fieldNamesA);
				% Print out missing fields
				for iField=1:numel(missingFields)
					jsonErrorReport = sprintf('%s\t\t\t\tMissing field: %s\n',jsonErrorReport,missingFields{iField});
				end
				
				extraFields = setdiff(fieldNamesA,fieldNamesB);
				% Print out missing fields
				for iField=1:numel(extraFields)
					jsonErrorReport = sprintf('%s\t\t\t\tExtra field: %s\n',jsonErrorReport,extraFields{iField});
				end
				
				% Now we can compare these fields like in the part above
				jsonErrorReport = [jsonErrorReport, compareFieldLists(jsonA,jsonB,sharedFieldsAB)];
				
				if ~isempty(jsonErrorReport)
                    if bPrintReport
                        fprintf('File:\t\t\t%s\n',allFiles{iFile});
                        fprintf('%s',jsonErrorReport);
                    end
                    
                    % Save difference
                    differences{dn,1} = ['Different file content: ', allFiles{iFile}, ' '];
                    dn = dn+1;
                end
                
			end
        elseif strcmp(extension,'.tsv') || strcmp(extension,'.txt') || strcmp(extension,'.csv')
            % Read files if they exist
            if (isfile(currentFileA) && isfile(currentFileB)) % xASL_exist somehow didn't work here (again)
                % Compare text files content directly
                currentFileTextA = fileread(currentFileA);
                currentFileTextB = fileread(currentFileB);
                if ~strcmp(currentFileTextA,currentFileTextB)
                    if bPrintReport
                        fprintf('%s:\t\t\n',allFiles{iFile});
                        fprintf('\t\t\t\tDifferent file content.\n');
                    end
                    identical = false;
                    
                    % Save difference
                    differences{dn,1} = ['Different file content: ', allFiles{iFile}, ' '];
                    dn = dn+1;
                    
                end
            end
        elseif strcmp(extension,'.nii') || contains(allFiles{iFile},'.nii.gz') % Warning: This unzips your NIFTIs right now!
            % Read files if they exist
            if (isfile(currentFileA) && isfile(currentFileB)) % xASL_exist somehow didn't work here (again)
                % Check file size (there were some 0KB images in the bids-examples)
                tmpFileA = dir(currentFileA);
                tmpFileB = dir(currentFileB);
                sizeA = tmpFileA.bytes;
                sizeB = tmpFileB.bytes;
                
                % Check file size (images with a file size lower than 1 byte are corrupt)
                if (sizeA>1 && sizeB>1)
                    if contains(allFiles{iFile},'.nii.gz')
                        % Unzip first
                        tmpPathImageA = xASL_adm_UnzipNifti(char(currentFileA));
                        tmpPathImageB = xASL_adm_UnzipNifti(char(currentFileB));
                        % Read in image
                        imageA = xASL_io_Nifti2Im(tmpPathImageA);
                        imageB = xASL_io_Nifti2Im(tmpPathImageB);
                        % GZIP Niftis afterwards
                        xASL_adm_GzipNifti(tmpPathImageA);
                        xASL_adm_GzipNifti(tmpPathImageB);
                    else
                        imageA = xASL_io_Nifti2Im(char(currentFileA));
                        imageB = xASL_io_Nifti2Im(char(currentFileB));
                    end
                    
                    % Check that matrix dimensions agree first
                    sizeA = size(imageA);
                    sizeB = size(imageB);
                    if length(sizeA)==length(sizeB)
                        if sum(sizeA==sizeB)==length(sizeA)
                            RMSE = sqrt(mean((imageA(:) - imageB(:)).^2))*2/sqrt(mean(abs(imageA(:)) + abs(imageB(:))).^2);
                            if (RMSE>threshRmseNii)
                                % Report function which prints to the console
                                if bPrintReport
                                    fprintf('File:\t\t\t%s\n',allFiles{iFile});
                                    fprintf('\t\t\t\tRMSE (%d) of NIFTIs above threshold.\n',RMSE);
                                end
                                identical = false;
                                
                                % Save difference
                                differences{dn,1} = ['RMSE of NIFTIs above threshold: ', allFiles{iFile}, ' '];
                                dn = dn+1;
                            end
                        else
                            if bPrintReport
                                fprintf('File:\t\t\t%s\n',allFiles{iFile});
                                fprintf('\t\t\t\tMatrix dimensions do not agree.\n');
                            end
                            identical = false;
                            
                            % Save difference
                            differences{dn,1} = ['Matrix dimensions do not agree: ', allFiles{iFile}, ' '];
                            dn = dn+1;
                        end
                    else
                        if bPrintReport
                            fprintf('File:\t\t\t%s\n',allFiles{iFile});
                            fprintf('\t\t\t\tMatrix dimensions do not agree.\n');
                        end
                        identical = false;
                        
                        % Save difference
                        differences{dn,1} = ['Matrix dimensions do not agree: ', allFiles{iFile}, ' '];
                        dn = dn+1;
                    end
                else
                    if bPrintReport
                        fprintf('%s:\t\t\n',allFiles{iFile});
                        fprintf('\t\t\t\tFile is too small to be a real image.\n');
                    end
                    % Save difference
                    differences{dn,1} = ['File is too small to be a real image: ', allFiles{iFile}, ' '];
                    dn = dn+1;
                end
            end
        end
    end
end

%% Compare field lists
function strError = compareFieldLists(jsonStructA,jsonStructB,fieldList)
    strError = '';
    
    % Threshold for the difference of numeric values
    threshNumeric = 1e-5;

    % Iterate over fields
    for iField=1:numel(fieldList)
        curFieldName = fieldList{iField};
        fieldContentA = jsonStructA.(fieldList{iField});
        fieldContentB = jsonStructB.(fieldList{iField});
        if isnumeric(fieldContentA) && isnumeric(fieldContentB)
            % Compare numbers
            if length(fieldContentA)==length(fieldContentB)
                if length(fieldContentA)==1
                    % Compare numbers (check absolute difference)
                    if abs(fieldContentA-fieldContentB)>threshNumeric
                        strError = sprintf('%s\t\t\t\tDifferent value: %s (%.6f vs %.6f)\n', strError,curFieldName,fieldContentA,fieldContentB);
                    end
                else
                    % Compare arrays (check sum of absolute differences)
                    sumDiff = sum(abs(fieldContentA-fieldContentB));
                    if sumDiff>threshNumeric
                        strError = sprintf('%s\t\t\t\tDifferent value: %s (check arrays)\n', strError,curFieldName);
                        % Set max number of elements to display
                        maxNumElements = 10;
                        if length(fieldContentA)<maxNumElements
                            maxNumElements = length(fieldContentA);
                        end
                        % Initialize array view
                        strError = sprintf('%s\t\t\t\t[',strError);
                        % Iterate over individual elements of array A
                        for elField=1:length(fieldContentA)
                            if elField<=maxNumElements
                                if elField<maxNumElements
                                    strError = sprintf('%s%.4f, ',strError,fieldContentA(elField));
                                elseif elField==length(fieldContentA)
                                    strError = sprintf('%s%.4f]\n',strError,fieldContentA(elField));
                                elseif elField==maxNumElements
                                    strError = sprintf('%s%.4f ...]\n',strError,fieldContentA(elField));
                                end
                            end
                        end
                        % Initialize array view
                        strError = sprintf('%s\t\t\t\t[',strError);
                        % Iterate over individual elements of array B
                        for elField=1:length(fieldContentB)
                            if elField<=maxNumElements
                                if elField<maxNumElements
                                    strError = sprintf('%s%.4f, ',strError,fieldContentB(elField));
                                elseif elField==length(fieldContentB)
                                    strError = sprintf('%s%.4f]\n',strError,fieldContentB(elField));
                                elseif elField==maxNumElements
                                    strError = sprintf('%s%.4f ...]\n',strError,fieldContentB(elField));
                                end
                            end
                        end
                    end
                end
            else
                strError = sprintf('%s\t\t\t\tDifferent dimension: %s\n', strError,curFieldName);
            end
        elseif ischar(fieldContentA) && ischar(fieldContentB)
            % Compare char arrays and strings
            if ~strcmp(fieldContentA,fieldContentB)
                strError = sprintf('%s\t\t\t\tDifferent value: %s (%s vs %s)\n', strError,curFieldName,fieldContentA,fieldContentB);
            end
        elseif iscell(fieldContentA) && iscell(fieldContentB)
            % Compare cell arrays
            if ~(isempty(setdiff(fieldContentA,fieldContentB)) && isempty(setdiff(fieldContentB,fieldContentA)))
                strError = sprintf('%s\t\t\t\tDifferent value: %s (array)\n', strError,curFieldName);
            end
        elseif isstruct(fieldContentA) && isstruct(fieldContentB)
			% Compare cell arrays
			if ~isequal(fieldContentA,fieldContentB)
				strError = sprintf('%s\t\t\t\tDifferent value: %s (array)\n', strError,curFieldName);
			end
		elseif islogical(fieldContentA) && islogical(fieldContentB)
            % Compare numbers
            if ~(fieldContentA==fieldContentB)
				if fieldContentA
					fieldContentA = 'true';
				else
					fieldContentA = 'false';
				end
				
				if fieldContentB
					fieldContentB = 'true';
				else
					fieldContentB = 'false';
				end
				
                strError = sprintf('%s\t\t\t\tDifferent value: %s (%s vs %s)\n', strError,curFieldName,fieldContentA,fieldContentB);
            end			
        else
            % Neither number nor text
			if ~isequal(fieldContentA,fieldContentB)
				strError = sprintf('%s\t\t\t\tDifferent value: %s (%s vs %s) - unknown or differing types\n', strError,curFieldName,fieldContentA,fieldContentB);
			end
        end
    end
    
end
