function [identical,results] = xASL_bids_CompareStructures(pathDatasetA,pathDatasetB,bPrintReport,threshRmseNii)
%xASL_bids_CompareStructures Function that compares two BIDS folders with several subfolders and studies and prints the differences.
%
% FORMAT: [identical,results] = xASL_bids_CompareStructures(pathDatasetA,pathDatasetB,[bPrintReport,threshRmseNii]);
%
% INPUT:
%        pathDatasetA       - path to first BIDS structure [char array] (REQUIRED)
%        pathDatasetB       - path to second BIDS structure [char array] (REQUIRED)
%        bPrintReport       - true or false to print console report (OPTIONAL, DEFAULT = true)
%        threshRmseNii      - normalized RMSE threshold for comparing NIFTI content (OPTIONAL, DEFAULT = 0.01)
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
%                   [identical,results] = xASL_bids_CompareStructures(pathDatasetA,pathDatasetB,true);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2020 ExploreASL


    %% Input Check

    % Check if both root folders are valid char arrays or strings
    if ~(ischar(pathDatasetA))
        error('The path of structure A is neither a char array not a string.');
    end
    if ~(ischar(pathDatasetB))
        error('The path of structure A is neither a char array not a string.');
    end

    % Check if both root folders exists
    if ~(xASL_exist(pathDatasetA)==7)
        error('The root folder of structure A does not exist: %s',pathDatasetA);
    end
    if ~(xASL_exist(pathDatasetB)==7)
        error('The root folder of structure B does not exist: %s',pathDatasetB);
    end
	
    % Default value for bPrintReport
	if nargin < 3 || isempty(bPrintReport)
		bPrintReport = true;
	end
    
    % Default value for RMSE threshold
    if nargin < 4
       threshRmseNii = 0.01; 
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
    if bPrintReport
        fprintf(strcat(repmat('=',100,1)','\n'));
        fprintf('Dataset:\t\t%s\n',datasetA)
        printList(results.(datasetA).missingFolders)
        printList(results.(datasetA).missingFiles)
		
		if isempty(results.(datasetA).missingFolders) && isempty(results.(datasetA).missingFiles)
			fprintf('\t\t\t%s\n','No missing files');
		end

        fprintf(strcat(repmat('=',100,1)','\n'));
        fprintf('Dataset:\t\t%s\n',datasetB)
        printList(results.(datasetB).missingFolders)
        printList(results.(datasetB).missingFiles)
		
		if isempty(results.(datasetB).missingFolders) && isempty(results.(datasetB).missingFiles)
			fprintf('\t\t\t%s\n','No missing files');
		end

        % End of report
        fprintf(strcat(repmat('=',100,1)','\n'));
    end
    
    % Compare file content
    identical = checkFileContents(fileListA,fileListB,pathDatasetA,pathDatasetB,identical,bPrintReport,threshRmseNii);

    
    

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
            fprintf('Missing:\t\t%s\n',currentList(iFile))
        end
    end
end

%% Check file contents
function identical = checkFileContents(filesDatasetA,filesDatasetB,pathDatasetA,pathDatasetB,identical,bPrintReport,threshRmseNii)
    
    % All files
    allFiles = unique([filesDatasetA',filesDatasetB']');
    
    % Iterate over list
    for iFile=1:length(allFiles)
        % Assign root directory of dataset A
        currentFileA = fullfile(pathDatasetA,allFiles(iFile));
        currentFileB = fullfile(pathDatasetB,allFiles(iFile));
        % Get extension
        [~,~,extension] = fileparts(allFiles(iFile));
        % Check extension
        if strcmp(extension,'.json')
			jsonErrorReport='';

			% Compare JSON files on field basis
			if (isfile(currentFileA) && isfile(currentFileB)) % xASL_exist somehow didn't work here (again)
				% Import JSON files
				jsonA = xASL_import_json(char(currentFileA));
				jsonB = xASL_import_json(char(currentFileB));
				
				% Get JSON field names
				fieldNamesA = fieldnames(jsonA);
				fieldNamesB = fieldnames(jsonB);
				
				% Check which fields are shared and which different
				sharedFieldsAB = intersect(fieldNamesB,fieldNamesA);

				% Fields that are in B, but missing in A
				missingFields = setdiff(fieldNamesB,fieldNamesA);
				% Print out missing fields
				for iField=1:numel(missingFields)
					jsonErrorReport = sprintf('%s\t\t\t\tMissing field: %s\n',jsonErrorReport,string(missingFields{iField}));
				end
				
				extraFields = setdiff(fieldNamesA,fieldNamesB);
				% Print out missing fields
				for iField=1:numel(extraFields)
					jsonErrorReport = sprintf('%s\t\t\t\tExtra field: %s\n',jsonErrorReport,string(extraFields{iField}));
				end
				
				% Now we can compare these fields like in the part above
				jsonErrorReport = [jsonErrorReport, compareFieldLists(jsonA,jsonB,sharedFieldsAB)];
				
				if bPrintReport && ~isempty(jsonErrorReport)
					fprintf('%s:\n',allFiles(iFile));
					fprintf('%s',jsonErrorReport);
				end
			end
        elseif strcmp(extension,'.tsv') || strcmp(extension,'.txt') || strcmp(extension,'.csv')
            % Read files if they exist
            if (isfile(currentFileA) && isfile(currentFileB)) % xASL_exist somehow didn't work here (again)
                % Compare text files content directly
                currentFileTextA = fileread(currentFileA);
                currentFileTextB = fileread(currentFileB);
                if bPrintReport
                    if ~strcmp(currentFileTextA,currentFileTextB)
						fprintf('%s:\t\t\n',allFiles(iFile));
                        fprintf('\t\t\t\tDifferent file content.\n');
                        identical = false;
                    end
                end
            end
        elseif strcmp(extension,'.nii') || contains(allFiles(iFile),'.nii.gz') % Warning: This unzips your NIFTIs right now!
            % Read files if they exist
            if (isfile(currentFileA) && isfile(currentFileB)) % xASL_exist somehow didn't work here (again)
                % Check file size (there were some 0KB images in the bids-examples)
                tmpFileA = dir(currentFileA);
                tmpFileB = dir(currentFileB);
                sizeA = tmpFileA.bytes;
                sizeB = tmpFileB.bytes;

                % Check file size (images with a file size lower than 1 byte are corrupt)
                if (sizeA>1 && sizeB>1)
                    if contains(allFiles(iFile),'.nii.gz')
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
                    % Report function which prints to the console
                    if bPrintReport
                        % differenceAB = sum(imageA-imageB,'all');
                        RMSE = sqrt(mean((imageA(:) - imageB(:)).^2))*2/sqrt(mean(abs(imageA(:)) + abs(imageB(:))).^2);
                        if (RMSE>threshRmseNii)
							fprintf('%s:\t\t\n',allFiles(iFile));
                            fprintf('\t\t\t\tRMSE of NIFTIs above threshold.\n');
                            identical = false;
                        end
                    end
                else
                    if bPrintReport
						fprintf('%s:\t\t\n',allFiles(iFile));
                        fprintf('\t\t\t\tFile is too small to be a real image.\n');
                    end
                end
            end
        end
    end
end

%% Compare field lists
function strError = compareFieldLists(jsonStructA,jsonStructB,fieldList)
    strError = '';

    % Iterate over fields
    for iField=1:numel(fieldList)
        curFieldName = string(fieldList(iField));
        fieldContentA = jsonStructA.(fieldList{iField});
        fieldContentB = jsonStructB.(fieldList{iField});
        if isnumeric(fieldContentA) && isnumeric(fieldContentB)
            % Compare numbers
            if ~isequal(fieldContentA,fieldContentB)
                strError = sprintf('%s\t\t\t\tDifferent value: %s (%f vs %f)\n', strError,curFieldName,fieldContentA,fieldContentB);
            end
        elseif ischar(fieldContentA) && ischar(fieldContentB)
            % Compare char arrays and strings
            if ~(strcmp(fieldContentA,fieldContentB))
                strError = sprintf('%s\t\t\t\tDifferent value: %s (%s vs %s)\n', strError,curFieldName,fieldContentA,fieldContentB);
            end
        elseif iscell(fieldContentA) && iscell(fieldContentB)
            % Compare cell arrays
            if ~(isempty(setdiff(fieldContentA,fieldContentB)) && isempty(setdiff(fieldContentB,fieldContentA)))
                strError = sprintf('%s\t\t\t\tDifferent value: %s (array)\n', strError,curFieldName);
            end
        elseif isstruct(fieldContentA) && isstruct(fieldContentB)
			% Compare cell arrays
			if ~(isequal(fieldContentA,fieldContentB))
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
