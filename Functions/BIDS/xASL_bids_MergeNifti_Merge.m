function pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths, indexSortedFile, nameMerged, bAlternatingControlLabel, priorityList)
%xASL_bids_MergeNifti_Merge Merge NiftiPaths & save to pathOut
%
% FORMAT: pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths, indexSortedFile, nameMerged, bAlternatingControlLabel[, priorityList])
% 
% INPUT:
%   NiftiPaths - Nifti paths
%   indexSortedFile - Index sorted file
%   nameMerged - Name merged
%   bAlternatingControlLabel - Alternating control label (BOOLEAN, REQUIRED)
%   priorityList - a vector of priorities indicating (higher value = higher priority) how to merge the JSON files (OPTIONAL, DEFAULT = equal priorities)
%
% OUTPUT:
%   pathOut    - output paths
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Merge NiftiPaths & save to pathOut.
%
% EXAMPLE:     ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright 2015-2022 ExploreASL


    %% Track if all went well
    bStatus = 1; 
    pathOut = '';
	
	if nargin < 4
		error('Require 4 input parameters.');
	end
	
	if nargin < 5 || isempty(priorityList)
		priorityList = ones(1, numel(NiftiPaths));
	end
	
    %% Start loading all the files
    for iFile=1:length(NiftiPaths)
        tempIM = xASL_io_Nifti2Im(NiftiPaths{indexSortedFile(iFile)});

        % Merging only 3D and 4D files
        if length(size(tempIM))>4
            error('Dimensionality incorrect for this ASL NIfTI file');
        end

        % Compare size of the new file and if similar to the previous than concatenate, otherwise report an error
        if iFile == 1
            sizeFirst = size(tempIM);
            sizeFirst = sizeFirst(1:3);
        else
            sizeNew = size(tempIM);
            sizeNew = sizeNew(1:3);
            if ~isequal(sizeNew, sizeFirst)
                bStatus = 0;
            end
        end

        if bAlternatingControlLabel
            % Always interlace the two following files
            % For the first file, create the interleaved first volume
            if iFile == 1
                lengthFirst = size(tempIM,4);
                IM = zeros([sizeFirst,lengthFirst*2]);
                IM(:,:,:,1:2:end) = tempIM;
            elseif mod(iFile,2)
                % For odd files, create a new interleaved addition
                lengthFirst = size(tempIM,4);
                IM(:,:,:,end+1:end+2*lengthFirst) = zeros([sizeFirst,lengthFirst*2]);
                IM(:,:,:,end+2:2:end+2*lengthFirst) = tempIM;
            else
                % For even files - fill in the interleave spaces
                IM(:,:,:,end-2*lengthFirst+2:2:end) = tempIM;
            end
        else
            % Simply merge files in the order in which they come
            if iFile==1
                % Get the size of the first file
                IM = tempIM;
            else
                if bStatus
                    IM(:,:,:,end+1:end+size(tempIM,4)) = tempIM;
                end
            end
        end
        
	end

	%% Merge all the JSONs according to their priority
	
	% Create an order in which the JSONs are going to be read
	[~, priorityListIndex] = sort(priorityList);
	outputJSON = struct();
	for iJSON = 1:length(NiftiPaths)
		% Go through the JSONs in increasing priority
		[Fpath, Ffile] = xASL_fileparts(NiftiPaths{priorityListIndex(iJSON)});
		if xASL_exist(fullfile(Fpath,[Ffile '.json']), 'file')
			currentJSON = spm_jsonread(fullfile(Fpath,[Ffile '.json']));
			
			% Go through all fields
			listFieldNames = fieldnames(currentJSON);
			for iFieldName = 1:length(listFieldNames)
				% Overwrite fields as we are reading JSONs with increasing priority
				outputJSON.(listFieldNames{iFieldName}) = currentJSON.(listFieldNames{iFieldName});
			end
		end
	end
	
	
    %% If at the end and all went well
    if bStatus
        fprintf('Warning: concatenating multiple NIfTIs & jsons as output from dcm2niiX\n');
        % Save the concatenated file to a given name
        pathOut = fullfile(Fpath,[nameMerged '.nii']);
        xASL_io_SaveNifti(NiftiPaths{indexSortedFile(1)}, pathOut, IM, [], 0);
        % Special treatment for Hadamard encoded files
        EchoTimes = cell(size(NiftiPaths,2),1);
        for iFileCheck = 1:size(NiftiPaths,2)
            % Get JSON
            [jsonPathX, jsonNameX] = xASL_fileparts(NiftiPaths{iFileCheck});
            if exist(fullfile(jsonPathX, [jsonNameX '.json']),'file')
                tmpCheckJSON = spm_jsonread(fullfile(jsonPathX, [jsonNameX '.json']));
                % Check EchoTimes
                if isfield(tmpCheckJSON,'SeriesDescription')
                    isHadamardFME = ~isempty(regexp(char(tmpCheckJSON.SeriesDescription),...
                                    '(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once')) ... % ASL format
                                 || ~isempty(regexp(char(tmpCheckJSON.SeriesDescription),...
                                    '(ss_TE)\d\d(_TI)\d\d\d\d', 'once'));                            % M0 format
                    if isHadamardFME
                        if isfield(tmpCheckJSON,'EchoTime')
                            EchoTimes{iFileCheck,1} = tmpCheckJSON.EchoTime;
                        end
                    end
                end
            end
        end
        % Add echo number array if it exists
        if sum(~cellfun(@isempty,EchoTimes))~=0
            fprintf('Merging the echo numbers of the Hadamard encoded sequence...\n');
            % Sort echo numbers
            if length(indexSortedFile)==length(EchoTimes)
                EchoTimesBackUp = EchoTimes;
                for iEchoNumber=1:length(EchoTimes)
                    EchoTimes(indexSortedFile(iEchoNumber)) = EchoTimesBackUp(iEchoNumber,1);
                end
            end
            if ~issorted(cell2mat(EchoTimes))
                fprintf('Warning: echo times do not increase, resorting will be applied...\n');
                try
                    EchoTimes = sortrows(EchoTimes,1);
                catch
                    fprintf('Sorting failed...\n');
                end
			end
			
            % Write changes to JSON
			outputJSON.EchoTime = EchoTimes;
		end
		spm_jsonwrite(fullfile(Fpath,[nameMerged '.json']), outputJSON);
    else
        fprintf('Warning: Cannot concatenate multiple NIfTIs & jsons as output from dcm2niiX\n');
    end

end
