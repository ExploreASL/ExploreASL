function pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths,indexSortedFile,nameMerged,bAlternatingControlLabel)
%xASL_bids_MergeNifti_Merge Merge NiftiPaths & save to pathOut
%
% FORMAT: pathOut = xASL_bids_MergeNifti_Merge(NiftiPaths,indexSortedFile,nameMerged,bAlternatingControlLabel)
% 
% INPUT:
%   NiftiPaths - Nifti paths
%   indexSortedFile - Index sorted file
%   nameMerged - Name merged
%   bAlternatingControlLabel - Alternating control label (BOOLEAN, REQUIRED)
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
% Copyright 2015-2021 ExploreASL


    %% Track if all went well
    bStatus = 1; 
    pathOut = '';
    firstJSON = '';

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

        % Check for the path to JSON if existing and keep only the first existing JSON
        if isempty(firstJSON)
            [Fpath, Ffile] = xASL_fileparts(NiftiPaths{indexSortedFile(iFile)});
            pathJSON = fullfile(Fpath,[Ffile '.json']);
            if exist(pathJSON, 'file')
                firstJSON = pathJSON;
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
                    isHadamardFME = ~isempty(regexp(char(tmpCheckJSON.SeriesDescription),'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once'));
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
                    EchoTimes = sortrows(EchoTimes',1)';
                catch
                    fprintf('Sorting failed...\n');
                end
            end
            % Write changes to JSON
            structFirstJSON = spm_jsonread(firstJSON);
            structFirstJSON.EchoTime = EchoTimes;
            spm_jsonwrite(firstJSON, structFirstJSON);
        end
        % Copy the first JSON to this name
        if ~isempty(firstJSON)
            xASL_Copy(firstJSON,fullfile(Fpath,[nameMerged '.json']),1);
        end
    else
        fprintf('Warning: Cannot concatenate multiple NIfTIs & jsons as output from dcm2niiX\n');
    end


end
