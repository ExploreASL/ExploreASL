function [x,nii_files, summary_line, globalCounts, ASLContext] = xASL_imp_DCM2NII_Subject_ShuffleTheDynamics(x,globalCounts, scanpath, scan_name, nii_files, iSubject, iSession, iScan)
%xASL_imp_DCM2NII_Subject_ShuffleTheDynamics Shuffle the dynamics.
%
% FORMAT: [nii_files, summary_line, globalCounts, ASLContext] = xASL_imp_DCM2NII_Subject_ShuffleTheDynamics(globalCounts, scanpath, scan_name, nii_files, iSubject, iSession, iScan)
% 
% INPUT:
%   x               - ExploreASL x struct (STRUCT, REQUIRED)
%   globalCounts    - Converted, skipped & missing scans (REQUIRED, STRUCT)
%   scanpath        - Scan path (CHAR ARRAY, PATH, REQUIRED)
%   scan_name       - Scan name (CHAR ARRAY, REQUIRED)
%   nii_files       - List of NIfTI files (CELL ARRAY, REQUIRED)
%   iSubject        - Subject ID (INTEGER, REQUIRED)
%   iSession        - Session ID (INTEGER, REQUIRED)
%   iScan           - Scan ID (INTEGER, REQUIRED)
%
% OUTPUT:
%   x               - ExploreASL x struct (STRUCT)
%   nii_files       - List of NIfTI files (CELL ARRAY)
%   summary_line    - Summary line
%   globalCounts    - Converted, skipped & missing scans (REQUIRED, STRUCT)
%   ASLContext      - ASL context text ('deltam,m0scan' e.g.)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Shuffle the dynamics.
%
% 1. Fallbacks
% 2. Fill NIfTI Table
% 3. Get ASL context if possible
% 4. Only try shuffling if you dont know the ASL context already
% 5. Merge NIfTIs if there are multiples for ASL or M0, merge multiple files
% 6. Extract relevant parameters from nifti header and append to summary file
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     n/a
%
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% 1. Fallbacks

    % Fallback
    ASLContext = '';
    
    % Fallback niiTable
    niiTable = cell(size(nii_files,2),size(nii_files,1)*4);
    
    %% 2. Fill NIfTI Table
    % The idea is to sort the NIfTIs based on the InstanceNumbers, to make sure that we have the same
    % ASL context on every OS. Before doing this part, there were some differences between Linux & Windows.

    % Iterate over NIfTI files
    for iNii = 1:size(niiTable,1)
    	% Get current NIfTI path
        filePath = nii_files{1,iNii};
        [rootName, fileName] = xASL_fileparts(filePath);
        niiTable{iNii,1} = fileName;
        % Open corresponding JSON file
        tmpJSON = spm_jsonread(fullfile(rootName,[fileName '.json']));
        % Try to extract the InstanceNumber field from the JSON file
        if isfield(tmpJSON,'InstanceNumber')
            niiTable{iNii,2} = xASL_str2num(tmpJSON.InstanceNumber);
        else
            niiTable{iNii,2} = 0;
        end
        % Check the Manufacturer
        if isfield(tmpJSON, 'Manufacturer')
        	% If we have GE scans, we use the xASL_bids_determineImageTypeGE script to determine the ImageType
            if ~isempty(strfind(tmpJSON.Manufacturer, 'GE'))
                niiTable{iNii,3} = xASL_bids_determineImageTypeGE(tmpJSON);
                if isempty(niiTable{iNii,3})
                    niiTable{iNii,3} = NaN;
                end
            else
                niiTable{iNii,3} = NaN;
            end
        else
            niiTable{iNii,3} = NaN;
        end
        % Add the corresponding file path to the Table
        niiTable{iNii,4} = filePath;
    end
    
    % Fallback
    if ~exist('tmpJSON','var')
        tmpJSON = struct;
    end
    
    %% 3. Get ASL context if possible
    
    % Check if we can sort by validating that we have InstanceNumbers
    if sum(cellfun(@isnumeric, niiTable(:,3)))==0 
    	% Sort table based on InstanceNumbers
        niiTable = sortrows(niiTable,2);
        ASLContext = '';
        % Build the ASL context string
        for iImageType = 1:length(niiTable(:,3)')
            ASLContext = [ASLContext ',' niiTable{iImageType,3}];
        end
        % Remove one comma
        ASLContext = ASLContext(2:end);
        fprintf('ASL context: %s\n', ASLContext);
        % Sort nii_files based on ASLContext/InstanceNumbers
        nii_files = niiTable(:,4)';
        fprintf('Sorted NIfTIs based on InstanceNumbers...\n');
    end
    
    %% 4. Only try shuffling if you dont know the ASL context already
    
    % Check if it was not possible to determine the ASL context before
    if isempty(ASLContext)
        [~,~,scanExtension] = xASL_fileparts(scanpath);
        if ~isempty(regexpi(scanExtension, '^\.(par|rec)$')) && length(nii_files)==1 && ~isempty(regexpi(scan_name, 'ASL'))
            % For a PAR/REC files that produces a single ASL4D NIFTI
            imASL = xASL_io_Nifti2Im(nii_files{1});
            % If multiple dynamics
            if size(imASL,4) > 1
                % Then reshuffle them
                imASLreordered = zeros(size(imASL));
                imASLreordered(:,:,:,1:2:end) = imASL(:,:,:,1:ceil(size(imASL,4)/2));
                imASLreordered(:,:,:,2:2:end) = imASL(:,:,:,ceil(size(imASL,4)/2)+1:end);
                xASL_io_SaveNifti(nii_files{1}, nii_files{1}, imASLreordered);
            end
        end
    end
    
    %% 5. Merge NIfTIs if there are multiples for ASL or M0, merge multiple files
    if length(nii_files)>1
        if ~isempty(strfind(scan_name,'ASL4D'))
            [nii_files,ASLContext] = xASL_bids_MergeNifti(nii_files, 'ASL');
        elseif  ~isempty(strfind(scan_name,'M0'))
            [nii_files,ASLContext] = xASL_bids_MergeNifti(nii_files, 'M0');
        end
    end
    
    %% Hadamard Check
    x.bHadamard = false;
    isHadamardFME = false;
    
    % Determine if we have a Hadamard sequence based on the parameters of the studyPar.json
    if xASL_exist(x.dir.studyPar,'file')
        studyPar = spm_jsonread(x.dir.studyPar);
        if isfield(studyPar,'HadamardMatrixType')
            x.bHadamard = true;
        end
    end
    
    
    %% Check for specific Hadamard sequences and reorder the NIfTI volumes if necessary
    if numel(nii_files)>=1 
        [resultPath, resultFile] = xASL_fileparts(nii_files{1});
        % Check if we have the corresponding JSON file
        if exist(fullfile(resultPath, [resultFile '.json']), 'file') 
        	% Load the JSON
            resultJSON = spm_jsonread(fullfile(resultPath, [resultFile '.json']));
            % Check if we have the SeriesDescription field
            if isfield(resultJSON,'SeriesDescription') || x.bHadamard
            	% Determine if we have the specific FME Hadamard sequence from Bremen
                isHadamardFME = ~isempty(regexp(char(resultJSON.SeriesDescription),'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once'));
                % If the FME sequence was detected we can always set the general bHadamard to true as well
                if isHadamardFME
                    x.bHadamard = true;
                end
            end
            % Check if we the current sequence is a Hadamard or not
            if x.bHadamard || isHadamardFME
                % Check image
                if xASL_exist(nii_files{1},'file')
                    % Determine the number of time points within each NIfTI
                    imASL = xASL_io_Nifti2Im(nii_files{1});
                    numberTEs = length(resultJSON.EchoTime);
                    numberPLDs = int32(size(imASL,4)/numberTEs);
                    
                    % Reorder TEs and PLDs - first cycle TE afterwards PLD
                    vectorOldOrder = zeros(size(imASL,4),1);
                    for iPLD = 1:(double(numberPLDs))
                        vectorOldOrder((1:numberTEs)+(iPLD-1)*numberTEs) = (iPLD-1)+1:numberPLDs:size(imASL,4);
                    end
                    imASL(:,:,:,1:end) = imASL(:,:,:,vectorOldOrder);
                    xASL_io_SaveNifti(nii_files{1},nii_files{1},imASL);
                    % Repeat Echo Times
                    resultJSON.EchoTime = repmat(resultJSON.EchoTime,numberPLDs,1);
                    % Save the JSON with the updated echo times
                    spm_jsonwrite(fullfile(resultPath, [resultFile '.json']),resultJSON);
                else
                    % Fallback (something went wrong)
                    warning('Hadamard sequence with 1 PLD only');
                end
            end
        end
    end
    
    
    %% 6. Extract relevant parameters from nifti header and append to summary file
    summary_line = xASL_imp_AppendNiftiParameters(nii_files);
    globalCounts.converted_scans(iSubject, iSession, iScan) = 1;
    

end