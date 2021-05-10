function [nii_files, summary_line, globalCounts, ASLContext] = xASL_imp_DCM2NII_Subject_ShuffleTheDynamics(globalCounts, scanpath, scan_name, nii_files, iSubject, iSession, iScan)
%xASL_imp_DCM2NII_Subject_ShuffleTheDynamics Shuffle the dynamics.
%
% FORMAT: [nii_files, summary_line, globalCounts, ASLContext] = xASL_imp_DCM2NII_Subject_ShuffleTheDynamics(globalCounts, scanpath, scan_name, nii_files, iSubject, iSession, iScan)
% 
% INPUT:
%   globalCounts    - Converted, skipped & missing scans (REQUIRED, STRUCT)
%   scanpath        - Scan path
%   scan_name       - Scan name
%   nii_files       - List of NIfTI files (CELL ARRAY, REQUIRED)
%   iSubject        - Subject ID (INTEGER, REQUIRED)
%   iSession        - Session ID (INTEGER, REQUIRED)
%   iScan           - Scan ID (INTEGER, REQUIRED)
%
% OUTPUT:
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
    
    for iNii = 1:size(niiTable,1)
        filePath = nii_files{1,iNii};
        [rootName, fileName, ~] = xASL_fileparts(filePath);
        niiTable{iNii,1} = fileName;
        % Open JSON file
        tmpJSON = spm_jsonread(fullfile(rootName,[fileName '.json']));
        if isfield(tmpJSON,'InstanceNumber')
            niiTable{iNii,2} = xASL_str2num(tmpJSON.InstanceNumber);
        else
            niiTable{iNii,2} = 0;
        end
        if isfield(tmpJSON, 'Manufacturer')
            if ~isempty(strfind(tmpJSON.Manufacturer, 'GE'))
                niiTable{iNii,3} = xASL_bids_determineImageTypeGE(tmpJSON);
            else
                niiTable{iNii,3} = NaN;
            end
        else
            niiTable{iNii,3} = NaN;
        end
        niiTable{iNii,4} = filePath;
    end
    
    % Fallback
    if ~exist('tmpJSON','var')
        tmpJSON = struct;
    end
    
    %% 3. Get ASL context if possible
    
    if sum(cellfun(@isnumeric, niiTable(:,3)))==0 
        niiTable = sortrows(niiTable,2);
        ASLContext = '';
        for iImageType = 1:length(niiTable(:,3)')
            ASLContext = [ASLContext ',' niiTable{iImageType,3}];
        end
        ASLContext = ASLContext(2:end);
        fprintf('ASL context: %s\n', ASLContext);
        % Sort nii_files based on ASLContext/InstanceNumbers
        nii_files = niiTable(:,4)';
        fprintf('Sorted NIfTIs based on InstanceNumbers...\n');
    end
    
    %% 4. Only try shuffling if you dont know the ASL context already
    
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
    
    %% FME (Hadamard) Check
    if numel(nii_files)>=1
        [resultPath, resultFile] = xASL_fileparts(nii_files{1});
        if exist(fullfile(resultPath, [resultFile '.json']), 'file')
            resultJSON = spm_jsonread(fullfile(resultPath, [resultFile '.json']));
            if isfield(resultJSON,'SeriesDescription')
                isHadamardFME = ~isempty(regexp(char(resultJSON.SeriesDescription),'(Encoded_Images_Had)\d\d(_)\d\d(_TIs_)\d\d(_TEs)', 'once'));
                if isHadamardFME
                    % Check image
                    if exist(nii_files{1} ,'file')
                        imASL = xASL_io_Nifti2Im(nii_files{1});
                        numOfTimePoints = size(imASL,4)/length(resultJSON.EchoTime);
                    else
                        numOfTimePoints = 1;
                    end
                    % Repeat Echo Times
                    resultJSON.EchoTime = repmat(resultJSON.EchoTime,numOfTimePoints,1);
                end
            end
        end
    end
    
    
    %% 6. Extract relevant parameters from nifti header and append to summary file
    summary_line = xASL_imp_AppendNiftiParameters(nii_files);
    globalCounts.converted_scans(iSubject, iSession, iScan) = 1;
    

end




