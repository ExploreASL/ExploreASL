function [x,nii_files, summary_line, globalCounts, ASLContext] = xASL_imp_DCM2NII_Subject_SortASLVolumes(x, globalCounts, scanpath, scan_name, nii_files, iSubject, iVisit, iSession, iScan)
%xASL_imp_DCM2NII_Subject_SortASLVolumes Sort ASL Volumes
%
% FORMAT: [x, nii_files, summary_line, globalCounts, ASLContext] = xASL_imp_DCM2NII_Subject_SortASLVolumes(x, globalCounts, scanpath, scan_name, nii_files, iSubject, iVisit, iSession, iScan)
% 
% INPUT:
%   x               - ExploreASL x struct (STRUCT, REQUIRED)
%   globalCounts    - Converted, skipped & missing scans (REQUIRED, STRUCT)
%   scanpath        - Scan path (CHAR ARRAY, PATH, REQUIRED)
%   scan_name       - Scan name (CHAR ARRAY, REQUIRED)
%   nii_files       - List of NIfTI files (CELL ARRAY, REQUIRED)
%   iSubject        - Subject ID (INTEGER, REQUIRED)
%   iVisit          - Visit ID (INTEGER, REQUIRED)
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
% DESCRIPTION: Sort ASL Volumes.
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
		if isfield(tmpJSON,'SeriesNumber')
            niiTable{iNii,5} = xASL_str2num(tmpJSON.SeriesNumber);
        else
            niiTable{iNii,5} = 0;
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
    
    if isempty(niiTable)
        warning('Empty table, sorting based on instance numbers will be skipped...');
    else
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
            [nii_files,ASLContext] = xASL_bids_MergeNifti(nii_files, 'ASL', niiTable);
        elseif  ~isempty(strfind(scan_name,'M0'))
            [nii_files,ASLContext] = xASL_bids_MergeNifti(nii_files, 'M0', niiTable);
        end
    end
    
    %% Checks for time encoded sequences (Hadamard etc.)
    bTimeEncoded = false;
    bTimeEncodedFME = false;
    
    % Determine if we have a Hadamard sequence based on the parameters of the studyPar.json
    [bTimeEncoded] = xASL_imp_DCM2NII_CheckIfTimeEncoded(x, bTimeEncoded, iSubject, iVisit, iSession);
    
    % Check if the current sequence is a FME (Fraunhofer Mevis) time encoded sequence
    [resultJSON, bTimeEncoded, ~] = xASL_imp_DCM2NII_CheckIfFME(nii_files, bTimeEncoded, bTimeEncodedFME);
    
    % Reorder TEs and PLDs accordingly for time encoded sequences
    xASL_imp_DCM2NII_ReorderTimeEncoded(nii_files, bTimeEncoded, resultJSON);
     
    
    %% 6. Extract relevant parameters from nifti header and append to summary file
    summary_line = xASL_imp_AppendNiftiParameters(nii_files);
    globalCounts.converted_scans(iSubject, iVisit, iSession, iScan) = 1;
    

end


%% Determine if we have a Hadamard sequence based on the parameters of the studyPar.json
function [bTimeEncoded] = xASL_imp_DCM2NII_CheckIfTimeEncoded(x, bTimeEncoded, iSubject, iVisit, iSession)

    if nargin<2 || isempty(bTimeEncoded)
        bTimeEncoded = false; % default
    end
    
    if isfield(x.dir, 'studyPar') && ~isempty(x.dir.studyPar)
        if xASL_exist(x.dir.studyPar, 'file')
            % Get the specific studyPar parameters
			studyParAll = xASL_io_ReadDataPar(x.dir.studyPar, true);
			structSubject = x.importOverview.(['subject_' num2str(iSubject,'%.3d')]);
			structVisit   = structSubject.(['visit_' num2str(iVisit,'%.3d')]);
			structRun     = structVisit.(['run_' num2str(iSession,'%.3d')]);
			studyParSpecificSubjVisitSess = xASL_imp_StudyParPriority(studyParAll, structSubject.name, structVisit.name, structRun.name(5:end));
			
			if isfield(studyParSpecificSubjVisitSess,'TimeEncodedMatrixSize') && ~isempty(studyParSpecificSubjVisitSess.TimeEncodedMatrixSize) || ... % Should be 4, 8 or 12
					isfield(studyParSpecificSubjVisitSess,'TimeEncodedMatrixType') % Natural or walsh
				bTimeEncoded = true;
			end
        end
    end

end



