function [bidsPar, TypeIs, pathOrig, pathDest] = xASL_bids_BIDS2Legacy_CompilePathsForCopying(bidsPar, TypeIs, ModalityIs, RunIs, iSubjSess, BIDS, TypeRunIndex, ModalityFields, pathLegacy_SubjectVisit)
%xASL_bids_BIDS2Legacy_CompilePathsForCopying Compile paths for BIDS to Legacy copying.
%
% FORMAT:     [bidsPar, TypeIs, pathOrig, pathDest] = xASL_bids_BIDS2Legacy_CompilePathsForCopying(bidsPar, TypeIs)
% 
% INPUT:      bidsPar        - BIDS par struct (STRUCT, REQUIRED)
%             TypeIs         - Type (REQUIRED)
%             ModalityIs     - Modality (REQUIRED)
%             RunIs          - Run (REQUIRED)
%             iSubjSess      - Subject session (INTEGER, REQUIRED)
%             BIDS           - BIDS struct (STRUCT, REQUIRED)
%             TypeRunIndex   - Type run index (INTEGER, REQUIRED)
%             ModalityFields - Modality fields (REQUIRED)
%             pathLegacy_SubjectVisit - Legacy path for subject visit (STRING, REQUIRED)
%   
% OUTPUT:     bidsPar  - BIDS par struct (STRUCT, REQUIRED)
%             TypeIs   - Type
%             pathOrig - Origin path (STRING)
%             pathDest - Destination path (STRING)
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Compile paths for BIDS to Legacy copying.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     [bidsPar, TypeIs, pathOrig, pathDest] = xASL_bids_BIDS2Legacy_CompilePathsForCopying(bidsPar, TypeIs);
% __________________________________
% Copyright 2015-2021 ExploreASL


    % ModalityIs = current modality (e.g. 'anat' 'perf')
    % TypeIs = current scantype, e.g. 'asl' 'm0' 't1w'
    % RunIs = current run (e.g. 1, 2, 3)
    % TypeRunIndex = index of current scantype & run inside the above created Reference Table

    % Define folder & filename
    ConfigIndex = cellfun(@(y) strcmp(TypeIs, y), bidsPar.BIDS2LegacyFolderConfiguration(2:end, 1)); % find scantype index
    ConfigIndex2 = cellfun(@(y) strcmp(ModalityIs, y), bidsPar.BIDS2LegacyFolderConfiguration(2:end, 2)); % find modality index
    ConfigIndex = find(ConfigIndex & ConfigIndex2); % Combine them & convert to index

    % xASL legacy subfolder name (empty means no subfolder)
    FolderIs = bidsPar.BIDS2LegacyFolderConfiguration{ConfigIndex+1,3};

    % xASL legacy filename
    FileIs = bidsPar.BIDS2LegacyFolderConfiguration{ConfigIndex+1,4};

    % Manage runs
    PrintRun = false;
    % if xASL legacy requires to specify the first run (e.g. T1.nii.gz vs ASL_1/ASL4D.nii.gz)
    if RunIs==1 && bidsPar.BIDS2LegacyFolderConfiguration{ConfigIndex+1, 6}==1
        % if we need to print first run
        PrintRun = true;
    elseif RunIs>1
        PrintRun = true;
    end

    if PrintRun==1
        % xASL legacy location of run specification (e.g. T1_2.nii.gz vs ASL_2/ASL4D.nii.gz for file vs folder location)
        if strcmp(bidsPar.BIDS2LegacyFolderConfiguration{ConfigIndex+1,5}, 'file')
            FileIs = [FileIs '_' xASL_num2str(RunIs)];
        elseif strcmp(bidsPar.BIDS2LegacyFolderConfiguration{ConfigIndex+1,5}, 'folder')
            FolderIs = [FolderIs '_' xASL_num2str(RunIs)];
        end
    end

    pathOrig = '';
    pathDest = '';

    pathOrig{1} = fullfile(BIDS.subjects(iSubjSess).path, ModalityIs, ModalityFields(TypeRunIndex).filename);
    [~, ~, Fext] = xASL_fileparts(ModalityFields(TypeRunIndex).filename);
    pathDest{1} = fullfile(pathLegacy_SubjectVisit, FolderIs, [FileIs Fext]);


end


