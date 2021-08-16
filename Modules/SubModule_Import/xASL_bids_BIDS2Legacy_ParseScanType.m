function xASL_bids_BIDS2Legacy_ParseScanType(modalityConfiguration, SubjectVisit, RunsUnique, RunsAre, bOverwrite)
%xASL_bids_BIDS2Legacy_ParseScanType Parse scan type during BIDS to Legacy conversion.
%
% FORMAT:     xASL_bids_BIDS2Legacy_ParseScanType(modalityConfiguration, SubjectVisit, RunsUnique, RunsAre, bOverwrite)
% 
% INPUT:      modalityConfiguration    - Modality Configuration (REQUIRED)
%             SubjectVisit             - Subject visits (REQUIRED)
%             RunsUnique               - Runs uniqure (REQUIRED)
%             RunsAre                  - Runs are (REQUIRED)
%             bOverwrite               - Overwrite (BOOLEAN, REQUIRED)
%   
% OUTPUT:     n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Parse scan type during BIDS to Legacy conversion.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_bids_BIDS2Legacy_ParseScanType(modalityConfiguration, SubjectVisit, RunsUnique, RunsAre, bOverwrite);
% __________________________________
% Copyright 2015-2021 ExploreASL


    for iType=1:size(modalityConfiguration, 1) % iterate scantypes in this Subject/Visit/Modality
        TypeIs = modalityConfiguration{iType,1};
        TypeIndices = cellfun(@(y) strcmp(y, TypeIs), Reference(2:end, 2)); % this are the indices for this ScanType

        if ~isempty(TypeIndices)
            for iRun=1:length(RunsUnique) % iterate runs in this Subject/Visit/Modality
                RunIs = RunsUnique(iRun);
                RunIndices = RunsAre==RunsUnique(iRun);
                TypeRunIndex = find(RunIndices & TypeIndices);
                if length(TypeRunIndex)>1
                    warning(['Multiple NIfTIs found for ' SubjectVisit '_run-' xASL_num2str(RunIs) '_' TypeIs ', using first only']);
                    TypeRunIndex = TypeRunIndex(1);
                end

                %% 4.2. Compile paths for copying
                if length(TypeRunIndex)==1 % if this scantype-run combination exists
                    [bidsPar, TypeIs, pathOrig, pathDest] = xASL_bids_BIDS2Legacy_CompilePathsForCopying(bidsPar, TypeIs);

                    %% 4.3. Manage sidecars to copy
                    % Sidecars definitions are loaded by xASL_bids_Config at function start

                    % Assuming that each .nii has a .json sidecar, do the same for .json (and for
                    % other sidecars only if they exist per bidsPar.sidecarRequired)
                    [bidsPar, pathOrig, pathDest, TypeIs] = xASL_bids_BIDS2Legacy_ManageSidecars(bidsPar, pathOrig, pathDest, TypeIs);


                    %% 4.4. Copy files
                    for iFile=1:length(pathOrig)
                        xASL_bids_BIDS2xASL_CopyFile(pathOrig{iFile}, pathDest{iFile}, bOverwrite);
                    end

                end
            end % iterate runs in this Subject/Visit/Modality
        end
    end


end


