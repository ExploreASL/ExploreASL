function xASL_bids_BIDS2Legacy_ParseModality(BIDS, bidsPar, SubjectVisit, iSubjSess, ModalitiesUnique, nModalities, bOverwrite, pathLegacy_SubjectVisit)
%xASL_bids_BIDS2Legacy_ParseModality Parse modality for BIDS to Legacy conversion.
%
% FORMAT:     xASL_bids_BIDS2Legacy_ParseModality(BIDS, bidsPar, SubjectVisit, iSubjSess, ModalitiesUnique, nModalities, bOverwrite)
% 
% INPUT:      BIDS             - BIDS struct (STRUCT, REQUIRED)
%             bidsPar          - BIDS par struct (STRUCT, REQUIRED)
%             SubjectVisit     - Subject visits (REQUIRED)
%             iSubjSess        - Current subject session (INTEGER, REQUIRED)
%             ModalitiesUnique - Unique modalities (REQUIRED)
%             nModalities      - Number of modalities (INTEGER, REQUIRED)
%             bOverwrite       - Overwrite (BOOLEAN, REQUIRED)
%             pathLegacy_SubjectVisit - Legacy path for subject visit (STRING, REQUIRED)
%   
% OUTPUT:     n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Parse modality for BIDS to Legacy conversion.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:     xASL_bids_BIDS2Legacy_ParseModality(BIDS, bidsPar, SubjectVisit, iSubjSess, ModalitiesUnique, nModalities, bOverwrite);
% __________________________________
% Copyright 2015-2021 ExploreASL

    for iModality=1:nModalities % iterate modalities in this Subject/Visit
        ModalityIs = ModalitiesUnique{iModality};
        if isfield(BIDS.subjects(iSubjSess), ModalityIs) && ~isempty(BIDS.subjects(iSubjSess).(ModalityIs))
            ModalityFields = BIDS.subjects(iSubjSess).(ModalityIs);
            nScans = length(ModalityFields);

            % Parse fields of this modality for combinations ScanType & Run, in a reference table
            Reference = {'index' 'ScanType' 'run'};
            for iScan=1:nScans % iterate NIfTIs in this Subject/Visit/Modality
                Reference{iScan+1, 1} = iScan; % index
                Reference{iScan+1, 2} = ModalityFields(iScan).type; % ScanType
                Reference{iScan+1, 3} = xASL_str2num(ModalityFields(iScan).run); % run
                if isempty(Reference{iScan+1, 3}) || isnan(Reference{iScan+1, 3})
                    Reference{iScan+1, 3} = 1; % default to 1st run
                end
            end
            Reference(2:end,:) = sortrows(Reference(2:end,:), [2, 3]); % first sort for ScanType then run

            RunsAre = cellfun(@(y) y, Reference(2:end, 3));
            RunsUnique = unique(RunsAre);


            %% 4.1. Parse scantype
            modalityIndices = find(strcmp(bidsPar.BIDS2LegacyFolderConfiguration(2:end,2), ModalityIs));
            modalityConfiguration = bidsPar.BIDS2LegacyFolderConfiguration(1+modalityIndices, :);
            xASL_bids_BIDS2Legacy_ParseScanType(modalityConfiguration, SubjectVisit, RunsUnique, RunsAre, bOverwrite, Reference, bidsPar, ModalityIs, iSubjSess, BIDS, ModalityFields, pathLegacy_SubjectVisit);


        end
    end


end


