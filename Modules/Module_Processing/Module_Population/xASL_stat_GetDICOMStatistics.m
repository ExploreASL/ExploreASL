function xASL_stat_GetDICOMStatistics(x, ScanType, HasSessions, bOverwrite)
%xASL_stat_GetDICOMStatistics Collect DICOM metadata in CSV file to check
%
% FORMAT: xASL_stat_GetDICOMStatistics(x, ScanType, HasSessions)
% 
% INPUT:
%   x           - struct containing pipeline environment parameters (REQUIRED)
%   ScanType    - type of data (e.g. T1, FLAIR, ASL) (REQUIRED)
%   HasSessions - true for multiple ASL sessions (OPTIONAL, DEFAULT=false)
%   bOverwrite  - boolean to overwrite existing TSV summary file (OPTIONAL, DEFAULT=true)
%
% OUTPUT: n/a
% ------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions prints DICOM metadata (e.g. parameters used
%              for quantification) and collects them in a single tsv (per BIDS).
%              Summarizes this for the total population.
%              Can be useful to detect software upgrades, where only slight
%              parameter changes can hint on quantification changes.
%              This function carries out the following steps:
%
%              1. Load & save individual parameter files
%              2. Print summary
%              3. Write TSV file
% ------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_stat_GetDICOMStatistics(x, 'ASL', true);
% __________________________________
% Copyright 2016-2020 ExploreASL


%% -----------------------------------------------------------------------------------------------
%% Admin
if nargin<4 || isempty(bOverwrite)
    bOverwrite = true;
end
if nargin<3 || isempty(HasSessions)
    HasSessions = false;
end

% matFields = {'RepetitionTime' 'EchoTime' 'NumberOfTemporalPositions' 'MRScaleSlope' 'RescaleSlopeOriginal' 'RescaleIntercept' 'Scanner' 'SliceReadoutTime'};
% nFields = length(matFields);
PathTSV = fullfile(x.D.DICOMparameterDir, ['QuantificationParameters_' ScanType '.tsv']);

% Print header
TSV = {'participant_id' 'session'};
% TSV(1,3:2+nFields) = matFields(1:nFields);

% Initialize empty table
% numElements = size(TSV,2);
% numSubjectsSessions = x.dataset.nSubjects*x.dataset.nSessions;
% for iSubjSess=1:numSubjectsSessions
%     TSV(1+iSubjSess,:) = repmat({' '},1, numElements);
% end

%% 1. Load & save individual parameter files
fprintf('%s\n', ['Loading ' ScanType '.json sidecars...  ']);

for iSubject=1:x.dataset.nSubjects
    for iSession=1:x.dataset.nSessions

        % Track progress
        iSubjSess = (iSubject-1)*x.dataset.nSessions+iSession;
        xASL_TrackProgress(iSubjSess, x.dataset.nSubjects * x.dataset.nSessions);
        
        % Default parameters
        TSV{1+iSubjSess,1} = x.SUBJECTS{iSubject};
        TSV{1+iSubjSess,2} = x.SESSIONS{iSession};
        
        % TSV(1+iSubjSess, 3:nFields+2) = repmat({NaN}, [1 nFields]);
        % x.S.par(iSubjSess,1:nFields) = NaN;
        
        % Define paths
		if HasSessions
            PathJSON = fullfile(x.dir.xASLDerivatives, x.SUBJECTS{iSubject}, x.SESSIONS{iSession}, [ScanType '.json']);
        else
            PathJSON = fullfile(x.dir.xASLDerivatives, x.SUBJECTS{iSubject}, [ScanType '.json']); 
		end
		if exist(PathJSON, 'file')
			Parms = xASL_io_ReadJson(PathJSON);
		else
			continue; % skip this iSujectSession as its missing
		end
        
        % print all fields for subject_session & load them into x.S.par variabele for later calculations
        matFields = fieldnames(Parms);
        nFields = length(matFields);
        for iField=1:nFields
            if ~isempty(Parms.(matFields{iField}))
                % Check if field already exists as column
                iColumn = find(strcmp(TSV(1,:), matFields{iField}));
                if isempty(iColumn)
                    TSV{1,end+1} = matFields{iField}; % create new column for this parameter
                    TSV{1+iSubjSess,end} = xASL_num2str(Parms.(matFields{iField}), [], [], [], true); % add the value to the new column, and run unique
                elseif numel(iColumn)>1
                    error('When making the DICOMparameters.TSV, columns with the same name');
                else
                    TSV{1+iSubjSess,iColumn} = xASL_num2str(Parms.(matFields{iField}), [], [], [], true); % add the value to the existing column, and run unique
                end
            end
        end
    end
end
fprintf('\n');


%% 2. Write TSV file

% % Ensure all elements of the TSV cell array have the correct format
% for iRow=1:size(TSV,1)
%     for iColumn=1:size(TSV,2)
%         if size(TSV{iRow,iColumn},1)>1
%             % We skip lists, as this was introduced with BIDS
%             % This function aims to do DICOM QC, BIDS lists are not necessary
%             % here, as we can now find them in the *asl.json sidecars
%             TSV{iRow,iColumn} = 'ListSkipped';
%         end
%         % Remove empty elements (Empty cells in the TSV can lead to reading errors)
%         if isempty(TSV{iRow,iColumn})
%             TSV{iRow,iColumn} = '_';
%         end
%         % Convert numeric NaNs to text n/a's (NaN can be exported as empty cells into the TSV, which can lead to reading errors. In addition n/a is BIDS compliant.)
%         if isnumeric(TSV{iRow,iColumn})
%             if isnan(TSV{iRow,iColumn})
%                 TSV{iRow,iColumn} = 'n/a';
%             else
%                 TSV{iRow,iColumn} = xASL_num2str(TSV{iRow,iColumn});
%             end
%         end
%     end
% end

% Write TSV
xASL_tsvWrite(TSV, PathTSV, bOverwrite);


end