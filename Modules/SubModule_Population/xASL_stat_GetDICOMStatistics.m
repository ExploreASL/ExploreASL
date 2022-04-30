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

% Mat fields that will be printed (some are in x.(field) and some in x.Q.(field))
matFields = {'RepetitionTime' 'EchoTime' 'NumberOfTemporalPositions' 'MRScaleSlope' 'RescaleSlopeOriginal' 'RescaleIntercept' 'Scanner' 'SliceReadoutTime'};
nFields = length(matFields);
PathTSV = fullfile(x.D.DICOMparameterDir, ['QuantificationParameters_' ScanType '.tsv']);

% Print header
TSV = {'participant_id' 'session'};
TSV(1,3:2+nFields) = matFields(1:nFields);

% Initialize empty table
numElements = size(TSV,2);
numSubjectsSessions = x.nSubjects*x.dataset.nSessions;
for iSubjSess=1:numSubjectsSessions
    TSV(1+iSubjSess,:) = repmat({' '},1, numElements);
end

%% 1) Load & save individual parameter files
fprintf('%s\n','Loading & saving individual parameter files...  ');

for iSubject=1:x.nSubjects
    for iSession=1:x.dataset.nSessions

        % Track progress
        iSubjSess = (iSubject-1)*x.dataset.nSessions+iSession;
        xASL_TrackProgress(iSubjSess,x.dataset.nSubjectsSessions);        
        
        % Default parameters
        TSV{1+iSubjSess,1} = x.SUBJECTS{iSubject};
        TSV{1+iSubjSess,2} = x.SESSIONS{iSession};
        
        TSV(1+iSubjSess, 3:nFields+2) = repmat({NaN}, [1 nFields]);
        x.S.par(iSubjSess,1:nFields) = NaN;
        
        % Define paths
		if HasSessions
            PathMAT = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, x.SESSIONS{iSession}, [ScanType '_Parms.mat']);% Legacy
            PathJSON = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, x.SESSIONS{iSession}, [ScanType '.json']);
        else
            PathMAT = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, [ScanType '_Parms.mat']); % Legacy
            PathJSON = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, [ScanType '.json']); 
		end
		if exist(PathMAT,'file') || exist(PathJSON,'file')
			Parms = xASL_adm_LoadParms(PathMAT,[],0);
			% Parms = spm_jsonread(PathJSON);
		else
			continue; % skip this iSujectSession as its missing
		end
        
        % print all fields for subject_session & load them into x.S.par variabele for later calculations
        for iField=1:nFields
            % Parms.(field)
            if isfield(Parms,matFields{iField}) && ~isempty(Parms.(matFields{iField}))
                % Get field from parameters
                TSV{1+iSubjSess, 2+iField} = xASL_num2str(Parms.(matFields{iField}));
                if isnumeric(Parms.(matFields{iField}))
                    TempData = Parms.(matFields{iField});
                    x.S.par(iSubjSess,iField) = min(TempData);
                else
                    x.S.par(iSubjSess,iField) = NaN; % fill empties with NaNs, not zeros!
                end
            end
            % Parms.Q.(field)
            if isfield(Parms,'Q')
                if isfield(Parms.Q,matFields{iField}) && ~isempty(Parms.Q.(matFields{iField}))
                    % Get field from parameters.Q
                    TSV{1+iSubjSess, 2+iField} = xASL_num2str(Parms.Q.(matFields{iField}));
                    if isnumeric(Parms.Q.(matFields{iField}))
                        TempData = Parms.Q.(matFields{iField});
                        x.S.par(iSubjSess,iField) = min(TempData);
                    else
                        x.S.par(iSubjSess,iField) = NaN; % fill empties with NaNs, not zeros!
                    end
                end
            end
        end
    end
end
fprintf('\n');


%% 2) Print summary
if ~isfield(x.S,'par') || isempty(x.S.par) || sum(sum(isfinite(x.S.par)))==0
    fprintf('%s\n',['Checking DICOM header values skipped for ' ScanType ' because no images existed']);
    return;
end

SummaryStats = {'mean' 'SD' '"CV (coeff var, %)"' 'up_threshold (mean+3SD)' 'lo_threshold (mean-3SD)'};

fprintf('%s\n',['Printing summary of ' ScanType ' DICOM values']);
TSV(end+1, :) = repmat({'_'},1, numElements);

for iField=1:size(x.S.par,2)
    temp_data = x.S.par(:,iField);
    statField{1}(iField) = xASL_stat_MeanNan( temp_data );
    statField{2}(iField) = xASL_stat_StdNan( temp_data);
    statField{3}(iField) = 100 * statField{2}(iField) / statField{1}(iField) ;
    statField{4}(iField) = statField{1}(iField) + (3 * statField{2}(iField));
    statField{5}(iField) = statField{1}(iField) - (3 * statField{2}(iField));
end

for iStat=1:length(SummaryStats)
    TSV{end+1, 1} = SummaryStats{iStat};
    for iField=1:size(x.S.par,2)
        TSV{end, 2+iField} = xASL_num2str(statField{iStat}(iField));
    end
end


%% 3) Write TSV file

% Ensure all elements of the TSV cell array have the correct format
for iRow=1:size(TSV,1)
    for iColumn=1:size(TSV,2)
        if size(TSV{iRow,iColumn},1)>1
            % We skip lists, as this was introduced with BIDS
            % This function aims to do DICOM QC, BIDS lists are not necessary
            % here, as we can now find them in the *asl.json sidecars
            TSV{iRow,iColumn} = 'ListSkipped';
        end
        % Remove empty elements (Empty cells in the TSV can lead to reading errors)
        if isempty(TSV{iRow,iColumn})
            TSV{iRow,iColumn} = '_';
        end
        % Convert numeric NaNs to text n/a's (NaN can be exported as empty cells into the TSV, which can lead to reading errors. In addition n/a is BIDS compliant.)
        if isnumeric(TSV{iRow,iColumn})
            if isnan(TSV{iRow,iColumn})
                TSV{iRow,iColumn} = 'n/a';
            else
                TSV{iRow,iColumn} = xASL_num2str(TSV{iRow,iColumn});
            end
        end
    end
end

% Write TSV
xASL_tsvWrite(TSV, PathTSV, bOverwrite);


end