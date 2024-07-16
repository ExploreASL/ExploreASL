function xASL_stat_GetDICOMStatistics(x, ScanType, HasSessions, bOverwrite)
%xASL_stat_GetDICOMStatistics Collect JSON sidecar data in TSV file to check
%
% FORMAT: xASL_stat_GetDICOMStatistics(x, ScanType, HasSessions)
% 
% INPUT:
%   x           - struct containing pipeline environment parameters (REQUIRED)
%   ScanType    - type of NIfTI (e.g. T1, FLAIR, ASL, M0) (REQUIRED)
%   HasSessions - true for e.g., multiple ASL sessions (OPTIONAL, DEFAULT=false)
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
%              1. Loading JSON sidecars
%              2. Write TSV file
%
%              Note that this function was revamped from its previous version to obtain stats from ASL4D_parms.mat files.
%              PM: when we fully move to BIDS derivatives, we can also rename the folder and this function DICOMparameters
% 
% ------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_stat_GetDICOMStatistics(x, 'ASL', true);
% __________________________________
% Copyright 2016-2024 ExploreASL


%% -----------------------------------------------------------------------------------------------
%% Admin
if nargin<4 || isempty(bOverwrite)
    bOverwrite = true;
end
if nargin<3 || isempty(HasSessions)
    HasSessions = false;
end

PathTSV = fullfile(x.D.DICOMparameterDir, ['QuantificationParameters_' ScanType '.tsv']);

% Print header
TSV = {'participant_id' 'session'};


%% 1. Loading JSON sidecars
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
xASL_tsvWrite(TSV, PathTSV, bOverwrite);


end