function xASL_stat_GetRegistrationStatistics(x)
%xASL_stat_GetRegistrationStatistics Collect the Tanimoto Coefficients from registration in TSV file to check
%
% FORMAT: xASL_stat_GetRegistrationStatistics(x)
% 
% INPUT:
%   x - struct containing pipeline environment parameters (REQUIRED)
%
% INPUT FILES:
%   The QC files for each subject
%   '/MyStudy/SUBJECTDIR/QC_collection_SUBJECTDIR.json'
%   Containing the calculated Tanimoto coefficients for registration between:
%   - ASL-T1w: TC_ASL2T1w_Perc
%
% OUTPUT FILES: 
%   Written to:
%   /MyStudy/Population/Stats/RegistrationTC.tsv
%                         
% ------------------------------------------------------------------------------------------------
% DESCRIPTION: Loads the data from the study given in the QC_collection*.json files. Goes through all subjects and 
%              sessions and prints the Tanimoto coefficients that define the quality of the registrations. Steps:
%
% 1. Load & extract parameters from individual parameter files
% 2. Write TSV file
%
% ------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_stat_GetRegistrationStatistics(x);
% __________________________________
% Copyright 2015-2020 ExploreASL


%% -----------------------------------------------------------------------------------------------
%% Admin
if nargin<1 || isempty(x)
    error('Missing the x-struct');
end

% Define a list of fields to write
jsonFields = {'TC_ASL2T1w_Perc'};
nFields = length(jsonFields);

% Create the file name for TSV file to save
PathTSV = fullfile(x.D.PopDir, 'Stats','RegistrationTC.tsv');

% Print header
TSV = {'participant_id' 'session'};
TSV(1,3:2+nFields) = jsonFields(1:nFields);

%% -----------------------------------------------------------------------------------------------
%% 1. Load & extract parameters from individual parameter files
fprintf('%s\n','Loading & saving individual parameter files...  ');

for iSubject=1:x.nSubjects
    for iSession=1:x.nSessions

        % Joint index for subject and session
        iSubjSess = (iSubject-1)*x.nSessions+iSession;
		
		% Track progress
        xASL_TrackProgress(iSubjSess,x.dataset.nSubjectsSessions);        
        
        % Write the subject and session name in the TSV table
        TSV{1+iSubjSess,1} = x.SUBJECTS{iSubject};
        TSV{1+iSubjSess,2} = x.SESSIONS{iSession};
        
		% Initialize with NaNs
        TSV(1+iSubjSess, 3:nFields+2) = repmat({NaN}, [1 nFields]);
        
        % Define path of the parameter file
		PathJSON = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, ['QC_collection_' x.SUBJECTS{iSubject} '.json']);

		if exist(PathJSON, 'file')
			% Load the file
            Parms = spm_jsonread(PathJSON);
            Parms = Parms.ASL;

			% print all fields for subject_session into the TSV array
			for iField=1:nFields
				if isfield(Parms,jsonFields{iField}) && ~isempty(Parms.(jsonFields{iField}))
					TSV{1+iSubjSess, 2+iField} = xASL_num2str(Parms.(jsonFields{iField}));
				end
			end
		end
    end
end

fprintf('\n');

%% -----------------------------------------------------------------------------------------------
%% 2. Write TSV file
xASL_tsvWrite(TSV, PathTSV, 1);

end