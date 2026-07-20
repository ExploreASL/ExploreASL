function xASL_stat_GetDCEstatistics(x)
%xASL_stat_GetDCEstatistics Summarize DCE QC values
%
% FORMAT: xASL_stat_GetDCEstatistics(x)
% 
% INPUT:
%   x - struct containing pipeline environment parameters (REQUIRED)
% INPUT FILES:
% All x files in the format:
%   '/MyStudy/Subject_visit/x.mat'
%
% Containing the following :
%   -   median, mean, mean absolute deviation (MAD) of position and motion,
%       where position is the actual position of the head and motion the
%       between-volume difference of position. Units are the net displacement vector
%       (NDV), which is the RMS of XYZ translations and XYZ rotations.
%   -   Percent exclusion - percentage of the control-label volumes that
%       was excluded because its motion was too high (spikes)
%   -   Pvalue Threshold Free Spike Exclusion - threshold used for
%       exclusion of motion spikes
%
% OUTPUT: n/a
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions collects motion stats, with the following steps:
%
% 1. Collect QC metadata
% 2. Add QC metadata to participants.tsv
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_stat_GetMotionStatistics(x);
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________



%% -----------------------------------------------------------------------------------------------
%% 1. Collect data

fprintf('%s\n','Collecting DCE QC metadata:  ');


x.SESSIONS = {'DCE_1'};

for iSubject=1:x.dataset.nSubjects
    for iSession=1:x.dataset.nSessions
        % Keeping track
        iSubjSess = ((iSubject-1)*x.dataset.nSessions)+iSession;
        xASL_TrackProgress(iSubjSess, x.dataset.nSubjects * x.dataset.nSessions);

        PathMAT = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, 'x.mat');
        loadX = load(PathMAT, '-mat');

        if isfield(loadX.x, 'Output') && isfield(loadX.x.Output, 'DCE')
            fieldsAre = fields(loadX.x.Output.DCE);

            for iField = 1:length(fieldsAre)
                DCE_QC.(fieldsAre{iField})(iSubjSess,1)  = loadX.x.Output.DCE.(fieldsAre{iField});
            end
        else
            % Add a NaN
            fieldsAre = fields(DCE_QC);
            for iField = 1:length(fieldsAre)
                DCE_QC.(fieldsAre{iField})(iSubjSess,1)  = NaN;
            end
        end

    end
end


%% -----------------------------------------------------------------------------------------------
%% 2) Add QC data to participants.tsv
fieldsAre = fields(DCE_QC);
for iField = 1:length(fieldsAre)
    ParticipantsColumn = x.SUBJECTS(:);
    ParticipantsColumn(:,2) = num2cell(DCE_QC.(fieldsAre{iField}));
    xASL_bids_Add2ParticipantsTSV(ParticipantsColumn, fieldsAre{iField}, x);
end




end
