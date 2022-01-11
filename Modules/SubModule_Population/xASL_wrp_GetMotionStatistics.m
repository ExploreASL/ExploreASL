function xASL_wrp_GetMotionStatistics(x)
%xASL_wrp_GetMotionStatistics Summarize motion values
%
% FORMAT: xASL_wrp_GetMotionStatistics(x)
% 
% INPUT:
%   x - struct containing pipeline environment parameters (REQUIRED)
% INPUT FILES:
% All files in the format:
%   '/MyStudy/Population/MotionASL/motion_correction_NDV_Sub-001_ASL_1.mat'
% Containing the following rigid-body motion estimation parameters:
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
% 1. Collect motion data
% 2. If no data, skip this function
% 3. Print motion vs exclusion overview
% 4. Add motion data to participants.tsv
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: xASL_wrp_GetMotionStatistics(x);
% __________________________________
% Copyright (c) 2016-2022 ExploreASL


%% -----------------------------------------------------------------------------------------------
%% 1) Collect motion data
x.S.NamesROI = {'median position' 'mean position' 'MAD position' 'median motion' 'mean motion' 'MAD motion' 'Percent Exclusion' 'tValue Threshold Free Spike Exclusion'};    
x.S.output_ID = 'motion';
x.S.unit = 'mm';

fprintf('%s\n',['Collecting motion metadata with ' x.S.output_ID ' statistics:  ']);

for iSubject=1:x.nSubjects
    for iSession=1:x.dataset.nSessions
        % Keeping track
        iSubjSess = ((iSubject-1)*x.dataset.nSessions)+iSession;
        xASL_TrackProgress(iSubjSess, x.dataset.nSubjectsSessions);

        PathMAT = fullfile(x.D.MotionDir, ['motion_correction_NDV_' x.SUBJECTS{iSubject} '_' x.SESSIONS{iSession} '.mat'] );

        % Define defaults
        MeanMotion{iSubjSess,1} = x.SUBJECTS{iSubject};
        MeanMotion{iSubjSess,2} = x.SESSIONS{iSession};            
        MeanMotion{iSubjSess,3} = NaN;
        
        motionMat(1:8,iSession,iSubject) = repmat(NaN, [1 8]);
        PercExclusion(iSubjSess,1) = NaN;
        MeanMotionNum(iSubjSess,1) = NaN;
        
        if exist(PathMAT,'file')
            tempMot = load(PathMAT, '-mat');

            % position
            motionMat(1,iSession,iSubject) = tempMot.median_NDV{1};      % median position
            motionMat(2,iSession,iSubject) = tempMot.mean_NDV{1};        % mean position
            motionMat(3,iSession,iSubject) = tempMot.MAD_NDV{1};         % MAD position

            % motion
            motionMat(4,iSession,iSubject) = tempMot.median_NDV{2};      % median motion
            motionMat(5,iSession,iSubject) = tempMot.mean_NDV{2};        % mean motion
            motionMat(6,iSession,iSubject) = tempMot.MAD_NDV{2};         % MAD motion

            % motion spikes exclusion
            motionMat(7,iSession,iSubject) = tempMot.PercExcl;           % MAD
            motionMat(8,iSession,iSubject) = tempMot.MinimumtValue;      % MAD

            % Save for stats
            MeanMotion{iSubjSess,3} = tempMot.mean_NDV{2};
            MeanMotionNum(iSubjSess,1) = tempMot.mean_NDV{2};
            PercExclusion(iSubjSess,1) = tempMot.PercExcl;
            % Mean motion is used for stats, since this includes high
            % motion spikes rather than excluding them (median)
        end
    end
end

fprintf('\n');

%% 2) If no data, skip this function
HasData = sum(isfinite(MeanMotionNum))>0;
if ~HasData
    fprintf('%s\n','No motion stats available, skipping');
    return;
end


%% 3) Print motion vs exclusion overview
if usejava('jvm')
    fig = figure('Visible','off');
    plot(MeanMotionNum, PercExclusion, 'b.')
    % axis([0 0.8 0 30])
    xlabel('Mean motion (Diff Net Displacement Vector, mm)');
    ylabel('Percentage excluded pairs (%)');
    title('Threshold-free motion spike exclusion to optimize statistical power');

    PathJPG = fullfile( x.D.MotionDir,'Overview_motion_pair-exclusion.jpg');
    fprintf('Saving motion plot to %s\n',PathJPG);
    saveas(fig,PathJPG,'jpg');
    close all;
else
    fprintf('Skipping motion vs exclusion overview, missing JVM\n');
end

%% -----------------------------------------------------------------------------------------------
%% 4) Add motion data to participants.tsv
xASL_bids_Add2ParticipantsTSV(MeanMotion, 'MeanMotion', x);




end
