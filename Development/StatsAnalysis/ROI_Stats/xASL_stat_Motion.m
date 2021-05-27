function [x] = xASL_stat_Motion(x)
%xASL_stat_Motion Permute over groups to compare motion stats

%% -----------------------------------------------------------------------------------------------
%% Run basic statistical comparisons

% Restructure DATA per ROI/measurement
fprintf('%s\n','Restructure DATA per ROI/measurement');
for iSubject=1:x.nSubjects
    for iSession=1:x.nSessions

        % Track progress
        iSubjectSession = (iSubject-1)* x.nSessions +iSession;
        x.S.SUBJECTID{iSubjectSession,1} = x.SUBJECTS{iSubject};
        xASL_TrackProgress(iSubjectSession,x.nSubjectsSessions);

        % Collect motion
        for iMeas = 1:length(x.S.NamesROI)
            if size(motionMat,1)<iMeas || size(motionMat,2)<iSession || size(motionMat,3)<iSubject
                x.S.DAT(iSubjectSession,iMeas) = NaN;
            else
                x.S.DAT(iSubjectSession,iMeas) = motionMat(iMeas,iSession,iSubject);
            end
        end
    end
end

fprintf('\n');

% Initiation parameters for xASL_wrp_PermuteSets1
x.S.StatsDir = x.S.StatsDir;
x.S.function2call = @xASL_stat_PrintBasicStats;
x.S.KISS = 1; % keep it simple    

xASL_wrp_PermuteSets1(x);

%     x.D.HistogramDir      = fullfile( x.D.HistogramDir,S.output_ID);
%     xASL_adm_CreateDir(x.D.HistogramDir);
%     x.S.oriDIR            = x.D.HistogramDir;
%     x.S.function2call     = @xASL_stat_CreateHistograms;
%     fprintf('%s\n',['Creating ' x.S.output_ID ' histograms']);
end


end

