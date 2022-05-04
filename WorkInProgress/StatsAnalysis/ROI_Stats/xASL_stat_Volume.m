function [x] = xASL_stat_Volume(x)
%xASL_stat_Volume Permute over groups to compare volume stats

% Restructure DATA per ROI/measurement
x.S.DAT = vol;

% Restructure x.S.SetsID into single session only
if  isfield(x.S,'SetsID')
    if  numel(x.S.SetsID)>0
        TempID  = x.S.SetsID;
        x.S = rmfield(x.S,'SetsID');
        for iSubject=1:x.nSubjects
             x.S.SetsID(iSubject,:) = TempID( ((iSubject-1)*x.dataset.nSessions)+1,:);
        end
    end
end

% Initiation parameters for xASL_wrp_PermuteSets1

x.dataset.nSessions = 1;
x.S.KISS = 1; % keeps it simple
x.S.StatsDir = x.S.StatsDir;
x.S.function2call = @xASL_stat_PrintBasicStats;

xASL_wrp_PermuteSets1(x);


end

