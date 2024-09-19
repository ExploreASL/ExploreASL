% Copyright 2015-2024 ExploreASL (Works In Progress code)
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

function [x] = xASL_stat_Volume(x)
%xASL_stat_Volume Permute over groups to compare volume stats

% Restructure DATA per ROI/measurement
x.S.DAT = vol;

% Restructure x.S.SetsID into single session only
if  isfield(x.S,'SetsID')
    if  numel(x.S.SetsID)>0
        TempID  = x.S.SetsID;
        x.S = rmfield(x.S,'SetsID');
        for iSubject=1:x.dataset.nSubjects
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

