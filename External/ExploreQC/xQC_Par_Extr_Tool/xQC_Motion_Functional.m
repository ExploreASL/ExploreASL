function [ Motion ] = xQC_Motion_Functional(PathToMotion, SubjID)

%PathToMotion is the path to the directory storing the mat files output of
%the motion correfction step



%% func motion
PathToSubjMotion = fullfile(PathToMotion,['motion_correction_NDV_' SubjID '_func.mat']);
if  exist(PathToSubjMotion,'file')
    MotionPars  = load(PathToSubjMotion);
    NDV = MotionPars.NDV;
    % motion stats
    for i = 1:length(NDV)
        Motion.MotionMean_mm(1,i) = xASL_stat_MeanNan(NDV{i});
        Motion.MotionSD_mm(1,i) = xASL_stat_StdNan(NDV{i});
        Motion.MotionMax_mm(1,i) = max(NDV{i});
    end
     % motion outliers
    Motion.MotionExcl_Perc     = xASL_round(MotionPars.PercExcl,3);
else
    Motion.MotionMean_mm       = [NaN, NaN];
    Motion.MotionExcl_Perc     = NaN;
    Motion.MotionMax_mm        = [NaN , NaN];
    Motion.MotionSD_mm         = [NaN , NaN];
end
