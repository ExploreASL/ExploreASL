function xASL_stat_GetMotionStatistics(x, GetStats)
%xASL_stat_GetMotionStatistics Summarize motion values
% INPUT
% x = from ExploreASL
% masks   = physically loaded masks for each subject
% ASL     = the data to be analyzed. Could be ASL, or e.g. SD or SNR masks
%
%
% By HJMM Mutsaerts, ExploreASL 2016
%


% NB: if you have registration mean control - label only,
% position will be [0 DIFF]
% so position==motion
% and mean(motion)==median(motion)==MAD(motion)
% So all numbers will be equal

if ~exist('GetStats','var')
    GetStats    = 0;
end


    
%% -----------------------------------------------------------------------------------------------    
%% Print motion & exclusion overview

FList   = {};
fprintf('Searching for motion files...   ');
for iL=1:x.nSubjects
    xASL_TrackProgress(iL,x.nSubjects);
    TempList   = xASL_adm_GetFileList(x.D.MotionDir, ['^motion_correction_NDV_' x.SUBJECTS{iL} '_ASL_\d*\.mat$'],'FPList',[0 Inf]);

    if ~isempty(TempList)
        FList(end+1:end+length(TempList),1)     = TempList;
    end
end
fprintf('\n');
    
if  length(FList)<=0

    fprintf('%s\n','No motion stats available');
else
    
    
    for iL  = 1:length(FList)
        tempMot                 = load(FList{iL});
        mean_motion(iL,1)       = tempMot.mean_NDV{2};
        PercExclusion(iL,1)     = tempMot.PercExcl;
    end

    fig = figure('Visible','off');
    plot(mean_motion,PercExclusion,'b.')
    % axis([0 0.8 0 30])
    xlabel('Mean motion (Diff Net Displacement Vector, mm)');
    ylabel('Percentage excluded pairs (%)');
    title('Threshold-free motion spike exclusion to optimize statistical power');

    jpgfile                     = fullfile( x.D.MotionDir,'Overview_motion_pair-exclusion.jpg');
    fprintf('Saving motion plot to %s\n',jpgfile);
    saveas(fig,jpgfile,'jpg');
    close all;

    
    
    
    %% -----------------------------------------------------------------------------------------------
    %% Load motion files

    x.S.NamesROI                              = {'median position' 'mean position' 'MAD position' 'median motion' 'mean motion' 'MAD motion' 'Percent Exclusion' 'Pvalue Threshold Free Spike Exclusion'};    
    x.S.output_ID                             = 'motion';
    x.S.unit                                  = 'mm';

    fprintf('%s\n',['Printing csv-files with ' x.S.output_ID ' statistics...  ']);

    for iSubject=1:x.nSubjects
        xASL_TrackProgress(iSubject,x.nSubjects);
        x.P.SubjectID      = x.SUBJECTS{iSubject};
        for iSession=1:x.nSessions
            x.P.SessionID  = x.SESSIONS{iSession};

            iSubjSess   = ((iSubject-1)*x.nSessions)+iSession;

            loadMot                         = fullfile( x.D.MotionDir, ['motion_correction_NDV_' x.SUBJECTS{iSubject} '_' x.SESSIONS{iSession} '.mat'] );
            
            % Save for stats
            MeanMotion{iSubjSess,1}      = x.P.SubjectID;
            MeanMotion{iSubjSess,2}      = x.P.SessionID;            
            
            
            if ~exist(loadMot,'file')
                fprintf('%s\n',['No motion data for ' x.SUBJECTS{iSubject} '_' x.SESSIONS{iSession}]);
                MeanMotion{iSubjSess,3}          = NaN;  % missing value
                
                
                
            else
                tempMot                         = load( loadMot );

                % position
                motionMat(1,iSession,iSubject)   = tempMot.median_NDV{1};      % median position
                motionMat(2,iSession,iSubject)   = tempMot.mean_NDV{1};        % mean position
                motionMat(3,iSession,iSubject)   = tempMot.MAD_NDV{1};         % MAD position

                % motion
                motionMat(4,iSession,iSubject)   = tempMot.median_NDV{2};      % median motion
                motionMat(5,iSession,iSubject)   = tempMot.mean_NDV{2};        % mean motion
                motionMat(6,iSession,iSubject)   = tempMot.MAD_NDV{2};         % MAD motion

    %             % motion spikes exclusion
                motionMat(7,iSession,iSubject)   = tempMot.PercExcl;           % MAD
                motionMat(8,iSession,iSubject)   = tempMot.MinimumtValue;      % MAD

                % Save for stats

                MeanMotion{iSubjSess,3}          = tempMot.mean_NDV{2};  % mean motion
                % Mean motion is used for stats, since this includes high
                % motion spikes rather than excluding them (median)
            end
        end
    end
    
    fprintf('\n');

    
    
    
    
    %% -----------------------------------------------------------------------------------------------
    %% Save mat-file for statistics later in pipeline
    matName         = fullfile( x.D.ROOT, 'MeanMotion.mat');
    save( matName , 'MeanMotion');

    if  GetStats

        % Restructure DATA per ROI/measurement
        fprintf('%s\n','Restructure DATA per ROI/measurement');
        for iSubject=1:x.nSubjects
            xASL_TrackProgress(iSubject,x.nSubjects);
            for iSession=1:x.nSessions
                % ID (which name, group etc), all for identification
                SUBJECT_SESSION                 = (iSubject-1)* x.nSessions +iSession;
                x.S.SUBJECTID{SUBJECT_SESSION,1}  = x.SUBJECTS{iSubject};

                % MOTION
                for iMeas = 1:length(x.S.NamesROI)
                    if size(motionMat,1)<iMeas || size(motionMat,2)<iSession || size(motionMat,3)<iSubject
                        x.S.DAT(SUBJECT_SESSION,iMeas) = NaN;
                    else
                        x.S.DAT(SUBJECT_SESSION,iMeas) = motionMat(iMeas,iSession,iSubject);
                    end
                end
            end
        end
        
        fprintf('\n');

        % Initiation parameters for xASL_wrp_PermuteSets1

        x.S.StatsDir          = x.S.StatsDir;
        x.S.function2call     = @xASL_stat_PrintBasicStats;
        x.S.KISS              = 1; % keeps it simple    

        xASL_wrp_PermuteSets1( x );

    %     x.S.HistogramDir      = fullfile( x.HistogramDir,S.output_ID);
    %     xASL_adm_CreateDir(x.S.HistogramDir);
    %     x.S.oriDIR            = x.S.HistogramDir;
    %     x.S.function2call     = @xASL_stat_CreateHistograms;
    %     fprintf('%s\n',['Creating ' x.S.output_ID ' histograms']);
    end




end
