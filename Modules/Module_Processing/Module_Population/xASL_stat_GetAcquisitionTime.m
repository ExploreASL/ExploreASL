function x = xASL_stat_GetAcquisitionTime(x)
%xASL_stat_GetAcquisitionTime Summarize the acquisition time of scans
%
% FORMAT: [x] = xASL_stat_GetAcquisitionTime(x)
% 
% INPUT:
%   x   - struct containing pipeline environment parameters (REQUIRED)
%
% INPUTFILES: can be any combination of *parms.mat (legacy) and *.json
%
% OUTPUT: /MyStudy/Population/Stats/AcquisitionTime.jpg - AcquisitionTime
%                                                         histogram for total study population
% 
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions collects the DICOM field AcquisitionTime from
% each json sidecar (& parms.mat for backward compatibility) and saves them
% in the participants.tsv. Additionally, it creates a AcquisitionTime histogram of the
% full study, which can be useful to check time of scanning -> can
% influence physiological CBF variability.
%
% 1. Collect times
% 2. Save times
% 3. Create time histogram
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: x = xASL_stat_GetAcquisitionTime(x);
% __________________________________
% Copyright 2016-2021 ExploreASL



    %% -----------------------------------------------------------------------------------------------
    %% 1) Collect times
    fprintf('%s\n','Loading AcquisitionTimes:  ');
    for iSubject=1:x.dataset.nSubjects
        for iSession=1:x.dataset.nSessions
            
            % Track progress
            iSubjSess = ((iSubject-1)*x.dataset.nSessions)+iSession;
            xASL_TrackProgress(iSubjSess,x.dataset.nSubjectsSessions);

            % Initialize defaults
            AcquisitionTime{iSubjSess,1} = x.SUBJECTS{iSubject};
            AcquisitionTime{iSubjSess,2} = x.SESSIONS{iSession};
            AcquisitionTime{iSubjSess,3} = NaN;
            MissingData(iSubjSess,1) = 1;

            % Define paths
            PathMAT = fullfile(x.dir.xASLDerivatives, x.SUBJECTS{iSubject}, x.SESSIONS{iSession}, 'ASL4D_parms.mat'); % legacy
            Parms = xASL_adm_LoadParms(PathMAT,[],0);
			
			%PathJSON = fullfile(x.dir.xASLDerivatives, x.SUBJECTS{iSubject}, x.SESSIONS{iSession}, 'ASL4D.json');
			%Parms = xASL_io_ReadJson(PathJSON);
		           
            if isfield(Parms,'AcquisitionTime')
                AcquisitionTime{iSubjSess,3} = min(Parms.AcquisitionTime);
                AcquisitionTimeN(iSubjSess,1) = min(Parms.AcquisitionTime); % for histogram below
                MissingData(iSubjSess,1) = 0;
            end
        end
    end
    fprintf('\n');
    
    %% -----------------------------------------------------------------------------------------------
    %% 2) Save times
    if sum(MissingData)<ceil(0.1*x.dataset.nSubjectsSessions) % allow 10 percent missing data
        xASL_bids_Add2ParticipantsTSV(AcquisitionTime, 'AcquisitionTime', x);
    end

    
    
    %% -----------------------------------------------------------------------------------------------    
    %% 3) Create time histogram
    if sum(MissingData)<ceil(0.1*x.dataset.nSubjectsSessions) % allow 10 percent missing data
        if usejava('jvm')
            % Only create figure when Java Virtual Machine is loaded
            PathFig = fullfile(x.S.StatsDir, 'AcquisitionTime.jpg');
            [N, X] = hist(round(AcquisitionTimeN./100)./100);
            N = N./sum(N);
            fig = figure('Visible','off');
            plot(X,N);
            title('Total population histogram of time of ASL scan');
            xlabel('Time (hours)');
            ylabel('Normalized frequency (%)');
            xASL_adm_CreateDir(x.S.StatsDir);
            print(gcf,'-djpeg','-r200', PathFig); % could be replaced by xASL_vis_Imwrite?
        else
            fprintf('Warning: skipping AcquisitiontTime.jpg, JVM missing\n');
        end
    end

end