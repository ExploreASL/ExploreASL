function PathFig = xASL_stat_GetAverageIntensity(x, subject)
%xASL_stat_GetAverageIntensity Summarize the acquisition time of scans
%
% FORMAT: [x] = xASL_stat_GetAverageIntensity(x)
% 
% INPUT:
%   x   - struct containing pipeline environment parameters (REQUIRED)
%
% INPUTFILES: can be any combination of *parms.mat (legacy) and *.json
%
% OUTPUT: /MyStudy/Population/Stats/AverageIntensity.jpg 
%                                                         
% 
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This functions collects the average intensity of the ASL image at any specific MultiPLD time
%
% 1. Collect PLD times
% 2. Calculate average intensity at specific PLD
% 3. Create a graph of intensity over time
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: x = xASL_stat_GetAverageIntensity(x);
% __________________________________
% Copyright 2016-2023 ExploreASL



    %% -----------------------------------------------------------------------------------------------
    %% 1) Collect times
    fprintf('%s\n','Loading AcquisitionTimes:  ');

    for iSession=1:x.dataset.nSessions
        
        % Initialize defaults
        AcquisitionTime{iSession,1} = subject;
        AcquisitionTime{iSession,2} = x.SESSIONS{iSession};
        AcquisitionTime{iSession,3} = NaN;
        MissingData(iSession,1) = 1;

        % Define paths
        PathMAT = fullfile(x.D.ROOT, subject, x.SESSIONS{iSession}, 'ASL4D_parms.mat'); % legacy
        Parms = xASL_adm_LoadParms(PathMAT,[],0);
        
        %PathJSON = fullfile(x.D.ROOT, x.SUBJECTS{iSubject}, x.SESSIONS{iSession}, 'ASL4D.json');
        %Parms = spm_jsonread(PathJSON);
                
        if isfield(Parms,'AcquisitionTime')
            AcquisitionTime{iSession,3} = min(Parms.AcquisitionTime);
            AcquisitionTimeN(iSession,1) = min(Parms.AcquisitionTime); % for histogram below
            MissingData(iSession,1) = 0;
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