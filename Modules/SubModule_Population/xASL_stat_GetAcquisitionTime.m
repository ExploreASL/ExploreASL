function x = xASL_stat_GetAcquisitionTime( x )
%xASL_stat_GetAcquisitionTime Summarizes the acquisition time of scans





%% -----------------------------------------------------------------------------------------------
%% Administration

SaveFile        = fullfile( x.D.ROOT, 'AcquisitionTime.mat');





%% -----------------------------------------------------------------------------------------------
%% Get times

if ~exist(SaveFile)
    fprintf('%s\n','Loading AcquisitionTimes...  ');
    for iF=1:x.nSubjects
        xASL_TrackProgress(iF,x.nSubjects);
        for iS=1:x.nSessions
            clear iSubjSess MatFile ParmsLoad
            
            iSubjSess   = ((iF-1)*x.nSessions)+iS;

            AcquisitionTime{iSubjSess,1}    = x.SUBJECTS{iF};
            AcquisitionTime{iSubjSess,2}    = x.SESSIONS{iS};
            
            MatFile                         = fullfile( x.D.ROOT, x.SUBJECTS{iF}, x.SESSIONS{iS}, 'ASL4D_parms.mat');
            
            if  exist(MatFile ,'file')
                ParmsLoad                       = load(MatFile);

                if      isfield(ParmsLoad.parms,'AcquisitionTime')
                        AcquisitionTime{iSubjSess,3}    = min(ParmsLoad.parms.AcquisitionTime);
                        AcquisitionTimeN(iSubjSess,1)   = min(ParmsLoad.parms.AcquisitionTime);
                        WasEmpty(iSubjSess,1)           = 0;
                else    WasEmpty(iSubjSess,1)           = 1;
                end
            else        WasEmpty(iSubjSess,1)           = 1;
            end
        end
    end
    
    fprintf('\n');
    
    if  sum(WasEmpty)< round(0.1*x.nSubjectsSessions) % allow 10 percent missing data
        save(SaveFile,'AcquisitionTime');
    end

    
    
    
    
%% -----------------------------------------------------------------------------------------------    
%% Create time histogram

    if  sum(WasEmpty)< round(0.1*x.nSubjectsSessions) % allow 10 percent missing data
        NameSave            = fullfile( x.S.StatsDir, 'AcquisitionTime.jpg' );
        [N X]               = hist(round(AcquisitionTimeN./100)./100);
        N                   = N./sum(N);
        fig                 = figure('Visible','off');
        plot(X,N);
        title('Total population histogram of time of ASL scan');
        xlabel('Time (hours)');
        ylabel('Normalized frequency (%)');
        print(gcf,'-djpeg','-r200', NameSave );
    end
end


end

%% -----------------------------------------------------------------------------------------------
% % % %% Restructure for multiple time points, check whether they are similar
% % % NextN   = 1;
% % % for ii=1:length(AcquisitionTime)
% % %     if  ii+2<length(AcquisitionTime)
% % %     
% % %         nameTP1     = AcquisitionTime{ii  ,1}(1:6);
% % %         nameTP2     = AcquisitionTime{ii+2,1}(1:6);
% % %         TP          = str2num(AcquisitionTime{ii  ,1}(end));
% % %         TP2         = str2num(AcquisitionTime{ii+2,1}(end));
% % %         if  strcmp(AcquisitionTime{ii,2},'ASL_1') && strcmp(AcquisitionTime{ii+2,2},'ASL_1') && TP==1 && TP2==2 && strcmp(nameTP1,nameTP2)    % only use first sessions for simplicity
% % %             NewStats(NextN,1)   = AcquisitionTime{ii,3};
% % %             NewStats(NextN,2)   = AcquisitionTime{ii+2,3};
% % %             NextN   = NextN + 1;
% % %         end
% % %     end
% % % end
% % % 
% % % NewStats    = round(NewStats./100)./100;
% % % 
% % % LowLim  = 9;
% % % UpLim   = 17;
% % % figure(1); plot(NewStats(:,1),NewStats(:,2),'.',[LowLim UpLim],[LowLim UpLim],'k');
% % % axis([LowLim UpLim LowLim UpLim])
% % % xlabel('Scan time baseline (hrs)');
% % % ylabel('Scan time follow up (hrs)');
% % % title('Comparison scan times baseline & follow-up')
% % %         
% % %         NextN(TP)   = NextN(TP)+1;
% % %     % Paste it where it should belong
% % %     for iL=1
% % %         NewStats(Index2,iN+(iL-1)*3+TP)   = AcquisitionTime(ii,1+iL);
% % %     end      
