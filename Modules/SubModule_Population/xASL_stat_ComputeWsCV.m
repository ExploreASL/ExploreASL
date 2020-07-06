function xASL_stat_ComputeWsCV(x)
%xASL_stat_ComputeWsCV Calculates the within and between-subject
%coefficient of variance (wsCV and bsCV respectively), to estimate the
%power to detect effects
%
% This requires 4D images that have been split
%

if usejava('jvm')
    fprintf('Skipping xASL_stat_ComputeWsCV, missing JVM\n');
    return;
end

CheckList1                   = xASL_adm_GetFileList(x.D.PopDir,'^PWI_part1_.*\.(nii|nii\.gz)$','FPListRec',[0 Inf]);
CheckList2                   = xASL_adm_GetFileList(x.D.PopDir,'^PWI_part2_.*\.(nii|nii\.gz)$','FPListRec',[0 Inf]);

if   isempty(CheckList1) || isempty(CheckList2)
     fprintf('%s\n','Skipping intra-scan reproducibility calculation, since there are no splitted ASL images available');

else
     fprintf('%s\n','Loading files to computing intra-scan reproducibility...  ');





    %% ------------------------------------------------------------------------------------------------------------
    %% Obtain within-subject coefficient of variance (wsCV)
    pGM                         = xASL_io_Nifti2Im(x.D.ResliceRef);
    pGM                         = pGM>0.7*(max(pGM(:)));

    AbsenceList     = '';

    for iS=1:x.nSubjects
        xASL_TrackProgress(iS,x.nSubjects);
        for iSess=1:x.nSessions
            iSubjSess   = (iS-1)*x.nSessions+iSess;

            FileN{1}   = xASL_adm_GetFileList(x.D.PopDir,['^PWI_part1_' x.SUBJECTS{iS} '_' x.SESSIONS{iSess} '\.(nii|nii\.gz)$'],'FPListRec',[0 Inf]);
            FileN{2}   = xASL_adm_GetFileList(x.D.PopDir,['^PWI_part2_' x.SUBJECTS{iS} '_' x.SESSIONS{iSess} '\.(nii|nii\.gz)$'],'FPListRec',[0 Inf]);

            if ~isempty(FileN{1}) && ~isempty(FileN{2})
                %% Calculate CBF & spatial CoV within GM for each ASL scan
                for ii=1:2
                    IM{ii}                      = xASL_io_Nifti2Im(FileN{ii}{1});
                end
                MaskIM                          = pGM & isfinite(IM{1}) & isfinite(IM{2});
                for ii=1:2
                    %GM{ii}                      = IM{ii}(MaskIM);
                    %ParaM{1}(iSubjSess,ii)      = mean(GM{ii}); % CBF
                    %ParaM{2}(iSubjSess,ii)      = 100*( ParaM{1}(iSubjSess,ii)/std(GM{ii}) ); % spatial CoV
					ParaM{1}(iSubjSess,ii)      = xASL_stat_ComputeMean(IM{ii},MaskIM); % CBF
                    ParaM{2}(iSubjSess,ii)      = 100*xASL_stat_ComputeSpatialCoV(IM{ii},MaskIM); % spatial CoV
                end
            else
                % Note if not present
                AbsenceList{end+1}  = [x.SUBJECTS{iS} '_' x.SESSIONS{iSess}];
            end
        end
    end
    fprintf('\n');

    if  length(AbsenceList)<x.nSubjects

        % Normalize CBF to 60
        ParaM{1}    = ParaM{1}.*60./xASL_stat_MeanNan(ParaM{1}(:));
        % Spatial CoV doesn't need scaling, already normalized




        %% ------------------------------------------------------------------------------------------------------------
        % Get Bland-Altman parameters & plot them
        for iP=1:2
            Fig(iP)         = figure('Visible','off');

            MeanPar{iP}     = (ParaM{iP}(:,1) + ParaM{iP}(:,2)) ./2;
            MeanAll(iP)     = xASL_stat_MeanNan(ParaM{iP}(:));
            SdPar(iP)       = xASL_stat_StdNan(ParaM{iP}(:));
            DiffPar{iP}     = ParaM{iP}(:,1) - ParaM{iP}(:,2);
            MeanDiff(iP)    = xASL_stat_MeanNan(DiffPar{iP});
            SdDiff(iP)      = xASL_stat_StdNan(DiffPar{iP});
            HiLOA(iP)       = MeanDiff(iP)+1.96*SdDiff(iP);
            LoLOA(iP)       = MeanDiff(iP)-1.96*SdDiff(iP);

            wsCV(iP)        = 100.*(SdDiff(iP) ./ MeanAll(iP));
            bsCV(iP)        = 100.*(SdPar(iP) ./ MeanAll(iP));

            minX            = 100*floor(min(MeanPar{iP})/100);
            maxX            = 100*ceil(max(MeanPar{iP})/100);
            minY            = 100*floor( min([LoLOA(iP)  min(DiffPar{iP})]) /100);
            maxY            = 100* ceil( max([HiLOA(iP)  max(DiffPar{iP})]) /100);
            plot(MeanPar{iP},DiffPar{iP},'ko');
            hold on
            plot([minX maxX],[0 0],'k:');
            hold on
            plot([minX maxX],[MeanDiff(iP) MeanDiff(iP)],'r-');
            hold on
            plot([minX maxX],[HiLOA(iP) HiLOA(iP)],'r--');
            hold on
            plot([minX maxX],[LoLOA(iP) LoLOA(iP)],'r--');
            
            axis([minX maxX+eps minY maxY+eps]); % eps is to avoid too small axes

            xlabel('mean 1st & 2nd half');
            ylabel('delta 1st & 2nd half');

            if  iP==1
                title(['Within-scan repeatability CBF (mL/100g/min): wsCV=' num2str(wsCV(iP),3) '%']);
                SavePath    = fullfile(x.S.StatsDir,'wsCV_CBF.png');
            else
                title(['Within-scan repeatability spatial CoV (%): wsCV=' num2str(wsCV(iP),3) '%']);
                SavePath    = fullfile(x.S.StatsDir,'wsCV_spatialCoV.png');
            end

            print(SavePath, '-dpng');
        end




        %% ------------------------------------------------------------------------------------------------------------
        % Calculation effect size
        SampleSize  = x.nSubjects;

        Z_beta              = 0.841621;              % Z(1-beta)  = Z(0.8)    = 0.84;
        Z_alpha             = 1.959964;              % Z(1-alpha) = Z(1-(0.05/2)) = Z(0.975) = 1.96;
        % Z-score calculator -> www.fourmilab.ch/rpkp/experiments/analysis/zCalc.html
        ConstantF           = 2*(Z_beta+Z_alpha)^2;

        UnitN               = {'mL/100g/min' '%'};

        MinEffSize          = SdDiff ./ (SampleSize ./ ConstantF).^0.5;
        MinEffSizeRel       = 100.*(MinEffSize./MeanAll);

        for iP=1:2
            fprintf('%s\n',['We can detect longitudinal/challenge effects larger than ' num2str(MinEffSize(iP),3) ' ' UnitN{iP} ' (=' num2str(MinEffSizeRel(iP),3) '%)']);
        end

        MinEffSize          = SdPar ./ (SampleSize ./ ConstantF).^0.5;
        MinEffSizeRel       = 100.*(MinEffSize./MeanAll);

        for iP=1:2
            fprintf('%s\n',['We can detect cross-sectional effects larger than ' num2str(MinEffSize(iP),3) ' ' UnitN{iP} ' (=' num2str(MinEffSizeRel(iP),3) '%)']);
        end

    end
    close all

end

end
