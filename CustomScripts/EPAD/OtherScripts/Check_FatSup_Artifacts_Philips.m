%% All sites together

for iS=1:x.nSubjects
    xASL_TrackProgress(iS,x.nSubjects);
    Path    = fullfile(x.D.PopDir,['SD_' x.SUBJECTS{iS} '_ASL_1.nii']);
    if  exist(Path, 'file')
        tIM             = xASL_io_Nifti2Im(Path);
        tIM(40:80,40:80,25:100)     = tIM(40:80,40:80,25:100).*2;
        sumSD(iS,1)     = xASL_stat_SumNan(tIM(:));
    end
end



Site    = x.S.SetsID(:,7);
nSites  = unique(Site);
clear sumSD

for iSite=1:nSites
    SiteN   = Site==iSite;

    for iS=1:x.nSubjects
        xASL_TrackProgress(iS,x.nSubjects);

        if  SiteN(iS)
            Path    = fullfile(x.D.PopDir,['SD_' x.SUBJECTS{iS} '_ASL_1.nii']);
            if  exist(Path, 'file')
                tIM                         = xASL_io_Nifti2Im(Path);
                tIM(40:80,40:80,25:100)     = tIM(40:80,40:80,25:100).*2;
                sumSD{iSite}(iS,1)          = xASL_stat_SumNan(tIM(:));
            end
        end
    end

    MasN            = isfinite(sumSD{iSite}) & sumSD{iSite}~=0;
    sumSD{iSite}   = sumSD{iSite}(MasN);
    CSFvol          = x.S.SetsID(:,3);
    CSFvol          = CSFvol(MasN);

figure(1);plot(CSFvol,sumSD{iSite},'.')
xlabel('CSF volume (L)');
ylabel('Sum temporalSD signal');
title([x.S.SetsOptions{7}{4} ' correlation fat suppression artifact ventricle size']);
[coef, pval] = corr(CSFvol(MaskN), sumSD(MaskN));
