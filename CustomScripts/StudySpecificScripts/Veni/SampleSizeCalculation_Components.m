ROOT    = 'C:\Backup\ASL\GENFI\GENFI_DF2_2_Long\analysis\dartel';

VendorName  = {'GE' 'PH_Achieva_Bsup' 'PH_Achieva_noBsup' 'SI_Trio'};

GMmask  = logical(xASL_io_Nifti2Im(fullfile(ROOT,'VBA_mask_final.nii')));


for iV=1:4
    Flist{iV}   = xASL_adm_GetFileList(ROOT, ['^qCBF_untreated_' VendorName{iV} '_.*_ASL_1\.(nii|nii\.gz)$']);
    fprintf('%s\n',['Loading images for ' VendorName{iV} '...  ']);
    for iL=1:40
        xASL_TrackProgress(iL,40);
        tIM             = xASL_io_Nifti2Im(Flist{iV}{iL});
        IM(:,iL,iV)     = tIM(GMmask);
    end
end
fprintf('\n');

% Within RMS
for iV=1:4
    for iL=1:20
        wsRMS(iL,iV)    = (xASL_stat_MeanNan((IM(:,iL,iV) - IM(:,iL+20,iV) ).^2))^0.5;
        MeanInt         = xASL_stat_MeanNan(xASL_stat_MeanNan(IM(:,[iL iL+20],iV)));
        wsRMS(iL,iV)    = wsRMS(iL,iV)/MeanInt;
    end
end

median(wsRMS(:))
mad(wsRMS(:))

% Between RMS
CompPatt  = {[1 2] [1 3]  [1 4] [2 3] [2 4] [3 4]};
for iC=1:length(CompPatt)
    for iL=1:20
        bsRMS(iL,iC)    = (xASL_stat_MeanNan((IM(:,iL,CompPatt{iC}(1)) - IM(:,iL,CompPatt{iC}(2)) ).^2))^0.5;
        MeanInt         = xASL_stat_MeanNan(xASL_stat_MeanNan(IM(:,iL,[CompPatt{iC}(1) CompPatt{iC}(2)])));
        bsRMS(iL,iC)    = bsRMS(iL,iC)/MeanInt;
    end
end

median(bsRMS(:))
mad(bsRMS(:))

for iW=1:4
    for iB=1:6
        for iT=1:20
            IndexN  = (iW-1)*6+iB;
            [p(iT,IndexN),h(iT,IndexN),stats] = ranksum(bsRMS(1:iT,iB),wsRMS(1:iT,iW));
        end
    end
end


figure(1);plot([2:20],median(p(2:end,:),2));
xlabel('Sample size (n scans/subjects');
ylabel('Median p-value');
title('Between-sequence normRMS vs. within-sequence normRMS');
