% This code created Figure with screenshots of CBF contrast, vascular contract & poor SNR categories

x = ExploreASL_Master('',0);

PathHarmy   = 'C:\BackupWork\ASL\Harmy\analysis_Harmy\Population';

NIIpath{1}{1}  = fullfile(PathHarmy,'qCBF_untreated_HD062_1_ASL_1.nii');
NIIpath{1}{2}  = fullfile(PathHarmy,'qCBF_untreated_HD098_1_ASL_1.nii');
% NIIpath{1}{3}  = fullfile(PathHarmy,'qCBF_untreated_HD017_2_ASL_1.nii');
% NIIpath{1}{4}  = fullfile(PathHarmy,'qCBF_untreated_HD150_1_ASL_1.nii');
% NIIpath{1}{5}  = fullfile(PathHarmy,'qCBF_untreated_HD026_1_ASL_1.nii');
% NIIpath{1}{6}  = fullfile(PathHarmy,'qCBF_untreated_HD036_2_ASL_1.nii');

% NIIpath{2}{1}  = fullfile(PathHarmy,'qCBF_untreated_HD462_1_ASL_1.nii');
% NIIpath{2}{2}  = fullfile(PathHarmy,'qCBF_untreated_HD481_1_ASL_1.nii');
NIIpath{2}{1}  = fullfile(PathHarmy,'qCBF_untreated_HD508_1_ASL_1.nii');
% NIIpath{2}{4}  = fullfile(PathHarmy,'qCBF_untreated_HD572_1_ASL_1.nii');
% NIIpath{2}{5}  = fullfile(PathHarmy,'qCBF_untreated_HD563_1_ASL_1.nii');
NIIpath{2}{2}  = fullfile(PathHarmy,'qCBF_untreated_HD578_1_ASL_1.nii');

% NIIpath{3}{1}  = fullfile(PathHarmy,'qCBF_untreated_HD360_1_ASL_1.nii');
% NIIpath{3}{2}  = fullfile(PathHarmy,'qCBF_untreated_HD408_1_ASL_1.nii');
% NIIpath{3}{3}  = fullfile(PathHarmy,'qCBF_untreated_HD304_1_ASL_1.nii');
NIIpath{3}{1}  = fullfile(PathHarmy,'qCBF_untreated_HD459_1_ASL_1.nii');
NIIpath{3}{2}  = fullfile(PathHarmy,'qCBF_untreated_HD360_1_ASL_1.nii');
% NIIpath{3}{6}  = fullfile(PathHarmy,'qCBF_untreated_HD403_1_ASL_1.nii');

SaveDir = 'C:\Users\kyrav\Desktop\Gdrive\XploreLab\ProjectsPending\ExploreASL manuscript\FiguresNew\S11_QC\Examples';

for iC=1:length(NIIpath)
    for iL=1:length(NIIpath{iC})
        % Load image
        TempIm = xASL_io_Nifti2Im(NIIpath{iC}{iL});
        % Scale image
        MinN = 0;
        finiteMask = isfinite(TempIm);
        MaxN = mean(TempIm(finiteMask)) + 3*std(TempIm(finiteMask));
        TempIm(TempIm<MinN) = MinN;
        TempIm(TempIm>MaxN) = MaxN;
        % reslice
        x.S.TraSlices = x.S.slicesLarge(3:end-1);
        TempIm = xASL_im_TransformData2View(TempIm, x);
        % save image
        SavePath = fullfile(SaveDir,[num2str(iC) '_' num2str(iL) '.png']);
        fig=figure(1);imshow(TempIm,[],'border','tight');
        print(fig, '-dpng', SavePath);
    end
end
