function SPMCorrClusTh2( x, p1, GoldPos )
%SPMCorrClusTh2 Summary of this function goes here
%   Detailed explanation goes here

    clear matlabbatch
    matlabbatch{1}.spm.stats.results.spmmat                 = { fullfile( x.StatsCmpDir, 'SPM.mat') };
    matlabbatch{1}.spm.stats.results.conspec.titlestr       = '';
    matlabbatch{1}.spm.stats.results.conspec.contrasts      = 1;
    matlabbatch{1}.spm.stats.results.conspec.threshdesc     = 'none';
    matlabbatch{1}.spm.stats.results.conspec.thresh         = 0.5;
    matlabbatch{1}.spm.stats.results.conspec.mask.none      = 1;
    matlabbatch{1}.spm.stats.results.units                  = 1;
    matlabbatch{1}.spm.stats.results.print                  = false;
%     matlabbatch{1}.spm.stats.results.write.binary.basename  = 'filter';
    matlabbatch{1}.spm.stats.results.write.tspm.basename = 'filter';
%     matlabbatch{1}.spm.stats.results.write.nary.basename = 'filter';

    %% Loop over cluster volume thresholds
    ItRes   = 100; % nSteps
    kMin    = 1; % minimal cluster size to search for
    kMax    = 1000000;

    for It=1:ItRes It=5
        StatCalc(It,1)                                      = round(((kMax-kMin)/ItRes)*It);
        matlabbatch{1}.spm.stats.results.conspec.extent     = StatCalc(It,1); % cluster-size threshold

        spm_jobman('run',matlabbatch);
        close
        close

        FilterFile          = fullfile( x.StatsCmpDir, 'spmT_0001_filter.nii');

        FilterTemp          = xASL_io_ReadNifti( FilterFile );
        FilterTemp          = FilterTemp.dat(:,:,:);

        FigureOut   = fullfile( x.StatsCmpDir, ['ROC_' num2str(iSub) '_subjects_vol_' num2str(ActivationVolume) '_CLUSTER.jpg']);
        [ StatCalc AUC(iVol,1)  ] = GetRocCurve( x, GoldPos, FigureOut, fullfile( x.StatsCmpDir, 'spmT_0001_filter.nii') );

        [StatCalc(It,2) StatCalc(It,3)]     = CalcROC( GoldPos, FilterTemp); % true positive rate, false positive rate
        delete( FilterFile );
    end



    % Housekeeping
    delete( fullfile( x.StatsCmpDir, 'beta_0001.nii') );
    delete( fullfile( x.StatsCmpDir, 'beta_0002.nii') );
    delete( fullfile( x.StatsCmpDir, 'con_0001.nii') );
    delete( fullfile( x.StatsCmpDir, 'mask.nii') );
    delete( fullfile( x.StatsCmpDir, 'ResMS.nii') );
    delete( fullfile( x.StatsCmpDir, 'RPV.nii') );
    delete( fullfile( x.StatsCmpDir, 'SPM.mat') );
    delete( fullfile( x.StatsCmpDir, 'spmT_0001.nii') );
    xASL_adm_DeleteFileList( x.StatsCmpDir, '^CBF2_\d+\.(nii|nii\.gz)$');



end
