function SPMDesignT2SnPMtrial( x )
%SPMDesignT2SnPMtrial = Wrapper for SPM factorial design 2-sample t-test
%PM: make into 1 or 2 sample t-test, dependent upon groups/conditions

    % These are mostly standard settings
    % Last 2 options (global mean normalization) are open to discussion
% addpath(genpath('c:\ASL_pipeline_HJ'))
    spm('defaults','fmri');

    %% 1    Design the permutation

% x.StatsCmpDir = 'E:\Backup\ASL_E\KCL_stats\INOX\analysis\dartel\StatsCompare\PutamenBLOB\SnPMTrial'

    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.DesignName               = '2 Groups: Two Sample T test; 1 scan per subject';
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.DesignFile               = 'snpm_bch_ui_TwoSampT';
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.dir                      = { x.StatsCmpDir };
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans1                   = xASL_adm_GetImageList3D( x.StatsCmpDir, '^CBF1_\d+\.(nii|nii\.gz)$');
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans2                   = xASL_adm_GetImageList3D( x.StatsCmpDir, '^CBF2_\d+\.(nii|nii\.gz)$');
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov                      = struct('c', {}, 'cname', {});
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.nPerm                    = 5000; % max n permutations
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.vFWHM                    = [0 0 0]; % variance smoothing
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.bVolm                    = 1;
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.ST.ST_none               = 0;
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.tm.tm_none       = 1;
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.im               = 1;    % implicit masking
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.em               = {''}; % explicit masking
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalc.g_omit           = 1;    % mean scaling
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.gmsca.gmsca_no   = 1;
    matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.glonorm          = 1;
    spm_jobman('run',matlabbatch);

    %% 2    Run the permutation
    clear matlabbatch
    matlabbatch{1}.spm.tools.snpm.cp.snpmcfg                            = { fullfile( x.StatsCmpDir, 'SnPMcfg.mat') };
    spm_jobman('run',matlabbatch);

    %% 3    Show results
    matlabbatch{1}.spm.tools.snpm.inference.SnPMmat                     = {'E:\Backup\ASL_E\KCL_stats\INOX\analysis\dartel\StatsCompare\PutamenBLOB\SnPMTrial\SnPM.mat'};
    matlabbatch{1}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FWEth        = 0.05; % this can be voxel/cluster thresholding
    matlabbatch{1}.spm.tools.snpm.inference.Tsign                       = 1; % positive effects
    matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name           = 'SnPM_filtered';
    matlabbatch{1}.spm.tools.snpm.inference.Report                      = 'MIPtable';





















    matlabbatch{1}.spm.stats.factorial_design.dir                   = {x.StatsCmpDir}; % output dir for design matrix
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1         = xASL_adm_GetImageList3D( x.StatsCmpDir, '^CBF1_\d+\.(nii|nii\.gz)$');
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2         = xASL_adm_GetImageList3D( x.StatsCmpDir, '^CBF2_\d+\.(nii|nii\.gz)$');
    % NB: because it is 2-sample (unpaired) test, sorting doesn't matter here. For other purposes, make sure
    % numbering is 01 02 etc.

    matlabbatch{1}.spm.stats.factorial_design.des.t2.dept           = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.variance       = 1; % assumes unequal variance
    matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca          = 0; % grand mean scaling
    matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova         = 0; % allows different subjects to have different relationships between local & global measures (not for 2nd level)
    matlabbatch{1}.spm.stats.factorial_design.cov                   = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov             = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none    = 1; % threshold masking
    matlabbatch{1}.spm.stats.factorial_design.masking.im            = 1; % implicit masking, masks out zeros/NaNs
    matlabbatch{1}.spm.stats.factorial_design.masking.em            = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit        = 1; % don't take mean into account
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no= 1; % grand mean scaling
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm       = 1; % global normalization
    spm_jobman('run',matlabbatch);

    % Results:
    % SPM.mat
    % Contains
    % SPM.xY.P  = nifti FileNames
    % SPM.xY.VP = nifti information (objects/memory mapping/matrix information etc.)
    % SPM.xX.X  = design
    % SPM.xX.I  = design (ScanNumber, GroupNumber, 1 1 (mean scaling/ancova probably)


    clear matlabbatch
    CurrentDir      = cd;
    cd( x.StatsCmpDir );

    SPMmat                                              = fullfile( x.StatsCmpDir, 'SPM.mat');

    matlabbatch{1}.spm.stats.review.spmmat              = { SPMmat };
    matlabbatch{1}.spm.stats.review.display.matrix      = 1;
    matlabbatch{1}.spm.stats.review.print               = 'pdf';
    spm_jobman('run',matlabbatch);

    cd( CurrentDir );
    close all

    %% 2    Model estimation

    clear matlabbatch
    matlabbatch{1}.spm.stats.fmri_est.spmmat            = { SPMmat };
    matlabbatch{1}.spm.stats.fmri_est.write_residuals   = 0; % don't write residual maps
    matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1; % restrict maximum likelihood
    spm_jobman('run',matlabbatch);

    % Creates nifti/object memory mapping Vbeta

    %% 3    Contrast creation

    ContrastName                                        = 'DiffTest';

    clear matlabbatch
    matlabbatch{1}.spm.stats.con.spmmat                 = { SPMmat };
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name   = ContrastName;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights= [-1 1];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep= 'none';
    matlabbatch{1}.spm.stats.con.delete                 = 0;
    spm_jobman('run',matlabbatch);

    % Performs the t-test & creates T-statistical map (i.e. SPM, Statistical Parametric Map)

    %% 4    Results print
    clear matlabbatch
    matlabbatch{1}.spm.stats.results.spmmat             = { SPMmat };
    matlabbatch{1}.spm.stats.results.conspec.titlestr   = '';
    matlabbatch{1}.spm.stats.results.conspec.contrasts  = Inf; % use all existing contrasts
    matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none'; % family-wise error
    matlabbatch{1}.spm.stats.results.conspec.thresh     = 0.001; % p-value threshold
    matlabbatch{1}.spm.stats.results.conspec.extent     = 0;
    matlabbatch{1}.spm.stats.results.conspec.mask.none  = 1;
    matlabbatch{1}.spm.stats.results.units              = 1;
    matlabbatch{1}.spm.stats.results.print              = 'pdf';
    matlabbatch{1}.spm.stats.results.write.none         = 1;
    spm_jobman('run',matlabbatch);
    close all

    %% 5    Write results image file

    xSPMinput.swd           = x.StatsCmpDir; % directory where to search for "SPM.mat"
    xSPMinput.title         = ContrastName;
    xSPMinput.Ex            = 'none';
    xSPMinput.u             = 0.001;
    xSPMinput.k             = 0;
    xSPMinput.thresDesc     = 'none';
    xSPMinput.Ic            = 1; % indice of contrasts within xCon
    xSPMinput.Im            = []; % indice of masking contrasts within xCon
    xSPMinput.pm            = []; % p-value for masking???


% Evaluated fields in xSPM (input)
%
% xSPM      - structure containing SPM, distribution & filtering details
% .swd      - SPM working directory - directory containing current SPM.mat
% .title    - title for comparison (string)
% .Ic       - indices of contrasts (in SPM.xCon)
% .n        - conjunction number <= number of contrasts
% .Im       - indices of masking contrasts (in xCon)
% .pm       - p-value for masking (uncorrected)
% .Ex       - flag for exclusive or inclusive masking
% .u        - height threshold
% .k        - extent threshold {voxels}
% .thresDesc - description of height threshold (string)

    [SPM xSPM]      = spm_getSPM( xSPMinput );
    FilterFile      = fullfile( x.StatsCmpDir, 'FilterSPM.nii');
    spm_write_filtered(xSPM.Z,xSPM.XYZ,xSPM.DIM,xSPM.M,'SPM-filtered',FilterFile);



end

