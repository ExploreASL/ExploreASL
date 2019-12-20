function SPMDesignT2Loop( x, nSubjects, CoVaryMean )
%SPMDesignT2Loop = Wrapper for SPM factorial design 2-sample t-test
%PM: make into 1 or 2 sample t-test, dependent upon groups/conditions

    % These are mostly standard settings
    % Last 2 options (global mean normalization) are open to discussion

    CurrentDir              = cd;
    cd( x.StatsCmpDir ); % just to be sure, this is what SPM wants for some of the below activities
    SPMmat                                              = fullfile( x.StatsCmpDir, 'SPM.mat');
    spm('defaults','fmri');

    xASL_adm_DeleteFileList(x.StatsCmpDir,'^beta_.*$');
    delete('con_0001.nii');
    delete('mask.nii');
    delete('ResMS.nii');
    delete('RPV.nii');
    delete('SPM.mat');
    delete('spmT_0001.nii');
    xASL_adm_DeleteFileList(x.StatsCmpDir,'^FilterSPM\.(nii|nii\.gz)$');


    %% 1    Factorial design


    matlabbatch{1}.spm.stats.factorial_design.dir                   = {x.StatsCmpDir}; % output dir for design matrix
    matlabbatch{1}.spm.stats.factorial_design.des.t2.dept           = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.variance       = 1; % assumes unequal variance
    matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca          = 0; % grand mean scaling
    matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova         = 0; % allows different subjects to have different relationships between local & global measures (not for 2nd level)
    matlabbatch{1}.spm.stats.factorial_design.cov                   = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov             = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none    = 1; % threshold masking
    matlabbatch{1}.spm.stats.factorial_design.masking.im            = 1; % implicit masking, masks out zeros/NaNs
    matlabbatch{1}.spm.stats.factorial_design.masking.em            = {''};

    % Subjects

    scans1                                                          = xASL_adm_GetImageList3D( x.StatsCmpDir, '^CBF1_\d+\.(nii|nii\.gz)$');
    scans2                                                          = xASL_adm_GetImageList3D( x.StatsCmpDir, '^CBF2_\d+\.(nii|nii\.gz)$');

    if  exist('nSubjects','var')
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1         = scans1(1:nSubjects);
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2         = scans2(1:nSubjects);
    else
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1         = scans1;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2         = scans2;
    end

    % NB: because it is 2-sample (unpaired) test, sorting doesn't matter here. For other purposes, make sure
    % numbering is 01 02 etc.



    % Co-varying for global mean
    if  exist('CoVaryMean','var')
        if CoVaryMean==1 % normalization global mean
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_mean                = 1;  % calculate individual whole-image mean
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no        = 1;  % don't scale, rather co-vary for global mean
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm               = 3;  % co-vary for global mean
        end

    else    % no normalization with global mean

        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit        = 1; % don't take mean into account
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no= 1; % grand mean scaling
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm       = 1; % global normalization
    end

    spm_jobman('run',matlabbatch);

    % Results:
    % SPM.mat
    % Contains
    % SPM.xY.P  = nifti FileNames
    % SPM.xY.VP = nifti information (objects/memory mapping/matrix information etc.)
    % SPM.xX.X  = design
    % SPM.xX.I  = design (ScanNumber, GroupNumber, 1 1 (mean scaling/ancova probably)




%     Skip model review
%     
%     clear matlabbatch
%     matlabbatch{1}.spm.stats.review.spmmat              = { SPMmat };
%     matlabbatch{1}.spm.stats.review.display.matrix      = 1;
%     matlabbatch{1}.spm.stats.review.print               = 'pdf';
%     spm_jobman('run',matlabbatch);
% 
%     close all

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
%     matlabbatch{1}.spm.stats.con.consess{2}.tcon.name   = ContrastName;
%     matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights= [1 -1];
%     matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep= 'none';
    matlabbatch{1}.spm.stats.con.delete                 = 0;
    spm_jobman('run',matlabbatch);

    % Performs the t-test & creates T-statistical map (i.e. SPM, Statistical Parametric Map)

%     %% 4    Results print
%     clear matlabbatch
%     matlabbatch{1}.spm.stats.results.spmmat             = { SPMmat };
%     matlabbatch{1}.spm.stats.results.conspec.titlestr   = '';
%     matlabbatch{1}.spm.stats.results.conspec.contrasts  = Inf; % use all existing contrasts
%     matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE'; % family-wise error
%     matlabbatch{1}.spm.stats.results.conspec.thresh     = 0.05; % p-value threshold
%     matlabbatch{1}.spm.stats.results.conspec.extent     = 0;
%     matlabbatch{1}.spm.stats.results.conspec.mask.none  = 1;
%     matlabbatch{1}.spm.stats.results.units              = 1;
%     matlabbatch{1}.spm.stats.results.print              = 'pdf';
%     matlabbatch{1}.spm.stats.results.write.none         = 1;
%     spm_jobman('run',matlabbatch);
%     close all

%     %% 5    Write results image file
% 
%     xSPMinput.swd           = x.StatsCmpDir; % directory where to search for "SPM.mat"
%     xSPMinput.title         = ContrastName;
%     xSPMinput.Ex            = 'none';
%     xSPMinput.u             = 0.05;
%     xSPMinput.k             = 0;
%     xSPMinput.thresDesc     = 'FWE';
%     xSPMinput.Ic            = 1; % indice of contrasts within xCon
%     xSPMinput.Im            = []; % indice of masking contrasts within xCon
%     xSPMinput.pm            = []; % p-value for masking???


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

%     [SPM xSPM]      = spm_getSPM( xSPMinput );
%     FilterFile      = fullfile( x.StatsCmpDir, 'FilterSPM.nii');
%     spm_write_filtered(xSPM.Z,xSPM.XYZ,xSPM.DIM,xSPM.M,'SPM-filtered',FilterFile);
    
    cd(CurrentDir);

end

