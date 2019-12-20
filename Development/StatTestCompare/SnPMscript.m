function SnPMscript( x, nSubjects, CoVaryMean, nPerm, PermutationRun, InferenceOption )
%SnPMscript = Wrapper for SnPM factorial design 2-sample t-test
%PM: make into 1 or 2 sample t-test, dependent upon groups/conditions
% InferenceOption should be one of: {'uncorrected' 'FWE' 'cluster1' 'cluster2' 'cluster3' 'cluster4'}

    % These are mostly standard settings
    % Last 2 options (global mean normalization) are open to discussion
% addpath(genpath('c:\ASL_pipeline_HJ'))
    spm('defaults','fmri');


    if PermutationRun

        %% 0    Delete existing files
        xASL_adm_DeleteFileList(x.StatsCmpDir,'^beta_.*$');
        xASL_adm_DeleteFileList(x.StatsCmpDir,'^lP.*$');
        xASL_adm_DeleteFileList(x.StatsCmpDir,'^ResMS.*$');
        xASL_adm_DeleteFileList(x.StatsCmpDir,'^SnPM.*$');
        xASL_adm_DeleteFileList(x.StatsCmpDir,'^snpm.*$');
        xASL_adm_DeleteFileList(x.StatsCmpDir,'^spm.*$');
        xASL_adm_DeleteFileList(x.StatsCmpDir,'^XYZ.*$');



        %% 1    Design the permutation

        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.DesignName               = '2 Groups: Two Sample T test; 1 scan per subject';
        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.DesignFile               = 'snpm_bch_ui_TwoSampT';
        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.dir                      = { x.StatsCmpDir };
        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov                      = struct('c', {}, 'cname', {});
        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.nPerm                    = nPerm; % max n permutations
        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.vFWHM                    = [0 0 0]; % variance smoothing
        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.bVolm                    = 1;
        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.ST.ST_later              = -1;
        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.tm.tm_none       = 1;
        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.im               = 1;    % implicit masking
        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking.em               = {''}; % explicit masking

        % Subjects
        scans1                                                          = xASL_adm_GetImageList3D( x.StatsCmpDir, '^CBF1_\d+\.(nii|nii\.gz)$');
        scans2                                                          = xASL_adm_GetImageList3D( x.StatsCmpDir, '^CBF2_\d+\.(nii|nii\.gz)$');

        if  exist('nSubjects','var')
            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans1         = scans1(1:nSubjects);
            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans2         = scans2(1:nSubjects);
        else
            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans1         = scans1;
            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans2         = scans2;
        end

        % Co-varying for global mean
        if  exist('CoVaryMean','var')
            if CoVaryMean==1 % normalization global mean
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalc.g_mean                = 1;  % calculate individual whole-image mean
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.gmsca.gmsca_no        = 1;  % don't scale, rather co-vary for global mean
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.glonorm               = 3;  % co-vary for global mean
            end

        else    % no normalization with global mean

            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalc.g_omit        = 1; % don't take mean into account
            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.gmsca.gmsca_no= 1; % grand mean scaling
            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm.glonorm       = 1; % global normalization
        end

        tic
        spm_jobman('run',matlabbatch);

        %% 2    Run the permutation
        clear matlabbatch
        matlabbatch{1}.spm.tools.snpm.cp.snpmcfg                            = { fullfile( x.StatsCmpDir, 'SnPMcfg.mat') };
        spm_jobman('run',matlabbatch);
        toc

    end



    %% 3    Show results

    if exist('InferenceOption','var')

        FilterFile      = fullfile( x.StatsCmpDir,'SnPM_filtered.nii');
        delete( FilterFile );
        xASL_adm_DeleteFileList(x.StatsCmpDir,'^STCS.*$');

        clear matlabbatch
        matlabbatch{1}.spm.tools.snpm.inference.SnPMmat                     = { fullfile( x.StatsCmpDir, 'SnPM.mat') };
        matlabbatch{1}.spm.tools.snpm.inference.Tsign                       = 1; % positive effects -> check sign!
        matlabbatch{1}.spm.tools.snpm.inference.WriteFiltImg.name           = 'SnPM_filtered';
        matlabbatch{1}.spm.tools.snpm.inference.Report                      = 'MIPtable';



        switch InferenceOption
        case 'uncorrected'
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Vox.VoxSig.Pth      = 0.001;
        case 'FWE'
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Vox.VoxSig.FWEth    = 0.05;
        case 'cluster1'
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = 0.01;
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.FWEthC = 0.05;
        case 'cluster2'
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = 0.001;
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.FWEthC = 0.05;
        case 'cluster3'
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusMass.PFilt = 0.05;
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusMass.PrimThresh = 0.01;
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusMass.Theta = 0.5;
        case 'cluster4'
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusMass.PFilt = 0.05;
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusMass.PrimThresh = 0.001;
            matlabbatch{1}.spm.tools.snpm.inference.Thr.Clus.ClusMass.Theta = 0.5;
        otherwise error('Invalid InferenceOption!!!');
        end

        spm_jobman('run',matlabbatch);
        close
        close

        % voxel output when writing SnPM_filtered file is equal to the mean GM mask >0.5 probabilities, i.e. 306626 voxels

        if ~exist( FilterFile ) % create empty dummy nifti if no voxels survived thresholding
            ExampleNii      = fullfile( x.D.PopDir, ['DARTEL_c1T1_' x.SUBJECTS{1} '.nii']);
            IM              = xASL_io_ReadNifti( ExampleNii );
            IM              = IM.dat(:,:,:);
            IM(:,:,:)       = 0;

            xASL_io_SaveNifti( ExampleNii, FilterFile, IM );
        end
    end



end
%
%                     % THIS GAVE AN ERROR

%                     % 5) FWE p=0.05 based on cluster-mass (p=0.01 initially, SPM standard)
%                     SnPMscript( x, iS       , CovMean-1, 10000, 0, 'cluster3' );
%                     StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'SnPM_filtered.nii'));
%                     StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
%                     H                             = StatResult>0;
%
%                     I(1)        = (CovMean*6)-1;
%                     [TPR(I(1),I(2),I(3),I(4),I(5)) FPR(I(1),I(2),I(3),I(4),I(5))]   = CalcROC( GoldPos, H);
%
%                     % 6) FWE p=0.05 based on cluster-mass (p=0.001 initially, SPM standard)
%                     SnPMscript( x, iS       , CovMean-1, 10000, 0, 'cluster4' );
%                     StatResult                    = xASL_io_ReadNifti( fullfile( x.StatsCmpDir, 'SnPM_filtered.nii'));
%                     StatResult                    = StatResult.dat(:,:,:); % Load SPM T-map
%                     H                             = StatResult>0;
%
%                     I(1)        = (CovMean*6)-0;
%                     [TPR(I(1),I(2),I(3),I(4),I(5)) FPR(I(1),I(2),I(3),I(4),I(5))]   = CalcROC( GoldPos, H);
%
% Running 'Inference'
%
% SnPM: snpm_combo_pp
% ========================================================================
% Failed  'Inference'
% Undefined function or variable 'pU_ST_Ut_filt'.
% In file "c:\werk\spm12\toolbox\snpm13\snpm_combo_pp.m" (???), function "snpm_combo_pp" at line 224.
% In file "c:\werk\spm12\toolbox\snpm13\config\snpm_run_pp.m" (v1716), function "snpm_run_pp" at line 22.
%
% The following modules did not run:
% Failed: Inference
%
% ??? Error using ==> MATLABbatch system
% Job execution failed. The full log of this run can be found in MATLAB
% command window, starting with the lines (look for the line showing the exact
% #job as displayed in this error message)
% ------------------
% Running job #2
% ------------------
%
