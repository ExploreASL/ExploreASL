function [x] = xASL_spm_GLMresults(x)
%xASL_spm_GLMresults ExploreASL wrapper for SPM GLM - results part
% Creates mask of significant clusters


    %% Load contrast & make absolute (+ve)
    if      x.S.nSets>2
            ContrastMapFile     = fullfile( x.S.SPMdir, 'spmF_0001.nii');
    else
            ContrastMapFile     = fullfile( x.S.SPMdir, 'spmT_0001.nii');
    end


    ContrastMap                 = abs(xASL_io_Nifti2Im( ContrastMapFile ));
    xASL_io_SaveNifti( ContrastMapFile, ContrastMapFile, ContrastMap ,[],0);


    %% -------------------------------------------------------------------------------------
    %% Delete previous runs
    xASL_adm_DeleteFileList(x.S.SPMdir,'^spm(T|F)_.*Masked\.(nii|nii\.gz)$');
    % Reload mat
    load(x.S.SPMmat3);
    save(x.S.SPMmat,'SPM');


    matlabbatch{1}.spm.stats.results.spmmat             = { x.S.SPMmat };
    matlabbatch{1}.spm.stats.results.conspec.titlestr   = '';
    matlabbatch{1}.spm.stats.results.conspec.contrasts  = Inf; % use all existing contrasts
    matlabbatch{1}.spm.stats.results.conspec.threshdesc = x.S.ThreshType; % 'FWE'; % family-wise error

    if  strcmp(x.S.MultiComparisonCorrType,'cluster')
        matlabbatch{1}.spm.stats.results.conspec.thresh = x.S.clusterPthr;  % height threshold
    else
        matlabbatch{1}.spm.stats.results.conspec.thresh = x.S.uncorrThresh; % height threshold
    end

    matlabbatch{1}.spm.stats.results.conspec.extent     = x.S.ClusterExtent;  % width threshold (k)
    matlabbatch{1}.spm.stats.results.conspec.mask.none  = 1; % masking is already done
    matlabbatch{1}.spm.stats.results.units              = 1;
    matlabbatch{1}.spm.stats.results.print              = 'pdf';

    if  x.S.PrintSPMOutput
        spm_get_defaults('cmdline',false);
    else
        spm_get_defaults('cmdline',true);
    end

    matlabbatch{1}.spm.stats.results.write.tspm.basename = 'Masked'; % save masked t-stats
    spm_jobman('run',matlabbatch);
    close all
    load(x.S.SPMmat);
    save(x.S.SPMmat4,'SPM');


end
