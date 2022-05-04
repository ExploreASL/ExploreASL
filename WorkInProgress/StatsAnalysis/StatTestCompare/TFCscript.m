function TFCscript( x, PermN )
%TFCscript performs TFC estimation

    delete( fullfile( x.StatsCmpDir),'T_0001.hdr');
    delete( fullfile( x.StatsCmpDir),'T_0001.img');
    delete( fullfile( x.StatsCmpDir),'T_0001.txt');
    delete( fullfile( x.StatsCmpDir),'T_log_p_0001.hdr');
    delete( fullfile( x.StatsCmpDir),'T_log_p_0001.img');
    delete( fullfile( x.StatsCmpDir),'T_log_pFDR_0001.hdr');
    delete( fullfile( x.StatsCmpDir),'T_log_pFDR_0001.img');
    delete( fullfile( x.StatsCmpDir),'T_log_pFWE_0001.hdr');
    delete( fullfile( x.StatsCmpDir),'T_log_pFWE_0001.img');
    delete( fullfile( x.StatsCmpDir),'TFCE_0001.hdr');
    delete( fullfile( x.StatsCmpDir),'TFCE_0001.img');
    delete( fullfile( x.StatsCmpDir),'TFCE_log_p_0001.hdr');
    delete( fullfile( x.StatsCmpDir),'TFCE_log_p_0001.img');
    delete( fullfile( x.StatsCmpDir),'TFCE_log_pFDR_0001.hdr');
    delete( fullfile( x.StatsCmpDir),'TFCE_log_pFDR_0001.img');
    delete( fullfile( x.StatsCmpDir),'TFCE_log_pFWE_0001.hdr');
    delete( fullfile( x.StatsCmpDir),'TFCE_log_pFWE_0001.img');


    SPMmatFile                                                  = fullfile( x.StatsCmpDir, 'SPM.mat');
    MaskFile                                                    = fullfile( x.StatsCmpDir, 'mask.nii');

    matlabbatch{1}.spm.tools.tfce_estimate.spmmat               = { SPMmatFile };
    matlabbatch{1}.spm.tools.tfce_estimate.mask                 = { MaskFile };
    matlabbatch{1}.spm.tools.tfce_estimate.conspec.titlestr     = '';
    matlabbatch{1}.spm.tools.tfce_estimate.conspec.contrasts    = 1;
    matlabbatch{1}.spm.tools.tfce_estimate.conspec.n_perm       = PermN;
    matlabbatch{1}.spm.tools.tfce_estimate.conspec.vFWHM        = 0;
    matlabbatch{1}.spm.tools.tfce_estimate.openmp               = 1;

    spm_jobman('run',matlabbatch);

end