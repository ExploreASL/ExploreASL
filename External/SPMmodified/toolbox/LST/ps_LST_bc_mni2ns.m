function out = ps_LST_bc_mni2ns(mni, type, Vf2, indx_mni, indx_brain, bp)
% Applies the inverse deformation field to the information that is stored
% in MNI space.

switch type
    case 'bp'
        
        % save brain position in MNI space to disk
        img = zeros(121, 145, 121);
        img(mni) = mni;
        
        V = spm_vol(fullfile(spm('dir'), 'tpm', 'TPM.nii'));
        V = V(1);
                
        V.fname = 'LST_bp_mni.nii';
        spm_write_vol(V, img);
        
        % Apply inverse deformation field
        clear job
        namFlair = Vf2.fname(2:end);
        job.comp{1}.def = {['iy_', namFlair]};        
        job.out{1}.pull.fnames = {'LST_bp_mni.nii'};
        job.out{1}.pull.savedir.savepwd = 1;
        job.out{1}.pull.interp = 0;
        job.out{1}.pull.mask = 1;
        job.out{1}.pull.fwhm = [0 0 0];
        spm_deformations(job);
        
        out = spm_read_vols(spm_vol('wLST_bp_mni.nii'));
        spm_unlink('LST_bp_mni.nii')
        spm_unlink('wLST_bp_mni.nii')
        
        out(:,:,1) = 0 .* out(:,:,1); out(:,:,end) = 0 .* out(:,:,end);
        out(:,1,:) = 0 .* out(:,1,:); out(:,end,:) = 0 .* out(:,end,:);
        out(1,:,:) = 0 .* out(1,:,:); out(end,:,:) = 0 .* out(end,:,:);
        
    case 'tissue'
        
        mni_img = zeros(121, 145, 121); mni_img(indx_mni) = mni;
        out = zeros(Vf2.dim(1), Vf2.dim(2), Vf2.dim(3));        
        out(indx_brain) = mni_img(bp(indx_brain));
    
end

end