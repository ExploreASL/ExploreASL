%% Convert analyze->nifti
%% PreDIVA turbo-QUASAR

ROOT    = 'E:\Backup\ASL_E\Pre_DIVA_transitt\Esben\Result';
dlist   = xASL_adm_GetFileList( ROOT, '^.*\.img$', 'FPListRec');

for iList=1:length(dlist)
    
    [path,file,ext]                                 = fileparts( dlist{iList} );
    
    matlabbatch{1}.spm.util.imcalc.input            = {[dlist{iList} ',1']};
    matlabbatch{1}.spm.util.imcalc.output           = [file '.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir           = {path};
    matlabbatch{1}.spm.util.imcalc.expression       = 'i1';
    matlabbatch{1}.spm.util.imcalc.var              = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx     = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask     = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp   = -4;
    matlabbatch{1}.spm.util.imcalc.options.dtype    = 16;
    
    spm_jobman('run',matlabbatch);
    clear path file ext
    delete( dlist{iList} );
end


dlist   = xASL_adm_GetFileList( ROOT, '^.*\.hdr$', 'FPListRec');

for iList=1:length(dlist)
    delete( dlist{iList} );
end

dlist   = xASL_adm_GetFileList( ROOT, '^ATT_11.(nii|nii\.gz)$', 'FPListRec');

for iList=1:length(dlist)
    [path,file,ext]         = fileparts( dlist{iList} );
    xASL_Rename( dlist{iList}, 'ATT.nii' );
end


dlist   = xASL_adm_GetFileList( ROOT, '^R1_11.(nii|nii\.gz)$', 'FPListRec');

for iList=1:length(dlist)
    [path,file,ext]         = fileparts( dlist{iList} );
    xASL_Rename( dlist{iList}, 'R1.nii' );
end

dlist   = xASL_adm_GetFileList( ROOT, '^Raw_11.(nii|nii\.gz)$', 'FPListRec');

for iList=1:length(dlist)
    [path,file,ext]         = fileparts( dlist{iList} );
    xASL_Rename( dlist{iList}, 'Raw.nii' );
end
