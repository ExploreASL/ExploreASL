function xASL_test_SPM(TestDirDest, bTestDataUsed)
    % Reset path
    path(pathdef);
    x = ExploreASL();
    % Remove ExploreASL paths
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(genpath(x.opts.MyPath));
    warning('on','MATLAB:rmpath:DirNotFound');
    
    % Add SPM path
    addpath(fullfile(x.opts.MyPath,'External','SPMmodified'));
    
    % Initialize SPM, but only SPM
    spm('defaults','FMRI');
    spm('asciiwelcome');
    spm_jobman('initcfg');
    spm_get_defaults('cmdline',true);
    
    x.settings.Quality = false;
    % Find the first directory and copy out the first T1 and ASL just for SPM testing
    
    Dlist = xASL_adm_GetFileList(TestDirDest,'^.*$','List',[0 Inf], true);
    
    if bTestDataUsed
        % TestDataSet detected
        xASL_Copy(fullfile(TestDirDest,Dlist{1},'rawdata','sub-Sub1','anat','sub-Sub1_T1w.nii'),fullfile(TestDirDest,'T1.nii'));
        xASL_Copy(fullfile(TestDirDest,Dlist{1},'rawdata','sub-Sub1','perf','sub-Sub1_asl.nii'),fullfile(TestDirDest,'ASL4D.nii'));
    else
        % Default
        xASL_Copy(fullfile(TestDirDest,Dlist{1},'derivatives','ExploreASL','001DM_1','T1.nii'),fullfile(TestDirDest,'T1.nii'));
        xASL_Copy(fullfile(TestDirDest,Dlist{1},'derivatives','ExploreASL','001DM_1','ASL_1','ASL4D.nii'),fullfile(TestDirDest,'ASL4D.nii'));
    end
    
    % Read Nifti files
    try
        xASL_io_ReadNifti(fullfile(TestDirDest,'ASL4D.nii'));
        xASL_io_ReadNifti(fullfile(TestDirDest,'T1.nii'));
    catch
        error('ASL and T1 data set import failed...')
    end
    
    
    
    % Test spm_realign on low quality
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.realign.write.data               = {fullfile(TestDirDest,'ASL4D.nii')};
    matlabbatch{1}.spm.spatial.realign.write.roptions.which     = [2 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.interp    = 0; % low quality
    matlabbatch{1}.spm.spatial.realign.write.roptions.wrap      = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.mask      = 1;
    matlabbatch{1}.spm.spatial.realign.write.roptions.prefix    = 'r';
    spm_jobman('run',matlabbatch);
    
    
    % Test spm_reslice on low quality
    xASL_spm_reslice(fullfile(TestDirDest,'ASL4D.nii'), fullfile(TestDirDest,'T1.nii'), [], [], 0);
    
    % Test spm_smooth on low quality (small kernel)
    xASL_spm_smooth(fullfile(TestDirDest,'rT1.nii'), [2 2 2], fullfile(TestDirDest,'rsT1.nii'));
    
    % Test coreg on low quality
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fullfile(TestDirDest,'ASL4D.nii,1')};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {fullfile(TestDirDest,'rT1.nii,1')};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'mi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = 9;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm_jobman('run',matlabbatch);
    
    % Test CAT12
    matlabbatch = [];
    SPMTemplateNII    = fullfile(x.opts.MyPath,'External','SPMmodified', 'tpm', 'TPM.nii');
    [~,catVer] = cat_version();
    if str2double(catVer) > 1500
        catTempDir = 'templates_volumes';
    else
        catTempDir = 'templates_1.50mm';
    end
    matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm              = {SPMTemplateNII};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP           = 0; % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr        = 0; % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr       = 0; % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox           = 6; % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr          = eps; % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.opts.samp             = 9;   % low quality
    matlabbatch{1}.spm.tools.cat.estwrite.data                  = {fullfile(TestDirDest,'T1.nii')}; % T1.nii
    matlabbatch{1}.spm.tools.cat.estwrite.nproc                 = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg           = 'mni'; % regularize affine registration for MNI European brains
    matlabbatch{1}.spm.tools.cat.estwrite.output.surface        = 0;   % don't do surface modeling
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.noROI  = struct([]); % don't do ROI estimations
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native      = 0;   % save c1T1 in native space
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod         = 0;   % don't save modulation
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel      = 0;   % don't save DARTEL space c1T1, this happens below in the reslice part
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native      = 0;   % save c2T1 in native space
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod         = 0;   % don't save modulation
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel      = 0;   % don't save DARTEL space c2T1, this happens below in the reslice part
    matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.warps          = [1 0]; % save warp to MNI
    matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped    = 0;   % don't save bias-corrected T1.nii
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed= [1 0.1]; % process everything on 1 mm fixed resolution (default)
    spm_jobman('run',matlabbatch);
    
    % Test deformations after CAT12
    xASL_Copy(fullfile(TestDirDest,'mri','y_T1.nii'), fullfile(TestDirDest,'y_T1.nii'));
    matlabbatch = [];
    matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(TestDirDest,'y_T1.nii')};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(TestDirDest,'rT1.nii')};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
    spm_jobman('run',matlabbatch);
    
    % Remove temporary derivatives
    xASL_adm_DeleteFileList(TestDirDest, '.*', false, [0 Inf]);
    DirsAre = {'err' 'ASL_1' 'mri' 'report'};
    for iDir=1:length(DirsAre)
        DirIs = fullfile(TestDirDest,DirsAre{iDir});
        xASL_adm_DeleteFileList(DirIs, '.*', true, [0 Inf]);
        xASL_delete(DirIs);
    end    
end