pathBasic = '/pet/projekte/asl/data/Craniosynostosis/';
pathTPM = '/pet/projekte/asl/data/Craniosynostosis/UNCInfant012Atlas_20140325';

pathTest{1} = 'infant-1yr';
pathTest{2} = 'infant-2yr';
pathTest{3} = 'infant-neo';

for ii=1:length(pathTest)

	% Delete old directories
	xASL_delete(fullfile(pathBasic,pathTest{ii}));

	% Copy new results
	xASL_Copy(fullfile(pathBasic,'orig'),fullfile(pathBasic,pathTest{ii}));

	% Get all the files there
	list = xASL_adm_GetFileList(fullfile(pathBasic,pathTest{ii}),'^.+\T1.nii');
	
	for jj = 1:length(list)
		matlabbatch = [];
		matlabbatch{1}.spm.spatial.preproc.channel.vols = {list{jj}};
		matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
		matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
		matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(pathTPM,[pathTest{ii} '-TPM.nii,1'])};
		matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
		matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(pathTPM,[pathTest{ii} '-TPM.nii,2'])};
		matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
		matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(pathTPM,[pathTest{ii} '-TPM.nii,3'])};
		matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
		matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(pathTPM,[pathTest{ii} '-TPM.nii,4'])};
		matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
		matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(pathTPM,[pathTest{ii} '-TPM.nii,5'])};
		matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
		matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(pathTPM,[pathTest{ii} '-TPM.nii,6'])};
		matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
		matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
		matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
		matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
		matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
		matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
		matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
		matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
		matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
		matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
		spm_jobman('run',matlabbatch);
	end
end