%% Init
rawDir = '/home/janpetr/tmp/PICTUREresegmented';

%% Reslice the second to the first
patientNameList = xASL_adm_GetFileList(fullfile(rawDir), '^P.*$', 'List', [], false);

for ii = 1:(length(patientNameList)/2)
	matlabbatch = [];
	fileA = patientNameList{2*ii-1};
	fileB = patientNameList{2*ii};
	matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(rawDir,fileA)};
	matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(rawDir,fileB)};
	matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
	matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
	matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
	matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';
	spm_jobman('run',matlabbatch);
end

%% Calculate the coefficients


patientNameList = xASL_adm_GetFileList(fullfile(rawDir), '^P.*$', 'List', [], false);

for ii = 1:(length(patientNameList)/2)
	fileA = patientNameList{2*ii-1};
	fileB = patientNameList{2*ii};
	imA = xASL_io_Nifti2Im(fullfile(rawDir,fileA));
	imB = xASL_io_Nifti2Im(fullfile(rawDir,['w' fileB]));
	resDice(ii) = xASL_im_ComputeDice(imA,imB);
	[resHaus(ii),resHausMod(ii)] = xASL_im_HausdorffDist(imA,imB);
	imABcross = (imA.*imB)>0;
	imABunit  = (imA+imB)>0;
	resGCI(ii) = sum(imABcross(:))./sum(imABunit(:));
	% Dice, Haussdorf, modified Haussdorf, GCI.
	fprintf('%30s : %30s : %6.3f : %6.3f : %6.3f: %6.3f\n',fileA,fileB,resDice(ii),resHaus(ii),resHausMod(ii),resGCI(ii));
	
end