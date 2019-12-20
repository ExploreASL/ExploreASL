pathTPM = '/pet/projekte/asl/data/Craniosynostosis/UNCInfant012Atlas_20140325';

nameTPM{1} = 'infant-1yr';
nameTPM{2} = 'infant-2yr';
nameTPM{3} = 'infant-neo';

for ii=1:length(nameTPM)
	imAll = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '.nii']));
	imWithSkull = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '-withSkull.nii']));
	imWithCerebellum = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '-withCerebellum.nii']));
	
	imWM = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '-seg-wm.nii']));
	imGM = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '-seg-gm.nii']));
	imCSF = xASL_io_Nifti2Im(fullfile(pathTPM,[nameTPM{ii} '-seg-csf.nii']));
	
	imWM = imWM/150;
	imGM = imGM/150;
	imCSF = imCSF/10;
	
	imAir = (imWithSkull < 25).*(imWM==0).*(imCSF==0).*(imGM==0);
	imCavitySkull = (imWithCerebellum < 25).*(imAir==0).*(imWM==0).*(imCSF==0).*(imGM==0);
	imCavity = (imWithSkull<120) .* imCavitySkull;
	imSkull = (imWithSkull>=120) .* imCavitySkull;
	
	imTPM = zeros([size(imAll), 6]);
	imTPM(:,:,:,1) = imGM;
	imTPM(:,:,:,2) = imWM;
	imTPM(:,:,:,3) = imCSF;
	imTPM(:,:,:,4) = imCavity;
	imTPM(:,:,:,5) = imSkull;
	imTPM(:,:,:,6) = imAir;
	
	xASL_io_SaveNifti(fullfile(pathTPM,[nameTPM{ii} '.nii']),fullfile(pathTPM,[nameTPM{ii} '-TPM.nii']),imTPM,16,0);
end
% TPM 1GM,2WM,3CSF,4cavity,5skull,6air