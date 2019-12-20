mainPath = '/pet/projekte/asl/data/ExploreASLreproAll';

studyName = {'EPAD','Novice','Sleep'};

versionName = {'NEW','OLD'};

matlabName = {'Lin18b','Win18b','Win16b','Win15a','WinLin16b','Mac14b','Mac16b'};
pathXasl = '/home/janpetr/code/ExploreASL';
% Mat1, try1, ver1 [, Mat2, try2, ver2]
processVec = [1,1,1;
	          1,2,1;
	          2,1,1;
			  1,1,2
			  2,1,2;
			  1,2,2;
		      3,1,2;
		      3,1,1;
		      4,1,2;
		      4,1,1;
			  5,1,1;
		      5,1,2;
			  7,1,1;
			  6,1,1;
		      7,1,2;
			  6,1,2];


compareVec = [1,1,1,1,2,1; % compare twice the same on linux NEW 18b
	          1,1,1,2,1,1; % Compare linux and windows 18b on NEW
			  1,1,2,2,1,2; % Compare linux and windows 18b on OLD
			  1,1,2,1,2,2;% compare twice the same on linux OLD 18b
		      2,1,2,3,1,2; % Compare windows 18bvs 16b on OLD
		      2,1,1,3,1,1; % Compare windows 18bvs 16b on NEW
		      2,1,2,4,1,2; % Compare windows 18bvs 15a on OLD
		      2,1,1,4,1,1; % Compare windows 18bvs 15a on NEW
			  5,1,1,2,1,1; % Compare WinLin16b with Win18b NEW
		      5,1,2,2,1,2; % Compare WinLin16b with Win18b OLD
			  5,1,1,1,1,1; % Compare WinLin16b with Lin18b NEW
			  5,1,2,1,1,2; % Compare WinLin16b with Lin18b OLD
			  7,1,1,1,1,1; % Compare Mac16b with Lin18b NEW
			  6,1,1,1,1,1; % Compare Mac14b with Lin18b NEW
		      7,1,2,1,1,2; % Compare Mac16b with Lin18b OLD
			  6,1,2,1,1,2; % Compare Mac14b with Lin18b OLD
		      7,1,1,2,1,1; % Compare Mac16b with Win18b NEW
		      7,1,2,2,1,2; % Compare Mac16b with Win18b OLD
			  4,1,1,1,1,1; % Compare Win15a with Lin18b NEW
		      4,1,2,1,1,2]; % Compare Win15a with Lin18b OLD

locDir = pwd;

% Calculate PVC
for iS = 1:length(studyName)
	for iC = 1:size(processVec,1)
		nameRef = [studyName{iS} '_' matlabName{processVec(iC,1)} '_' versionName{processVec(iC,3)} '_' num2str(processVec(iC,2))];

		% Downsample the last GM and WM to ASL space

		% Presmooth
		xASL_im_PreSmooth(fullfile(mainPath,nameRef,'Sub-001','ASL_1','CBF.nii'),fullfile(mainPath,nameRef,'Sub-001','c1T1.nii'),fullfile(mainPath,nameRef,'Sub-001','smc1T1.nii'));
		xASL_im_PreSmooth(fullfile(mainPath,nameRef,'Sub-001','ASL_1','CBF.nii'),fullfile(mainPath,nameRef,'Sub-001','c2T1.nii'),fullfile(mainPath,nameRef,'Sub-001','smc2T1.nii'));

		% Interpolate
		xASL_spm_reslice(fullfile(mainPath,nameRef,'Sub-001','ASL_1','CBF.nii'), fullfile(mainPath,nameRef,'Sub-001','smc1T1.nii'),[],[],1,fullfile(mainPath,nameRef,'Sub-001','ASL_1','pGM.nii'));
		xASL_spm_reslice(fullfile(mainPath,nameRef,'Sub-001','ASL_1','CBF.nii'), fullfile(mainPath,nameRef,'Sub-001','smc2T1.nii'),[],[],1,fullfile(mainPath,nameRef,'Sub-001','ASL_1','pWM.nii'));

		% Delete the pre-smoothed images
		xASL_delete(fullfile(mainPath,nameRef,'Sub-001','smc1T1.nii'));
		xASL_delete(fullfile(mainPath,nameRef,'Sub-001','smc2T1.nii'));

		% Load GM, WM, CBF
		pGM = xASL_io_Nifti2Im(fullfile(mainPath,nameRef,'Sub-001','ASL_1','pGM.nii'));
		pWM = xASL_io_Nifti2Im(fullfile(mainPath,nameRef,'Sub-001','ASL_1','pWM.nii'));
		imPV = pGM;
		imPV(:,:,:,2) = pWM;

		imCBF = xASL_io_Nifti2Im(fullfile(mainPath,nameRef,'Sub-001','ASL_1','CBF.nii'));

		% Run PVEC
		[imPVEC,imCBFrec,imResidual] = xASL_im_PVCkernel(double(imCBF), double(imPV),[5,5,1],'asllani');

		% Save PVEC results
		xASL_io_SaveNifti(fullfile(mainPath,nameRef,'Sub-001','ASL_1','CBF.nii'),fullfile(mainPath,nameRef,'Sub-001','ASL_1','CBF-GM.nii'),imPVEC(:,:,:,1));
		xASL_io_SaveNifti(fullfile(mainPath,nameRef,'Sub-001','ASL_1','CBF.nii'),fullfile(mainPath,nameRef,'Sub-001','ASL_1','CBF-WM.nii'),imPVEC(:,:,:,2));
	end
end


% Run the PVC comparison
for iS = 1:length(studyName)
	for iC = 1:size(compareVec,1)
		nameRef = [studyName{iS} '_' matlabName{compareVec(iC,1)} '_' versionName{compareVec(iC,3)} '_' num2str(compareVec(iC,2))];
		nameSrc = [studyName{iS} '_' matlabName{compareVec(iC,4)} '_' versionName{compareVec(iC,6)} '_' num2str(compareVec(iC,5))];

		% Load the GM, WM, GMCBF, WMCBF
		pGMref = xASL_io_Nifti2Im(fullfile(mainPath,nameRef,'Sub-001','ASL_1','pGM.nii'));
		pGMsrc = xASL_io_Nifti2Im(fullfile(mainPath,nameSrc,'Sub-001','ASL_1','pGM.nii'));
		imGMCBFref = xASL_io_Nifti2Im(fullfile(mainPath,nameRef,'Sub-001','ASL_1','CBF-GM.nii'));
		imGMCBFsrc = xASL_io_Nifti2Im(fullfile(mainPath,nameSrc,'Sub-001','ASL_1','CBF-GM.nii'));


		% Calculate 70% GM Mask
		imMask = (pGMref > 0.7).*(pGMsrc > 0.7);

		minV = -50;-100;
		maxV = 300;400;
		imGMCBFref(imGMCBFref<minV) = minV;
		imGMCBFref(imGMCBFref> maxV) = maxV;
		imGMCBFsrc(imGMCBFsrc<minV) = minV;
		imGMCBFsrc(imGMCBFsrc>maxV) = maxV;

		% Calculate mean GM on that and mean voxel-wise difference
		meanGMCBF = sum(imMask(:).*imGMCBFref(:))/sum(imMask(:));
		meanDifGMCBF = sum(imMask(:).*abs(imGMCBFref(:)-imGMCBFsrc(:)))/sum(imMask(:));

		% Return the percentage of that
		display([nameRef '-' nameSrc ' diff ' num2str(meanDifGMCBF/meanGMCBF*100,',%.2f')]);

	end
end
