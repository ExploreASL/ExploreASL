% Copyright 2015-2025 ExploreASL (Works In Progress code)
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

%% Please add here the path to the atlas files and to ExploreASL

pathTPM = '/home/janpetr/ExploreASL/Atlas';
pathTPMAtlasDir = 'mni_icbm152_nlin_asym_09c';
pathTPMAtlasFile = 'mni_icbm152_gm_tal_nlin_asym_09c.nii';
pathTPMAtlasFileWM = 'mni_icbm152_wm_tal_nlin_asym_09c.nii';
pathTPMAtlasFileGMWM = 'mni_icbm152_gmwm_tal_nlin_asym_09c.nii';

% Do this also for the symetrical one
pathTPMSymAtlasDir = 'mni_icbm152_nlin_sym_09c';
pathTPMSymAtlasFile = 'mni_icbm152_gm_tal_nlin_sym_09c.nii';
pathTPMSymAtlasFileWM = 'mni_icbm152_wm_tal_nlin_sym_09c.nii';
pathTPMSymAtlasFileGMWM = 'mni_icbm152_gmwm_tal_nlin_sym_09c.nii';

pathColinDir = 'mni_colin27_2008';
pathColinFileCLS = 'colin27_cls_tal_hires.nii';
pathColinFileGM = 'colin27_gm.nii';
pathColinFileWM = 'colin27_wm.nii';
pathColinFileGMWM = 'colin27_gmwm.nii';

pathFSL152Dir = 'fsl_mni152';
pathFSL152File = 'FSL_MNI152_FreeSurferConformed_1mm.nii';

pathSchaefer = 'Schaefer';

pathExploreASL = '/home/janpetr/ExploreASL/ExploreASL';
pathSPM = fullfile(pathExploreASL,'External', 'SPMmodified');
pathCATtemplates = fullfile(pathSPM, 'toolbox', 'cat12', 'templates_volumes');

%% Load the ICMB-152 atlas 
% We used this one http://www.bic.mni.mcgill.ca/~vfonov/icbm/2009/mni_icbm152_nlin_asym_09a_nifti.zip
% This is the reference template that we assume the atlases are registered to

% We create a combined image
imGM = xASL_io_Nifti2Im(fullfile(pathTPM, pathTPMAtlasDir, pathTPMAtlasFile));
imWM = xASL_io_Nifti2Im(fullfile(pathTPM, pathTPMAtlasDir, pathTPMAtlasFileWM));
imCombined = imGM;
imCombined(:,:,:,2) = imWM;

xASL_io_SaveNifti(fullfile(pathTPM,pathTPMAtlasDir, pathTPMAtlasFile),fullfile(pathTPM,pathTPMAtlasDir, pathTPMAtlasFileGMWM),imCombined);

% We create a combined image for the symetrical atlas as well
imGM = xASL_io_Nifti2Im(fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFile));
imWM = xASL_io_Nifti2Im(fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFileWM));
imCombined = imGM;
imCombined(:,:,:,2) = imWM;

xASL_io_SaveNifti(fullfile(pathTPM,pathTPMSymAtlasDir, pathTPMSymAtlasFile),fullfile(pathTPM,pathTPMSymAtlasDir, pathTPMSymAtlasFileGMWM),imCombined);


% Do the same for the Colin-27 atlas
% We create a combined image
imGM = xASL_io_Nifti2Im(fullfile(pathTPM, pathColinDir, pathColinFileCLS));
imGM = round(imGM);
imCombined = imGM == 2;
imCombined(:,:,:,2) = imGM == 3;

xASL_io_SaveNifti(fullfile(pathTPM,pathColinDir, pathColinFileCLS),fullfile(pathTPM,pathColinDir, pathColinFileGM),imCombined(:,:,:,1));
xASL_io_SaveNifti(fullfile(pathTPM,pathColinDir, pathColinFileCLS),fullfile(pathTPM,pathColinDir, pathColinFileWM),imCombined(:,:,:,2));
xASL_io_SaveNifti(fullfile(pathTPM,pathColinDir, pathColinFileCLS),fullfile(pathTPM,pathColinDir, pathColinFileGMWM),imCombined);

%% Presmooth the atlas and move it to the IXI resolution
xASL_im_PreSmooth(fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMAtlasDir, pathTPMAtlasFileGMWM));
xASL_im_PreSmooth(fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMAtlasDir, pathTPMAtlasFile));
xASL_im_PreSmooth(fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMAtlasDir, pathTPMAtlasFileWM));

xASL_im_PreSmooth(fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFileGMWM));
xASL_im_PreSmooth(fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFile));
xASL_im_PreSmooth(fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFileWM));

xASL_im_PreSmooth(fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathColinDir, pathColinFileGM));
xASL_im_PreSmooth(fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathColinDir, pathColinFileWM));
xASL_im_PreSmooth(fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathColinDir, pathColinFileGMWM));

% And resample to 1.5mm DARTEL space
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM, pathTPMAtlasDir, ['s' pathTPMAtlasFileGMWM])
	fullfile(pathTPM, pathTPMAtlasDir, ['s' pathTPMAtlasFile])
	fullfile(pathTPM, pathTPMAtlasDir, ['s' pathTPMAtlasFileWM])
	fullfile(pathTPM, pathTPMSymAtlasDir, ['s' pathTPMSymAtlasFileGMWM])
	fullfile(pathTPM, pathTPMSymAtlasDir, ['s' pathTPMSymAtlasFile])
	fullfile(pathTPM, pathTPMSymAtlasDir, ['s' pathTPMSymAtlasFileWM])
	fullfile(pathTPM, pathColinDir, ['s' pathColinFileGM])
	fullfile(pathTPM, pathColinDir, ['s' pathColinFileWM])
	fullfile(pathTPM, pathColinDir, ['s' pathColinFileGMWM])};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
spm_jobman('run',matlabbatch);
%% Align the reference MNI2009c template to the IXI512 template

matlabbatch = [];
matlabbatch{1}.spm.tools.dartel.warp1.images = {
                                                {fullfile(pathTPM, pathTPMAtlasDir, ['ws' pathTPMAtlasFile])}
												{fullfile(pathTPM, pathTPMAtlasDir, ['ws' pathTPMAtlasFileWM])}
                                                }';
matlabbatch{1}.spm.tools.dartel.warp1.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).template = {fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).template = {fullfile(pathCATtemplates ,'Template_2_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).template = {fullfile(pathCATtemplates ,'Template_3_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).template = {fullfile(pathCATtemplates ,'Template_4_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).template = {fullfile(pathCATtemplates ,'Template_5_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).template = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;

spm_jobman('run',matlabbatch);
%% Transform the template to the correct template, just as a test
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {fullfile(pathTPM, pathTPMAtlasDir, ['u_ws' pathTPMAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM, pathTPMAtlasDir, ['ws' pathTPMAtlasFile])};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {pathTPM};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'o';

spm_jobman('run',matlabbatch);
%% Align the reference MNI2009c template to the IXI512 template

matlabbatch = [];
matlabbatch{1}.spm.tools.dartel.warp1.images = {
                                                {fullfile(pathTPM, pathTPMSymAtlasDir, ['ws' pathTPMSymAtlasFile])}
												{fullfile(pathTPM, pathTPMSymAtlasDir, ['ws' pathTPMSymAtlasFileWM])}
                                                }';
matlabbatch{1}.spm.tools.dartel.warp1.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).template = {fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).template = {fullfile(pathCATtemplates ,'Template_2_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).template = {fullfile(pathCATtemplates ,'Template_3_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).template = {fullfile(pathCATtemplates ,'Template_4_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).template = {fullfile(pathCATtemplates ,'Template_5_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).template = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;

spm_jobman('run',matlabbatch);
%% Transform the template to the correct template, just as a test
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {fullfile(pathTPM, pathTPMSymAtlasDir, ['u_ws' pathTPMSymAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM, pathTPMSymAtlasDir, ['ws' pathTPMSymAtlasFile])};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {pathTPM};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'o';

spm_jobman('run',matlabbatch);
%% Align the reference Colin27 template to the IXI512 template

matlabbatch = [];
matlabbatch{1}.spm.tools.dartel.warp1.images = {
                                                {fullfile(pathTPM, pathColinDir, ['ws' pathColinFileGM])}
												{fullfile(pathTPM, pathColinDir, ['ws' pathColinFileWM])}
                                                }';
matlabbatch{1}.spm.tools.dartel.warp1.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).template = {fullfile(pathCATtemplates ,'Template_1_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).template = {fullfile(pathCATtemplates ,'Template_2_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).template = {fullfile(pathCATtemplates ,'Template_3_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).template = {fullfile(pathCATtemplates ,'Template_4_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).template = {fullfile(pathCATtemplates ,'Template_5_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).template = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;

spm_jobman('run',matlabbatch);
%% Transform the Colin27 template to the correct template, just as a test
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {fullfile(pathTPM, pathColinDir, ['u_ws' pathColinFileGM])};
matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM, pathColinDir, ['ws' pathColinFileGM])};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {pathTPM};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'o';

spm_jobman('run',matlabbatch);
%% Transform the Brainstem atlas from sym-MNI2009c to IXI512
imAtlasOrig = xASL_io_Nifti2Im(fullfile(pathTPM, 'FreeSurferSubfields', 'BrainstemProbs.MNIsymSpace.nii.gz')); % Load atlas
[~,imAtlasMax] = max(imAtlasOrig, [], 4); %Convert probabilities to labels

imAtlas = zeros(size(imAtlasMax(:,:,:,1))); % Create an empty atlas
TSV = {};
imAtlas(imAtlasMax==2) = 1;TSV{1,1} = 'Medulla';
imAtlas(imAtlasMax==3) = 2;TSV{2,1} = 'Pons';
imAtlas(imAtlasMax==4) = 3;TSV{3,1} = 'SCP';
imAtlas(imAtlasMax==5) = 4;TSV{4,1} = 'Midbrain';
xASL_io_SaveNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'BrainstemProbs.MNIsymSpace.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferBrainstem.nii'), imAtlas);

% Nearest neighbor transformation
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {fullfile(pathTPM, pathTPMSymAtlasDir, ['u_ws' pathTPMSymAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFile)
												   fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferBrainstem.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(pathTPM, 'FreeSurferSubfields')};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'o';

spm_jobman('run',matlabbatch);

% Create MAT file
xASL_tsvWrite(TSV, fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferBrainstem.tsv'), 1);
IM = xASL_io_Nifti2Im(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferBrainstem.nii'));
IM = uint8(IM);
save(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferBrainstem.mat'),'IM');
xASL_adm_GzipNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferBrainstem.nii'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferBrainstem.nii'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferBrainstem.nii.mat'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferBrainstem.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferBrainstem.nii.gz'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferBrainstem.mat'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferBrainstem.nii.mat'));
%% Transform the Brainstem atlas from sym-MNI2009c to IXI512
imAtlasOrig = xASL_io_Nifti2Im(fullfile(pathTPM, 'FreeSurferSubfields', 'ThalamusProbs.MNIsymSpace.nii.gz')); % Load atlas
[~,imAtlasMax] = max(imAtlasOrig, [], 4); %Convert probabilities to labels

imAtlas = zeros(size(imAtlasMax(:,:,:,1))); % Create an empty atlas
TSV = {};
imAtlas(imAtlasMax==2) = 1;TSV{1,1} = 'LGN';
imAtlas(imAtlasMax==3) = 1;
imAtlas(imAtlasMax==4) = 2;TSV{2,1} = 'MGN';
imAtlas(imAtlasMax==5) = 2;
imAtlas(imAtlasMax==6) = 3;TSV{3,1} = 'PuI';
imAtlas(imAtlasMax==29) = 3;
imAtlas(imAtlasMax==7) = 4;TSV{4,1} = 'PuM';
imAtlas(imAtlasMax==30) = 4;
imAtlas(imAtlasMax==8) = 5;TSV{5,1} = 'L-Sg';
imAtlas(imAtlasMax==31) = 5;
imAtlas(imAtlasMax==9) = 6;TSV{6,1} = 'VPL';
imAtlas(imAtlasMax==32) = 6;
imAtlas(imAtlasMax==10) = 7;TSV{7,1} = 'CM';
imAtlas(imAtlasMax==33) = 7;
imAtlas(imAtlasMax==11) = 8;TSV{8,1} = 'VLa';
imAtlas(imAtlasMax==34) = 8;
imAtlas(imAtlasMax==12) = 9;TSV{9,1} = 'PuA';
imAtlas(imAtlasMax==35) = 9;
imAtlas(imAtlasMax==13) = 10;TSV{10,1} = 'MDm';
imAtlas(imAtlasMax==36) = 10;
imAtlas(imAtlasMax==14) = 11;TSV{11,1} = 'Pf';
imAtlas(imAtlasMax==37) = 11;
imAtlas(imAtlasMax==15) = 12;TSV{12,1} = 'VAmc';
imAtlas(imAtlasMax==38) = 12;
imAtlas(imAtlasMax==16) = 13;TSV{13,1} = 'MD1';
imAtlas(imAtlasMax==39) = 13;
imAtlas(imAtlasMax==17) = 14;TSV{14,1} = 'CeM';
imAtlas(imAtlasMax==40) = 14;
imAtlas(imAtlasMax==18) = 15;TSV{15,1} = 'VA';
imAtlas(imAtlasMax==41) = 15;
imAtlas(imAtlasMax==19) = 16;TSV{16,1} = 'MV(Re)';
imAtlas(imAtlasMax==42) = 16;
imAtlas(imAtlasMax==20) = 17;TSV{17,1} = 'R';
imAtlas(imAtlasMax==43) = 17;
imAtlas(imAtlasMax==21) = 18;TSV{18,1} = 'VM';
imAtlas(imAtlasMax==44) = 18;
imAtlas(imAtlasMax==22) = 19;TSV{19,1} = 'CL';
imAtlas(imAtlasMax==46) = 19;
imAtlas(imAtlasMax==23) = 20;TSV{20,1} = 'PuL';
imAtlas(imAtlasMax==45) = 20;
imAtlas(imAtlasMax==24) = 21;TSV{21,1} = 'Pt';
imAtlas(imAtlasMax==49) = 21;
imAtlas(imAtlasMax==25) = 22;TSV{22,1} = 'AV';
imAtlas(imAtlasMax==50) = 22;
imAtlas(imAtlasMax==26) = 23;TSV{23,1} = 'Pc';
imAtlas(imAtlasMax==48) = 23;
imAtlas(imAtlasMax==27) = 24;TSV{24,1} = 'VLp';
imAtlas(imAtlasMax==47) = 24;
imAtlas(imAtlasMax==28) = 25;TSV{25,1} = 'LP';
imAtlas(imAtlasMax==51) = 25;
imAtlas(imAtlasMax==52) = 26;TSV{26,1} = 'LD';
imAtlas(imAtlasMax==53) = 26;

% Merged Histology-based
imAtlasLargeHist = zeros(size(imAtlasMax(:,:,:,1))); % Create an empty atlas with merged regions histology based
TSVLargeHist = {};
imAtlasLargeHist(imAtlasMax==25) = 1;TSVLargeHist{1,1} = 'Anterior Thalamus';%AV
imAtlasLargeHist(imAtlasMax==50) = 1;
imAtlasLargeHist(imAtlasMax==24) = 1;%Pt
imAtlasLargeHist(imAtlasMax==49) = 1;
imAtlasLargeHist(imAtlasMax==16) = 2;TSVLargeHist{2,1} = 'Mediodorsal Thalamus';%MD1
imAtlasLargeHist(imAtlasMax==39) = 2;
imAtlasLargeHist(imAtlasMax==13) = 2;%MDm
imAtlasLargeHist(imAtlasMax==36) = 2;
imAtlasLargeHist(imAtlasMax==18) = 3;TSVLargeHist{3,1} = 'Ventral Motor Complex';%VA
imAtlasLargeHist(imAtlasMax==41) = 3;
imAtlasLargeHist(imAtlasMax==15) = 3;%VAmc
imAtlasLargeHist(imAtlasMax==38) = 3;
imAtlasLargeHist(imAtlasMax==11) = 3;%VLa
imAtlasLargeHist(imAtlasMax==34) = 3;
imAtlasLargeHist(imAtlasMax==21) = 3;%VM
imAtlasLargeHist(imAtlasMax==44) = 3;
imAtlasLargeHist(imAtlasMax==27) = 3;%VLp
imAtlasLargeHist(imAtlasMax==47) = 3;
imAtlasLargeHist(imAtlasMax==9) = 4;TSVLargeHist{4,1} = 'Ventral Sensory Complex';%VPL
imAtlasLargeHist(imAtlasMax==32) = 4;
imAtlasLargeHist(imAtlasMax==17) = 5;TSVLargeHist{5,1} = 'Intralaminar midline Complex';%CeM
imAtlasLargeHist(imAtlasMax==40) = 5;
imAtlasLargeHist(imAtlasMax==10) = 5;%CM
imAtlasLargeHist(imAtlasMax==33) = 5;
imAtlasLargeHist(imAtlasMax==22) = 5;%CL
imAtlasLargeHist(imAtlasMax==46) = 5;
imAtlasLargeHist(imAtlasMax==14) = 5;%Pf
imAtlasLargeHist(imAtlasMax==37) = 5;
imAtlasLargeHist(imAtlasMax==19) = 5;%MV(Re)
imAtlasLargeHist(imAtlasMax==42) = 5;
imAtlasLargeHist(imAtlasMax==26) = 5;%Pc
imAtlasLargeHist(imAtlasMax==48) = 5;
imAtlasLargeHist(imAtlasMax==12) = 6;TSVLargeHist{6,1} = 'Pulvinar Complex';%PuA
imAtlasLargeHist(imAtlasMax==35) = 6;
imAtlasLargeHist(imAtlasMax==6) = 6;%PuI
imAtlasLargeHist(imAtlasMax==29) = 6;
imAtlasLargeHist(imAtlasMax==7) = 6;%PuM
imAtlasLargeHist(imAtlasMax==30) = 6;
imAtlasLargeHist(imAtlasMax==23) = 6;%PuL
imAtlasLargeHist(imAtlasMax==45) = 6;
imAtlasLargeHist(imAtlasMax==2) = 7;TSVLargeHist{7,1} = 'Posterior Sensory';%LGN
imAtlasLargeHist(imAtlasMax==3) = 7;
imAtlasLargeHist(imAtlasMax==4) = 7;%MGN
imAtlasLargeHist(imAtlasMax==5) = 7;
imAtlasLargeHist(imAtlasMax==8) = 7;%L-Sg
imAtlasLargeHist(imAtlasMax==31) = 7;

%imAtlasLargeHist(imAtlasMax==20) = 17;TSV{17,1} = 'R';
%imAtlasLargeHist(imAtlasMax==43) = 17;
%imAtlasLargeHist(imAtlasMax==28) = 25;TSV{25,1} = 'LP';
%imAtlasLargeHist(imAtlasMax==51) = 25;
%imAtlasLargeHist(imAtlasMax==52) = 26;TSV{26,1} = 'LD';
%imAtlasLargeHist(imAtlasMax==53) = 26;

% Merged Functionally-based
imAtlasLargeFunc = zeros(size(imAtlasMax(:,:,:,1))); % Create an empty atlas with merged regions functionally based
TSVLargeFunc = {};
imAtlasLargeFunc(imAtlasMax==25) = 1;TSVLargeFunc{1,1} = 'Limbic Prefrontal Thalamus';%AV
imAtlasLargeFunc(imAtlasMax==50) = 1;
imAtlasLargeFunc(imAtlasMax==13) = 1;%MDm
imAtlasLargeFunc(imAtlasMax==36) = 1;
imAtlasLargeFunc(imAtlasMax==17) = 1;%CeM
imAtlasLargeFunc(imAtlasMax==40) = 1;
imAtlasLargeFunc(imAtlasMax==16) = 1;%MD1
imAtlasLargeFunc(imAtlasMax==39) = 1;
imAtlasLargeFunc(imAtlasMax==19) = 1;%MV(Re)
imAtlasLargeFunc(imAtlasMax==42) = 1;
imAtlasLargeFunc(imAtlasMax==24) = 1;%Pt
imAtlasLargeFunc(imAtlasMax==49) = 1;
imAtlasLargeFunc(imAtlasMax==26) = 1;%Pc
imAtlasLargeFunc(imAtlasMax==48) = 1;
imAtlasLargeFunc(imAtlasMax==18) = 2;TSVLargeFunc{2,1} = 'Sensorimotor Intralaminar Reticular Thalamus';%VA
imAtlasLargeFunc(imAtlasMax==41) = 2;
imAtlasLargeFunc(imAtlasMax==15) = 2;%VAmc
imAtlasLargeFunc(imAtlasMax==38) = 2;
imAtlasLargeFunc(imAtlasMax==11) = 2;%VLa
imAtlasLargeFunc(imAtlasMax==34) = 2;
imAtlasLargeFunc(imAtlasMax==27) = 2;%VLp
imAtlasLargeFunc(imAtlasMax==47) = 2;
imAtlasLargeFunc(imAtlasMax==9) = 2;%VPL
imAtlasLargeFunc(imAtlasMax==32) = 2;
imAtlasLargeFunc(imAtlasMax==21) = 2;%VM
imAtlasLargeFunc(imAtlasMax==44) = 2;
imAtlasLargeFunc(imAtlasMax==10) = 2;%CM
imAtlasLargeFunc(imAtlasMax==33) = 2;
imAtlasLargeFunc(imAtlasMax==22) = 2;%CL
imAtlasLargeFunc(imAtlasMax==46) = 2;
imAtlasLargeFunc(imAtlasMax==14) = 2;%Pf
imAtlasLargeFunc(imAtlasMax==37) = 2;
imAtlasLargeFunc(imAtlasMax==20) = 2;%R
imAtlasLargeFunc(imAtlasMax==43) = 2;
imAtlasLargeFunc(imAtlasMax==12) = 3;TSVLargeFunc{3,1} = 'Posterior Association Pulvinar Thalamus';%PuA
imAtlasLargeFunc(imAtlasMax==35) = 3;
imAtlasLargeFunc(imAtlasMax==7) = 3;%PuM
imAtlasLargeFunc(imAtlasMax==30) = 3;
imAtlasLargeFunc(imAtlasMax==23) = 3;%PuL
imAtlasLargeFunc(imAtlasMax==45) = 3;
imAtlasLargeFunc(imAtlasMax==6) = 3;%PuI
imAtlasLargeFunc(imAtlasMax==29) = 3;
imAtlasLargeFunc(imAtlasMax==52) = 3;%LD
imAtlasLargeFunc(imAtlasMax==53) = 3;
imAtlasLargeFunc(imAtlasMax==28) = 3;%LP
imAtlasLargeFunc(imAtlasMax==51) = 3;
imAtlasLargeFunc(imAtlasMax==8) = 3;%L-Sg
imAtlasLargeFunc(imAtlasMax==31) = 3;
imAtlasLargeFunc(imAtlasMax==2) = 4;TSVLargeFunc{4,1} = 'Primary Sensory Relay';%LGN
imAtlasLargeFunc(imAtlasMax==3) = 4;
imAtlasLargeFunc(imAtlasMax==4) = 4;%MGN
imAtlasLargeFunc(imAtlasMax==5) = 4;

xASL_io_SaveNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'ThalamusProbs.MNIsymSpace.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii'), imAtlas);
xASL_io_SaveNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'ThalamusProbs.MNIsymSpace.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHistological.nii'), imAtlasLargeHist);
xASL_io_SaveNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'ThalamusProbs.MNIsymSpace.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunctional.nii'), imAtlasLargeFunc);

% Nearest neighbor transformation
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {fullfile(pathTPM, pathTPMSymAtlasDir, ['u_ws' pathTPMSymAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFile)
												   fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii')
												   fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHistological.nii')
												   fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunctional.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(pathTPM, 'FreeSurferSubfields')};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'o';

spm_jobman('run',matlabbatch);

% Create MAT file
xASL_tsvWrite(TSV, fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.tsv'), 1);
xASL_tsvWrite(TSVLargeHist, fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHistological.tsv'), 1);
xASL_tsvWrite(TSVLargeFunc, fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunctional.tsv'), 1);
IM = xASL_io_Nifti2Im(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamus.nii'));
IMHist = xASL_io_Nifti2Im(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusHistological.nii'));
IMFunc = xASL_io_Nifti2Im(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusFunctional.nii'));
IM = uint8(IM);
IMHist = uint8(IMHist);
IMFunc = uint8(IMFunc);
save(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.mat'),'IM');
save(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHistological.mat'),'IMHist');
save(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunctional.mat'),'IMFunc');
xASL_adm_GzipNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamus.nii'));
xASL_adm_GzipNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusHistological.nii'));
xASL_adm_GzipNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusFunctional.nii'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHistological.nii'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunctional.nii'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii.mat'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHistological.nii.mat'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunctional.nii.mat'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamus.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii.gz'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusHistological.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHistological.nii.gz'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusFunctional.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunctional.nii.gz'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.mat'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii.mat'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHistological.mat'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHistological.nii.mat'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunctional.mat'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunctional.nii.mat'));

%% Transform the AAN atlas from sym-MNI2009c to IXI512
imAtlasOrig = xASL_io_Nifti2Im(fullfile(pathTPM, 'AAN', 'AAN_Brainstem_MNI152_1mm_v2p0.nii')); % Load atlas

imAtlas = zeros(size(imAtlasOrig(:,:,:,1))); % Create an empty atlas
TSV = {};
imAtlas(imAtlasOrig==7201) = 1;TSV{1,1} = 'DR';
imAtlas(imAtlasOrig==7202) = 2;TSV{2,1} = 'MnR';
imAtlas(imAtlasOrig==7203) = 3;TSV{3,1} = 'PAG';
imAtlas(imAtlasOrig==7204) = 4;TSV{4,1} = 'VTA';
imAtlas(imAtlasOrig==7205) = 5;TSV{5,1} = 'LC';
imAtlas(imAtlasOrig==7301) = 5;
imAtlas(imAtlasOrig==7401) = 5;
imAtlas(imAtlasOrig==7206) = 6;TSV{6,1} = 'LDTg';
imAtlas(imAtlasOrig==7302) = 6;
imAtlas(imAtlasOrig==7402) = 6;
imAtlas(imAtlasOrig==7207) = 7;TSV{7,1} = 'mRt';
imAtlas(imAtlasOrig==7303) = 7;
imAtlas(imAtlasOrig==7403) = 7;
imAtlas(imAtlasOrig==7208) = 8;TSV{8,1} = 'PBC';
imAtlas(imAtlasOrig==7304) = 8;
imAtlas(imAtlasOrig==7404) = 8;
imAtlas(imAtlasOrig==7209) = 9;TSV{9,1} = 'PnO';
imAtlas(imAtlasOrig==7305) = 9;
imAtlas(imAtlasOrig==7405) = 9;
imAtlas(imAtlasOrig==7210) = 10;TSV{10,1} = 'PTg';
imAtlas(imAtlasOrig==7306) = 10;
imAtlas(imAtlasOrig==7406) = 10;

xASL_io_SaveNifti(fullfile(pathTPM, 'AAN', 'AAN_Brainstem_MNI152_1mm_v2p0.nii'),fullfile(pathTPM, 'AAN', 'AAN_Brainstem.nii'), imAtlas);

% Nearest neighbor transformation
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {fullfile(pathTPM, pathTPMSymAtlasDir, ['u_ws' pathTPMSymAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFile)
												   fullfile(pathTPM, 'AAN', 'AAN_Brainstem.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(pathTPM, 'AAN')};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'o';

spm_jobman('run',matlabbatch);

% Create MAT file
xASL_tsvWrite(TSV, fullfile(pathTPM, 'AAN', 'AAN_Brainstem.tsv'), 1);
IM = xASL_io_Nifti2Im(fullfile(pathTPM, 'AAN', 'oAAN_Brainstem.nii'));
IM = uint8(IM);
save(fullfile(pathTPM, 'AAN', 'AAN_Brainstem.mat'),'IM');
xASL_adm_GzipNifti(fullfile(pathTPM, 'AAN', 'oAAN_Brainstem.nii'));
xASL_delete(fullfile(pathTPM, 'AAN', 'AAN_Brainstem.nii'));
xASL_delete(fullfile(pathTPM, 'AAN', 'AAN_Brainstem.nii.mat'));
xASL_Move(fullfile(pathTPM, 'AAN', 'oAAN_Brainstem.nii.gz'),fullfile(pathTPM, 'AAN', 'AAN_Brainstem.nii.gz'));
xASL_Move(fullfile(pathTPM, 'AAN', 'AAN_Brainstem.mat'),fullfile(pathTPM, 'AAN', 'AAN_Brainstem.nii.mat'));

%% Transform the atlases from MNI2009c to IXI512
% This has to be done label per label

% We load atlas and count all labels
imAtlas = xASL_io_Nifti2Im(fullfile(pathTPM, 'WMPM', 'WMPM_Type_III.nii'));
max(imAtlas(:))

% Nearest neighbor transformation

matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {fullfile(pathTPM, pathTPMAtlasDir, ['u_ws' pathTPMAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM, pathTPMAtlasDir, pathTPMAtlasFile)
												   fullfile(pathTPM, 'WMPM', 'WMPM_Type_III.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {pathTPM};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'o';

spm_jobman('run',matlabbatch);

% Create MAT file
xASL_adm_AtlasConvert_LeftRight2Bilateral(fullfile(pathTPM, 'oWMPM_Type_III.nii'), fullfile(pathTPM, 'WMPM', 'WMPM.txt'), fullfile(pathTPM, 'oWMPM_Type_III.tsv'), 'WMPM_Type_III');
%% Transform the atlases from Colin27 to IXI512
% This has to be done label per label

% We load atlas and count all labels
imAtlas = xASL_io_Nifti2Im(fullfile(pathTPM, 'AAL3', 'AAL3v1_1mm.nii'));
max(imAtlas(:))

% Nearest neighbor transformation

matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {fullfile(pathTPM, pathColinDir, ['u_ws' pathColinFileGM])};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {fullfile(pathCATtemplates ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM, 'AAL3', 'AAL3v1_1mm.nii')
	                                               fullfile(pathTPM, pathColinDir, pathColinFileGM)};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {pathTPM};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'o';

spm_jobman('run',matlabbatch);

% Create MAT file
xASL_Copy(fullfile(pathTPM, 'AAL3', 'AAL3v1.nii.txt'), fullfile(pathTPM, 'oAAL3v1_1mm.tsv'));
xASL_adm_AtlasConvert_LeftRight2Bilateral(fullfile(pathTPM, 'oAAL3v1_1mm.nii'), fullfile(pathTPM, 'oAAL3v1_1mm.tsv'), fullfile(pathTPM, 'oAAL3v1_1mm.tsv'), 'AAL3v1');


%IM = xASL_io_Nifti2Im(fullfile(pathTPM, 'oAAL3v1_1mm.nii'));
%IM = uint8(ceil(IM/2.0));
%save(fullfile(pathTPM, 'oAAL3v1_1mm.nii.mat'),'IM');
%xASL_io_SaveNifti(fullfile(pathTPM, 'oAAL3v1_1mm.nii'), fullfile(pathTPM, 'oAAL3v1_1mm.nii'), IM);
%xASL_adm_GzipNifti(fullfile(pathTPM, 'oAAL3v1_1mm.nii'));


%% Spatially normalize and segment the FSL_MNI152 atlas
% The Basic atlas is from here (note that it is the same dataset, but different registration than the standard MNI atlases like MNI2009c etc
% https://neurovault.org/images/14245/

% First align to IXI512 and segment

% Unzip the atlases
xASL_adm_UnzipNifti(fullfile(pathTPM, pathFSL152Dir, pathFSL152File));

% Set CAT12 template & registration settings
matlabbatch = [];
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm                         = {fullfile(pathSPM, 'tpm', 'TPM.nii')};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.xasl_savesteps           = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.xasl_quality             = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.xasl_disabledartel       = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {fullfile(pathCATtemplates, 'Template_0_IXI555_MNI152_GS.nii')}; % Runs Geodesic Shooting to this n=555 subjects template  %pathCATtemplates
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP           = 1070; % full cleanup. 1070 light cleanup
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr        = 0.5; % 0.5; % strength local adaptive segmentation
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr       = 2; % using SPM approach -> 0.5 GCUT may be more robust, to avoid stripping GM at brain poles
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox           = 1.5; % voxelsize on which registration is run (1.5 == default)
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr          = 0.5; % SPM bias-correction strength
matlabbatch{1}.spm.tools.cat.estwrite.opts.samp             = 3;   % spm sampling distance
% Add the FSL152MNI reference
matlabbatch{1}.spm.tools.cat.estwrite.data                  = {fullfile(pathTPM, pathFSL152Dir, pathFSL152File)}; % T1.nii
matlabbatch{1}.spm.tools.cat.estwrite.nproc                 = 0; % don't split the segmentation in multiple processes.
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg           = 'mni'; % regularize affine registration for MNI European brains
matlabbatch{1}.spm.tools.cat.estwrite.output.surface        = 0;   % don't do surface modeling
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native      = 1;   % save c1T1 in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod         = 0;   % don't save modulation
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel      = 0;   % don't save DARTEL space c1T1, this happens below in the reslice part
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native      = 1;   % save c2T1 in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod         = 0;   % don't save modulation
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel      = 0;   % don't save DARTEL space c2T1, this happens below in the reslice part
matlabbatch{1}.spm.tools.cat.estwrite.output.warps          = [1 0]; % save warp to MNI
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped    = 0;   % don't save bias-corrected T1.nii
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.noROI  = struct([]); % don't do ROI estimations
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed= [1 0.1]; % process everything on 1 mm fixed resolution (default)
% Run CAT12 segmentation
spm_jobman('run',matlabbatch); % Run CAT12

xASL_adm_UnzipNifti(fullfile(pathTPM, pathSchaefer, 'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz'));
xASL_adm_UnzipNifti(fullfile(pathTPM, pathSchaefer, 'Schaefer2018_100Parcels_17Networks_order_FSLMNI152_1mm.nii.gz'));

% Transform to IXI512 template
matlabbatch = [];
matlabbatch{1}.spm.tools.cat.tools.defs.field1 = {[fullfile(pathTPM, pathFSL152Dir, 'mri', 'y_FSL_MNI152_FreeSurferConformed_1mm.nii') ',1']};
matlabbatch{1}.spm.tools.cat.tools.defs.images = {[fullfile(pathTPM, pathSchaefer, 'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii') ',1']};
matlabbatch{1}.spm.tools.cat.tools.defs.interp = 0;
matlabbatch{1}.spm.tools.cat.tools.defs.modulate = 0;
spm_jobman('run',matlabbatch);
matlabbatch{1}.spm.tools.cat.tools.defs.images = {[fullfile(pathTPM, pathSchaefer, 'Schaefer2018_100Parcels_17Networks_order_FSLMNI152_1mm.nii') ',1']};
spm_jobman('run',matlabbatch);

% Convert atlas to a bilateral version for 7 Networks
imAtlasOrig = xASL_io_Nifti2Im(fullfile(pathTPM, pathSchaefer, 'wSchaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii'));
imAtlas = zeros(size(imAtlasOrig(:,:,:,1))); % Create an empty atlas
TSV = {};

% Relabel the atlas
imAtlas(imAtlasOrig==51) = 1;TSV{1,1} = '7Networks_Vis_1';
imAtlas(imAtlasOrig==52) = 2;TSV{2,1} = '7Networks_Vis_2';
imAtlas(imAtlasOrig==53) = 3;TSV{3,1} = '7Networks_Vis_3';
imAtlas(imAtlasOrig==54) = 4;TSV{4,1} = '7Networks_Vis_4';
imAtlas(imAtlasOrig==55) = 5;TSV{5,1} = '7Networks_Vis_5';
imAtlas(imAtlasOrig==56) = 6;TSV{6,1} = '7Networks_Vis_6';
imAtlas(imAtlasOrig==57) = 7;TSV{7,1} = '7Networks_Vis_7';
imAtlas(imAtlasOrig==58) = 8;TSV{8,1} = '7Networks_Vis_8';
TSV{9,1} = '7Networks_Vis_9';% RH_Vis_9 isn't defined
imAtlas(imAtlasOrig==59) = 10;TSV{10,1} = '7Networks_SomMot_1';
imAtlas(imAtlasOrig==60) = 11;TSV{11,1} = '7Networks_SomMot_2';
imAtlas(imAtlasOrig==61) = 12;TSV{12,1} = '7Networks_SomMot_3';
imAtlas(imAtlasOrig==62) = 13;TSV{13,1} = '7Networks_SomMot_4';
imAtlas(imAtlasOrig==63) = 14;TSV{14,1} = '7Networks_SomMot_5';
imAtlas(imAtlasOrig==64) = 15;TSV{15,1} = '7Networks_SomMot_6';
imAtlas(imAtlasOrig==67) = 16;TSV{16,1} = '7Networks_DorsAttn_Post_1';
imAtlas(imAtlasOrig==68) = 17;TSV{17,1} = '7Networks_DorsAttn_Post_2';
imAtlas(imAtlasOrig==69) = 18;TSV{18,1} = '7Networks_DorsAttn_Post_3';
imAtlas(imAtlasOrig==70) = 19;TSV{19,1} = '7Networks_DorsAttn_Post_4';
imAtlas(imAtlasOrig==71) = 20;TSV{20,1} = '7Networks_DorsAttn_Post_5';
TSV{21,1} = '7Networks_DorsAttn_Post_6';
imAtlas(imAtlasOrig==72) = 22;TSV{22,1} = '7Networks_DorsAttn_PrCv_1';
imAtlas(imAtlasOrig==73) = 23;TSV{23,1} = '7Networks_DorsAttn_FEF_1';
TSV{24,1} = '7Networks_SalVentAttn_ParOper_';
imAtlas(imAtlasOrig==76) = 25;TSV{25,1} = '7Networks_SalVentAttn_FrOperIns_1';
TSV{26,1} = '7Networks_SalVentAttn_FrOperIns_2';
TSV{27,1} = '7Networks_SalVentAttn_PFCl_1';
imAtlas(imAtlasOrig==77) = 28;TSV{28,1} = '7Networks_SalVentAttn_Med_1';
imAtlas(imAtlasOrig==78) = 29;TSV{29,1} = '7Networks_SalVentAttn_Med_2';
TSV{30,1} = '7Networks_SalVentAttn_Med_3';
imAtlas(imAtlasOrig==79) = 31;TSV{31,1} = '7Networks_Limbic_OFC_1';
imAtlas(imAtlasOrig==80) = 32;TSV{32,1} = '7Networks_Limbic_TempPole_1';
TSV{33,1} = '7Networks_Limbic_TempPole_2';
imAtlas(imAtlasOrig==81) = 34;TSV{34,1} = '7Networks_Cont_Par_1';
imAtlas(imAtlasOrig==83) = 35;TSV{35,1} = '7Networks_Cont_PFCl_1';
imAtlas(imAtlasOrig==89) = 36;TSV{36,1} = '7Networks_Cont_pCun_1';
imAtlas(imAtlasOrig==87) = 37;TSV{37,1} = '7Networks_Cont_Cing_1';
imAtlas(imAtlasOrig==91) = 38;TSV{38,1} = '7Networks_Default_Temp_1';
imAtlas(imAtlasOrig==92) = 39;TSV{39,1} = '7Networks_Default_Temp_2';
imAtlas(imAtlasOrig==90) = 40;TSV{40,1} = '7Networks_Default_Par_1';
TSV{41,1} = '7Networks_Default_Par_2';
TSV{42,1} = '7Networks_Default_PFC_1';
TSV{43,1} = '7Networks_Default_PFC_2';
TSV{44,1} = '7Networks_Default_PFC_3';
TSV{45,1} = '7Networks_Default_PFC_4';
TSV{46,1} = '7Networks_Default_PFC_5';
TSV{47,1} = '7Networks_Default_PFC_6';
TSV{48,1} = '7Networks_Default_PFC_7';
imAtlas(imAtlasOrig==99) = 49;TSV{49,1} = '7Networks_Default_pCunPCC_1';
imAtlas(imAtlasOrig==100) = 50;TSV{50,1} = '7Networks_Default_pCunPCC_2';
imAtlas(imAtlasOrig==65) = 51;TSV{51,1} = '7Networks_SomMot_7';
imAtlas(imAtlasOrig==66) = 52;TSV{52,1} = '7Networks_SomMot_8';
imAtlas(imAtlasOrig==74) = 53;TSV{53,1} = '7Networks_SalVentAttn_TempOccPar_1';
imAtlas(imAtlasOrig==75) = 54;TSV{54,1} = '7Networks_SalVentAttn_TempOccPar_2';
imAtlas(imAtlasOrig==82) = 55;TSV{55,1} = '7Networks_Cont_Par_2';
imAtlas(imAtlasOrig==84) = 56;TSV{56,1} = '7Networks_Cont_PFCl_2';
imAtlas(imAtlasOrig==85) = 57;TSV{57,1} = '7Networks_Cont_PFCl_3';
imAtlas(imAtlasOrig==86) = 58;TSV{58,1} = '7Networks_Cont_PFCl_4';
imAtlas(imAtlasOrig==88) = 59;TSV{59,1} = '7Networks_Cont_PFCmp_1';
imAtlas(imAtlasOrig==93) = 60;TSV{60,1} = '7Networks_Default_Temp_3';
imAtlas(imAtlasOrig==94) = 61;TSV{61,1} = '7Networks_Default_PFCv_1';
imAtlas(imAtlasOrig==95) = 62;TSV{62,1} = '7Networks_Default_PFCv_2';
imAtlas(imAtlasOrig==96) = 63;TSV{63,1} = '7Networks_Default_PFCdPFCm_1';
imAtlas(imAtlasOrig==97) = 64;TSV{64,1} = '7Networks_Default_PFCdPFCm_2';
imAtlas(imAtlasOrig==98) = 65;TSV{65,1} = '7Networks_Default_PFCdPFCm_3';

% Save the relabeled atlas and TSV for 7Networks
xASL_io_SaveNifti(fullfile(pathTPM, pathSchaefer, 'wSchaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii'),fullfile(pathTPM, pathSchaefer, 'wSchaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii'), imAtlas);% Save the new version
xASL_tsvWrite(TSV, fullfile(pathTPM, pathSchaefer, 'Schaefer_100Parcels_7Networks.tsv'), 1);

% Convert atlas to a bilateral version for 17 Networks
imAtlasOrig = xASL_io_Nifti2Im(fullfile(pathTPM, pathSchaefer, 'wSchaefer2018_100Parcels_17Networks_order_FSLMNI152_1mm.nii'));
imAtlas = zeros(size(imAtlasOrig(:,:,:,1))); % Create an empty atlas
TSV = {};

% Relabel the atlas
imAtlas(imAtlasOrig==51) = 1;TSV{1,1} = '17Networks_VisCent_ExStr_1';
imAtlas(imAtlasOrig==52) = 2;TSV{2,1} = '17Networks_VisCent_ExStr_2';
TSV{3,1} = '17Networks_VisCent_Striate_1';
imAtlas(imAtlasOrig==53) = 4;TSV{4,1} = '17Networks_VisCent_ExStr_3';
imAtlas(imAtlasOrig==55) = 5;TSV{5,1} = '17Networks_VisPeri_ExStrInf_1';
imAtlas(imAtlasOrig==54) = 6;TSV{6,1} = '17Networks_VisPeri_StriCal_1';
imAtlas(imAtlasOrig==56) = 7;TSV{7,1} = '17Networks_VisPeri_ExStrSup_1';
imAtlas(imAtlasOrig==57) = 8;TSV{8,1} = '17Networks_SomMotA_1';
imAtlas(imAtlasOrig==58) = 9;TSV{9,1} = '17Networks_SomMotA_2';
imAtlas(imAtlasOrig==61) = 10;TSV{10,1} = '17Networks_SomMotB_Aud_1';
imAtlas(imAtlasOrig==62) = 11;TSV{11,1} = '17Networks_SomMotB_S2_1';
imAtlas(imAtlasOrig==63) = 12;TSV{12,1} = '17Networks_SomMotB_S2_2';
imAtlas(imAtlasOrig==64) = 13;TSV{13,1} = '17Networks_SomMotB_Cent_1';
imAtlas(imAtlasOrig==65) = 14;TSV{14,1} = '17Networks_DorsAttnA_TempOcc_1';
imAtlas(imAtlasOrig==66) = 15;TSV{15,1} = '17Networks_DorsAttnA_ParOcc_1';
imAtlas(imAtlasOrig==67) = 16;TSV{16,1} = '17Networks_DorsAttnA_SPL_1';
imAtlas(imAtlasOrig==68) = 17;TSV{17,1} = '17Networks_DorsAttnB_PostC_1';
imAtlas(imAtlasOrig==69) = 18;TSV{18,1} = '17Networks_DorsAttnB_PostC_2';
TSV{19,1} = '17Networks_DorsAttnB_PostC_3';
imAtlas(imAtlasOrig==70) = 20;TSV{20,1} = '17Networks_DorsAttnB_FEF_1';
imAtlas(imAtlasOrig==71) = 21;TSV{21,1} = '17Networks_SalVentAttnA_ParOper_1';
imAtlas(imAtlasOrig==72) = 22;TSV{22,1} = '17Networks_SalVentAttnA_Ins_1';
TSV{23,1} = '17Networks_SalVentAttnA_Ins_2';
imAtlas(imAtlasOrig==73) = 24;TSV{24,1} = '17Networks_SalVentAttnA_ParMed_1';
imAtlas(imAtlasOrig==74) = 25;TSV{25,1} = '17Networks_SalVentAttnA_FrMed_1';
imAtlas(imAtlasOrig==76) = 26;TSV{26,1} = '17Networks_SalVentAttnB_PFCl_1';
imAtlas(imAtlasOrig==77) = 27;TSV{27,1} = '17Networks_SalVentAttnB_PFCmp_1';
imAtlas(imAtlasOrig==78) = 28;TSV{28,1} = '17Networks_LimbicB_OFC_1';
imAtlas(imAtlasOrig==79) = 29;TSV{29,1} = '17Networks_LimbicA_TempPole_1';
TSV{30,1} = '17Networks_LimbicA_TempPole_2';
imAtlas(imAtlasOrig==80) = 31;TSV{31,1} = '17Networks_ContA_IPS_1';
imAtlas(imAtlasOrig==81) = 32;TSV{32,1} = '17Networks_ContA_PFCl_1';
imAtlas(imAtlasOrig==82) = 33;TSV{33,1} = '17Networks_ContA_PFCl_2';
imAtlas(imAtlasOrig==86) = 34;TSV{34,1} = '17Networks_ContB_PFClv_1';
imAtlas(imAtlasOrig==88) = 35;TSV{35,1} = '17Networks_ContC_pCun_1';
TSV{36,1} = '17Networks_ContC_pCun_2';
imAtlas(imAtlasOrig==87) = 37;TSV{37,1} = '17Networks_ContC_Cingp_1';
imAtlas(imAtlasOrig==90) = 38;TSV{38,1} = '17Networks_DefaultA_PFCd_1';
imAtlas(imAtlasOrig==91) = 39;TSV{39,1} = '17Networks_DefaultA_pCunPCC_1';
imAtlas(imAtlasOrig==92) = 40;TSV{40,1} = '17Networks_DefaultA_PFCm_1';
TSV{41,1} = '17Networks_DefaultB_Temp_1';
TSV{42,1} = '17Networks_DefaultB_Temp_2';
TSV{43,1} = '17Networks_DefaultB_IPL_1';
imAtlas(imAtlasOrig==93) = 44;TSV{44,1} = '17Networks_DefaultB_PFCd_1';
TSV{45,1} = '17Networks_DefaultB_PFCl_1';
imAtlas(imAtlasOrig==94) = 46;TSV{46,1} = '17Networks_DefaultB_PFCv_1';
imAtlas(imAtlasOrig==95) = 47;TSV{47,1} = '17Networks_DefaultB_PFCv_2';
imAtlas(imAtlasOrig==96) = 48;TSV{48,1} = '17Networks_DefaultC_Rsp_1';
imAtlas(imAtlasOrig==97) = 49;TSV{49,1} = '17Networks_DefaultC_PHC_1';
imAtlas(imAtlasOrig==98) = 50;TSV{50,1} = '17Networks_TempPar_1';
imAtlas(imAtlasOrig==59) = 51;TSV{51,1} = '17Networks_SomMotA_3';
imAtlas(imAtlasOrig==60) = 52;TSV{52,1} = '17Networks_SomMotA_4';
imAtlas(imAtlasOrig==75) = 53;TSV{53,1} = '17Networks_SalVentAttnB_IPL_1';
imAtlas(imAtlasOrig==83) = 54;TSV{54,1} = '17Networks_ContB_Temp_1';
imAtlas(imAtlasOrig==84) = 55;TSV{55,1} = '17Networks_ContB_IPL_1';
imAtlas(imAtlasOrig==85) = 56;TSV{56,1} = '17Networks_ContB_PFCld_1';
imAtlas(imAtlasOrig==89) = 57;TSV{57,1} = '17Networks_DefaultA_IPL_1';
imAtlas(imAtlasOrig==99) = 58;TSV{58,1} = '17Networks_TempPar_2';
imAtlas(imAtlasOrig==100) = 59;TSV{59,1} = '17Networks_TempPar_3';

% Save the relabeled atlas and TSV for 17Networks
xASL_io_SaveNifti(fullfile(pathTPM, pathSchaefer, 'wSchaefer2018_100Parcels_17Networks_order_FSLMNI152_1mm.nii'),fullfile(pathTPM, pathSchaefer, 'wSchaefer2018_100Parcels_17Networks_order_FSLMNI152_1mm.nii'), imAtlas);% Save the new version
xASL_tsvWrite(TSV, fullfile(pathTPM, pathSchaefer, 'Schaefer_100Parcels_17Networks.tsv'), 1);

% Save the atlas
for iNet = [7,17]
	IM = xASL_io_Nifti2Im(fullfile(pathTPM, pathSchaefer, ['wSchaefer2018_100Parcels_' num2str(iNet) 'Networks_order_FSLMNI152_1mm.nii']));
	IM = uint8(IM);
	save(fullfile(pathTPM, pathSchaefer, ['Schaefer_100Parcels_' num2str(iNet) 'Networks.nii.mat']),'IM');
	xASL_adm_GzipNifti(fullfile(pathTPM, pathSchaefer, ['wSchaefer2018_100Parcels_' num2str(iNet) 'Networks_order_FSLMNI152_1mm.nii']));
	xASL_delete(fullfile(pathTPM, pathSchaefer, ['Schaefer_100Parcels_' num2str(iNet) 'Networks.nii']));
	xASL_Move(fullfile(pathTPM, pathSchaefer, ['wSchaefer2018_100Parcels_' num2str(iNet) 'Networks_order_FSLMNI152_1mm.nii.gz']),fullfile(pathTPM, pathSchaefer, ['Schaefer_100Parcels_' num2str(iNet) 'Networks.nii.gz']));
end
