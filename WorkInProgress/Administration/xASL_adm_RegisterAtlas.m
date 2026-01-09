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
pathExploreASL = '/home/janpetr/ExploreASL/ExploreASL';

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
xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMAtlasDir, pathTPMAtlasFileGMWM));
xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMAtlasDir, pathTPMAtlasFile));
xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMAtlasDir, pathTPMAtlasFileWM));

xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFileGMWM));
xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFile));
xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFileWM));

xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathColinDir, pathColinFileGM));
xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathColinDir, pathColinFileWM));
xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathColinDir, pathColinFileGMWM));

% And resample to 1.5mm DARTEL space
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii')};
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
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_2_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_3_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_4_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_5_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;

spm_jobman('run',matlabbatch);
%% Transform the template to the correct template, just as a test
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {fullfile(pathTPM, pathTPMAtlasDir, ['u_ws' pathTPMAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
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
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_2_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_3_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_4_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_5_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;

spm_jobman('run',matlabbatch);
%% Transform the template to the correct template, just as a test
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {fullfile(pathTPM, pathTPMSymAtlasDir, ['u_ws' pathTPMSymAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
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
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_2_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_3_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_4_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_5_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;

spm_jobman('run',matlabbatch);
%% Transform the Colin27 template to the correct template, just as a test
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {fullfile(pathTPM, pathColinDir, ['u_ws' pathColinFileGM])};
matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
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
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {fullfile(pathTPM, pathTPMSymAtlasDir, ['u_ws' pathTPMSymAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
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
imAtlasLargeFunc(imAtlasMax==2) = 4;TSVLargeFunc{1,1} = 'Primary Sensory Relay';%LGN
imAtlasLargeFunc(imAtlasMax==3) = 4;
imAtlasLargeFunc(imAtlasMax==4) = 4;%MGN
imAtlasLargeFunc(imAtlasMax==5) = 4;

xASL_io_SaveNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'ThalamusProbs.MNIsymSpace.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii'), imAtlas);
xASL_io_SaveNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'ThalamusProbs.MNIsymSpace.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHist.nii'), imAtlasLargeHist);
xASL_io_SaveNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'ThalamusProbs.MNIsymSpace.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunc.nii'), imAtlasLargeFunc);

% Nearest neighbor transformation
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {fullfile(pathTPM, pathTPMSymAtlasDir, ['u_ws' pathTPMSymAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM, pathTPMSymAtlasDir, pathTPMSymAtlasFile)
												   fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii')
												   fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHist.nii')
												   fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunc.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {fullfile(pathTPM, 'FreeSurferSubfields')};
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'o';

spm_jobman('run',matlabbatch);

% Create MAT file
xASL_tsvWrite(TSV, fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.tsv'), 1);
xASL_tsvWrite(TSVLargeHist, fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHist.tsv'), 1);
xASL_tsvWrite(TSVLargeFunc, fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunc.tsv'), 1);
IM = xASL_io_Nifti2Im(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamus.nii'));
IMHist = xASL_io_Nifti2Im(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusHist.nii'));
IMFunc = xASL_io_Nifti2Im(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusFunc.nii'));
IM = uint8(IM);
IMHist = uint8(IMHist);
IMFunc = uint8(IMFunc);
save(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.mat'),'IM');
save(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHist.mat'),'IMHist');
save(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunc.mat'),'IMFunc');
xASL_adm_GzipNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamus.nii'));
xASL_adm_GzipNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusHist.nii'));
xASL_adm_GzipNifti(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusFunc.nii'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHist.nii'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunc.nii'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii.mat'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHist.nii.mat'));
xASL_delete(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunc.nii.mat'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamus.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii.gz'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusHist.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHist.nii.gz'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'oFreesurferThalamusFunc.nii.gz'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunc.nii.gz'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.mat'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamus.nii.mat'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHist.mat'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusHist.nii.mat'));
xASL_Move(fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunc.mat'),fullfile(pathTPM, 'FreeSurferSubfields', 'FreesurferThalamusFunc.nii.mat'));

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
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {fullfile(pathTPM, pathTPMSymAtlasDir, ['u_ws' pathTPMSymAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
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
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {fullfile(pathTPM, pathTPMAtlasDir, ['u_ws' pathTPMAtlasFile])};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
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
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.flowfield = {fullfile(pathTPM, pathColinDir, ['u_ws' pathColinFileGM])};
matlabbatch{1}.spm.util.defs.comp{2}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{2}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{2}.dartel.template = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_6_IXI555_MNI152.nii')};
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
