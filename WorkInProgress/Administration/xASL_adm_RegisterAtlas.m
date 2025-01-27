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

xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathColinDir, pathColinFileGM));
xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathColinDir, pathColinFileWM));
xASL_im_PreSmooth(fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii'),fullfile(pathTPM, pathColinDir, pathColinFileGMWM));

% And resample to 1.5mm DARTEL space
matlabbatch = [];
matlabbatch{1}.spm.util.defs.comp{1}.id.space = {fullfile(pathExploreASL, 'External', 'SPMmodified', 'toolbox', 'cat12', 'templates_volumes' ,'Template_1_IXI555_MNI152.nii')};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {fullfile(pathTPM, pathTPMAtlasDir, ['s' pathTPMAtlasFileGMWM])
	fullfile(pathTPM, pathTPMAtlasDir, ['s' pathTPMAtlasFile])
	fullfile(pathTPM, pathTPMAtlasDir, ['s' pathTPMAtlasFileWM])
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
