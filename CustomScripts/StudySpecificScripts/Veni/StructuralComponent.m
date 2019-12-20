%% Create synthetic CBF image

pGM     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\Cmp_rc1T1.nii');
pWM     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\Cmp_rc2T1.nii');
CBF     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\Philips_2DEPI_Bsup_CBF_Template.nii');

IM1     = pGM.*60+pWM.*10;

IM1     = IM1./max(IM1(:));

%Smoothed    = dip_array(smooth(IM1,1.5));
Smoothed    = xASL_im_ndnanfilter(IM1,'gauss',[1.5 1.5 1.5]*2.335,0);
Smoothed    = Smoothed./max(Smoothed(:));
CBF         = CBF./max(CBF(:));
Smoothed(isnan(Smoothed))   = 0;
CBF(isnan(CBF))             = 0;
CBF(CBF<0.1)                = 0;
Smoothed(Smoothed<0.05)     = 0;

%dip_image([xASL_im_rotate(IM1(:,:,53),90)./max(IM1(:)) xASL_im_rotate(Smoothed(:,:,53),90) xASL_im_rotate(CBF(:,:,53),90)])

xASL_io_SaveNifti('C:\ExploreASL\Maps\Templates\Cmp_rc1T1.nii','C:\ExploreASL\Maps\Templates\Veni_DARTEL\Smoothed_structInfo.nii', Smoothed );
xASL_io_SaveNifti('C:\ExploreASL\Maps\Templates\Cmp_rc1T1.nii','C:\ExploreASL\Maps\Templates\Veni_DARTEL\CBF.nii', CBF );

Diff1   = abs(CBF - IM1);
Diff2   = abs(CBF - Smoothed);

%dip_image([xASL_im_rotate(Diff1(:,:,53),90) xASL_im_rotate(Diff2(:,:,53),90)]  )

%% Estimate deformations

clear matlabbatch

matlabbatch{1}.spm.tools.dartel.warp.images = {
                                               {
                                               'C:\ExploreASL\Maps\Templates\Veni_DARTEL\CBF.nii,1'
                                               'C:\ExploreASL\Maps\Templates\Veni_DARTEL\Smoothed_structInfo.nii,1'
                                               }
                                               }';
matlabbatch{1}.spm.tools.dartel.warp.settings.template = 'Template';
matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;

spm_jobman('run',matlabbatch);

%% Apply deformation

clear matlabbatch

matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {'C:\ExploreASL\Maps\Templates\Veni_DARTEL\u_Smoothed_structInfo_Template.nii'};
matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {''};
matlabbatch{1}.spm.util.defs.comp{2}.inv.comp{1}.dartel.flowfield = {'C:\ExploreASL\Maps\Templates\Veni_DARTEL\u_CBF_Template.nii'};
matlabbatch{1}.spm.util.defs.comp{2}.inv.comp{1}.dartel.times = [1 0];
matlabbatch{1}.spm.util.defs.comp{2}.inv.comp{1}.dartel.K = 6;
matlabbatch{1}.spm.util.defs.comp{2}.inv.comp{1}.dartel.template = {''};
matlabbatch{1}.spm.util.defs.comp{2}.inv.space = {'C:\ExploreASL\Maps\rbrainmask.nii'};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {'C:\ExploreASL\Maps\Templates\Veni_DARTEL\Smoothed_structInfo.nii'};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savesrc = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';

spm_jobman('run',matlabbatch);


IM3     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\Veni_DARTEL\wSmoothed_structInfo.nii');
IM3     = IM3./max(IM3(:));

Diff3   = abs(CBF - IM3);
%dip_image([xASL_im_rotate(Diff1(:,:,53),90) xASL_im_rotate(Diff2(:,:,53),90) xASL_im_rotate(Diff3(:,:,53),90)])
%dip_image(xASL_im_rotate(IM3(:,:,53),90))
%dip_image([xASL_im_rotate(pGM,90) xASL_im_rotate(pWM,90)])
