%% PVC example Figure 1
% This creates the example Figure in the pipeline flow chart figure of the ExploreASL paper

ExploreASL_Master('',0);

pGM     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\rc1T1_ASL_res.nii');
pWM     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\rc2T1_ASL_res.nii');
CBF     = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\Philips_2DEPI_Bsup_CBF_Template.nii');

WBmask  = (pGM+pWM)>0.1 & isfinite(CBF);

gwpv                       = pGM(WBmask);
gwpv(:,2)                  = pWM(WBmask);
CBFin                      = CBF(WBmask);

gwcbf                      = (CBFin')*pinv(gwpv');

pGM_mask                   = single(pGM>0.5);
pWM_mask                   = single(pWM>0.5);

meanGMcbf                  = mean(CBF(pGM_mask & isfinite(CBF)));
meanWMcbf                  = mean(CBF(pWM_mask & isfinite(CBF)));
ScaleGM                    = gwcbf(1)/meanGMcbf;
ScaleWM                    = gwcbf(2)/meanWMcbf;

CBF_GM                     = CBF.*pGM_mask.*ScaleGM;
CBF_WM                     = CBF.*pWM_mask.*ScaleWM;

CBF_GM(CBF_GM==0)          = NaN;
CBF_WM(CBF_WM==0)          = NaN;

CBF_GM                     = xASL_im_ndnanfilter(CBF_GM,'gauss',[8 8 8],1);
CBF_WM                     = xASL_im_ndnanfilter(CBF_WM,'gauss',[8 8 8],1);

CBF_GM(isnan(CBF_GM))          = 0;
CBF_WM(isnan(CBF_WM))          = 0;

CBF_PVC                    = CBF_GM+CBF_WM;

jet256      = jet(256);
jet256(1,:) = 0;

figure(1);imshow(imrotate(CBF_PVC(:,:,53),90),[0 135],'colormap',jet256,'border','tight')
