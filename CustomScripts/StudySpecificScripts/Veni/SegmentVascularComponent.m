clear IM

ddir    = 'C:\Backup\ASL\GENFI\GENFI_DF2_2_Long\analysis\dartel';
Flist   = xASL_adm_GetFileList(ddir, '^qCBF_untreated_PH_Achieva_Bsup_GRN\d{3}_1_ASL_1\.nii$');

for iL=1:10
    IM(:,:,:,iL)    = xASL_io_Nifti2Im(Flist{iL});
end

% 2D EPI timeseries -> provide more ways to acquire a vascular template

% 3D single images do not allow for that,
% but still, they allow to do the same between subjects
% also between subjects, the same holds true that the structural component
% is very comparable, but the vascular component is variable, except for
% perhaps the major arteries

% they have in common that they have extreme values, variable, % do not
% follow the slow frequencies, nor the WM/GM differentiation

% diff with template
% within-subject or between-subject variance
% difference with smoothed CBF -> noise + structure (vessels)
% off from structural/CBF relation (PVEc, pGM/CBF, but requires multiple
% images, or iteration)
% can do the same by creating a distance map from the deep WM skeleton,
% which should on average keep the same sigmoidal relation (from WM coming
% close to WM border, then steep increase at the border, then GM

% but this relation is not the same for subcortical structures as for
% GM cortex


IM              = xASL_im_rotate(IM,90);

IM=IM(:,:,:,1);

%dip_image([xASL_stat_MeanNan(IM,4) dip_array(smooth(xASL_stat_StdNan(IM,[],4),2))])

%smoothIM        = dip_array(smooth(IM,2));
smoothIM        = xASL_im_ndnanfilter(IM,'gauss',[2 2 2]*2.335,0);
diffIM          = smoothIM-IM;

%dip_image([IM smoothIM 10.*abs(dip_array(smooth(diffIM,2)))])
