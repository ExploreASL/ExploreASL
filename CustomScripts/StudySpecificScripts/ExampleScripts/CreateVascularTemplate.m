%% Create vascular template from BioCog ASL data

% VascularList    = {''};
load('C:\Backup\ASL\BioCog\VascularList.mat');
for iS=1:length(VascularList)
    IM(:,:,:,iS)    = xASL_io_Nifti2Im(fullfile('C:\Backup\ASL\BioCog\DataFreeze1\dartel',['PWI_' VascularList{iS} '_ASL_1.nii']) );
end

% Average vessel image
VascTemplate        = ClipVesselImage(xASL_stat_MeanNan(IM,4));
VascTemplate2       = ClipVesselImage(xASL_stat_MedianNan(IM,4));
VascTemplate3       = abs((VascTemplate+VascTemplate2.*2)./2);
%VascTemplateSmooth  = dip_array(smooth(VascTemplate3 ,2/2.355));
VascTemplateSmooth  = xASL_im_ndnanfilter(VascTemplate3,'gauss',[2 2 2],0);
VascTemplateSmooth  = VascTemplateSmooth-min(VascTemplateSmooth(:));


% Maximal vessel image
MaxTemplate     = abs(max(IM,[],4));
MaxTemplate     = ClipVesselImage(MaxTemplate);
%MaxTemplate     = dip_array(smooth(MaxTemplate,2/2.355));
MaxTemplate     = xASL_im_ndnanfilter(MaxTemplate,'gauss',[2 2 2],0);
MaxTemplate     = MaxTemplate-min(MaxTemplate(:));
MaxTemplate(MaxTemplate<150)     = 0;
%MaxTemplate     = dip_array(smooth(MaxTemplate,2/2.355));
MaxTemplate     = xASL_im_ndnanfilter(MaxTemplate,'gauss',[2 2 2],0);

% Mirror maximal vessel image left & right -> make symmetrical
Lmax            = MaxTemplate;
Rmax            = MaxTemplate;
Lmax(1:60,:,:)  = 0;
Rmax(61:end,:,:)= 0;

Lmax             = Lmax + xASL_im_Flip(Lmax,1);
Rmax             = Rmax + xASL_im_Flip(Rmax,1);

Mtemp(:,:,:,1)  = Lmax;
Mtemp(:,:,:,2)  = Rmax;

MaxTemplate     = max(Mtemp,[],4);
%MaxTemplate     = dip_array(smooth(MaxTemplate,2/2.355));
MaxTemplate     = xASL_im_ndnanfilter(MaxTemplate,'gauss',[2 2 2],0);

xASL_io_SaveNifti('C:\ExploreASL\Maps\Templates\MeanVesselTemplate.nii', 'C:\ExploreASL\Maps\Templates\MaxVesselTemplate.nii', MaxTemplate);

xASL_io_SaveNifti('c:\ExploreASL\Maps\rgrey.nii','C:\ExploreASL\Maps\Templates\MaxVesselTemplate.nii',MaxTemplate);
xASL_io_SaveNifti('c:\ExploreASL\Maps\rgrey.nii','C:\ExploreASL\Maps\Templates\MeanVesselTemplate.nii',VascTemplateSmooth);




%dip_image([xASL_im_rotate(VascTemplate,90) xASL_im_rotate(VascTemplate2,90).*2 xASL_im_rotate(VascTemplateSmooth,90) xASL_im_rotate(VascTemplate2Smooth,90).*2 xASL_im_rotate(VascTemplate3Smooth,90) abs(xASL_im_rotate(VascTemplate3Smooth,90))])

%dip_image([xASL_im_rotate(VascTemplateSmooth,90).*10 xASL_im_rotate(MaxTemplate,90)])
