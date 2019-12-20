%% Create vascular probability map

MeanVess                = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\MeanVesselTemplate.nii');
MeanVess                = xASL_im_rotate(MeanVess(:,:,53),90);
MeanVess(MeanVess<5)    = 5;
MeanVess(MeanVess>30)   = 30;
MeanVess                = MeanVess-5;
MeanVess                = MeanVess./max(MeanVess(:));

MaxVess                = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\MaxVesselTemplate.nii');
MaxVess                = xASL_im_rotate(MaxVess(:,:,53),90);
MaxVess(MaxVess<150)   = 150;
MaxVess(MaxVess>300)   = 300;
MaxVess                = MaxVess-150;
MaxVess                = MaxVess./max(MaxVess(:));

IndMax                  = [0:0.25:1].^0.5;
IndMean                 = [1:-0.25:0].^2;

Bmask                   = xASL_io_Nifti2Im('C:\ExploreASL\Maps\Templates\brainmask.nii');
Bmask                   = xASL_im_rotate(Bmask(:,:,53),90);
%Bmask                   = dip_array(smooth(Bmask,3));
Bmask                   = xASL_im_ndnanfilter(Bmask,'gauss',[3 3 3]*2.335,0);
Bmask                   = Bmask.^3;

for iV=1:5
    ProbVess(:,:,iV)   = IndMean(iV).*MeanVess + IndMax(iV).* MaxVess;
    %ProbVess(:,:,iV)   = dip_array(smooth(ProbVess(:,:,iV),1));
	ProbVess(:,:,iV)   = xASL_im_ndnanfilter(ProbVess(:,:,iV),'gauss',[1 1 1]*2.335,0);
    ProbVess(:,:,iV)   = ProbVess(:,:,iV) .* Bmask;
end

jet256         = jet(256);
jet256(1,:)    = 0;




figure(1);imshow(singlesequencesort(ProbVess),[],'InitialMagnification',250,'border', 'tight','Colormap',jet256);



for iV=1:5
    figure(iV);imshow(ProbVess(:,:,iV),[],'InitialMagnification',250,'border', 'tight','Colormap',jet256);
end
