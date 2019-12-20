%% Creation EPI

ROOT    = 'C:\Users\henkj\Google Drive\AMC\Presentations and other pieces\2017 Neuroscience Master lecture ASL\ExampleDataLecture';

EPI     = xASL_io_Nifti2Im( fullfile(ROOT, 'mMean_2DEPI_Siemens.nii'));
EPI(:,:,[1:12 107:end])   = 0;
% dip_image(EPI)

CBF     = xASL_io_Nifti2Im( fullfile(ROOT, 'Mean_2D_EPI_CBF_Siemens.nii' ));
max(EPI(:))

EPI     = EPI + 2.5.*CBF;

%% Creation EPI with Bsup
SliceG                          = ones(size(EPI));
for iSl=1:size(EPI,3)
    SliceG(:,:,iSl)             = iSl;
end

SliceG                          = SliceG./max(SliceG(:));

EPI_Bsup                        = EPI.*SliceG;

dip_image([EPI EPI_Bsup])

xASL_io_SaveNifti('C:\ExploreASL\Maps\Templates\Siemens_3DGRASE_Template.nii','C:\ExploreASL\Maps\Templates\Siemens_2DEPI_noBsup.nii',EPI);
xASL_io_SaveNifti('C:\ExploreASL\Maps\Templates\Siemens_3DGRASE_Template.nii','C:\ExploreASL\Maps\Templates\Siemens_2DEPI_Bsup.nii',EPI_Bsup);

%% No motion without Bsup
UpMoveN     = 10;

EPI1        = EPI + CBF.*5;
EPI2        = EPI - CBF.*5;


%% Motion without Bsup
UpMoveN     = 10;

EPI1        = EPI + CBF.*5;
EPI2        = EPI - CBF.*5;
for ii=1:UpMoveN
    EPI1(:,:,1+ii:end)     = (EPI2(:,:,1:end-ii) + EPI1(:,:,1+ii:end) )./2;
end



%% Motion with Bsup
UpMoveN     = 10;

EPI1        = EPI_Bsup + CBF.*5;
EPI2        = EPI_Bsup - CBF.*5;
for ii=1:UpMoveN
    EPI1(:,:,1+ii:end)     = (EPI2(:,:,1:end-ii) + EPI1(:,:,1+ii:end) )./2;
end




deltaM  = abs((EPI1-EPI2) + CBF.*10);

dip_image([EPI2(:,:,53) EPI1(:,:,53) deltaM(:,:,53).*1.5])






% CBF signal ~ 5%, so if CBF = 50 mL/100g/min, then raw signal = 1000
