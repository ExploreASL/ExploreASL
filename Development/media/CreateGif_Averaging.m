tIM     = xASL_io_Nifti2Im('C:\Backup\ASL\BioCog\DataFreeze1\BCU005_1\ASL_1\ASL4D.nii');

Pairs   = tIM(:,:,:,[1:2:end]) - tIM(:,:,:,[2:2:end]);
Pairs   = squeeze(Pairs(:,:,8,:));
Pairs   = xASL_im_rotate(Pairs,90);
% Pairs   = Pairs(:,:,1:30);

Pairs   = xASL_im_CropParmsApply(Pairs,10,70,18,62);

clear MeanIM
for iP=1:size(Pairs,3)
    MeanIM(:,:,iP)  = mean(Pairs(:,:,1:iP),3);
end

% [xmin xmax ymin ymax] = xASL_im_CropParmsAcquire(MeanIM(:,:,end));

SelectIM                = [1:10 12 14 16 18 20 23 26 30];

RootDir     = 'C:\Users\henkj\Google Drive\AMC\Presentations and other pieces\2017 Neuroscience Master lecture ASL\ExampleDataLecture';
SaveFile    = fullfile(RootDir, 'AveragingASL.gif');
SaveNII     = fullfile(RootDir, 'AveragingASL.nii');
SaveNIIpairs= fullfile(RootDir, 'SaveNIIpairs.nii');
SaveNII_Ori = fullfile(RootDir, 'CBF_2DEPI.nii');
xASL_io_SaveNifti(SaveNII_Ori, SaveNII, MeanIM);
xASL_io_SaveNifti(SaveNII_Ori, SaveNIIpairs, Pairs);
UpsampleNII( SaveNII, SaveNII, [0.35 0.35 1]);
UpsampleNII( SaveNIIpairs, SaveNIIpairs, [0.35 0.35 1.5]);

MeanIM      = xASL_io_Nifti2Im(SaveNII);
Pairs       = xASL_io_Nifti2Im(SaveNIIpairs);
MeanIM      = MeanIM(:,:,2:end);
% Convert to "colorimage"
MeanIM      = MeanIM+30;
MeanIM(MeanIM>100)   = 100;
MeanIM      = MeanIM./max(MeanIM(:));
delete(SaveNII);

for iP=1:size(MeanIM,3)
    [ImIndex,ColorMap]   = rgb2ind(repmat(MeanIM(:,:,iP),[1 1 3]),256);
    if iP==1
        imwrite(ImIndex,ColorMap,SaveFile,'gif', 'Loopcount',inf, 'TransparentColor',double(min(ImIndex(:))),'DelayTime',0); %,'Screensize',[80 80]
    else
        imwrite(ImIndex,ColorMap,SaveFile,'gif', 'WriteMode','append', 'TransparentColor',double(min(ImIndex(:))),'DelayTime',0);
    end
end

TotalView               = singlesequencesort(Pairs(:,:,SelectIM),6);

% jet_256     = jet(256);
% jet_256(1,:)= 0;

figure(1);imshow(TotalView,[-10 70],'InitialMagnification',250,'Border','Tight')
