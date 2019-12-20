%% Create pGM atlas

List    = xASL_adm_GetFileList('C:\Backup\ASL','^rc1T1_.*\.(nii|nii\.gz)$','FPListRec');

fprintf('%s\n','Loading...  ')
Nn  = 1;
for iL=1:length(List)
    xASL_TrackProgress(iL,length(List));

    try

        tIM                 = xASL_io_ReadNifti(List{iL});
        tIM                 = tIM.dat(:,:,53);

        pGMsort(Nn,1)       = xASL_stat_SumNan(tIM(:));
        IM(:,:,Nn)          = tIM;
        Nn=Nn+1;

    end
end

%% Smooth images
for iS=1:size(IM,3)
    xASL_TrackProgress(iS,size(IM,3));
    %IMsmooth(:,:,iS) = dip_array(smooth(IM(:,:,iS),[3 3]));
	IMsmooth(:,:,iS) = xASL_im_ndnanfilter(IM(:,:,iS),'gauss',[3 3]*2.335,0);
end

%% Reshape the pGM

IMsmooth     = reshape(IMsmooth,[121*145 size(IMsmooth,3)]);
IM2          = reshape(IM      ,[121*145 size(IM      ,3)]);
for iS=1:size(IM2,1)
    clear tIM
    xASL_TrackProgress(iS,size(IM2,1));
    tIM         = IMsmooth(iS,:)';
    tIM(:,2)    = [1:1:length(tIM)];
    tIM         = sortrows(tIM,1);
    IM2(iS,:)   = IM2(iS,tIM(:,2));
end

IM3     = xASL_im_rotate(reshape(IM2,[121 145 size(IM,3)]),90);
%% Sort pGM image


%
%
%
% [X N]   = hist(pGMsort);
%
% figure(1);plot(N,X)
%
% pGMsort(:,2)    = [1:1:length(pGMsort)];
% pGMsort         = sortrows(pGMsort,1);
%
% IMsorted        = IM(:,:,pGMsort(:,2)');
% IMsorted        = xASL_im_rotate(IMsorted,90);


IMmean(:,:,1)   = mean(IM3(:,:,1:100),3);
IMmean(:,:,2)   = mean(IM3(:,:,101:500),3);
IMmean(:,:,3)   = mean(IM3(:,:,501:1000),3);
IMmean(:,:,4)   = mean(IM3(:,:,1001:2000),3);
IMmean(:,:,5)   = mean(IM3(:,:,2001:end),3);

%dip_image(singlesequencesort(IMmean,5))
