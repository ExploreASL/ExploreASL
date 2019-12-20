%% Create pGM atlas

List    = xASL_adm_GetFileList('C:\Backup\ASL','^c1T1\.(nii|nii\.gz)$','FPListRec');

fprintf('%s\n','Loading...  ')
Nn  = 1;
for iL=1:100 % length(List)
    xASL_TrackProgress(iL,100);

    unzipped_files  = xASL_adm_UnzipNifti( List{iL} );

    try
        [Fpath Ffile Fext]  = fileparts(unzipped_files{1});
        rFile               = fullfile(Fpath,['r' Ffile Fext]);
        T1File              = fullfile(Fpath,[Ffile(3:end) Fext]);
        rT1File             = fullfile(Fpath,['r' Ffile(3:end) Fext]);
        xASL_adm_UnzipNifti(T1File);

        clear matlabbatch
        matlabbatch{1}.spm.spatial.coreg.write.ref              = {'C:\ExploreASL\Maps\rbrainmask.nii,1'};
        matlabbatch{1}.spm.spatial.coreg.write.source           = {[unzipped_files{1} ',1']};
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp  = 1;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap    = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask    = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix  = 'r';

        spm_jobman('run',matlabbatch);

        clear matlabbatch
        matlabbatch{1}.spm.spatial.coreg.write.ref              = {'C:\ExploreASL\Maps\rbrainmask.nii,1'};
        matlabbatch{1}.spm.spatial.coreg.write.source           = {[T1File ',1']};
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp  = 1;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap    = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask    = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix  = 'r';

        spm_jobman('run',matlabbatch);





        tIM                 = xASL_io_ReadNifti(rFile);
        tIM                 = tIM.dat(:,:,56);

        pGMsort(Nn,1)       = xASL_stat_SumNan(tIM(:));
        IM(:,:,Nn)          = tIM;

        tIM                 = xASL_io_ReadNifti(rT1File);
        IM_T1(:,:,Nn)       = tIM.dat(:,:,56);

        Nn=Nn+1;

        delete(rFile);
        delete(rT1File);
    end
end





[X N]   = hist(pGMsort);

figure(1);plot(X, N)

pGMsort(:,2)    = [1:1:length(pGMsort)];
pGMsort         = sortrows(pGMsort,1);

IMsorted        = IM(:,:,pGMsort(:,2)');
IMsorted        = xASL_im_rotate(IMsorted,90);


IMmean(:,:,1)   = mean(IMsorted(:,:,1:10),3);
IMmean(:,:,2)   = mean(IMsorted(:,:,11:50),3);
IMmean(:,:,3)   = mean(IMsorted(:,:,51:1000),3);
IMmean(:,:,4)   = mean(IMsorted(:,:,1001:2000),3);
IMmean(:,:,5)   = mean(IMsorted(:,:,2001:end),3);

dip_image(singlesequencesort(IMmean,5))
