%% CreateFiguresHongerWinter Cohorts


IM{1}   = squeeze(ASL_untreated.Data.data(x.S.SetsID(:,16)==1,:,:,:));
IM{2}   = squeeze(ASL_untreated.Data.data(x.S.SetsID(:,16)==2,:,:,:));

IM{1}(isnan(IM{1}))     = 0;
IM{2}(isnan(IM{2}))     = 0;

SaveDir     = fullfile(x.D.PopDir,'CreateImageVesselsStats');

for iC=1:length(IM)

    for iI=1:size(IM{iC},1)
        xASL_TrackProgress(iI,size(IM{iC},1));
        clear matlabbatch tFile tIM

        tFile   = fullfile(SaveDir,[num2str(iI) '.nii']);

        xASL_io_SaveNifti(x.D.ResliceRef,tFile,squeeze(IM{iC}(iI,:,:,:)));

        %% Edge-preserved smoothing
        matlabbatch{1}.spm.tools.cat.tools.sanlm.data       = {tFile};
        matlabbatch{1}.spm.tools.cat.tools.sanlm.prefix     = 'smooth';
        matlabbatch{1}.spm.tools.cat.tools.sanlm.NCstr      = 1;
        matlabbatch{1}.spm.tools.cat.tools.sanlm.rician     = 1;

        spm_jobman('run',matlabbatch);
        delete(tFile);

        %% Remove bias-field
        tFile   = fullfile(SaveDir,['smooth' num2str(iI) '.nii']);

        tIM             = xASL_io_Nifti2Im(tFile);

        %tBiasField      = dip_array(smooth(tIM,8));
		tBiasField      = xASL_im_ndnanfilter(tIM,'gauss',[8 8 8]*2.335,0);
        IM{iC}(iI,:,:,:) = tIM./tBiasField;
        delete(tFile);
    end
end

%% Smooth again
for iC=1:length(IM)
    for iI=1:size(IM{iC},1)
        xASL_TrackProgress(iI,size(IM{iC},1));
        %IM{iC}(iI,:,:,:) = dip_array(smooth(IM{iC}(iI,:,:,:),1));
		IM{iC}(iI,:,:,:) = xASL_im_ndnanfilter(IM{iC}(iI,:,:,:),'gauss',[1 1 1]*2.335,0);
    end
end

%dip_image([xASL_im_rotate(squeeze(IM{2}(end,:,:,:)),90) dip_array(smooth(xASL_im_rotate(squeeze(IM{2}(end,:,:,:)),90),1))])

%dip_image(smooth([xASL_im_rotate(squeeze(std(IM{1},[],1)),90) xASL_im_rotate(squeeze(std(IM{2},[],1)),90)],1.5))


%% Individual examples

spatialCoV          = x.S.SetsID(:, 4);
spatialCoV(:,2)     = [1:1:length(spatialCoV)];
spatialCoV(:,3)     = x.S.SetsID(:,16);
spatialCoV          = sortrows(spatialCoV,1);

[N1 X1]             = hist(spatialCoV(x.S.SetsID(:,16)==1));
[N2 X2]             = hist(spatialCoV(x.S.SetsID(:,16)==2));
N1                  = N1./sum(N1);
N2                  = N2./sum(N2);

figure(1);plot(X1,N1,'r-',X2,N2,'b-')

x.SUBJECTS{113}
