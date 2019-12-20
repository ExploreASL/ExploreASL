%% Delete tempfiles

ExploreASL_Master('',0);

ROOT{1}         = 'C:\Backup\ASL';

% %%
% File2Del    = {'RawTemplate' 'Mean_CBF_Template' 'ATT_BiasField' 'Mask_Template' 'VascularArtifact_Template'};
%
% for iR=1:length(ROOT)
%     for iD=1:length(File2Del)
%         xASL_adm_DeleteFileList(ROOT{iR},['^' File2Del{iD} '.(nii|nii\.gz)$'],1,[0 Inf]);
%     end
% end

% %% Unzip if necessary
%
% for iR=1:length(ROOT)
%     PathList    = xASL_adm_GetFileList(ROOT{iR},['^.*\.nii.gz$'],'FPListRec',[0 Inf]);
%         for iP=1:length(PathList)
%             try
%                 xASL_adm_UnzipNifti(PathList{iP});
%             end
%         end
% end

% %% Convert to 16 bit, if not already
% NiftiNames  = {'slice_gradient' 'WMH_SEGM' 'T1' 'FLAIR' 'SNR_control' 'SNR_control_beforeMoCo' 'SNR_PWI_beforeMoCo' 'SNR' 'SD' 'SD_PWI_beforeMoCo' 'SD_control_beforeMoCo' 'SD_control' 'SD' 'mean_control'};
%
% for iR=1:length(ROOT)
%     for iN=1:length(NiftiNames)
%         PathList    = xASL_adm_GetFileList(ROOT{iR},['^(r*)' NiftiNames{iN} '.*\.(nii|nii\.gz)$'],'FPListRec',[0 Inf]);
%
%         for iP=1:length(PathList)
%             try
%                 clear Fname tIM
%                 Fname   = PathList{iP};
%                 tIM     = xASL_io_ReadNifti(Fname);
%                 if  tIM.hdr.bitpix>16
%                     xASL_io_SaveNifti(Fname,Fname,tIM.dat(:,:,:,:,:),16);
%                 end
%             catch
%             end
%         end
%     end
% end


%% Zip all nii files
for iR=1:length(ROOT)
    PathList    = xASL_adm_GetFileList(ROOT{iR},'^.*\.nii$','FPListRec',[0 Inf]);
    fprintf('%s',[ROOT{iR} '...   ']);
    for iP=1:length(PathList)
        xASL_TrackProgress(iP,length(PathList));
        NiiFile     = PathList{iP} % show which file we're working on, to avoid corrupt gz file
        GzFile      = [PathList{iP} '.gz'];

        if  exist(NiiFile) && ~exist(GzFile)
            gzip(NiiFile); % stopping here can create a corrupt gz file
        end

        if  exist(NiiFile) && exist(GzFile) % stopping here is fine
            delete(NiiFile);
        end
    end
    fprintf('\n');
end
