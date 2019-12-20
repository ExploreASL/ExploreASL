function xASL_wrp_VisualCheckLesionRemoval(x, Lesion_T1_list, Lesion_FLAIR_list)
% xASL_wrp_VisualCheckLesionRemoval Creates for each subject  a JPEG image containing
% the segmented image with removed lesion, T1+lesion ROI, segmentation before removal

    if ~isempty(Lesion_T1_list) || ~isempty(Lesion_FLAIR_list)

        % c1T1 & c2T1 before lesion masking
        if  xASL_exist(x.P.Path_c1T1_ORI,'file') && xASL_exist(x.P.Path_c2T1_ORI,'file') % if T1w was lesion-filled or bias-field corrected
            xASL_io_ReadNifti(x.P.Path_c1T1_ORI); % unzip if needed
            xASL_io_ReadNifti(x.P.Path_c2T1_ORI); % unzip if needed
            xASL_spm_deformations(x,{x.P.Path_c1T1_ORI;x.P.Path_c2T1_ORI},{x.P.Pop_Path_rc1T1_ORI;x.P.Pop_Path_rc2T1_ORI}); % no DARTEL yet, no LongReg yet
        end

        LesionIM    = zeros(121,145,121); % assuming 1.5 mm MNI

		% Call this for Lesions only and not ROIs
        [INname,OUTname]     = xASL_wrp_LesionResliceList(x,1,1,0,0);

        if ~isempty(INname) && ~isempty(OUTname)
            % First dilate ROIs, if they were e.g. used for annotation (single voxel only)
            % Do linear interpolation to avoid negative edge effects
            for iO=1:length(OUTname)
                if ~xASL_exist(OUTname{iO})
                    for iL=1:length(INname)
                        xASL_im_dilateROI(INname{iL});
                    end
                    xASL_spm_deformations(x,INname,OUTname,1);
                end
            end
        end

        for iL=1:length(Lesion_T1_list)
            LesionFile      = fullfile( x.D.PopDir, ['rLesion_' x.P.STRUCT '_' num2str(iL) '_' x.P.SubjectID '.nii']);
            TempIM          = xASL_io_Nifti2Im(LesionFile);
            LesionIM        = min(1,LesionIM+xASL_im_ConvertMap2Mask(TempIM(:,:,:,1)));
        end
        for iL=1:length(Lesion_FLAIR_list)
            LesionFile      = fullfile( x.D.PopDir, ['rLesion_' x.P.FLAIR '_' num2str(iL) '_' x.P.SubjectID '.nii']);
            TempIM          = xASL_io_Nifti2Im(LesionFile);
            LesionIM        = min(1,LesionIM+xASL_im_ConvertMap2Mask(TempIM(:,:,:,1)));
		end

		ImOut1 = xASL_im_CreateVisualFig( x, {x.P.Pop_Path_rc1T1_ORI x.P.Pop_Path_rc2T1_ORI}, [], [1 0.8], [], []);
		ImOut2 = xASL_im_CreateVisualFig( x, {x.P.Pop_Path_rT1 LesionIM},[], [1 0.8], [], []);
		ImOut3 = xASL_im_CreateVisualFig( x, {x.P.Pop_Path_rc1T1 x.P.Pop_Path_rc2T1},[], [1 0.8], [], []);

        IM = [ImOut1,ImOut2,ImOut3];
        xASL_adm_CreateDir(x.D.LesionCheckDir);
        xASL_imwrite((IM+eps)./max(IM(:)), fullfile(x.D.LesionCheckDir , ['Lesion_corr_' x.P.SubjectID '.jpg']));

    end
end
