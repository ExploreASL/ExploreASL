function xASL_im_VisualizeROIs(x, ROI_T1_list, ROI_FLAIR_list)
% xASL_im_VisualizeROIs Creates for each subject  a JPEG image containing
% the original T1w, WMH_SEGM and T1w after lesion-filling

    if ~isempty(ROI_T1_list)
        T1File = fullfile(x.D.PopDir, ['r' x.P.STRUCT '_' x.P.SubjectID '.nii']);

        for iL=1:length(ROI_T1_list)
            LesionFile = fullfile(x.D.PopDir, ['rROI_' x.P.STRUCT '_' num2str(iL) '_' x.P.SubjectID '.nii']);
            TempIM = xASL_io_Nifti2Im(LesionFile);
            LesionIM = xASL_im_ConvertMap2Mask(TempIM(:,:,:,1));

            %OutIm1      = visual_registration_MNI_1_5_INPUT(xASL_io_Nifti2Im(T1File), LesionIM, 0.75,3,'green');
			OutIm1 = xASL_im_CreateVisualFig(x, {T1File LesionIM}, [], [0.75 0.35], [], {x.S.gray x.S.green});

            fprintf('%s\n', ['Printing check ROI image ' num2str(iL)]);
			xASL_adm_CreateDir(x.D.ROICheckDir);
            xASL_imwrite((OutIm1+eps)./max(OutIm1(:)), fullfile(x.D.ROICheckDir ,['ROI_' x.P.SubjectID '_' num2str(iL) '.jpg']));
        end
    end


end
