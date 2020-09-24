function xASL_vis_VisualizeROIs(x, ROI_list)
% xASL_vis_VisualizeROIs Creates for each subject a JPEG image containing
% the original T1w, WMH_SEGM and T1w after lesion-filling
%
% FORMAT:       xASL_vis_VisualizeROIs(x, ROI_list)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Creates for each subject a JPEG image containing
%               the original T1w, WMH_SEGM and T1w after lesion-filling.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

    if ~isempty(ROI_list)
        T1File = fullfile(x.D.PopDir, ['r' x.P.STRUCT '_' x.P.SubjectID '.nii']);

        for iL=1:length(ROI_list)
			[~, Ffile] = xASL_fileparts(ROI_list{iL});
            LesionFile = fullfile(x.D.PopDir, ['r' Ffile '_' x.P.SubjectID '.nii']);

            TempIM = xASL_io_Nifti2Im(LesionFile);
            LesionIM = xASL_im_ConvertMap2Mask(TempIM(:,:,:,1));

			OutIm1 = xASL_vis_CreateVisualFig(x, {T1File LesionIM}, [], [0.75 0.35], [], {x.S.gray x.S.green});

            fprintf('%s\n', ['Printing check ROI image ' num2str(iL)]);
			xASL_adm_CreateDir(x.D.ROICheckDir);
            xASL_vis_Imwrite((OutIm1+eps)./max(OutIm1(:)), fullfile(x.D.ROICheckDir ,['ROI_' x.P.SubjectID '_' num2str(iL) '.jpg']));
        end
    end

end
