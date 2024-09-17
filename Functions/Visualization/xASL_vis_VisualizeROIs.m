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
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    if ~isempty(ROI_list)
        T1File = fullfile(x.D.PopDir, ['r' x.P.STRUCT '_' x.P.SubjectID '.nii']);
        for iROI=1:length(ROI_list)
			[~, Ffile] = xASL_fileparts(ROI_list{iROI});
            ROIFile = fullfile(x.D.PopDir, ['r' Ffile '_' x.P.SubjectID '.nii']);
            TempIM = xASL_io_Nifti2Im(ROIFile);
            ROIIM = xASL_im_ConvertMap2Mask(TempIM(:,:,:,1));
			OutIm1 = xASL_vis_CreateVisualFig(x, {T1File ROIIM}, [], [0.75 0.35], [], {x.S.gray x.S.green});
            fprintf('%s\n', ['Printing check ROI image ' num2str(iROI)]);
            xASL_vis_Imwrite((OutIm1+eps)./max(OutIm1(:)), fullfile(x.D.ROICheckDir ,['ROI_' x.P.SubjectID '_' num2str(iROI) '.jpg']));
        end
    end
end
