function xASL_vis_OverlapT1_ASL( x, ASL)
%xASL_vis_OverlapT1_ASL Part of ExploreASL.
% Shows spatial agreement ASL and probability maps.
%
% FORMAT:       xASL_vis_OverlapT1_ASL( x, ASL)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Part of ExploreASL.
%               Shows spatial agreement ASL and probability maps.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL



    SaveDir                     = fullfile(x.D.PopDir, 'check_overlap_T1_ASL');
    xASL_adm_CreateDir(SaveDir);
    SaveFile                   = fullfile(SaveDir,'check_overlap_T1_ASL.jpg');

    if ~exist(SaveFile,'file')

        fprintf('%s\n','Printing spatial agreement mean ASL with mean GM probability map');

        % Load GM probability maps
		fprintf('%s\n','Loading GM probability maps...  ')
        for ii=1:x.nSubjects
            xASL_TrackProgress(ii,x.nSubjects);

            x.P.SubjectID              = xASL_adm_GetFileList(x.D.PopDir,['^' x.c_PreFix{1} '_' x.SUBJECTS{ii} '.(nii|nii\.gz)$']);

            if ~isempty(x.P.SubjectID)
                prob_map{1}(:,:,:,ii)   = xASL_io_Nifti2Im( x.P.SubjectID{1} );
            end

        end

        fprintf('\n');

        ASL(~isfinite(ASL))                                         = 0;
        prob_map{1}(~isfinite(prob_map{1}))                         = 0;

        scale                       = mean(ASL(:))./mean(prob_map{1}(:));

        T1_temp                     = xASL_vis_CropParmsApply(xASL_im_rotate( mean(prob_map{1}(:,:,x.slicesLarge,:),4),90),x.S.TransCrop(1),x.S.TransCrop(2),x.S.TransCrop(3),x.S.TransCrop(4) );
        ASL_temp                    = xASL_vis_CropParmsApply(xASL_im_rotate( squeeze(mean(ASL(:,:,:,x.slicesLarge),1)),90),x.S.TransCrop(1),x.S.TransCrop(2),x.S.TransCrop(3),x.S.TransCrop(4) );

        T1_color                    = xASL_vis_TileImages(T1_temp , 4);
        ASL_color                   = xASL_vis_TileImages(ASL_temp, 4);
        T1_color                    = ind2rgb(round( T1_color  .*255),x.red);
        ASL_color                   = ind2rgb(round( ASL_color .*255./scale),x.yellow);
        check_nii                   = T1_color + ASL_color;


        xASL_vis_Imwrite(check_nii, SaveFile);
    end

end
