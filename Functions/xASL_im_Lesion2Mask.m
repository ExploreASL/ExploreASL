function LesionIM = xASL_im_Lesion2Mask(LesionPath, T1path, pGMpath, pWMpath, x)
% xASL_im_Lesion2Mask For a standard space lesion mask (or map), this stores
% the lesion mask, and in additional its perimask (15 mm) and contralateral
% mask, as 2nd and 3rd volumes
% It plots the masks on a T1 image, and masks the new masks with the
% subjects' brainmask (pGM+pWM)
%
% FORMAT:       LesionIM = xASL_im_Lesion2Mask(LesionPath, T1path, pGMpath, pWMpath, x)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  For a standard space lesion mask (or map), this stores
%               the lesion mask, and in additional its perimask (15 mm) and contralateral
%               mask, as 2nd and 3rd volumes.
%               It plots the masks on a T1 image, and masks the new masks with the
%               subjects' brainmask (pGM+pWM).
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL


% Distinguish between lesion & ROI masks
[Fpath, Ffile, ~]  = xASL_fileparts(LesionPath);
if ~isempty(regexp(Ffile,'ROI'))
	fprintf('%s\n','Creating ROI masks');
	if exist('x', 'var')
		if isfield(x.D, 'ROICheckDir')
			ImageSaveDir = x.D.ROICheckDir;
		end
	end

elseif ~isempty(regexp(Ffile,'Lesion'))
	fprintf('%s\n','Creating Lesion masks');
	if  exist('x', 'var')
		if isfield(x.D,'LesionCheckDir')
			ImageSaveDir  = x.D.LesionCheckDir;
		end
	end

else
	warning('Wrong lesion/ROI mask name');
	ImageSaveDir = Fpath;
end


%% -------------------------------------------------------------------------
if ~xASL_exist(LesionPath,'file')
    fprintf('%s\n',['Skipped because mask didnt exist: ' LesionPath]);
else

    LesionIM    = xASL_io_Nifti2Im(LesionPath);

    % If lesion is empty, skip this & delete the file
    if  sum(LesionIM(:))<0.00001
        clear LesionIM
        xASL_delete(LesionPath);
        fprintf('%s\n',[LesionPath ' removed because no mask was provided (empty mask)']);
    else
        fprintf('%s\n','Printing lesion QC image');
        LesionIM    = LesionIM(:,:,:,1);
        LesionIM    = xASL_im_ConvertMap2Mask(LesionIM);

        DistanceMap = xASL_im_DistanceTransform(LesionIM);
        DistanceMap = 1.5 .* DistanceMap; % convert voxels to 1.5 mm
        PeriMask    = DistanceMap>0 & DistanceMap<=16.6; % 25 mm

        % BrainMasking
        pGM         = xASL_io_Nifti2Im(x.P.Pop_Path_rc1T1);
        pWM         = xASL_io_Nifti2Im(x.P.Pop_Path_rc2T1);
        BrainMask   = (pGM+pWM)>0.5;
        PeriMask    = PeriMask.*BrainMask;
        AddMask     = LesionIM+PeriMask;
        ContraMask  = xASL_im_Flip(AddMask,1).*BrainMask;

        % Create hemispheres
        LeftHemisphere  = ones(size(BrainMask));
        RightHemisphere = ones(size(BrainMask));
        LeftHemisphere(1:round(size(LeftHemisphere,1)/2),:,:)     = 0;
        RightHemisphere(round(size(RightHemisphere,1)/2):end,:,:) = 0;

        if      sum(sum(sum(LeftHemisphere(logical(AddMask)))))>sum(sum(sum(RightHemisphere(logical(AddMask)))))
                Hemisphere      = (LeftHemisphere .* BrainMask) & ~AddMask;
                ContraLateral   = RightHemisphere .* BrainMask;
        else    Hemisphere      = (RightHemisphere .* BrainMask) & ~AddMask;
                ContraLateral   = LeftHemisphere .* BrainMask;
        end

		OutIm1 = xASL_vis_CreateVisualFig( x, {x.P.Pop_Path_rT1 LesionIM}, [], [0.75 0.35], [], {x.S.gray x.S.red});
		OutIm2 = xASL_vis_CreateVisualFig( x, {x.P.Pop_Path_rT1 PeriMask}, [], [0.75 0.35], [], {x.S.gray x.S.green});
		OutIm3 = xASL_vis_CreateVisualFig( x, {x.P.Pop_Path_rT1 ContraMask}, [], [0.75 0.35], [], {x.S.gray x.S.blue});

		% Check masks, where to combine the images
        OutIm1Mask  = OutIm1(:,:,1)~=OutIm1(:,:,2) | OutIm1(:,:,2)~=OutIm1(:,:,3)  | OutIm1(:,:,1)~=OutIm1(:,:,3);
        OutIm1Mask  = repmat(OutIm1Mask,[1 1 3]);

        OutIm2Mask  = OutIm2(:,:,1)~=OutIm2(:,:,2) | OutIm2(:,:,2)~=OutIm2(:,:,3)  | OutIm2(:,:,1)~=OutIm2(:,:,3);
        OutIm2Mask  = repmat(OutIm2Mask,[1 1 3]);

        OutIm3Mask  = OutIm3(:,:,1)~=OutIm3(:,:,2) | OutIm3(:,:,2)~=OutIm3(:,:,3)  | OutIm3(:,:,1)~=OutIm3(:,:,3);
        OutIm3Mask  = repmat(OutIm3Mask,[1 1 3]);

        % Combine OutIm1 2 & 3
        OutIm4                  = zeros(size(OutIm1));
        OutIm4(OutIm1==OutIm2)  = OutIm1(OutIm1==OutIm2);
        OutIm4(OutIm1Mask)      = OutIm1(OutIm1Mask);
        OutIm4(OutIm2Mask)      = OutIm2(OutIm2Mask);
        OutIm4(OutIm3Mask)      = OutIm3(OutIm3Mask);

        jpgFile = fullfile(ImageSaveDir, ['Regions_' Ffile '.jpg']);
		xASL_adm_CreateDir(ImageSaveDir);
        xASL_vis_Imwrite(OutIm4, jpgFile);

        % Save segmentations
        % intra, pGM+pWM
        LesionIM(:,:,:,2) = PeriMask; % peri, pGM+pWM
        LesionIM(:,:,:,3) = AddMask;  % intra+peri, pGM+pWM
        LesionIM(:,:,:,4) = ContraMask; % contralateral, pGM+pWM
        LesionIM(:,:,:,5) = Hemisphere; % hemisphere - AddMask
        LesionIM(:,:,:,6) = ContraLateral; % contralateral hemisphere

        LesionIM = logical(LesionIM);
        xASL_io_SaveNifti(LesionPath,LesionPath,LesionIM,8);
    end
end

end
