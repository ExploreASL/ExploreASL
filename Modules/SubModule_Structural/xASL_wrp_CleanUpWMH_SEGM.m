function xASL_wrp_CleanUpWMH_SEGM(x)
%xASL_wrp_CleanUpWMH_SEGM Submodule of ExploreASL Structural Module, that cleans up under- and over-segmentations of WMH SEGM
%
% FORMAT: xASL_wrp_CleanUpWMH_SEGM(x)
%
% INPUT:
%   x 	    - structure containing fields with all information required to run this submodule (REQUIRED)
%   x.P     - paths with NIfTIs for which this function should be applied to (REQUIRED)
%
% OUTPUT: n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule aims to clean up WMH under- or oversegmentations in a conservatively & robust way,
% i.e. erring on the side of caution. It uses input from the tissue class segmentation (e.g. CAT12) to repair the
% WMH segmentation (e.g. LST LPA/LGA or externally provided). Note that before running the tissue segmentation, the T1w was
% (conservatively) filled for WMH lesions. % This function is not tested a lot, so mainly conservatively set up to improve the
% WMH volumetrics, rather than improve the registration.
%
% This submodule contains the following steps:
% 0) Administration
% 1) Correct pGM islands inside pWM
%    WMH can have an intensity similar to GM on the T1w, which erroneously classifies them as GM instead of WM(H).
%    The rule used here, is to define GM islands within the WM as clusters of pGM>0.05 for which 3 layers (dilations)
%    have at least 95% pWM. For these islands, pGM is given 100% to pWM. 50% of pWM is given to pWMH (the pWMH/pNAWM distinction is
%    made later in the pipeline, here still pWM=pWMH+pNAWM). The reason is that not all low T1w intensities within the WM
%    are WMH, we still expect some lacunes, perivascular (Virchow-Robin) spaces, which could be considered pNAWM rather than pWMH.
% 2) Perform brainmasking & join masks
% 3) Correct any WMH inside GM or CSF -> here we assume that CAT12 did a good segmentation job. If pGM is larger than pWM & larger
%    than pWMH, we consider a voxel to be pGM and remove the pWMH. This effectively removes pWMH segmentation noise in the GM or CSF,
%    it doesn't correct any significant misclassification of WMH in the GM or CSF. If the WMH segmentation does a significant misclassification
%    (e.g. setting the pWMH inside GM or CSF to a probability higher than GM or CSF is by tissue segmentation), this is lesion filled
%    after the WMH segmentation, on the T1w, hence the tissue segmentation won't have a chance to correct this. Fortunately, most
%    oversegmentations in the GM/CSF have low pWMH, as WMH segmentation algorithms already perform a light tissue prior-based
%    clean up themselves.
% 4) Saving & file management
% 5) Prepare visuals for visual QC & file management
%
%
% EXAMPLE: xASL_wrp_CleanUpWMH_SEGM(x);
% __________________________________
% Copyright 2015-2019 ExploreASL
%
% 2019-05-02 HJM



%% -------------------------------------------------------------------------------
%% 0) File & paths management
CleanUpFile = [x.P.Path_WMH_SEGM(1:end-4) '_CleanUp.nii'];
PreMaskFile = fullfile(x.SUBJECTDIR, 'Ones.nii');
MaskFile = fullfile(x.SUBJECTDIR, 'mask_Ones.nii');
PathMNIMask = fullfile(x.D.MapsSPMmodifiedDir, 'brainmask_supratentorial.nii');

if xASL_exist(x.P.Path_rWMH_SEGM, 'file')
    WMH_File = x.P.Path_rWMH_SEGM;
elseif xASL_exist(x.P.Path_WMH_SEGM, 'file')
    WMH_File = x.P.Path_WMH_SEGM;
else
    fprintf('%s\n', 'No WMH_SEGM file, skipping WMH_SEGM correction');
    return;
end

pGMim = xASL_io_Nifti2Im(x.P.Path_c1T1);
pWMim = xASL_io_Nifti2Im(x.P.Path_c2T1);
WMHim = xASL_io_Nifti2Im(WMH_File);

if xASL_stat_SumNan(pGMim(:))==0
    error('Invalid pGM image');
elseif xASL_stat_SumNan(pWMim(:))==0
    error('Invalid pWM image');
end

if xASL_exist(x.P.Path_c3T1, 'file')
    pCSFim = xASL_io_Nifti2Im(x.P.Path_c3T1);
else
    pCSFim = max(0, 1-pGMim-pWMim);
end

fprintf('%s', 'Cleaning up WMH_SEGM: 0%');



%% -------------------------------------------------------------------------------
%% 1) Correct pGM islands inside pWM
% Create pGMrois
WMHnew = zeros(size(WMHim)); % create dummy WMHmask
GMroi = pGMim>0.05; % find any GM that is an island
GMLabels = spm_bwlabel(double(GMroi)); % create labels for all connected GM regions

if max(size(WMHim)~=size(pGMim)) % if FLAIR & T1w werent same voxelsize
    error('FLAIR wasnt resampled to T1w space');
end

if xASL_stat_SumNan(GMLabels(:))==0
    warning('No GM islands found, something wrong with GM image?');
else
    for iVol=1:max(GMLabels(:)) % get volumes for all regions
        GMvolume(iVol,1) = iVol;
        GMvolume(iVol,2) = sum(GMLabels(:)==iVol);
    end

    ROIs = find(GMvolume(:,2)<max(GMvolume(:,2))); % select GM regions that are smaller than largest GM region

    for iROI=1:length(ROIs) % loop over all regions
        xASL_TrackProgress(iROI, length(ROIs));
        CurrROI = GMLabels==GMvolume(ROIs(iROI), 1); % current region (skipping the largest)
        DilatedROI = xASL_im_DilateErodeFull(CurrROI, 'dilate', xASL_im_DilateErodeSphere(2)); % dilate region
        % PM: here we could do multiple erosions, for multiple layers with different conditions
        Rim = max(0,DilatedROI - CurrROI); % obtain the rim of the region (i.e. the dilation part)
        if  sum(Rim(:))>0 % bugfix
            if sum(sum(sum(pWMim(logical(Rim))))) > 0.95*sum(Rim(:))
                % this checks whether the GMroi is truely an island
                % (the 0.95 allowing for a little bit of GM near the edge, this number can go up)
                % now we have any GM islands within the WM.

                % we could also try multiple layers with different conditions, e.g. assuming that the first layer
                % should contain 100% WM, second layer >0.75, third layer >0.50. But for now keep it simple.

                % These non-WM islands are partly WMH, partly WM
                % 1) in this ROI, divide the pGM & pCSF over WMHnew & pWM
                pWMim(CurrROI) = min(1,pGMim(CurrROI) + pCSFim(CurrROI) + pWMim(CurrROI)); % give fully to pWM
                WMHnew(CurrROI) = min(1,0.5.*(pGMim(CurrROI) + pCSFim(CurrROI)) + WMHnew(CurrROI)); % give half to pWMH
                pGMim(CurrROI) = 0; % we have given all pGM to pWM
                pCSFim(CurrROI) = 0; % we have given all pCSF to pWM
            end
        end
    end
    fprintf('\n');
end



%% -------------------------------------------------------------------------------
%% 2) Remove anything outside the brainmask
xASL_Copy(x.P.Path_c1T1, PreMaskFile, 1);
OnesIm = xASL_io_Nifti2Im(PreMaskFile);
OnesIm(:) = 1;
xASL_io_SaveNifti(PreMaskFile, PreMaskFile, OnesIm, [], 0);

xASL_im_SkullStrip(PreMaskFile, PathMNIMask, x);
Bmask = logical(xASL_io_Nifti2Im(MaskFile));
WMHnew(~Bmask) = 0;

% Combine existing & new WMH masks
WMHim = min(1, WMHim + WMHnew);
% dip_image([FLAIRim+WMHim.*2 pGMim FLAIRim+WMHnew.*2])




%% -------------------------------------------------------------------------------
%% 3) Correct oversegmentation of WMH inside the pGM or inside the pCSF
% We assume that CAT12 does a very good job segmenting, so we trust pGM (CAT12) & pCSF (CAT12) more than pWMH (LST)
% But only when we are really certain that this is GM, and WMH will usually
% have very low probabilities in these regions
WMHim(WMHim>0 & (pGMim.*0.8)>(WMHim+pWMim)) = 0;
WMHim(WMHim>0 & (pCSFim.*0.8)>(WMHim+pWMim)) = 0;
% we don't have to change the pGM & pCSF here, they were segmented by a different program (CAT12) to total 1.

% dip_image([pGMim+3.*(WMHim>0 & pGMim<WMHim) pGMim+3.*(WMHim>0 & pGMim>WMHim)])
% dip_image([pGMim+3.*(WMHim>0 & pCSFim<WMHim) pGMim+3.*(WMHim>0 & pCSFim>WMHim)])



%% -------------------------------------------------------------------------------
%% 4) Save segmentations & file management
xASL_io_SaveNifti(WMH_File, CleanUpFile, WMHim, [], false);
xASL_io_SaveNifti(x.P.Path_c1T1, x.P.Path_c1T1, pGMim, [], false);
xASL_io_SaveNifti(x.P.Path_c2T1, x.P.Path_c2T1, pWMim, [], false);

if xASL_exist(x.P.Path_c3T1, 'file')
    xASL_io_SaveNifti(x.P.Path_c3T1, x.P.Path_c3T1, pCSFim, [], false);
end

xASL_delete(MaskFile);
xASL_delete(PreMaskFile);




%% 5 File management & visual QC
CleanUpFile = [x.P.Path_WMH_SEGM(1:end-4) '_CleanUp.nii'];
CleanUpFile_Pop = [x.P.Pop_Path_WMH_SEGM(1:end-4) '_CleanUp.nii'];
INname  = {x.P.Path_FLAIR x.P.Path_WMH_SEGM CleanUpFile};
OUTname = {x.P.Pop_Path_rFLAIR x.P.Pop_Path_rWMH_SEGM CleanUpFile_Pop};
xASL_spm_deformations(x, INname, OUTname);

OutIM1 = xASL_im_CreateVisualFig(x, {x.P.Pop_Path_rFLAIR x.P.Pop_Path_rWMH_SEGM}, [], [1 1.5], [], {x.S.gray x.S.red}, false);
OutIM2 = xASL_im_CreateVisualFig(x, {x.P.Pop_Path_rFLAIR CleanUpFile_Pop}, [], [1 1.5], [], {x.S.gray x.S.red}, false);
xASL_imwrite([OutIM1,OutIM2], fullfile(x.D.FLAIR_CheckDir, ['CleanUp_WMH_SEGM_' x.P.SubjectID '.jpg']));

xASL_delete(x.P.Pop_Path_rWMH_SEGM);
xASL_delete(CleanUpFile_Pop);
xASL_Move(CleanUpFile, x.P.Path_WMH_SEGM, true);


end
