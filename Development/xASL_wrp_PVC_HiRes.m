function [ x ] = xASL_wrp_PVC_HiRes( x )
% xASL_wrp_PVC_HiRes This function does a standard partial volume error correction,
% but on high resolution using 3D B-splines, which is faster


if  length(xASL_adm_GetFileList( x.D.PopDir, ['^q' x.P.CBF '_(GM|WM)_PVEC_.*_.*.nii']))==x.dataset.nSubjectsSessions*2
    % All images are already created, skip this
    fprintf('%s\n','PV corrected CBF images were already created, skipping...')
else
    fprintf('%s\n','Creating PV-corrected CBF images using B-splines...')

    for iS=1:x.nSubjects % iterate over subjects
        for iSess=1:x.dataset.nSessions % iterate over sessions
            clear CBFname pGMname CBFim pGMim CBFim pGMim CBFgmName CBFwmName

            CBFgmName   = fullfile( x.D.PopDir, ['q' x.P.CBF '_GM_PVEC_' x.SUBJECTS{iS} '_' x.SESSIONS{iSess} '.nii']);
            CBFwmName   = fullfile( x.D.PopDir, ['q' x.P.CBF '_WM_PVEC_' x.SUBJECTS{iS} '_' x.SESSIONS{iSess} '.nii']);

            CBFname     = fullfile( x.D.PopDir, ['q' x.P.CBF '_' x.SUBJECTS{iS} '_' x.SESSIONS{iSess} '.nii']);
            pGMname     = fullfile( x.D.PopDir, ['PV_pGM_' x.SUBJECTS{iS} '.nii']);
            pWMname     = fullfile( x.D.PopDir, ['PV_pWM_' x.SUBJECTS{iS} '.nii']);

            CBFim       = xASL_io_ReadNifti(CBFname);
            pGMim       = xASL_io_ReadNifti(pGMname);
            pWMim       = xASL_io_ReadNifti(pWMname);

            CBFim       = CBFim.dat(:,:,:);
            pGMim       = pGMim.dat(:,:,:);
            pWMim       = pWMim.dat(:,:,:);



            %% Perform PVEc
            % number of B-splines doesn't matter, FoV is same for all
            % sequences, since it is in MNI space
            PVEc_kernel                             = [14 16 14]; % 8x8x8 mm @ 3 mm resolution [20 20 20]; % [5 5 5]; alternative = 16 16 8 or 20 20 8
            PVim(:,:,:,1)                           = pGMim;
            PVim(:,:,:,2)                           = pWMim;

            %% Mask CBF image
%             for iSlice=1:size(pGMim,3)
%                 MaskpGM(:,:,iSlice)                 = pGMim(:,:,iSlice)>0.1*max(max((pGMim(:,:,iSlice))));
%             end
            MaskpGM                                 = x.GMSliceMask;
            MaskpWM                                 = pWMim>0.8;
            MaskIM                                  = (pGMim+pWMim)>0.1;

            CBFim                                   = CBFim.*MaskIM;

            %% All NaNs & non-finite should be set to 0
            CBFim(~isfinite(CBFim))                 = 0;
            CBFim                                   = double(CBFim);

            tic
            [imPVEC,imCBFrec,imResidual,FWHM]       = xASL_im_PVCbspline(CBFim,PVim,PVEc_kernel);
            toc

            CBFgm                                   = imPVEC(:,:,:,1).*MaskpGM;
            CBFwm                                   = imPVEC(:,:,:,2).*MaskpWM;

%             % Clip images
%             LowClip                                 = mean(nonzeros(CBFgm)) - (3 .* std(nonzeros(CBFgm)));
%             HighClip                                = mean(nonzeros(CBFgm)) + (3 .* std(nonzeros(CBFgm)));
%
%             CBFgm(CBFgm>HighClip)                   = HighClip;
%             CBFgm(CBFgm<LowClip)                    = LowClip;
%
%             LowClip                                 = mean(nonzeros(CBFwm)) - (3 .* std(nonzeros(CBFwm)));
%             HighClip                                = mean(nonzeros(CBFwm)) + (3 .* std(nonzeros(CBFwm)));
%
%             CBFwm(CBFwm>HighClip)                   = HighClip;
%             CBFwm(CBFwm<LowClip)                    = LowClip;

            xASL_io_SaveNifti( x.D.ResliceRef, CBFgmName, CBFgm);
            xASL_io_SaveNifti( x.D.ResliceRef, CBFwmName, CBFwm);


    %         mean(mean(mean(CBFgm(logical(MaskpGM)))))
    %         mean(mean(mean(CBFwm(logical(MaskpWM)))))
    %
    %         figure(1);imshow(xASL_im_rotate(double(imPVEC(:,:,53,1))./3.19.*double(x.S.masks.skull(:,:,53)),90),[0 400/3.19],'colormap',jet,'InitialMagnification',250)
    %         figure(2);imshow(xASL_im_rotate(double(imPVEC(:,:,53,2))./3.19.*double(x.S.masks.skull(:,:,53)),90),[0 60],'colormap',jet,'InitialMagnification',250)
    %
    %         figure(1);imshow(xASL_im_rotate(double(CBFgm(:,:,53))./3.19.*double(x.S.masks.skull(:,:,53)),90),[0 120],'colormap',jet,'InitialMagnification',250)
    %         figure(2);imshow(xASL_im_rotate(double(CBFwm(:,:,53))./3.19.*double(x.S.masks.skull(:,:,53)),90),[0 50],'colormap',jet,'InitialMagnification',250)

        end
    end

end

end

















%             %% Estimate effective resolution
%             % Obtain native image resolution, to estimate FWHM
%             NativeName  = fullfile( x.D.ROOT, x.SUBJECTS{iS}, x.SESSIONS{iSess}, [x.P.ASL4D  '.nii']);
%             NativeIM    = xASL_io_ReadNifti(NativeName);
%             NativeRes   = NativeIM.hdr.pixdim(2:4);
%
%             % Determine estimated FWHM to convert native image resolution to
%             % effective resolution
%             if      strcmpi(x.Sequence,'2D_EPI')
%                     EstimatedEffectiveResolution    = [(3.7+3.6)/2 (3.7+3.6)/2 7.6]; % Paper Jan Petr, 3.7 & 3.6 where X Y, which should be similar
%                     BasedOnResolution               = [3 3 7];
%                     Estimated_FWHM                  = EstimatedEffectiveResolution./BasedOnResolution;
%
%             elseif  strcmpi(x.Sequence,'3D_spiral')
%                     EstimatedEffectiveResolution    = [(4.9+5.1)/2 (4.9+5.1)/2 9.5]; % Paper Jan Petr, 4.9 & 5.1 where X Y, which should be similar
%                     BasedOnResolution               = [1.875 1.875 4];
%                     Estimated_FWHM                  = EstimatedEffectiveResolution./BasedOnResolution;
%
%             elseif  strcmpi(x.Sequence,'3D_GRASE')
%                     Estimated_FWHM                  = ([1.2167 1.2167 1.0857] + [2.6667 2.6667 2.3750]) ./2;
%                     % This is a rough educated guess, FWHM of 3D GRASE should
%                     % be in between 3D spiral & 2D EPI.
%             end
%
%             EffectiveResolution                     = NativeRes.*Estimated_FWHM;
%
%
%             %% Smooth pGM to same resolution
%             Smoothing_FWHM                          = (EffectiveResolution.^2 - [1.5 1.5 1.5].^2).^0.5;
%             VoxelSize                               = 1.5;
%             Smoothing_FWHM                          = Smoothing_FWHM./VoxelSize;
%             FwHm2SD                                 = (2*(2*reallog(2))^0.5);
%             Smoothing_SD                            = (Smoothing_FWHM./FwHm2SD);
