%% 1.5  Get tissue contrast for control & PWI images, if both exist
% Rationale here was to paste this in the registration function, to check
% which image has better contrast available for registration
if xASL_exist(x.P.Path_mean_control, 'file') && xASL_exist(x.P.Path_mean_PWI_Clipped, 'file')
    
    % 1) Get tissue masks
    EstimatedResolution = xASL_init_DefaultEffectiveResolution(x.P.Path_mean_control, x);
    
    xASL_im_PreSmooth(x.P.Path_mean_control, x.P.Path_c1T1, x.P.Path_PVgm, EstimatedResolution);
    xASL_im_PreSmooth(x.P.Path_mean_control, x.P.Path_c2T1, x.P.Path_PVwm, EstimatedResolution);
    xASL_im_PreSmooth(x.P.Path_mean_control, x.P.Path_c3T1, x.P.Path_PVcsf, EstimatedResolution);
    xASL_spm_reslice(x.P.Path_mean_control, x.P.Path_PVgm, [], [], x.Quality, x.P.Path_PVgm, 1);
    xASL_spm_reslice(x.P.Path_mean_control, x.P.Path_PVwm, [], [], x.Quality, x.P.Path_PVwm, 1);
    xASL_spm_reslice(x.P.Path_mean_control, x.P.Path_PVcsf, [], [], x.Quality, x.P.Path_PVcsf, 1);

    pGM = xASL_io_Nifti2Im(x.P.Path_PVgm);
    pWM = xASL_io_Nifti2Im(x.P.Path_PVwm);
    pCSF = xASL_io_Nifti2Im(x.P.Path_PVcsf);
    
    WBmask = (pGM+pWM+pCSF)>0.5;
    WBmaskeroded = xASL_im_DilateErodeFull(WBmask, 'erode', xASL_im_DilateErodeSphere(2));
    for iErode=1:3
        WBmaskeroded = xASL_im_DilateErodeFull(WBmaskeroded, 'erode', xASL_im_DilateErodeSphere(2));
    end
    
    GMmask = pGM>(pWM+pCSF);
    WMmask = pWM>(pGM+pCSF);
    CSFmask = pCSF>(pGM+pWM);
    CSFmask(~WBmaskeroded) = 0;
    
    % 2) Get median values per tissue type
    ControlIm = xASL_io_Nifti2Im(x.P.Path_mean_control);
    ControlGM = median(ControlIm(GMmask & isfinite(ControlIm)));
    ControlWM = median(ControlIm(WMmask & isfinite(ControlIm)));
    ControlParenchyma = (ControlGM+ControlWM)/2;
    ControlCSF = median(ControlIm(CSFmask & isfinite(ControlIm)));
    
    CSFcontrast = ControlCSF/ControlParenchyma;
    
    % 3) Get (robust?) spatial CoV of Control
    % ControlSmooth = xASL_im_ndnanfilter(ControlIm, 'gauss', [2 2 2]);
    ControlMedian = median(ControlIm(WBmask & isfinite(ControlIm)));
    ControlMAD = xASL_stat_MadNan(ControlIm(WBmask & isfinite(ControlIm)));
    Robust_sCoV = ControlMAD/ControlMedian;
    ControlMean = mean(ControlIm(WBmask & isfinite(ControlIm)))
    ControlSD = std(ControlIm(WBmask & isfinite(ControlIm)))
    ControlsCoV = ControlSD/ControlMean;
    
    % 4) Get spatial CoV of PWI
    PWI = xASL_io_Nifti2Im(x.P.Path_mean_PWI_Clipped);
    PWIMean = mean(PWI(WBmask & isfinite(PWI)));
    PWISD = std(PWI(WBmask & isfinite(PWI)));
    PWIsCoV = PWISD/PWIMean;
    
    dip_image(PWI)
    
    
end