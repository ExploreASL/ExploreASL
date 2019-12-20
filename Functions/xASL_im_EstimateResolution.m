function [ resFWHM, resSigma,resErr,imSmo,imMask ] = xASL_im_EstimateResolution(imCBF,imGM,imWM,imMaskOrig,PSFtype,maxIter)

% NB: everything in this code is in voxels, not in mm

%%%% Remove too high values
imCBF(imCBF<0) = 0;
imCBF(isnan(imCBF)) = 0;
imGM(imGM<0) = 0;
imGM(imGM>1) = 1;
imWM(imWM<0) = 0;
imWM(imWM>1) = 1;
imGM(isnan(imGM)) = 0;
imWM(isnan(imWM)) = 0;

imPV(:,:,:,1) = imGM;
imPV(:,:,:,2) = imWM;

%%%% Filter settings
tx = 10;
for tp=1:length(PSFtype)
    switch(PSFtype{tp})
        case 'gaussian'
            startSigmaVec = [1.2 1.2 3];
            txVec = 0.01;
        case 'lorentzian'
            startSigmaVec = [1.2 1.2 3];
            txVec = 0.01;
        case 'sinc'
            startSigmaVec = [0.4 0.4 0.25];
            txVec = 0.005;
            %    case 'gaussexp'
            %        startSigma = [1 1.5 1.5 1 1.5 1.5 1 3 3];
            %        tx = 0.01;
        case 'flat'
            startSigmaVec = [1.2 1.2 3];
            txVec = 0.01;
    end;
    tx = min(tx,txVec);
    optimSigma(tp) = startSigmaVec(tp);
    %opt = optimset('MaxIter',maxIter,'TolX',tx,'Display','iter');
    opt = optimset('MaxIter',maxIter,'TolX',tx,'Display','off');
end;

%%%% Create a mask where the CBF values are reasonable
% Restrict mask to GM and WM
imMask = ones(size(imGM));
imMask = (imMask.*imMaskOrig)>0;

% Run PV-correction with smooth kernel
[imPVEC,~,imResidual,~] = xASL_im_PVCbspline(imCBF,imPV,[11 11 11]);

% Remove from mask all voxels with too high residuals after PVEc
imErr = imResidual(imMask);
meanErr = mean(imErr);
stdErr = std(imErr);
imMask = imMask.*(imResidual < (meanErr + 2*stdErr));

%%%% Run several iterations of pseudoCBF intensity check and smoothing
for i=1:3
    % Re-run Bspline PVEC with presmoothed imPV according to the last
    % estimate
    [imGMSmo,~,~,~] = xASL_im_Smooth3D(optimSigma,imGM,PSFtype);
    [imWMSmo,~,~,~] = xASL_im_Smooth3D(optimSigma,imWM,PSFtype);
    imPVSmo(:,:,:,1) = imGMSmo;
    imPVSmo(:,:,:,2) = imWMSmo;
    [imPVEC,imCBFrec,imResidual,~] = xASL_im_PVCbspline(imCBF,imPVSmo,[9 9 9]);
    
    % Reset the pseudoCBF image according to locally adjusted GM/WM-CBF
    % contrast obtained from PVEc
    imPVEC(imPVEC<-20) = -20;
    imPVEC(imPVEC>200) = 200;
    imPseudoCBF = imGM.*imPVEC(:,:,:,1) + imWM.*imPVEC(:,:,:,2);
    
    imPseudoCBF(imPseudoCBF<0) = 0;
    imPseudoCBF(isnan(imPseudoCBF)) = 0;
    
    % The pseudoCBF image with the currently estimated GM/WM-CBF is
    % smoothed to match the CBF image
    minFuncLoc = @(sigma)minFunc(sigma,imPseudoCBF,imCBF,PSFtype,imMask);
    [optimSigma,optimErr] = fminsearch(minFuncLoc,optimSigma,opt);
end;

%%%% Calculate resolution, smoothed image and so on
[imSmo,imGauss{1},imGauss{2},imGauss{3}] = xASL_im_Smooth3D(optimSigma,imPseudoCBF,PSFtype);

imGauss{1} = imGauss{1}(:);
imGauss{2} = imGauss{2}(:);
imGauss{3} = imGauss{3}(:);
for tp = 1:length(PSFtype)
    switch(PSFtype{tp})
        case 'gaussian'
            resFWHM(tp) = sqrt((1./2.355).^2+(optimSigma(tp).^2))*2.355;
        case 'lorentzian'
            resFWHM(tp) = sum(imGauss{tp}>(0.5*max(imGauss{tp})));
        case 'sinc'
            resFWHM(tp) = sum(imGauss{tp}>(0.5*max(imGauss{tp})));
        case 'flat'
            resFWHM(tp) = sqrt((1./2.355).^2+(optimSigma(tp).^2))*2.355;
    end;
end;
resSigma = sqrt((1./2.355).^2+(optimSigma.^2));
resErr  = optimErr;

return;

function err = minFunc(sigma,imPseudoCBF,imCBF,PSFtype,imMask)
        
% Do the smoothing
[imSmo,~,~,~] = xASL_im_Smooth3D(sigma,imPseudoCBF,PSFtype);

% Restrict only to voxels on the mask
imSmoMasked = imSmo(imMask>0);
imCBFMasked = imCBF(imMask>0);

% Adjust the contrast for the best fitting using linear regression globally
% for all values
X = imSmoMasked(:);
X = [ones(length(X),1),X];
Y = imCBFMasked(:);
sol = inv(X'*X)*X'*Y;

% Calculate the error
err = mean(mean(mean((imSmoMasked*sol(2)+sol(1)-imCBFMasked).^2)));

return;
