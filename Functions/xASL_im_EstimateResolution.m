function [resFWHM, resSigma, resErr, imSmo, imMask] = xASL_im_EstimateResolution(imCBF, imGM, imWM, imMaskOrig, PSFtype, maxIter)
%xASL_im_EstimateResolution Estimates resolution of a CBF image by iterative matching with a pseudo-CBF based on GM and WM masks
%
% FORMAT:       [resFWHM, resSigma, resErr, imSmo, imMask] = xASL_im_EstimateResolution(imCBF, imGM, imWM[, imMaskOrig, PSFtype, maxIter])
% 
% INPUT: imCBF      - Input CBF image for which it estimates the resolution (REQUIRED)
%        imGM       - Input GM segmentation map with the same matrix and voxel size, but at higher resolution (REQUIRED)
%        imWM       - Input WM segmentation map with the same matrix and voxel size, but at higher resolution (REQUIRED)
%        imMaskOrig - A binary mask that restricts the calculation of the matching between pseudoCBF and CBF (OPTIONAL, DEFAULT GM+WM>0.1)
%        PSFtype    - The shape of the point-spread-function used for the resolution estimation. Cell of three char arrays is expected,
%                     the options are 'gaussian','lorentzian','sinc','flat' (OPTIONAL, DEFAULT = {'gaussian' 'gaussian' 'gaussian'})
%        maxIter    - Maximum number of iterations run (OPTIONAL, DEFAULT = 5)
%
% OUTPUT: resFWHM   - Full-width at half maximum of the estimated PSF of the imCBF - in voxels
%         resSigma  - SD of the estimated PSF of the imCBF - in voxels
%         resErr    - The minimal error reacher between the CBF and smoothed pseudoCBF
%         imSmo     - PseudoCBF image smoothed to the CBF resolution
%         imMask    - The final mask upon which the mismatch was evaluated
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Creates a high-resolution pseudo-CBF image based on segmented GM and WM maps and iteratively adjusts its resolution 
%               by smoothing until reaching a perfect fit with the CBF image thus obtaining the resolution difference between the GM and CBF image and
%               uses this to calculate the estimated effective resolution of hte CBF. Note that all the calculations are done using voxels as measures
%               and not mm, so the output resolution is also in voxels and has to be transfered to mm by using the knowledge about the voxel size.
%               It is assumed that for imGM and imWM, the voxel size equals the resolution, and the imCBF is upsampled to the smaller voxels of imGM.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [resFWHM, resSigma, resErr, imSmo, imMask] = xASL_im_EstimateResolution(imCBF, imGM, imWM, [], {'gaussian' 'gaussian' 'gaussian'}, 5);
% __________________________________
% Copyright 2015-2021 ExploreASL

%% Parameter admin
if nargin < 3 || isempty(imCBF) || isempty(imGM) || isempty(imWM)
	error('Need to provide imCBF, imGM, and imWM');
end

if nargin < 4 || isempty(imMaskOrig)
	imMaskOrig = (imGM+imWM)>0.1;
end

if nargin < 5 || isempty(PSFtype)
	PSFtype = {'gaussian' 'gaussian' 'gaussian'};
end

if nargin < 6 || isempty(maxIter)
	maxIter = 5;
end

%% Clean data 
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

%% Set the shape of the PSF function
tx = 10;
for tp=1:length(PSFtype)
	switch (PSFtype{tp})
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
	end
    tx = min(tx,txVec);
    optimSigma(tp) = startSigmaVec(tp);
    %opt = optimset('MaxIter',maxIter,'TolX',tx,'Display','iter');
    opt = optimset('MaxIter',maxIter,'TolX',tx,'Display','off');
end

%% Prepare the estimation
% Create a mask where the CBF values are reasonable
imMask = ones(size(imGM));
imMask = (imMask.*imMaskOrig)>0;

% Run PV-correction with smooth kernel
[~,~,imResidual,~] = xASL_im_PVCbspline(imCBF,imPV,[11 11 11]);

% Remove from mask all voxels with too high residuals after PVEc
imErr = imResidual(imMask);
meanErr = mean(imErr);
stdErr = std(imErr);
imMask = imMask.*(imResidual < (meanErr + 2*stdErr));

%% Run several iterations of pseudoCBF intensity check and smoothing
for i=1:3
    % Re-run Bspline PVEC with presmoothed imPV according to the last
    % estimate
    [imGMSmo,~,~,~] = xASL_im_Smooth3D(optimSigma,imGM,PSFtype);
    [imWMSmo,~,~,~] = xASL_im_Smooth3D(optimSigma,imWM,PSFtype);
    imPVSmo(:,:,:,1) = imGMSmo;
    imPVSmo(:,:,:,2) = imWMSmo;
    [imPVEC,~,~,~] = xASL_im_PVCbspline(imCBF,imPVSmo,[9 9 9]);
    
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
end

%% Calculate resolution, smoothed image and so on
[imSmo,imGauss{1},imGauss{2},imGauss{3}] = xASL_im_Smooth3D(optimSigma,imPseudoCBF,PSFtype);

imGauss{1} = imGauss{1}(:);
imGauss{2} = imGauss{2}(:);
imGauss{3} = imGauss{3}(:);

resFWHM = zeros(1,length(PSFType));
for tp = 1:length(PSFtype)
	switch (PSFtype{tp})
        case {'gaussian','flat'}
            resFWHM(tp) = sqrt((1./2.355).^2+(optimSigma(tp).^2))*2.355;
        case 'lorentzian'
            resFWHM(tp) = sum(imGauss{tp}>(0.5*max(imGauss{tp})));
        case 'sinc'
            resFWHM(tp) = sum(imGauss{tp}>(0.5*max(imGauss{tp})));
	end
end
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
