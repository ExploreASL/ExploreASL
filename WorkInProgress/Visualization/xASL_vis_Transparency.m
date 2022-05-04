function [ImageOut] = xASL_vis_Transparency(ForegroundImage, BackgroundImage, Transparency, bWhite)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Make foreground transparent over background
% If no background image provided, this will be a black image
% Optional -> transparancy to white


if nargin<4 || isempty(bWhite)
    bWhite= false; % black background by default
end
if nargin<3 || isempty(Transparency)
    Transparency = 0.5;
end
if nargin<2 || isempty(BackgroundImage)
    BackgroundImage = zeros(size(ForegroundImage));
end

if ~isequal(size(ForegroundImage),size(BackgroundImage))
    warning('Input images need to have the same size, skipping');
    return;
elseif max(ForegroundImage(:))>1 || min(ForegroundImage(:))<0
    warning('Foreground image outside of colorscale range, skipping');
    return;
elseif max(BackgroundImage(:))>1 || min(BackgroundImage(:))<0
    warning('Background image outside of colorscale range, skipping');
    return;
end
    

% 
% ImageOut = BackgroundIm;

DeltaImage = ForegroundImage - BackgroundImage;

% ImageOut(~Mask) = BackgroundIm.*Transparency + DeltaImage.*Transparency;
ImageOut = BackgroundImage + DeltaImage.*Transparency;
% EmptyIm = zeros(size(DeltaImage));
% FullIm = ones(size(DeltaImage));
% ImageOut = max(ImageOut, EmptyIm);
% ImageOut = min(ImageOut, FullIm);

% figure(1);imshow(ImageOut, 'InitialMagnification',500)


end

