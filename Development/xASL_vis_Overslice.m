function [ImageOut] = xASL_vis_Overslice(ForegroundIm, BackgroundIm, Mask, Overslice, TransparencyBackground, bWhite)
%xASL_vis_Overslice Overlap images for nice visualization

% Overslice from -1 to +1 (completely aside)
% [horizontal vertical] in fraction

% Overslice 0 will not move the images
% +1 will move the image up and to the right
% -1 down and to the left

%% -----------------------------------------------------------
%% Admin

% ForegroundIm = ThisIM;
% BackgroundIm = ThisIM;

% Mask = max(ForegroundIm==0 | isnan(ForegroundIm),3);
% Mask2 = ForegroundIm(:,:,1)==ForegroundIm(:,:,2) & ForegroundIm(:,:,2)==ForegroundIm(:,:,3);
% Mask = Mask & Mask2;

% TransparencyBackground = 0.75;
% Overslice = [0.5 0];

% if no mask, create one

ForegroundIm = xASL_io_Nifti2Im(ForegroundIm); % allow inserting path or image matrix
BackgroundIm = xASL_io_Nifti2Im(BackgroundIm); % allow inserting path or image matrix
ImageOut = NaN;

if ~isequal(size(ForegroundIm),size(BackgroundIm))
    warning('Input images need to have the same size, skipping');
    return;
elseif max(ForegroundIm(:))>1 || min(ForegroundIm(:))<0
    warning('Foreground image outside of colorscale range, skipping');
    return;
elseif max(BackgroundIm(:))>1 || min(BackgroundIm(:))<0
    warning('Background image outside of colorscale range, skipping');
    return;    
elseif ~isequal(size(ForegroundIm),size(Mask))
    warning('Mask needs to have the same size as images, skipping');
    return;
end
if nargin<6 || isempty(bWhite)
    bWhite = 0; % black background
end


%% -----------------------------------------------------------
%% If no mask specified, create it ourselves
if nargin>3 && ~isempty(Mask)
    if ndims(Mask)~=3
        warning('Incorrect Mask, skipping');
        return;
    elseif length(unique(Mask))~=2
        warning('Mask wasnt logical, skipping');
    else
        Mask = logical(Mask);
    end
else
    Mask = ForegroundIm(:,:,1)==ForegroundIm(:,:,2) & ForegroundIm(:,:,1)==ForegroundIm(:,:,3);
    if bWhite
        Mask2 = max(ForegroundIm==1, [], 3);
    else
        Mask2 = max(ForegroundIm==0, [], 3);
    end
    Mask = logical(repmat(Mask & Mask2, [1 1 3]));
end
        

%% -----------------------------------------------------------
%% Make BackgroundIm transparent
BackgroundIm = xASL_vis_Transparency(BackgroundIm, [], TransparencyBackground);

Overslice = flip(Overslice); % switch LR/UD
Overslice(2) = -Overslice(2); % switch LeftRight

%% Translate Overslide percentage to nVoxels
Dim = size(BackgroundIm);


%% Split in positive (up/right) & negative (down/left) components
OverslicePos = Overslice.*double(Overslice>0);
OversliceNeg = abs(Overslice).*double(Overslice<0);

PadVoxelsPos = round(Dim(1:2).*OverslicePos);
PadVoxelsNeg = round(Dim(1:2).*OversliceNeg);


%% Create new image dimensions
NewDim = Dim(1:2)+PadVoxelsPos+PadVoxelsNeg;
NewForegroundIm = zeros([NewDim 3]);
NewMask = logical(ones([NewDim 3]));
NewBackgroundIm = zeros([NewDim 3]);


%% Move foreground and background images
ForeIndex = [PadVoxelsNeg(1)+1 PadVoxelsNeg(1)+Dim(1) PadVoxelsNeg(2)+1 PadVoxelsNeg(2)+Dim(2)];
NewForegroundIm(ForeIndex(1):ForeIndex(2),ForeIndex(3):ForeIndex(4),:) = ForegroundIm;

NewMask(ForeIndex(1):ForeIndex(2),ForeIndex(3):ForeIndex(4),:) = Mask;

NewBackgroundIm(PadVoxelsPos(1)+1:end-PadVoxelsNeg(1), PadVoxelsPos(2)+1:end-PadVoxelsNeg(2),:) = BackgroundIm;


%% Merge images
ImageOut = NewForegroundIm;
ImageOut(NewMask) = NewBackgroundIm(NewMask);


end