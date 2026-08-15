function ImageThrough = xASL_vis_AddColorBorder(ImageThrough, BorderColor, BorderWidth, bBlackOuterLine)
%xASL_vis_AddColorBorder Add a colored border to an RGB image
%
% FORMAT: ImageThrough = xASL_vis_AddColorBorder(ImageThrough, BorderColor[, BorderWidth, bBlackOuterLine])
%
% INPUT:
%   ImageThrough         - Input 3D image (2D image with 3rd dimension RGB) (REQUIRED)
%   BorderColor          - RGB triplet specifying the border color, with values
%                          between 0 and 1 (OPTIONAL, DEFAULT = [0 1 0], green)
%                          see examples below
%   BorderWidth          - Border width in pixels (OPTIONAL, DEFAULT = 2% of smallest dimension, minimal 2 lines)
%   bBlackOuterLine      - boolean for keeping the outermost line black
%                          which could facilitate visualizing concatenated images with
%                          differently colored borders (OPTIONAL, DEFAULT = true)
%
% OUTPUT:
%   ImOut                - Output image with the colored border added
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function adds a colored border to the inside of a 2D or RGB image. A
% grayscale input image is converted to RGB before adding the border.
%
% The RGB values supplied in BorderColor should be between 0 and 1. They
% are converted automatically to the numeric range and data type of ImIn.
% The dimensions of the output image are identical to those of the input image.
%
% This function can be used to indicate the status of visual quality
% control results. For example, a green border can indicate an accepted
% processing step and a red border can indicate a rejected processing
% step.
%
% EXAMPLE:
% % Add a four-pixel green border:
% ImOut = xASL_vis_AddColorBorder(ImIn, [0 1 0], 4);
%
% % Add a six-pixel red border:
% ImOut = xASL_vis_AddColorBorder(ImIn, [1 0 0], 6);
%
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________


if nargin < 4 || isempty(bBlackOuterLine)
    bBlackOuterLine = true;
end

if nargin < 3 || isempty(BorderWidth)
    minSize = min(size(ImageThrough(:,:,1)));
    BorderWidth = ceil(0.02*minSize); % take 2% of the smallest dimension
    BorderWidth = max(BorderWidth, 2); % take at least 2 lines
end

% Ensure that the image has an RGB dimension
if ndims(ImageThrough)==2
    ImageThrough = repmat(ImageThrough, [1 1 3]);
elseif ndims(ImageThrough)~=3 || size(ImageThrough, 3) ~= 3
    error('Image has to be two-dimensional or with RGB 3rd dimension');
end

% Numeric range: use either 0-1 range or 0-255 range
if isinteger(ImageThrough)
    BorderColor = cast(BorderColor .* double(intmax(class(ImageThrough))), ...
                       class(ImageThrough));
else
    ImageMaximum = max(ImageThrough(:));
    if ImageMaximum > 1
        BorderColor = BorderColor .* 255;
    end
end


% Draw each side of the border
for iColor = 1:3
    ImageThrough(1:BorderWidth, :, iColor) = BorderColor(iColor);
    ImageThrough(end-BorderWidth+1:end, :, iColor) = BorderColor(iColor);
    ImageThrough(:, 1:BorderWidth, iColor) = BorderColor(iColor);
    ImageThrough(:, end-BorderWidth+1:end, iColor) = BorderColor(iColor);
end

% Make the outermost line black
ImageThrough([1, end], :, :) = 0;
ImageThrough(:, [1, end], :) = 0;


end