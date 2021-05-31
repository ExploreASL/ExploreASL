function [ImageOut] = xASL_im_Column2IM(ColumnIn, BrainMask)
%xASL_im_IM2Column QC Converts an single-dimensional column back to an
%image matrix
%
% FORMAT: [ImageOut] = xASL_im_Column2IM(ColumnIn, BrainMask)
%
% INPUT:
%   ColumnIn        - image matrix (REQUIRED)
%   BrainMask       - image matrix containing mask (REQUIRED)
%   
% OUTPUT:
%    ImageOut       - column containing "compressed" image matrix
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function "decompresses" an image matrix (or multiple matrices)
%              from a single-dimensional column, by reconstructing the image matrix 
%              from the voxel positions within the BrainMask. 
%              NB: Important to use the same BrainMask as used for converting the
%              image matrix to the column!
%              See also: xASL_im_IM2Column.m
% 
%              The mask mostly used for xASL_im_IM2Column is x.S.masks.WBmask, which completely
%              engulfes pGM, pWM & pCSF.
%
% EXAMPLE: ImageOut = xASL_im_Column2IM(ColumnIn, BrainMask);
% __________________________________
% Copyright (C) 2015-2020 ExploreASL


BrainMask = logical(BrainMask);
ClassIM = class(ColumnIn);

if islogical(ColumnIn) % this is required for early versions of Matlab
    ClassIM = 'uint8';
    Convert2Logical = true;
else
    Convert2Logical = false;
end


%%  ------------------------------------------------------------
%% Shift dims, if multiple subjects
if size(ColumnIn,2)/size(ColumnIn,1)>50
    ApplyShiftDim = true;
    ColumnIn = shiftdim(ColumnIn,1);
else
    ApplyShiftDim = false;
end


%%  ------------------------------------------------------------
%%  Apply the dimension change (i.e. decompression)
if size(ColumnIn,1)==sum(BrainMask(:)) % otherwise, we could have the incorrect mask
    Sz = size(BrainMask);
    ImageOut = zeros([Sz size(ColumnIn,2) size(ColumnIn,3)],ClassIM);
    for iS=1:size(ColumnIn,2)
        for iY=1:size(ColumnIn,3)
            IMtemp = zeros(Sz,ClassIM);
            IMtemp(BrainMask) = ColumnIn(:,iS,iY);
            ImageOut(:,:,:,iS,iY) = IMtemp;
        end
    end
else
    ImageOut = ColumnIn;
end

if ApplyShiftDim && length(size(ImageOut))~=3
    ImageOut = shiftdim(ImageOut,3); % this is always 3, because of the image size 3
end

if Convert2Logical % this is required for early versions of Matlab
    ImageOut = logical(ImageOut);
end

end

