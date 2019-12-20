function [ColumnOut] = xASL_im_IM2Column(ImageIn, BrainMask, ApplyShiftDim)
%xASL_im_IM2Column QC Converts an image matrix to a single-dimensional
%column to save memory space & computation time
%
% FORMAT: [ColumnOut] = xASL_im_IM2Column(ImageIn, BrainMask[, ApplyShiftDim])
%
% INPUT:
%   ImageIn         - image matrix (REQUIRED)
%   BrainMask       - image matrix containing mask (REQUIRED)
%   ApplyShiftDim   - false to avoid dimension shifting (OPTIONAL, DEFAULT=true)
%   
% OUTPUT:
%    ColumnOut      - column containing "compressed" image matrix
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function "compresses" an image matrix (or multiple matrices)
%              for optimization of memory and CPU resources. The output column only includes 
%              voxels that lie within the BrainMask. This excludes extracranial
%              zero-information voxels from computations and memory use.
%              NB: Important to use the same BrainMask for converting the
%              column back to an image matrix!
%              See also: xASL_im_Column2IM.m
% 
%              The mask mostly used for xASL_im_IM2Column is x.WBmask, which completely
%              engulfes pGM, pWM & pCSF
%
% EXAMPLE: ColumnOut = xASL_im_IM2Column(ImageIn, BrainMask);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL


%%  ------------------------------------------------------------
%% Shift dims, if multiple subjects
if nargin<3 || isempty(ApplyShiftDim)
    ApplyShiftDim = true;
end

if ApplyShiftDim && length(size(ImageIn))>3
    ImageIn = shiftdim(ImageIn,length(size(ImageIn))-3);
else
    ApplyShiftDim = false; % to avoid re-shifting below
end

SizeIM = size(ImageIn);

%%  ------------------------------------------------------------
%%  Apply dimension change (compression, IM2column)
if size(ImageIn,3)>1 && min(SizeIM(1:3)==[121 145 121]) % here we apply iteration, to avoid using too much memory
    for iS=1:size(ImageIn,4)
        for iY=1:size(ImageIn,5)
            IMtemp = ImageIn(:,:,:,iS,iY);
            ColumnOut(:,iS,iY) = IMtemp(logical(BrainMask));
        end
    end
else
    ColumnOut = ImageIn;
end

if ApplyShiftDim
    ColumnOut = shiftdim(ColumnOut,1);
end


end