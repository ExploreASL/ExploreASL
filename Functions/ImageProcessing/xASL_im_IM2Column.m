function [ColumnOut] = xASL_im_IM2Column(ImageIn, BrainMask, ApplyShiftDim)
%xASL_im_IM2Column QC Converts an image matrix to a single-dimensional
%column to save memory space & computation time
%
% FORMAT: [ColumnOut] = xASL_im_IM2Column(ImageIn, BrainMask[, ApplyShiftDim])
%
% INPUT:
%   ImageIn         - image matrix (MATRIX, REQUIRED)
%   BrainMask       - image matrix containing mask (MATRIX, REQUIRED)
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
%              See also: `xASL_im_Column2IM.m`
% 
%              The mask mostly used for xASL_im_IM2Column is `x.S.masks.WBmask`, which completely
%              engulfes pGM, pWM & pCSF
%
%              3D example: an [121 145 121] image becomes [2122945 1] (a
%              single column with 121x145x121 values), but then compressed
%              by the 3D mask.
%
%              4D example: ten images [121 145 121 10] become [2122945 10],
%              but then compressed by the 3D mask.
%
%              So e.g., if the BrainMask has 500000 voxels to include, the
%              above examples would become [500000 1] for the 3D example
%              and [500000 10] for the 4D example
%
% EXAMPLE: ColumnOut = xASL_im_IM2Column(ImageIn, BrainMask);
% __________________________________
% Copyright (C) 2015-2021 ExploreASL


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
    try
        ColumnOut = ImageIn(logical(BrainMask));
    catch ME
        warning(['Image with matrix size [' xASL_num2str(size(ImageIn)) '] instead of [121 145 121]']);
        error('%s\n', ME.message);
    end
end

if ApplyShiftDim
    ColumnOut = shiftdim(ColumnOut,1);
end


end