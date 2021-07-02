function [DiceCoeff] = xASL_im_ComputeDice(imA, imB)
%xASL_im_ComputeDice Calculate Dice coefficient of image overlap
%
% FORMAT:       [DiceCoeff] = xASL_im_ComputeDice(imA, imB)
% 
% INPUT:
%   imA         - Image matrix A (numerical or logical, REQUIRED)
%   imB         - Image matrix B (numerical or logical, REQUIRED)
%
% OUTPUT:
%    DiceCoeff  - Dice coefficient, 0 = no image overlap, 1 = perfect image overlap.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  This function calculates the Dice coefficient of image overlap.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: DiceCoeff = xASL_im_ComputeDice([1 0 1], [1 1 0]); % mask has
%          two voxels, of which 1 voxel overlaps, so DiceCoeff=0.5
% __________________________________
% Copyright 2015-2021 ExploreASL


%% 0. Admin
if nargin<2
    error('Input image(s) missing');
end    
if ~islogical(imA) && isnumeric(imA)
    warning('Input image A was not dichotomous, dichotomizing now');
    imA = imA>0;
elseif ~islogical(imA) && ~isnumeric(imA)
    error('Input image A has the wrong format');
end
if ~islogical(imB)  && isnumeric(imB)
    warning('Input image B was not dichotomous, dichotomizing now');
    imB = imB>0;
elseif ~islogical(imB) && ~isnumeric(imB)
    error('Input image B has the wrong format');
end

if ~isequal(size(imA), size(imB))
    error('Input image matrices A & B must have the same size');
end

%% 1. Calculate Dice coefficient
if sum(imA(:))==0 && sum(imB(:))==0
    warning('Two empty maps, DICE overlap 1');
    DiceCoeff = 1;
else
    DiceCoeff = (2 .* sum(sum(sum(imA & imB)))) ./ (sum(imA(:)) + sum(imB(:)) );
end


end

