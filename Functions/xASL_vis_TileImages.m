function [ImOut] = xASL_vis_TileImages(ImIn, nColumns)
%xASL_vis_TileImages Merges selected slices (3D) into one single 2D picture.
%
% FORMAT:       ...
% 
% INPUT:        ImIn = input 3D volume, as stack of 2D slices
%               [nColumns = desired number of columns]
%
% OUTPUT:       ImOut = 2D output image, with all slices tiled
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Merges selected slices (3D) into one single 2D picture.
%               Plots all slices in one figure with specified rows and
%               columns, aiming for a square tile.
%
%               PM: can be extended to multiple slices
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

%% Admin
if nargin<1 || isempty(ImIn)
    error('No input image provided');
else
    Dim = size(ImIn);
    if length(Dim)<3 % no need to process
        ImOut = ImIn;
        return;
    elseif length(Dim)>3
        warning('Image has more dimensions than will be tiled');
        ImOut = ImIn;
        return;        
    end
end
    
if  nargin<2 || isempty(nColumns)
    % Try to make the image as square as possible
    nColumns = floor(Dim(3)^0.5);
end

nColumns = round(nColumns);  % make sure to have integers as column count                                                                                                
nRows = ceil(Dim(3)/nColumns); % number of rows=slices/columns



%% Tile slices
ImOut = zeros([Dim(1)*nRows Dim(2)*nColumns]);
iSlice = 1;
for iI=1:nRows
     for iJ=1:nColumns
         if iSlice<=Dim(3)
             ImOut(1+(iI-1)*Dim(1):iI*Dim(1), 1+(iJ-1)*Dim(2):iJ*Dim(2)) = ImIn(:,:,iSlice);
             iSlice = iSlice+1;
         end
     end
end


end