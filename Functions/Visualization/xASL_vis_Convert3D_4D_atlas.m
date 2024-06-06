function [AtlasOut] = xASL_vis_Convert3D_4D_atlas(AtlasIn)
%xASL_vis_Convert3D_4D_atlas Converts 3D atlas with ordinal integers to 4D binary image
%
% FORMAT: [AtlasOut] = xASL_vis_Convert3D_4D_atlas(AtlasIn)
%
% INPUT:
%   AtlasIn   - 3D image matrix with ordinal integers where each integer defines an ROI mask (REQUIRED)
%
% OUTPUT:
%   AtlasOut  - binary 4D matrix, concatenation of 3D binary masks, where
%               fourth dim is equal to the mask integer of AtlasIn
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function converts 3D atlas with ordinal integers to 4D binary image
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: AtlasOut = xASL_vis_Convert3D_4D_atlas(AtlasIn);
% __________________________________
% Copyright 2017-2024 ExploreASL
    
    AtlasOut = zeros([size(AtlasIn(:,:,:,1,1,1)) max(AtlasIn(:))],'uint8');
    for iL=1:max(AtlasIn(:))
        tempIM = zeros(size(AtlasIn(:,:,:,1,1,1)));
        tempIM(AtlasIn==iL) = 1;
        AtlasOut(:,:,:,iL) = tempIM;
    end
    
end