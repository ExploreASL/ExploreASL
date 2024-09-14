function [AtlasOut] = xASL_vis_Convert4D_3D_atlas(AtlasIn)
    %xASL_vis_Convert4D_3D_atlas Converts 4D binary image to 3D atlas with ordinal integers
    %
    % FORMAT: [AtlasOut] = xASL_vis_Convert4D_3D_atlas(AtlasIn)
    %
    % INPUT:
    %   AtlasIn   - binary 4D matrix, concatenation of 3D binary masks, where
    %               fourth dim is equal to the mask number (REQUIRED)
    %
    % OUTPUT:
    %   AtlasOut  - 3D image matrix with ordinal integers where each integer
    %               defines an ROI mask, equal to the number of the fourth dim of the input
    %               atlas
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % DESCRIPTION: This function converts 4D binary image to 3D atlas with ordinal integers
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % EXAMPLE: AtlasOut = xASL_vis_Convert4D_3D_atlas(AtlasIn);
    % __________________________________
    % Copyright 2017-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    
    AtlasOut = zeros(size(AtlasIn(:,:,:,1,1,1)),'uint8');
    for iL=1:size(AtlasIn,4)
        AtlasOut(logical(AtlasIn(:,:,:,iL))) = iL;
    end
end