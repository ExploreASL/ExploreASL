% Copyright 2015-2024 ExploreASL (Works In Progress code)
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

function xASL_stat_Plot2Dsurface( IM2D, FilePath )
% Plot2Dsurface Visualizes input 2D image in 3D surface plot

    %% Check 2D
    if length(size(IM2D))~=2
        error('Dimensions of input image should be 2!!!');
    end

    IM2D    = double(IM2D);

    %% View in meshplot & print

    [X Y]           = meshgrid([-(size(IM2D,2)-1)/2:(size(IM2D,2)-1)/2],[-(size(IM2D,1)-1)/2:(size(IM2D,1)-1)/2]);
    Z               = IM2D+eps; % add very small number to avoid top clipping

    if      exist('FilePath','var') % don't show when printing
            fig = figure('Visible','off');
    else    fig             = figure;
    end

    jet_1024        = jet(1024);
%     jet_1024(1,:)   = 0;

    colormap(jet_1024)
    surf(X,Y,Z,'EdgeColor','none');
    camlight left
    lighting phong

    if exist('FilePath','var')
        saveas( fig, FilePath, 'jpg');
    end

end

