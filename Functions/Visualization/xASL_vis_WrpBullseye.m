function  x = xASL_vis_WrpBullseye(x)
    %xASL_vis_WrpBullseye uses xASL_vis_Bullseye to make a bullseye plot in subject directory.
    %
    % FORMAT:       xASL_vis_Bullseye(x);
    % 
    % INPUT:        x - x struct (REQUIRED)
    %
    % OUTPUT:       x - x struct with the path of x.P.Path_Bullseye_CBF added.
    % 
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % DESCRIPTION:  Generates a bullseye plot.
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % EXAMPLE:      ...
    % __________________________________
    % Copyright 2015-2023 ExploreASL
    
    if nargin <1 || isempty(x)
        warning('x Struct not provided in Bullseye wrapper')
    end

    segmentCount = 9;
    layerCount = 4;
    data = rand(segmentCount, layerCount);

    if ~isfield(x.P, 'Path_Bullseye_CBF')
        x.P.Path_Bullseye_CBF = fullfile(x.D.ROOT, x.SUBJECT, 'Bullseye.png')
    else
        warning('Bullseye plot already exists, overwriting')
    end

    xASL_vis_Bullseye(x.P.Path_Bullseye_CBF, data, [], [], x.S.hot)

end