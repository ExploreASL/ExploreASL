function  xASL_vis_Bullseye(x, savePath)
    %xASL_vis_Bullseye uses plotting to make a bullseye plot.
    %
    % FORMAT:       xASL_vis_Bullseye(x,);
    % 
    % INPUT:        x - x struct (REQUIRED)
    %               savePath - char path to save file (OPTIONAL)
    %
    % OUTPUT:       ...
    % 
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % DESCRIPTION:  Generates a bullseye plot.
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % EXAMPLE:      ...
    % __________________________________
    % Copyright 2015-2023 ExploreASL
    
%% Admin
if nargin<2
    savePath = 'Bullseye';
end
% Test Values

% Define the number of segments
segmentCount = 9;
layerCount = 4;

% Generate n random values between 0 and 1
values = rand(segmentCount, layerCount);

% Actual function
% Define the segment labels
labels = {'Frontal', 'Parietal', 'Temporal', 'Occipital', 'SubCortical', 'Parietal', 'Temporal', 'Occipital', 'Frontal'};

% Define the segment boundaries in degrees
ANGLE_OFFSET = 90;
segment_boundaries = linspace(0, 360, segmentCount+1) + ANGLE_OFFSET;

% Initialize Figure
fig = figure;
colormap('autumn');
title('Cardiac Bullseye Plot');
colorbar('Ticks', 1:length(values), 'TickLabels', labels);
axis off;

% Plot the segments with different colors
for iLayer = 1:layerCount
    for iSegment = 1:segmentCount
        angle1 = segment_boundaries(iSegment);
        angle2 = segment_boundaries(iSegment+1);
        depth2 = iLayer/(layerCount+1);
        depth1 = (iLayer+1)/(layerCount+1);
        vertexX = [depth1*[cosd(angle2), cosd(angle1)], depth2*[cosd(angle1), cosd(angle2)]];  
        vertexY = [depth1*[sind(angle2), sind(angle1)], depth2*[sind(angle1), sind(angle2)]];
        patch(vertexX, vertexY, values(iSegment,iLayer), 'EdgeColor', 'k');
    end
end

% Add the segment labels
for iSegment = 1:segmentCount
    angle = mean([segment_boundaries(iSegment), segment_boundaries(iSegment+1)]);
    text(cosd(angle), sind(angle), labels{iSegment}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

text(0, 0, 'VALUES ARE GENERATED RANDOMLY, NOT ACTUAL VALUES', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
print(fig, savePath, '-dpng');
close all
clear
    
    
end