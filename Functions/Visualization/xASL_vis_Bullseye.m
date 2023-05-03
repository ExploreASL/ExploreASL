function  xASL_vis_Bullseye(savePath, dataMatrix, dataLabels)
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
if nargin<1 || isempty(savePath)
    savePath = 'Bullseye';
end

if nargin<2 || isempty(dataMatrix)
    segmentCount = 9;
    layerCount = 4;
    dataMatrix = rand(segmentCount, layerCount);
end

if nargin<3 || isempty(dataLabels)
    dataLabels = {'Frontal', 'Parietal', 'Temporal', 'Occipital', 'SubCortical', 'Parietal', 'Temporal', 'Occipital', 'Frontal'};
end

[segmentCount, layerCount, ~ ] = size(dataMatrix);

% Define the segment boundaries in degrees
ANGLE_OFFSET = 90;
segmentBoundaryAngles = linspace(0, 360, segmentCount+1) + ANGLE_OFFSET;

% Initialize Figure
fig = figure;
colormap('autumn');
title('Coronal Bullseye Plot');
colorbar('Ticks', 1:length(dataMatrix), 'TickLabels', dataLabels);
axis off;

% Plot the segments with different colors
for iLayer = 1:layerCount
    for iSegment = 1:segmentCount
        angle1 = segmentBoundaryAngles(iSegment);
        angle2 = segmentBoundaryAngles(iSegment+1);
        depth2 = iLayer/(layerCount+1);
        depth1 = (iLayer+1)/(layerCount+1);
        vertexX = [depth1*[cosd(angle2), cosd(angle1)], depth2*[cosd(angle1), cosd(angle2)]];  
        vertexY = [depth1*[sind(angle2), sind(angle1)], depth2*[sind(angle1), sind(angle2)]];
        patch(vertexX, vertexY, dataMatrix(iSegment,iLayer), 'EdgeColor', 'k');
    end
end

% Add the segment labels
for iSegment = 1:segmentCount
    angle = mean([segmentBoundaryAngles(iSegment), segmentBoundaryAngles(iSegment+1)]);
    
    text(cosd(angle), sind(angle), dataLabels{iSegment}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'rotation', xASL_sub_textAngle(angle));
end

text(0, 0, 'VALUES ARE GENERATED RANDOMLY, NOT ACTUAL VALUES', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
print(fig, savePath, '-dpng');
close all
clear
    
    
end

function output = xASL_sub_textAngle(inputAngle)
    % Quick function to ensure the text in the plot is always facedown.
    output = inputAngle + 90;
    while output > 90 
        output = output - 180;
    end
end
    