function  xASL_vis_Bullseye(savePath, dataMatrix, dataLabels, dataTitle, dataColorMap)
    %xASL_vis_Bullseye uses plotting to make a bullseye plot.
    %
    % FORMAT:       xASL_vis_Bullseye(savepath, dataMatrix[, dataLabels, DataTitle, dataColorMap)];
    % 
    % INPUT:        savePath - char path to save file (REQUIRED)
    %               dataMatrix - 2D matrix with values to plot (REQUIRED)
    %               dataLabels - cell array with labels for each segment (OPTIONAL - default is 9 segments of brain regions if a 9 by N matrix is provided)
    %               dataTitle - char title for the plot (OPTIONAL - default is 'Coronal Bullseye Plot')
    %               dataColorMap - colormap to use for the plot (OPTIONAL - default is autumn colormap)
    %
    % OUTPUT:       None (saves a file at savePath)
    % 
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % DESCRIPTION:  Generates a bullseye plot from the dataMatrix. The dataMatrix should be 2D matrix with the first dimension corresponding 
    %               to the segments and the second dimension corresponding to the layers. The plot will be saved at savePath.
    %               The dataLabels is a cell array with labels for each segment.
    %               The plot will be saved at savePath.
    %
    % -----------------------------------------------------------------------------------------------------------------------------------------------------
    % EXAMPLE:      xASL_vis_Bullseye('Bullseye.png', meshgrid(0.1:0.1:0.4, 0.1:0.1:0.9), {'Frontal', 'Parietal', 'Temporal', 'Occipital', 'SubCortical', 'Occipital', 'Temporal', 'Parietal', 'Frontal'}, 'This is an Expample of 9 segments with 4 layers', flipud(autumn));
    %               xASL_vis_Bullseye('Bullseye.png', meshgrid(0.1:0.1:0.4, 0.1:0.1:0.9));
    %               xASL_vis_Bullseye('Bullseye.png', meshgrid(0.25:0.25:1, 0.25:0.25:1), {'segment1', 'segment2', 'segment3', 'segment4'}, 'This is an Expample of 4 segments with 4 layers', flipud(autumn));
    %               
    % __________________________________
    % Copyright 2015-2024 ExploreASL
    
%% Admin
if nargin<1 || isempty(savePath)
% Default savedirectory if not provided is in the current directory
    savePath = 'Bullseye.png';
end

if nargin<2 || isempty(dataMatrix)
    warning('No Values provided, Bullseye plot not performed. Please provide a 2D dataMatrix.')
    return
end

[segmentCount, layerCount, ~ ] = size(dataMatrix);

if nargin<3 || isempty(dataLabels)
    % Default labels if dataMatrix has 9 segments
    if segmentCount == 9
        dataLabels = {'Frontal', 'Parietal', 'Temporal', 'Occipital', 'SubCortical', 'Occipital', 'Temporal', 'Parietal', 'Frontal'};
    else
        dataLabels = {};
    end
end

if nargin<4 || isempty(dataTitle)
    dataTitle = 'Bullseye Plot';
end

if nargin<4 || isempty(dataColorMap)
    dataColorMap= flipud(autumn);
end



% Define the segment boundaries in degrees
ANGLE_OFFSET = 90;
segmentBoundaryAngles = linspace(0, 360, segmentCount+1) + ANGLE_OFFSET;

% Initialize Figure
fig = figure('visible', 'off');
colormap(dataColorMap);
title(dataTitle);
colorbar();
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
    if ~isempty(dataLabels)
        text(cosd(angle), sind(angle), dataLabels{iSegment}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'rotation', xASL_sub_textAngle(angle));
    end
end

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
    