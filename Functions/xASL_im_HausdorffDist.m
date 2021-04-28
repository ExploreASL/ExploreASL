function [hausDist,modHausDist] = xASL_im_HausdorffDist(imIn1,imIn2)
%xASL_im_HausdorffDist(imIn1,imIn2) Calculate Hausdorff distance between two binary ROIs
%
% FORMAT: xASL_im_HausdorffDist(imIn1,imIn2)
% 
% INPUT:
%   imIn1       - First ROI (REQUIRED)
%   imIn2       - Second ROI (REQUIRED)
%
% OUTPUT:
%   hausDist    - Hausdorff distance
%   modHausDist - Modified Hasudorff distance
%                         
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Calculate Hausdorff and modified Hausdorff distance between two ROIs in volumes imIn1, imIn2.
%              Input images are binarized as 0 and non-0
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: 
% [hausDist,modHausDist] = xASL_im_HausdorffDist(imIn1,imIn2)
%
% __________________________________
% Copyright 2015-2021 ExploreASL

%% First does an erosion to obtain the border of the two ROIs
% Create binary
imIn1 = double(imIn1>0);
imIn2 = double(imIn2>0);

% The smallest sphere element
el = xASL_im_DilateErodeSphere(1);

% Extract the inner border of the ROIs by subtraction of the ROI with an eroded ROI
imIn1er = double(xASL_im_DilateErodeFull(imIn1,'erode',el)>0);
imIn1er = imIn1-imIn1er;

imIn2er = double(xASL_im_DilateErodeFull(imIn2,'erode',el)>0);
imIn2er = imIn2-imIn2er;

%% Crop the images to increase efficiency of the calculation
[x1,y1,z1] = ind2sub(size(imIn1er),find(imIn1er>0));
[x2,y2,z2] = ind2sub(size(imIn2er),find(imIn2er>0));

% Find the min and max coordinates for cropping
xMin = max(1,min([x1;x2]));
yMin = max(1,min([y1;y2]));
zMin = max(1,min([z1;z2]));

xMax = min(size(imIn1er,1),max([x1;x2]));
yMax = min(size(imIn1er,2),max([y1;y2]));
zMax = min(size(imIn1er,3),max([z1;z2]));

% Perform the cropping
imIn1er = imIn1er(xMin:xMax,yMin:yMax,zMin:zMax);
imIn2er = imIn2er(xMin:xMax,yMin:yMax,zMin:zMax);

%% Calculate the Hausdorff distance using the Euclidean distance from the ROI

% Calculate the Euclidean distance map from both borders
[distMap1,~,~,~] = xASL_im_DistanceTransform(imIn1er);
[distMap2,~,~,~] = xASL_im_DistanceTransform(imIn2er);

distList1 = distMap1(imIn2er>0);
distList2 = distMap2(imIn1er>0);

% Calculate Hausdorff and modified Hausdorff distance
hausDist = max(max(distList1),max(distList2));
modHausDist = max(mean(distList1),mean(distList2));

end
