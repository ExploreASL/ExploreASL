function ImageOut = xASL_vis_CropParmsApply(ImageIn,CropParameters,Xmax,Ymin,Ymax)
%xASL_vis_CropParmsApply Crops 2D image data
%
% FORMAT: ImageOut = xASL_vis_CropParmsApply(ImageIn,CropParameters)
%
% INPUT:
%   ImageIn        - single or series of 2D image matrix/matrices that need to be cropped (REQUIRED)
%                    series of images not tested
%   CropParameters - vector of [Xmin Xmax Ymin Ymax] (REQUIRED)
%
% OUTPUT:
%   ImageOut       - cropped image matrix/matrices
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function crops 2D image matrices
%
% EXAMPLE to crop a [121 145] MNI 1.5 mm slice to a [121 121] square:
%         ImageOut = xASL_vis_CropParmsApply(ImageIn,[12 133 0 121])
% __________________________________
% Copyright 2015-2019 ExploreASL


%% Manage crop parameters

Xmin = CropParameters(1);
if nargin<3 && length(CropParameters)==4
    Xmax = CropParameters(2);
    Ymin = CropParameters(3);
    Ymax = CropParameters(4);
end    
    
    
%% Get image size
dim = size(ImageIn);

%% Add space if requested dimensions are outside dimensions of original image
Add_Xmin    = max(0,1-Xmin); % default is 0, if Xmin is smaller than 1, add this
Xmin        = max(1,Xmin);   % set Xmin to 1, if it  is smaller than 1
Add_Ymin    = max(0,1-Ymin); 
Ymin        = max(1,Ymin);   

Add_Xmax    = max(0,Xmax-dim(1)); % if Xmax is higher than dim(1), add this. default=0
Xmax        = min(Xmax,dim(1));   % if Xmax is higher than dim(1), set to dim(1), otherwise leave at Xmax
Add_Ymax    = max(0,Ymax-dim(2));
Ymax        = min(Ymax,dim(2));

dim(3)      = max(1,size(ImageIn,3)); % if no 3D dim, 1 by default. don't replace this by dim(3), this will crash when length(size())<3
dim(4)      = max(1,size(ImageIn,4)); % if no 4D dim, 1 by default

% Recreate dims
Dim1        = dim(1)+Add_Xmin+Add_Xmax;
Dim2        = dim(2)+Add_Ymin+Add_Ymax;

%% Process the image
ImageOut = zeros(Dim1,Dim2,dim(3),size(ImageIn,4),size(ImageIn,5),size(ImageIn,6),size(ImageIn,7)); % create new image with zeros, including added borders
ImageOut(1+Add_Xmin:dim(1)+Add_Xmin,1+Add_Ymin:dim(2)+Add_Ymin,1:dim(3),:,:,:,:) = ImageIn; % put image into new image, on new position within border
ImageOut = ImageOut(Xmin:Xmax+Add_Xmin+Add_Xmax,Ymin:Ymax+Add_Ymin+Add_Ymax,1:dim(3),:,:,:,:); % crop new image




end

