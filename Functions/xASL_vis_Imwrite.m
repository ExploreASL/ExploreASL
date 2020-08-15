function [ImOut] = xASL_vis_Imwrite(ImIn, PathOut, ColorMap, bRescale)
% Saves an image to JPG while skipping the visual output to the screen.
%
% FORMAT:   [ImOut] = xASL_vis_Imwrite(ImIn, PathOut[, ColorMap, bRescale])
%
% INPUT:
%   ImIn    - input image, has to be a 2D (NX x NY) or a 3D matrix (NX x NY x 3),
%             otherwise it takes the first 2D image and issues a warning(REQUIRED)
%   PathOut - file path and name where the image will be saved (REQUIRED)
%   ColorMap- N by 3 matrix, specifying the RGB colors the image should
%             be saved in. N>=64, preferably N=256, N=2^M (OPTIONAL, default = gray colorscale)
%   bRescale - true to scale to maximum intensity (OPTIONAL, DEFAULT=true)
%             
% OUTPUT:
%   ImOut   - output image matrix, in Matlab colorscale (3 images for RGB) (OPTIONAL)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%   
% DESCRIPTION: This functions takes an input image matrix, interpolates it
%              to HD resolution (1920x1080) for visibility, and saves the image as jpg.
%              This function avoids the graphic interface of Matlab, for running from CLI
%              Careful: this function overwrites any existing PathOut.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%              
% EXAMPLE: xASL_vis_Imwrite([80x80x3 matrix], '/MyOutputFolder/MyFileName.jpg');
%          xASL_vis_Imwrite([80x80   matrix], '/MyOutputFolder/MyFileName.jpg', [], false); 
%
% __________________________________
% Copyright 2015-2020 ExploreASL
% 
% 2015-01-01 HJ

% Admin
    if nargin < 2
		error('xASL_vis_Imwrite: Needs at least two input parameters.');
	end
	
	% Only work for 2D files or files with XxYx3, otherwise reduce to 2D and issue a warning
	if ndims(ImIn) > 3
		warning('xASL_vis_Imwrite: The number of dimension must be <=3. Ignoring all above 3D.');
		ImIn = ImIn(:,:,:,1);
	end
	
	if ndims(ImIn) < 2
		error('xASL_vis_Imwrite: The number of dimension must be 2 or 3.');
	end
	
	if (ndims(ImIn) == 3) && (size(ImIn,3)~=3)
		warning('xASL_vis_Imwrite: The 3rd dimension can only have 1 or 3 components. Taking the first 2D image only.');
		ImIn = ImIn(:,:,1);
    end
    
    % manage extension, force jpg
    [Fpath, Ffile] = xASL_fileparts(PathOut);
    PathOut = fullfile(Fpath, [Ffile '.jpg']);

	if nargin<4 || isempty(bRescale)
		bRescale = true;
	end
	
    dimLow = size(ImIn);
    dimHD = [1080 1920];
    % Rescale to Full HD size screen (1920x1080)
    ScaleF = max(dimLow(1:2)./dimHD);

	iNewX = 0:ScaleF:(dimLow(2)-1);
	iNewY = 0:ScaleF:(dimLow(1)-1);
	% Allocate the matrix
	ImOut = zeros([length(iNewY),length(iNewX), size(ImIn,3)]);
	
	% apply the rescaling for multiple dimensions
	for iK=1:size(ImIn,3)
		tIM = ImIn(:,:,iK);
				
		[Xold, Yold] = meshgrid(0:     1:(dimLow(2)-1),0:     1:(dimLow(1)-1));
		[Xnew, Ynew] = meshgrid(iNewX,iNewY);
		ImOut(:,:,iK) = interp2(Xold,Yold,tIM,Xnew,Ynew,'cubic');
	end
	
	ImOut = ImOut+eps;
    
    if bRescale
        ImOut = ImOut./max(ImOut(:));
    end
    
    xASL_delete(PathOut); % delete if the output file already exists
    if nargin<3 || isempty(ColorMap)
        imwrite(ImOut, PathOut); % overwrite is implied, JPEG is implied
    else
        imwrite(ImOut, ColorMap, PathOut); % overwrite is implied, JPEG is implied
    end

end
