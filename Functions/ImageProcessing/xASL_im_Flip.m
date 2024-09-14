function [MatrixOut] = xASL_im_Flip(MatrixIn, varargin)
%xASL_im_Flip Backwards compatibility for flipping left-right in standard
%space (NB: this can be different than in native space!)
%
% FORMAT:       [MatrixOut] = xASL_im_Flip(MatrixIn, varargin)
% 
% INPUT:        MatrixIn    - Input matrix (MATRIX, REQUIRED)
%               varargin    - Second argument can be the dimension (OPTIONAL)
%
% OUTPUT:       MatrixOut   - Output matrix (MATRIX)
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Backwards compatibility for flipping left-right in standard
%               space (NB: this can be different than in native space!).
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      [MatrixOut] = xASL_im_Flip(MatrixIn, [dim])
% __________________________________
% Copyright 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

if nargin < 2
	dim = 1;
end
if nargin == 2
	dim = varargin{1};
end
if nargin < 1
	error('xASL_im_Flip requires at least on parameter.');
end
if nargin > 2
	error('xASL_im_Flip requires at maximum two parameters.');
end
dim = round(dim);
if (dim < 1) || (dim > 3)
	error('Dim must be between 1 and 3');
end
switch (dim)
	case 1
		MatrixOut = MatrixIn(end:-1:1,:,:,:,:,:,:);
	case 2
		MatrixOut = MatrixIn(:,end:-1:1,:,:,:,:,:);
	case 3
		MatrixOut = MatrixIn(:,:,end:-1:1,:,:,:,:);
end
% MatrixOut = flip(MatrixIn,1);
end
