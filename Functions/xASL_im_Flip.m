function [ MatrixOut ] = xASL_im_Flip( MatrixIn, varargin )
%xASL_im_Flip Backwards compatibility for flipping left-right in standard
%space (NB: this can be different than in native space!)

% [ MatrixOut ] = xASL_im_Flip( MatrixIn, [dim] )


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

