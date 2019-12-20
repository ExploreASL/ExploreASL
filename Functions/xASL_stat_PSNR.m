function PSNR=xASL_stat_PSNR(imRef,imSrc)
% Calculates peak signal-to-noise ratio.
%
% FORMAT: PSNR=xASL_stat_PSNR(imRef,imSrc)
%
% INPUT:
%   imRef   - the noiseless reference image, can be 3D
%   imSrc   - noisy source image of the same size
%
% OUTPUT:
%   PSNR    - calculated PSNR
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Calculates the PSNR, needs two input arguments - 3D images of the same size. 
%              Uses 95% percentile instead of MAX.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2017-00-00 HJ

    % Admin
	if nargin ~= 2
		error('xASL_stat_PSNR: Needs two input arguments.');
	end
	
	if ~isequal(size(imRef),size(imSrc))
		error('xASL_stat_PSNR: Both inputs need to have the same size.');
	end
	
    % Create finite mask
    FiniteMask  = isfinite(imRef) & isfinite(imSrc);

	% Calculate 5 and 95 percentile of the reference image to be robust to outliers
    SortRef     = sort(imRef(FiniteMask),'ascend');
    RobustMin   = SortRef(max(0,ceil(0.05*length(SortRef))));
    RobustMax   = SortRef(ceil(0.95*length(SortRef)));

    PeakValue   = RobustMax - RobustMin; % this is 2^16 for a 16 bit color image, 
    % provides proxy for dynamic range. with ASL there can also be negative values
    

    SqDiff      = double(imRef-imSrc).^2; % squared difference image
    MeanSqDiff  = xASL_stat_MeanNan(SqDiff(:));
    PSNR        = 10.0 * (log10((PeakValue^2) / MeanSqDiff));

end