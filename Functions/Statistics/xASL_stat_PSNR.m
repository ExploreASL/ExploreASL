function PSNR=xASL_stat_PSNR(imRef,imSrc)
%xASL_stat_PSNR Calculates peak signal-to-noise ratio.
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
% EXAMPLE:      ...
% __________________________________
% Copyright � 2015-2020 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

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