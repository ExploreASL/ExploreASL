function [Y,fil] = xASL_im_ndnanfilter(X,filterType,F,WNAN)
% xASL_im_ndnanfilter   3-dimensional convolution, ignoring NaNs.

% By Jan Petr/ HJ Mutsaerts:
%   Syntax:
%         Y = xASL_im_ndnanfilter(X,filterType,F);
%         Y = xASL_im_ndnanfilter(X,filterType,F,WNAN);
%     [Y,W] = xASL_im_ndnanfilter(...);
%
%   Input:
%     X          - Data to be filtered with/without NaNs.
%     filterType - Smoothing kernel type. 
%                  Can be either:
%                    'rect' for flat kernel 
%                    'gauss' for Gaussian kernel
%     F      - A vector specifying the semi-width of the window or FWHM of Gaussian
%              Length 1-3, depending on which dimension are smoothed
%              Give zeros for dimension that should not be smoothed at all
%                'rect' - a flat kernel of width 2*F+1. 
%                'gauss' - Gaussian kernel of FWHM F in voxels
%     WNAN   - Integer indicating NaNs treatment and program behaviour!:
%              0: Filters data and interpolates NaNs with values from non-NaNs (default). 
%              1: Filters data but leaves NaNs as NaNs
%              2: Leaves non-NaNs how they were and interpolate NaNs with values form non-NaNs
%              See the NOTEs below
%
%   Output:
%     Y      - Filtered X data (same size as X!).
%     W      - Cell 1x3 with the separable filters in each dimension.
%
%   Description:
%     This function applies a 3-dimensional convolution of X with given kernel.
%     NaNs elements are taken into account (ignored).
%
%     By default, edges are not padded and one-sided filter is used at the image edges.
%
%   Notes: 
%     * Accepts empty value for any input. When X is empty, the program can
%       be used as a N-dimensional window generator.
%     * NaNs elements surrounded by no-NaNs elements (which will depend on
%       window width) are the ones that will be interpolated. The others
%       are leaved untouched.
%     * When WNAN=2, the programs acts like an NAN-interpolat/GAP-filling,
%       leaving untouched the no-NaNs elements but the filtering is
%       perfomed anyway. I recommend the default behaviour (WNAN=0) in order
%       to keep the filtered data in the workspace, and then use the code
%       at the end of this function to get/remove the interpolated NaNs
%     * To achieve similar results as ndnanfilter previously, use same F
%       as with the 'rect' filter.
%     * Note that the FWHM of Gaussian is given in VOXELS, not in mm
%     * For the Gaussian filter, use (previous N, new FWHM)
					 % N= 1 ~ FWHM 0.94
					 % N= 2 ~ FWHM 1.885 
					 % N= 4 ~ FWHM 3.76  
					 % N= 6 ~ FWHM 5.652
					 % N= 8 ~ FWHM 7.536 
					 % N=10 ~ FWHM 9.42  
					 % N=12 ~ FWHM 11.3
					 % N=16 ~ FWHM 15.07
					 % N=20 ~ FWHM 18.84
					 % N=10/2.355 ~ FWHM 4
					 % Basically divide by 1.06


% Check inputs and sets defaults of main input arguments:
if nargin<3 || nargin>4
 error('xASL_im_ndnanfilter:IncorrectNumberOfInputs',...
  'At least 3 inputs are needed and less than 5.')
end

if ~exist('filterType','var')
    filterType = 'rect';
elseif isempty(filterType)
    filterType = 'rect';
end

if ~exist('F','var')
    F = 3;
elseif isempty(F)
    F = 3;
end

if ~exist('WNAN','var')
    WNAN = 0;
elseif isempty(WNAN)
    WNAN = 0;
end

% Filter has maximum three dimensions, all initialized to no filtering
fil{1} = 0;
fil{2} = 0;
fil{3} = 0;

% For all three dimensions
for k=1:3
	% If the image has this dimension and the FWHM of the filter is defined
	if (size(X,k) > 1) && (length(F)>=k)
		% And as non-zeros
		if F(k)>0
			switch(filterType)
				% Then either do a 2*F+1 rectangular filter
				case {'rect'}
					F = round(F);
					fil{k} = ones(2*F(k) + 1,1)/(2*F(k)+1);
				% or a Gaussian filter
				case {'gauss'}
					 % Calculating sigma in voxels
					 sigma = F(k)/2.355;
					 
					 % Set the semi-width of the window to 2.5 sigma
					 window = round(2.5*sigma);
					 
					 % Calculate the Gaussian
					 fil{k} = exp((-((-window:window)').^2)/(2.0*sigma*sigma));
					 
					 % Normalize sum to 1
					 fil{k} = fil{k}/sum(fil{k});
				otherwise
					error(['xASL_im_ndnanfilter: Unknown filter type : ' filterType]);
			end
		end
	end
end

% If no input data but two outputs then generates the window only:
if isempty(X)
	Y = [];
	return
end

% Check for NaN's:
inan = isnan(X);
ynan = any(inan(:));                       % Bug fixed 30/jun/2008
if ynan
	X(inan) = 0;
end

%%% No padding at all
Y = X;
 
% Convolves both arrays:
Y = xASL_mex_conv3Dsep(double(Y),fil{1},fil{2},fil{3});

% Estimates the averages when NaNs are present:
if ynan
	factor = ~inan;
	factor      = xASL_mex_conv3Dsep(double(factor),fil{1},fil{2},fil{3});
	%factor (factor < 0.0001) = 1;
	Y = Y./factor;
	Y(factor < 0.0001) = NaN;
end

% What about NaNs?:
if     WNAN == 1       % Leave NaNs elements untouched!
 Y(inan) = NaN;
elseif WNAN == 2       % Leave no-NaNs elements untouched!!!
 X(inan) = Y(inan);
 Y = X;
end  

