function TC = xASL_stat_TanimotoCoeff(im1, im2, imMask, type)
% Calculates the Tanimoto similarity coefficient for a binary/fuzzy set/continuous images.

% FORMAT: TC = xASL_stat_TanimotoCoeff(im1, im2[, imMask, type])
%
% INPUT:
%   im1    - First input image
%   im2    - Second input image
%   imMask - Mask over which the similarity is calculated (DEFAULT all)
%   type   - 1 binary input image
%            2 fuzzy set input (values in [0,1])
%            3 continuous TC with any real values
%            DEFAULT 1/2/3 - based on the data
% OUTPUT:
%   TC - Tanimoto similarity coefficient
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Compares images IM1 and IM2 over the mask IMMASK. TYPE specifies the input data type.
% EXAMPLE: TC = xASL_stat_TanimotoCoeff([0,1],[0,0],[1,1],1)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:
%    Tanimoto, Taffee T. (17 Nov 1958). "An Elementary Mathematical theory of Classification and Prediction". 
%        Internal IBM Technical Report. 1957
%
%    W. R. Crum, O. Camara and D. L. G. Hill, "Generalized Overlap Measures for Evaluation and Validation in 
%        Medical Image Analysis," in IEEE Transactions on Medical Imaging, vol. 25, no. 11, pp. 1451-1461, Nov. 2006.
% __________________________________
% Copyright © 2015-2019 ExploreASL
%
% 2019-04-17 JP

    % Parameter admin
	if nargin < 2
		error('xASL_stat_TanimotoCoeff: Need to have at least two inputs.');
	end
	
	if ~isequal(size(im1),size(im2))
		error('xASL_stat_TanimotoCoeff: The two input images need to have the same size.');
	end
	
	if nargin < 3 || isempty(imMask)
		imMask = ones(size(im1));
	end
	
	if nargin < 4 || isempty(type)
		if islogical(im1) && islogical(im2)
			% Binary input
			type = 1;
		else
			if max(im1(:))<=1 && max(im2(:))<=1 && min(im1(:))>=0 && min(im2(:))>=0
				% Fuzzy set input
				type = 2;
			else
				% Any value input
				type = 3;
			end
		end
	end
	
	imMask = imMask.*(1-isnan(im1)).*(1-isnan(im2));
	imMask = imMask > 0;
	im1 = im1(imMask);
	im2 = im2(imMask);
	
	switch (type)
		case 1
			% Needs to have a logical input
			if ~islogical(im1)
				im1 = im1>0;
			end
			
			if ~islogical(im2)
				im2 = im2>0;
			end
			TC = sum(im1&im2)/sum(im1|im2);
			
		case 2
			TC = sum(min([im1,im2],[],2))/sum(max([im1,im2],[],2));
		case 3
			TC = sum(im1.*im2) / (sum(im1.^2) + sum(im2.^2) - sum(im1.*im2));
		otherwise
			error('xASL_stat_TanimotoCoeff: Unknown type.');
	end
end
