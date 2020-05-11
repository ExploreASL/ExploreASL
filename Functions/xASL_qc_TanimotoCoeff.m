function TC = xASL_qc_TanimotoCoeff(Image1, Image2, imMask, type, bClip, bSmooth)
%xASL_qc_TanimotoCoeff Calculates the Tanimoto similarity coefficient for a binary/fuzzy set/continuous images.
% FORMAT: TC = xASL_qc_TanimotoCoeff(Image1, Image2[, imMask, type])
%
% INPUT:
%   Image1 - First input image, can also be a path
%   Image2 - Second input image, can also be a path
%   imMask - Mask over which the similarity is calculated (DEFAULT all)
%   type   - 1 binary input image
%            2 fuzzy set input (values in [0,1])
%            3 continuous TC with any real values
%            DEFAULT 1/2/3 - based on the data
%   bClip  - specifies whether we clip the images below 0 and
%            above a specified threshold (OPTIONAL, DEFAULT = no clipping)
%            Any non-zero bClip value will clip below 0. When bClip is
%            provided as percentile (e.g. bClip==0.975), so between 0 and
%            1, it will be used to clip peak values. And the images will be
%            scaled as well.
%            So for e.g. bClip==1, there will be clipped below 0 but not at
%            the peak values.
%   bSmooth - specifies how to smooth images with an isotropic kernel.
%             A single FWHM kernel value needs to be specified for each image in voxels, e.g. [4 4] would smooth both images with [4 4 4]
%             but [3 0] would only smooth the first image with a [3 3 3] kernel. (OPTIONAL, DEFAULT = [0 0],
%             which is no smoothing
%             Note that these are FWHM in voxels, not in mm!
%
% OUTPUT:
%   TC - Tanimoto similarity coefficient
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Compares images Image1 and Image2 within the mask imMask. TYPE specifies the input data type.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% RATIONALE:   Note that the Tanimoto Coefficient is a measure of image
%              overlap/intersection, similar to the Dice coefficient. With
%              the option type 3, this is a fuzzy coefficient, which
%              doesn't require to convert the two images to a binary mask.
%              The TC can be interpreted as a stringent Kappa, ranging from
%              0 (completely dissimilar) to 100% (identical images).
%              Assuming that perfect registration should not lead to
%              identical images but still retain physiological differences,
%              TC>70% can be regarded as excellent image agreement. The TC
%              will be overestimated when smoothing, but this may lead to more stable artifact detection.
%
% GENERAL EXAMPLE: TC = xASL_qc_TanimotoCoeff([0,1],[0,0],[1,1],1);
% ExploreQC EXAMPLE: TC = xASL_qc_TanimotoCoeff(CBFpath, TemplatePath, x.WBmask, 3, 0.975, [4 0]);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES:
%    Tanimoto, Taffee T. (17 Nov 1958). "An Elementary Mathematical theory of Classification and Prediction".
%        Internal IBM Technical Report. 1957
%
%    W. R. Crum, O. Camara and D. L. G. Hill, "Generalized Overlap Measures for Evaluation and Validation in
%        Medical Image Analysis," in IEEE Transactions on Medical Imaging, vol. 25, no. 11, pp. 1451-1461, Nov. 2006.
%    David C. Anastasiu, George Karypis, "Efficient identification of Tanimoto nearest neighbors",
%        International Journal of Data Science and Analytics, 2017.
% __________________________________
% Copyright 2015-2020 ExploreASL


    %% Administration
    if nargin<6 || isempty(bSmooth)
        bSmooth = [0 0];
    end
    if nargin<5 || isempty(bClip)
        bClip = 0;
    end
	if nargin < 2 || isempty(Image1) || isempty(Image2)
		error('Need to have at least two inputs');
    end

    Image1 = xASL_io_Nifti2Im(Image1); % this allows for image matrix or NIfTI path input
    Image2 = xASL_io_Nifti2Im(Image2); % this allows for image matrix or NIfTI path input

	if ~isequal(size(Image1),size(Image2))
		error('The two input images need to have the same size');
	end

	if nargin < 3 || isempty(imMask)
		imMask = ones(size(Image1));
	end

	if nargin < 4 || isempty(type)
        fprintf('Warning, didnt know image type, automatically selecting one\n');
		if islogical(Image1) && islogical(Image2)
			% Binary input
			type = 1;
		else
			if max(Image1(:))<=1 && max(Image2(:))<=1 && min(Image1(:))>=0 && min(Image2(:))>=0
				% Fuzzy set input
				type = 2;
			else
				% Any value input
				type = 3;
			end
		end
    end

    %% Deal with clipping
    if bClip~=0
        % Clip below zero
        Image1(Image1<0) = 0;
        Image2(Image2<0) = 0;
    end
    if bClip>0 && bClip<1
        SortedIntensity1 = sort(Image1(isfinite(Image1)));
        SortedIntensity2 = sort(Image2(isfinite(Image2)));
        Threshold1 = SortedIntensity1(round(bClip*length(SortedIntensity1)));
        Threshold2 = SortedIntensity2(round(bClip*length(SortedIntensity2)));
        Image1(Image1>Threshold1) = Threshold1;
        Image2(Image2>Threshold2) = Threshold2;
        Image1 = Image1./Threshold1;
        Image2 = Image2./Threshold2;
    end

    %% Apply smoothing
    if bSmooth(1)~=0
        Image1 = xASL_im_ndnanfilter(Image1,'gauss',[bSmooth(1) bSmooth(1) bSmooth(1)]);
    end
    if bSmooth(2)~=0
        Image2 = xASL_im_ndnanfilter(Image2,'gauss',[bSmooth(2) bSmooth(2) bSmooth(2)]);
    end

    %% Deal with mask
	imMask = imMask.*(1-isnan(Image1)).*(1-isnan(Image2));
	imMask = imMask > 0;
	Image1 = Image1(imMask);
	Image2 = Image2(imMask);

    %% Compute the TC
	switch (type)
		case 1
			% Needs to have a logical input
			if ~islogical(Image1)
				Image1 = Image1>0;
			end

			if ~islogical(Image2)
				Image2 = Image2>0;
			end
			TC = sum(Image1&Image2)/sum(Image1|Image2);

		case 2
			TC = sum(min([Image1,Image2],[],2))/sum(max([Image1,Image2],[],2));
		case 3
			TC = sum(Image1.*Image2) / (sum(Image1.^2) + sum(Image2.^2) - sum(Image1.*Image2));
		otherwise
			error('Unknown Tanimoto Coefficient type to compute');
    end


end
