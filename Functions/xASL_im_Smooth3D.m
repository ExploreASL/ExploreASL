function [imSmo,imGaussX,imGaussY,imGaussZ] = xASL_im_Smooth3D(sigma,imIn,PSFtype)
%xASL_im_Smooth3D ...
%
% FORMAT:       [imSmo,imGaussX,imGaussY,imGaussZ] = xASL_im_Smooth3D(sigma,imIn,PSFtype)
% 
% INPUT:        sigma   - ...
%               imIn    - ...
%               PSFtype - ...
%
% OUTPUT:       imSmo    - ...
%               imGaussX - ...
%               imGaussY - ...
%               imGaussZ - ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2021 ExploreASL

for tp=1:length(PSFtype)
    switch(PSFtype{tp})
        case 'gaussian'
            windowSize = round(2.5*sigma);
            imGauss{tp} = exp( -((-windowSize(tp):windowSize(tp)).^2)./2./sigma(tp)./sigma(tp));
        case 'lorentzian'
            windowSize = [10 10 10];
            imGauss{tp} = 1./(1+ ((-windowSize(tp):windowSize(tp)).^2)/sigma(tp)/sigma(tp));
        case 'sinc'
            windowSize = [12 12 12];
            imGauss{tp} = sin(sigma(tp)*pi*(-windowSize(tp):windowSize(tp)))./(sigma(tp)*pi*(-windowSize(tp):windowSize(tp)));imGauss{tp}(windowSize(tp)+1) = 1;
            %    case 'gaussexp'
            %        % a * exp(-r*r/(2*b*b)) + (1-a)*exp(-r/c)
            %        imGaussX = sigma(1) * exp( -((-windowSize(1):windowSize(1)).^2)./2./sigma(2)./sigma(2)) + (1-sigma(1)) *  exp( -(abs(-windowSize(1):windowSize(1)))./sigma(3));
            %        imGaussY = sigma(4) * exp( -((-windowSize(1):windowSize(1)).^2)./2./sigma(5)./sigma(5)) + (1-sigma(4)) *  exp( -(abs(-windowSize(1):windowSize(1)))./sigma(6));
            %        imGaussZ = sigma(7) * exp( -((-windowSize(1):windowSize(1)).^2)./2./sigma(8)./sigma(8)) + (1-sigma(7)) *  exp( -(abs(-windowSize(1):windowSize(1)))./sigma(9));
        case 'flat'
            %imGauss{tp} = ones(1,windowSize(tp)*2+1);
            windowSize = [10 10 10];
            fwhm = 2.355*sigma(tp);
            imGauss{tp} = zeros(1,windowSize(tp)*2+1);
            if fwhm>=1 
                imGauss{tp}(windowSize(tp)+1) = 1;
			end
            marg = floor((fwhm-1)/2);
            imGauss{tp}(windowSize(tp)+1+(-marg:marg)) = 1;
            
            pos = marg + 1;
            fwhm = fwhm-2*marg-1;
            if (windowSize(tp)+1-pos)>0
                imGauss{tp}(windowSize(tp)+1-pos) = fwhm/2;
			end
            if (windowSize(tp)+1+pos)<=length(imGauss{tp})
                imGauss{tp}(windowSize(tp)+1+pos) = fwhm/2;
			end
	end
end
imGaussX = imGauss{1};
imGaussY = imGauss{2};
imGaussZ = imGauss{3};

imGaussX = imGaussX(:)/sum(imGaussX(:));
%locGMX = convn(imIn,imGaussX,'same');

imGaussY = imGaussY/sum(imGaussY(:));
%locGMY = convn(locGMX,imGaussY,'same');

imGaussZ = repmat(imGaussZ/sum(imGaussZ(:)),[1 1 1]);
imGaussZ = shiftdim(imGaussZ,-1);
%imSmo = convn(locGMY,imGaussZ,'same');

imSmo = xASL_im_conv3Dsep(imIn,imGaussX,imGaussY,imGaussZ);

return
