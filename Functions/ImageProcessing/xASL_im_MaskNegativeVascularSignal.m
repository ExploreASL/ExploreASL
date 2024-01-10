function [NegativeMask, TreatedPWI] = xASL_im_MaskNegativeVascularSignal(x, IsSpace)
%xASL_im_MaskNegativeVascularSignal Segment ASL image clusters with
%significantly negative signal
%
% FORMAT: [NegativeMask, TreatedPWI] = xASL_quant_DetectNegativeVascularSignal(x)
%
% INPUT:
%   x                   - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging (REQUIRED)
%   x.P.Path_PWI        - path to PWI image
% & x.P.Pop_Path_PWI
%   x.P.Path_PVgm       - path to GM segmentation map
% & x.P.Pop_Path_PV_pGM
%   IsSpace             - 1 for native ASL space, 2 for standard space (REQUIRED)
%
% OUTPUT:
%   NegativeMask    - Binary mask with 1 for negative signal
%   TreatedPWI      - Original PWI image, with negative vascular clusters interpolated
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function segments clusters with significant negative
%              ASL signal. This can be tricky as there is also the negative tail of Gaussian noise
%              from the ASL subtraction. The image feature we use here, is that negative
%              vascular signal will be a relatively large region with
%              significant median negative value, whereas noise will be
%              regions with relatively small negative signal.
%              Negative signal from wrong background suppression timing
%              (e.g. in the first slice with **2D EPI**) can be masked out with
%              this as well.
%              The procedure works as follows:
%
%              1. Obtain mask of negative voxels within `pGM>0.5` mask
%              2. Obtain distribution of subzero clusters
%              3. Define the negative threshold
%              4. Create mask by thresholding whole image
%
%              Note that the definition of the threshold is obtained within
%              the GM only, but that this threshold is applied to the full image.
%
%              Note that instead of the PWI path input, a CBF image should
%              work equally well, as we don't expect a smooth M0 biasfield
%              to change the distribution of negative clusters
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: NegativeMask = xASL_im_MaskNegativeVascularSignal(x);
% __________________________________
% Copyright (c) 2015-2024 ExploreASL
 
    if nargin<2 || isempty(IsSpace)
        warning('Didnt know space to calculate negative mask in, skipping');
        return;
    end
 
    fprintf('%s\n','Detecting vascular clusters with negative ASL signal:...');
 
    if IsSpace==1 % native space
        PWIpath = x.P.Path_PWI;
        GMpath = x.P.Path_PVgm;
    elseif IsSpace==2 % standard space
        PWIpath = x.P.Pop_Path_PWI;
        GMpath = x.P.Pop_Path_PV_pGM;
        % no need to reslice
    end
    
    if ~xASL_exist(PWIpath, 'file')
        warning([PWIpath ' missing, skipping']);
        return;
    elseif ~xASL_exist(GMpath, 'file')
        warning([GMpath ' missing, skipping']);
        return;
    end
        
 
    %% 1. Obtain mask of negative voxels within pGM>0.5 mask
    % Create GMmask in ASL space from SPM
    PWIim = xASL_io_Nifti2Im(PWIpath);
    if size(PWIim,4)>1
        PWIim = xASL_stat_MeanNan(PWIim,4);
    end
    GMmask = xASL_io_Nifti2Im(GMpath)>0.5;
    
    if ~isequal(size(PWIim),size(GMmask))
        warning('Sizes of PWI and GM images differed');
    end
    
    NegativeMask = PWIim<0 & GMmask;
 
    if IsSpace==2
        % dilate to connect neighboring labels
        NegativeLabels = xASL_im_DilateErodeFull(NegativeMask,'dilate',xASL_im_DilateErodeSphere(2));
    else
        NegativeLabels = NegativeMask;
    end
 
    if sum(NegativeMask(:))>0 % if there are subzero voxels
 
        %% 2. Obtain distribution of subzero clusters
        % Label (i.e. assign a number) to each subzero clusters
        labeled = single(spm_bwlabel(double(NegativeLabels),26));
%         labeled = xASL_im_DilateErodeFull(labeled,'dilate',xASL_im_DilateErodeSphere(1));
        % 26 = corner criterion for connectivity (being strict in the definition of "clusters")
 
        % mask labels with NegativeMask again
        labeled(~NegativeMask) = 0;
 
        % Get the mean value of each subzero clusters
        for iP=1:max(labeled(:))
            xASL_TrackProgress(iP, max(labeled(:)));
            MeanValue(iP) = mean(PWIim(labeled==iP));
        end
 
        %% 3. Define the negative threshold
        medianValue = xASL_stat_MedianNan(MeanValue);
        madValue = xASL_stat_MadNan(MeanValue,1);
        ClipThr = medianValue - (4*madValue);
 
%         [N, X] = hist(MeanValue);
%         N = N./sum(N(:));
%         figure(1);plot(X, N)
 
    else
        ClipThr = -max(PWIim(:)); % don't clip
    end
 
    fprintf('\b');
    xASL_TrackProgress(0.5,1);
 
    %% 4. Create mask by thresholding whole image
    TreatedPWI = PWIim;
    NegativeMask = TreatedPWI<ClipThr;
%     NegativeMask = xASL_im_DilateErodeFull(NegativeMask,'dilate',xASL_im_DilateErodeSphere(1));
    TreatedPWI(NegativeMask) = NaN;
 
    if nargout>1
        xASL_TrackProgress(0.75,1);
        TreatedPWI = xASL_im_ndnanfilter(TreatedPWI,'gauss', [8 8 8], 2); % extrapolation only
    end
 
    xASL_TrackProgress(1,1);
    fprintf('\n');
    
end