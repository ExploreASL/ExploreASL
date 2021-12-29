function [ovol,varargout]=bilateralFilterASL_HJM_LV(volIM,VoxelSize,varargin)
% function [ovol,[filter]]=bilateralFilterASL_HJM_LV(vol,voxelsize,[label],[mask])
%
% Previous version from 31.05.2018 - kept for testing
%
% Apply bilateral filter to unlabeled/labeled ASL data separately.
% This filters out temporal noise induced by too strong STIR fat suppression.
% Apply this filter after motion correction.
%
% DIPimage must be installed (http://www.diplib.org)
%
% Inputs:
%     vol:       4D-matrix with ASL-dynamics
%     voxelsize: voxel size in mm (3x1)
%     label    : binary vector, 0=unlabeled, 1=labeled 
%                (optional, default=[0 1 0 1 ...])
%     mask     : 3D brain mask (optional, computed if not given)
% Outputs:
%     ovol     : filtered data
%     filter   : applied correction (for debugging purposes)
%
% Matthan Caan, AMC 2016

% HJM:
% Masking disabled in script
% xASL_stat_MeanNan & xASL_stat_MedianNan instead of mean & median
% avoiding NaN-smoothing
% clearing to free memory

%% variables
tonSigma=0.1; % tonal sigma, relative to median signal intensity
spatSigma=4; % smoothing in mm (isotropic)
labelIDs=[0 1]; % unlabeled/labeled IDs in label input

%% process and sanity check
try
  dip_image;
catch
  error('DIPimage not installed.')
end

sz=size(volIM);
if length(sz)~=4
  error('No 4D input')
end

% voxel size
if length(VoxelSize)~=3
  error('Voxel size must contain 3 elements.')  
end
if min(VoxelSize)<1
  warning('Small voxel size for ASL-data (<1mm), check if correct.')  
end
spatSigma=spatSigma./double(VoxelSize); % apply spacing 
disp(['Spatial sigma in voxels: ' num2str(spatSigma)])

% label vector
if nargin<3 || isempty(varargin{1})
  disp('Assuming alternating unlabeled/labeled dynamics')
  label=[zeros(1,sz(4)/2);ones(1,sz(4)/2)];
  label=label(:);
else
  disp('User specified labels.')
  label=varargin{1};
end
if length(label)~=sz(4)
  error('Label vector does not match nr of dynamics')
end
if ~all(sort(unique(label(:)))==[0;1])
  error('Label vector should contain only 0s and 1s.')
end
if sum(label==0)~=sum(label==1)
  error('Number of unlabeled and labeled dynamics does not match.')
end
disp(['Nr of unlabeled/labeled pairs: ' num2str(sum(label==0))])

% mask
if nargin<4 || isempty(varargin{2});
  disp('Calculating mask:')
  mask=zeros(sz(1:3));
  % 2D-thresholding, erosion and closing
  for ii=1:sz(3)
    tmp = xASL_stat_MeanNan(volIM(:,:,ii,:),4);
	if any(tmp(:))
		tmp  		 = xASL_im_IsoDataThreshold(tmp);
		tmp  		 = xASL_im_DilateErodeSeparable(double(tmp),'erode',[1 1 1],[1 1 1],1);
		tmp  		 = xASL_im_DilateErodeFull(tmp,'dilate',[0 1 0;1 1 1;0 1 0]);
		tmp  		 = xASL_im_DilateErodeSeparable(tmp,'dilate',[1 1 1],[1 1 1],1);
		tmp  		 = xASL_im_DilateErodeFull(tmp,'erode',[0 1 0;1 1 1;0 1 0]);
		mask(:,:,ii) = xASL_im_DilateErodeSeparable(tmp,'erode',[1 1 1],[1 1 1],1);
	end
  end
  disp([num2str(sum(mask(:))) ' voxels.'])
else
  disp('Using user-specified mask.')
  mask=varargin{2};
end
if ~all(size(mask)==sz(1:3))
  error('Mask and volume size do not match.')
end
if ~all(sort(unique(mask(:)))==[0;1])
  error('Mask is not binary.')
end
mask=mask>0;

%% adjust tonal sigma to median intensity
tmp=reshape(volIM,[],sz(4));
tmpMask=tmp(mask>0,:);
scaleFactor=xASL_stat_MedianNan(tmpMask(:));
tonSigma=tonSigma*scaleFactor;
disp(['Tonal sigma bilateral filter: ' num2str(tonSigma)]);

%% filter data, unlabeled/labeled separately
ovol=zeros(size(volIM)); ofilter=ovol;
for iLabel=1:2
      clear idx tmp sz mtime mspace std_before NaNmask gtmp otmp MeanV

      disp(['Label ' num2str(labelIDs(iLabel))]);
      idx=label==labelIDs(iLabel);
      disp(['Selected ' num2str(sum(idx)) ' volumes.'])
      tmp=squeeze(volIM(:,:,:,idx)); %odd/even only
      sz=size(tmp);
      
      % First remove temporal low-frequency changes
      for iV=1:sz(4)
          TempTmp           = tmp(:,:,:,iV);
          TempTmp           = TempTmp(mask);
          MeanV(iV,1)       = xASL_stat_MeanNan(TempTmp(:));
      end
      MeanV                 = xASL_im_ndnanfilter(MeanV,'gauss',15.07);
      MeanV                 = MeanV./mean(MeanV);
      for iV=1:sz(4)
        tmp(:,:,:,iV)       = tmp(:,:,:,iV)./MeanV(iV,1);
      end

      % Now compute & remove the smooth residuals (i.e. the artifacts)
      tmp               = reshape(tmp,[],sz(end));
      mtime             = xASL_stat_MeanNan(tmp,2)*ones(1,sz(end));
      tmp               = tmp-mtime; % remove mean over time
      mspace            = xASL_stat_MeanNan(tmp(mask>0,:),1); % calc mean
      tmp               = bsxfun(@minus,tmp,mspace); % remove mean over space
    %   tmp(~mask,:)=0;
      std_before        = std(tmp(mask,:),[],2);
      tmp               = reshape(tmp,sz);

      %% Avoid NaN-smoothing
      NaNmask           = isnan(tmp);
      tmp(NaNmask)      = 0;

      gtmp=zeros(size(tmp));otmp=gtmp;
      fprintf('%s\n','Smoothing dynamics for filter, some patience please...  ');
      for ii=1:size(tmp,4)
        xASL_TrackProgress(ii,size(tmp,4));
        gtmp(:,:,:,ii)=dip_array(medif(bilateralf(tmp(:,:,:,ii),spatSigma,tonSigma),2)); %bilateral filter 
        otmp(:,:,:,ii)=tmp(:,:,:,ii)-gtmp(:,:,:,ii);       % subtract
      end

      otmp(NaNmask)     = NaN;

      otmp=reshape(otmp,[],sz(4));
      std_after=std(otmp(mask,:),[],2);
      otmp=bsxfun(@plus,otmp,mspace); % add mean back
      otmp=otmp+mtime; % add mtime back
      otmp=reshape(otmp,sz);
    %   otmp=otmp.*repmat(mask,[1 1 1 sz(4)]); % re-apply mask

     % Re-apply the temporal low-frequency changes
      for iV=1:sz(4)
        otmp(:,:,:,iV)       = otmp(:,:,:,iV).*MeanV(iV,1);
      end
      
      ovol(:,:,:,idx)=otmp;
      ofilter(:,:,:,idx)=gtmp;

      disp(['Std from ' num2str(median(std_before)) ' to ' num2str(median(std_after))]);
end

if nargout==2
  varargout{1}=ofilter;
end
