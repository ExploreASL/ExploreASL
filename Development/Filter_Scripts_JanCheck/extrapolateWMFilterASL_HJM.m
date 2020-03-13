function varargout=extrapolateWMFilterASL_HJM(varargin)
% function [corr_cbf,corr]=extrapolateWMFilterASL_HJM(niifile|matrix,[voxelsize],[label],[mask])
%
% Extrapolate WM-signal to correct GM-signal in ASL.
% NaNs are treated as zeros during filtering, and then set back.
%
% Data need to be motion-corrected, e.g. run first:
%     spm_realign(spm_vol('ASL4D.nii'))
%     changeResolutionNifti('ASL4D.nii',[3 3 3])
%     spm_reslice({'wASL4D.nii','ASL4D.nii'})
%
% DIPimage must be installed (http://www.diplib.org)
%
% Inputs:
%     niifile|matrix  : nifti filename or matrix
%     (voxelsize): voxel size
%     (label)  : binary vector, 0=unlabeled, 1=labeled
%                (optional, default=[0 1 0 1 ...])
%     (mask)   : 3D brain mask (optional, computed if not given)
% Outputs:
%     corr_cbf : corrected difference images (labeled-unlabeled)
%     corr     : correction
% No output arguments: saves a few files to disk
%
% Matthan Caan, AMC 2016

%% variables
spatSigma=[16 16 16]; %7; % smoothing in N (isotropic)
labelIDs=[0 1]; % unlabeled/labeled IDs in label input

%% CHANGED by HJM
% doPlot was 1, set to 0, since older Matlab versions do not support the function "arr"

doPlot=1; % do not plot images



%% process and sanity check
% try
%   dip_image;
% catch
%   error('DIPimage not installed.')
% end

if ischar(varargin{1})
  niifile=varargin{1};
  if ~xASL_exist(niifile)
    error('Does not exist.')
  end
  disp(niifile)
  n=xASL_io_ReadNifti(niifile);
  sz=n.dat.dim;
  vol=n.dat(:,:,:,:);
  vol(isnan(vol))=0;
  voxelsize=n.hdr.pixdim(2:4);
elseif isnumeric(varargin{1})
  vol=varargin{1};
  voxelsize=varargin{2};
  sz=size(vol);
else
  error('Unexpected input.')
end

if length(sz)~=4
  error('No 4D input')
end
% voxel size
if length(voxelsize)~=3
  error('Voxel size must contain 3 elements.')
end
if min(voxelsize)<1
  warning('Small voxel size for ASL-data (<1mm), check if correct.')
end
% spatSigma=spatSigma./double(voxelsize); % apply spacing
FwHm2SD     = (2*(2*reallog(2))^0.5);
% FWHM = (2*spatSigma+1)/5*voxelsize*FwHm2SD;
% disp(['Applying smoothing kernel ' num2str(spatSigma,3)])

% label vector
if nargin<3 || isempty(varargin{3})
  disp('Assuming alternating unlabeled/labeled dynamics')
  label=[zeros(1,sz(4)/2);ones(1,sz(4)/2)];
  label=label(:);
else
  disp('User specified labels.')
  label=varargin{3};
end
if length(label)~=sz(4)
  error('Label vector does not match nr of dynamics')
end
if ~all(sort(unique(label(:)))==[0;1])
  warning('Label vector does not contain both 0s and 1s.');
end
if sum(label==0)~=sum(label==1)
  warning('Number of unlabeled and labeled dynamics does not match.');
end
disp(['Nr of unlabeled pairs: ' num2str(sum(label==0))]);
disp(['Nr of labeled pairs: ' num2str(sum(label==1))]);

% create mask, if no mask was given
if nargin<4 || isempty(varargin{4});
  disp('Calculating mask:');
  mask=zeros(sz(1:3));
  % 2D-thresholding, erosion and closing
  for ii=1:sz(3)
    tmp=mean(vol(:,:,ii,:),4);
    if any(tmp(:))
      %mask(:,:,ii)=dip_array(closing(erosion(threshold(tmp),3),5));
      tmp = xASL_im_IsoDataThreshold(tmp);
      tmp = xASL_im_DilateErodeSeparable(double(tmp),'erode',[1 1 1],[1 1 1],1);
      tmp = xASL_im_DilateErodeFull(tmp,'dilate',[0 1 0;1 1 1;0 1 0]);
      tmp = xASL_im_DilateErodeSeparable(tmp,'dilate',[1 1 1],[1 1 1],1);
      tmp = xASL_im_DilateErodeFull(tmp,'erode',[0 1 0;1 1 1;0 1 0]);
      mask(:,:,ii) = xASL_im_DilateErodeSeparable(tmp,'erode',[1 1 1],[1 1 1],1);
    end
  end
  disp([num2str(sum(mask(:))) ' voxels.'])
else
  disp('Using user-specified mask.');
  mask=varargin{4};
end
if ~all(size(mask)==sz(1:3))
  error('Mask and volume size do not match.');
end
if  length(unique(mask(:)))>2
    error('Mask is not binary.');
end
if ~all(sort(unique(mask(:)))==[0;1])
   error('Mask is not binary.');
end
mask=mask>0;

%% prepare
% deal with NaNs
nanmask=any(isnan(vol),4);
vol=reshape(vol,[],sz(4));
vol(nanmask,:)=0;
vol=reshape(vol,sz);

idx1=label==labelIDs(1);
idx2=label==labelIDs(2);
d=vol(:,:,:,idx1)-vol(:,:,:,idx2); % difference
cbf=median(d,4); % cbf as median - robust to outliers

% % revert labels if needed -> PM: taken care of by the filter wrapper of ExploreASL already
% if mean(cbf(mask(:)))<0
%   idx1=label==labelIDs(2);
%   idx2=label==labelIDs(1);
%   d=vol(:,:,:,idx1)-vol(:,:,:,idx2); % difference
%   cbf=median(d,4); % cbf
% end

md=bsxfun(@minus,d,cbf); % demean

% gm=xASL_im_IsoDataThreshold(cbf); % GM-mask

sortInt     = sort(cbf(:));
GMthr       = sortInt(round(0.85.*length(sortInt)));
gm          = cbf>GMthr;

wm=~gm;           % non-GM-mask (incl background)

%warning('dd')
%wm=ishow('wm_rtemp_ASL4D.nii')>0;
%% do correction
%extrapolate dCBF from WM-mask
sz=size(md);

% make robust for
%[xx,yy,zz]=meshgrid(1:sz(1),1:sz(2),1:sz(3));
[xx,yy,zz]=ndgrid(1:sz(1),1:sz(2),1:sz(3));
% xx=permute(xx,[2 1 3]);
% yy=permute(yy,[2 1 3]);
% zz=permute(zz,[2 1 3]);

% loop over dynamics
fprintf('%s\n','Smoothing dynamics for filter, some patience please...    ');
for ii=1:sz(4)
    xASL_TrackProgress(ii,sz(4));
%   disp(num2str(ii))
  %tmp=dip_array(gaussf(md(:,:,:,ii),spatSigma)); % blur
  tmp=md(:,:,:,ii); % do NOT blur to avoid WM/GM signal contamination
  tmp2=tmp;
  tmp2(~wm)=0;
  %im=vdt(~wm); % vector-distance transform
  %x2=xx+dip_array(im{1});
  %y2=yy+dip_array(im{2});
  %z2=zz+dip_array(im{3});
  [~,vdtx,vdty,vdtz] = xASL_mex_chamfers3D(double(wm));
  x2 = xx+vdtx;
  y2 = yy+vdty;
  z2 = zz+vdtz;
  %o=interp3(xx,yy,zz,tmp,x2,y2,z2); % NN-interpolation
  o=interp3(yy,xx,zz,tmp,y2,x2,z2); % NN-interpolation
  tmp2(~wm)=o(~wm);
  % remove baseline WM signal, i.e. do not correct for 0-order WM signal
  % (more stringent filtering, leaves more noise when turned on)
  %tmp2=tmp2-mean(tmp2(wm>0));
  tmp2(~mask)=0;
  %md(:,:,:,ii)=dip_array(gaussf(tmp2,spatSigma)); % blur a little bit
  md(:,:,:,ii)=xASL_im_ndnanfilter(tmp2,'gauss',spatSigma*2.335,0); % blur a little bit

    tmpFilter   = tmp2;
%     tmpFilter   = xASL_im_ndnanfilter(tmpFilter,'rect',repmat(spatSigma,[1 3])); % blur a little bit
%     for ii=1:16
%         tmpFilter   = xASL_im_ndnanfilter(tmpFilter,'rect',[2 2 2]);
%     end

   md(:,:,:,ii)     = tmpFilter;

end

% calc corrected cbf
cbf_corr=d-md;

%% stats and plot images
tmp=reshape(d,[],sz(4));
std_before=std(tmp(mask,:),[],2);
tmp=reshape(cbf_corr,[],sz(4));
std_after=std(tmp(mask,:),[],2);
disp(['Std from ' num2str(median(std_before)) ' to ' num2str(median(std_after))])

% get mean GM/WM orig&corrected timeseries
tmp=reshape(d,[],sz(4));
gmorig=mean(tmp(gm(:)>0,:),1);
tmp=reshape(cbf_corr,[],sz(4));
gmcorr=mean(tmp(gm(:)>0,:),1);
tmp=reshape(d,[],sz(4));
wmorig=mean(tmp(wm(:)>0,:),1);
tmp=reshape(cbf_corr,[],sz(4));
wmcorr=mean(tmp(wm(:)>0,:),1);
% diff first/last quartile
disp('CoV first / last quartile:')
idx1=1:round(sz(4)/4);
idx2=sz(4)-round(sz(4)/4)+1:sz(4);
disp(['Before: ' num2str(mean(gmorig(idx1)-gmorig(idx2))/std(gmorig))])
disp(['After:  ' num2str(mean(gmcorr(idx1)-gmcorr(idx2))/std(gmcorr))]);

% if doPlot
%   dipshow(arr(mean(cbf_corr(:,:,:,1:30),4)),[0 75])
%   try,arr(wm),end
%
%   figure
%   subplot(2,1,1)
%   plot([gmorig;gmcorr]')
%   legend({'orig','corr'})
%   title('mean GM')
%   subplot(2,1,2)
%   plot([wmorig;wmcorr]')
%   legend({'orig','corr'})
%   title('mean non-GM')
% end

%% set back NaNs
nanmask=double(nanmask);
nanmask(nanmask>0)=NaN;
nanmask(~isnan(nanmask))=1;
cbf_corr=bsxfun(@times,cbf_corr,nanmask);
md=bsxfun(@times,md,nanmask);

%% write output (specific for large series pharmaASL)
if nargout
  varargout{1}=cbf_corr;
  varargout{2}=md;
else
  disp('Writing output...')
  [path,name,ext]=fileparts(niifile);
  ofile=fullfile(path,['wmfilt_diff_' name ext]);
  n.dat.fname=ofile;
  n.dat.dim=size(cbf_corr);
  n.dat.scl_slope=max(abs(cbf_corr(:)))/1e4;
  create(n);
  n.dat(:,:,:,:)=cbf_corr;

  [path,name,ext]=fileparts(niifile);
  ofile=fullfile(path,['orig_diff_' name ext]);
  n.dat.fname=ofile;
  n.dat.dim=size(d);
  n.dat.scl_slope=max(abs(d(:)))/1e4;
  create(n);
  n.dat(:,:,:,:)=d;

  ofile=fullfile(path,['cbf_first30_' name ext]);
  n.dat.fname=ofile;
  n.dat.dim=size(cbf);
  tmp=mean(cbf_corr(:,:,:,1:30),4);
  n.dat.scl_slope=max(abs(tmp(:)))/1e4;
  create(n);
  n.dat(:,:,:,:)=tmp;

  ofile=fullfile(path,['cbf_last30_' name ext]);
  n.dat.fname=ofile;
  n.dat.dim=size(cbf);
  tmp=mean(cbf_corr(:,:,:,end-29:end),4);
  n.dat.scl_slope=max(abs(tmp(:)))/1e4;
  create(n);
  n.dat(:,:,:,:)=tmp;

  ofile=fullfile(path,['cbf_first30_orig_' name ext]);
  n.dat.fname=ofile;
  n.dat.dim=size(cbf);
  tmp=mean(d(:,:,:,1:30),4);
  n.dat.scl_slope=max(abs(tmp(:)))/1e4;
  create(n);
  n.dat(:,:,:,:)=tmp;

  ofile=fullfile(path,['cbf_last30_orig_' name ext]);
  n.dat.fname=ofile;
  n.dat.dim=size(cbf);
  tmp=mean(d(:,:,:,end-29:end),4);
  n.dat.scl_slope=max(abs(tmp(:)))/1e4;
  create(n);
  n.dat(:,:,:,:)=tmp;

  ofile=fullfile(path,['wm_' name ext]);
  n.dat.fname=ofile;
  n.dat.dim=size(wm);
  n.dat.scl_slope=max(abs(wm(:)))/1e4;
  create(n);
  n.dat(:,:,:,:)=wm;

end
disp('Done.')
