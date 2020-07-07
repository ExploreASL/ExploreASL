function xASL_im_CreateSliceGradient(x)
% xASL_im_CreateSliceGradient
%
% FORMAT:       xASL_im_CreateSliceGradient(x)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  1. Create slice gradient in same space as input file
%               2. Reslice slice gradient to MNI (using existing ASL matrix changes from e.g. registration to MNI, motion correction, registration to GM)
%               3. Creating average slice gradient
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
%
% BACKGROUND INFORMATION:
%
% When a 2D readout is used with ASL, post-label delay and hence T1 decay
% will be dependent on slice timing.
% Therefore, quantification part needs slice reference to quantify per
% slice and correct for effective post-label delay differences.
%
% This function uses exact same ASL matrix changes that occurred due to
% registration to MNI, motion correction and registration to T1.
%
% Script dependencies: SPM12
%
% HJ Mutsaerts, ExploreASL 2016
%
% __________________________________
% Copyright 2015-2020 ExploreASL


%% Administration

xASL_delete(x.P.Pop_Path_SliceGradient);
xASL_delete(x.P.Path_SliceGradient_extrapolated);


%% ------------------------------------------------------------------------------------
%% 1    Create slice gradient in same space as input file
tNII = xASL_io_ReadNifti(x.P.Path_ASL4D);
dim = size(tNII.dat(:,:,:,:,:,:,:,:,:));
dim = dim(1:3);
SGim = zeros(dim);

for iSlice=1:dim(3)
    SGim(:,:,iSlice) = iSlice;
end

% We skip motion correction here, to avoid edge artifacts
xASL_io_SaveNifti(x.P.Path_ASL4D,x.P.Path_SliceGradient,SGim,8,0);
xASL_delete(x.P.Path_SliceGradient_mat);
xASL_Copy(x.P.Path_SliceGradient, x.P.Path_SliceGradient_extrapolated,1);

%ResampleRes     = repmat(min(tNII.hdr.pixdim(2:4))-0.3,1,3);
PaddingN = 10;
PaddingN = round(PaddingN/2)*2; % Make PaddingN even

%% ------------------------------------------------------------------------------------
%% Here we add padding near the edges
xASL_im_Upsample( x.P.Path_SliceGradient_extrapolated, x.P.Path_SliceGradient_extrapolated, tNII.hdr.pixdim(2:4),0,[PaddingN PaddingN PaddingN]);
SGim = xASL_io_Nifti2Im(x.P.Path_SliceGradient_extrapolated);

PN = PaddingN/2;
SGim(:,:,1:PN) = 1; % Fill the bottom padded slices
SGim(:,:,dim(3)+PN+1:end) = dim(3); % Fill the top padded slices (top slice is also corrected here, was set to 0)
for iSlice=1:dim(3)
    SGim(:,:,iSlice+PN) = iSlice;
end

xASL_io_SaveNifti(x.P.Path_SliceGradient_extrapolated,x.P.Path_SliceGradient_extrapolated,SGim,8,0);

%% ------------------------------------------------------------------------------------
%%     Reslice slice gradient to MNI (using existing ASL matrix changes from e.g. registration to MNI, motion correction, registration to GM)

if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') % BACKWARDS COMPATIBILITY, CAN BE REMOVED
	AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
else
	AffineTransfPath = [];
end

xASL_spm_deformations(x, x.P.Path_SliceGradient, x.P.Pop_Path_SliceGradient, 0, [], AffineTransfPath, x.P.Path_y_ASL); % nearest neighbor

if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') % BACKWARDS COMPATIBILITY, CAN BE REMOVED
	AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
else
	AffineTransfPath = [];
end

xASL_spm_deformations(x, x.P.Path_SliceGradient_extrapolated, x.P.Pop_Path_SliceGradient_extrapolated,1,[], AffineTransfPath, x.P.Path_y_ASL);
% Since smoothing is desirable for this slice gradient, and spline edge artifacts are not, we use trilinear here as interpolation method

%% Add some more extrapolation to be sure
SGim = xASL_io_Nifti2Im(x.P.Pop_Path_SliceGradient_extrapolated);
SGim(SGim==0) = NaN;
countNanLast = numel(SGim)+1;
while sum(isnan(SGim(:)))<countNanLast
	countNanLast = sum(isnan(SGim(:)));
	SGim = xASL_im_ndnanfilter(SGim,'rect',[8 8 8],2);
end

xASL_io_SaveNifti(x.P.Pop_Path_SliceGradient_extrapolated,x.P.Pop_Path_SliceGradient_extrapolated,SGim,[],0);

end
