function PathFig = xASL_im_CreateSlabImage(x)
% xASL_im_CreateSlabImage
%
% FORMAT:       xASL_im_CreateSlabImage(x)
% 
% INPUT:        x - ExploreASL x structure (STRUCT, REQUIRED)
%
% OUTPUT:       Path to the image that is generated.
% OUTPUT:       /MyStudy/Population/T1Check/Slab.jpg 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  1. Get slab location from ASL sidecar
%               2. Check if anat image exists, if so plot slab and readout location to t1w
%               3. Save Slab image to disk
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
%
% BACKGROUND INFORMATION:
%
% __________________________________
% Copyright 2015-2023 ExploreASL


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
%================================================================================================
if ~xASL_exist(x.D.ROOT); return; end
Dlist = xASL_adm_GetFileList(x.D.ROOT,'^.*$','List',[0 Inf], true);
Flist = xASL_adm_GetFileList(x.D.ROOT,'^.nii*$','List',[0 Inf], false);
bAnat = false;
NiftiList = {};

%Load T1W image
if contains(fieldnames(Flist), 'T1*')
	bAnat = true;
	NiftiList{1,1} = fullfile(x.D.ROOT, 'T1.nii.gz');
end

if contains(fieldnames(Flist), 'T2*')
	bAnat = true;
	NiftiList{end,1} = fullfile(x.D.ROOT, 'T2.nii.gz');
end

if contains(fieldnames(Flist), 'FLAIR*')
	bAnat = true;
	NiftiList{end,1} = fullfile(x.D.ROOT, 'Flair.nii.gz');
end

if ~(bAnat)
	warning('No anatomical scans found, cannot plot orientation.');
end

%Load ASL image
if contains(fieldnames(Dlist), 'ASL')
	
end

xASL_qc_PrintOrientation(NiftiList, x.D.ROOT, '');
SlabPlaced = [];
% =====================================================================================================
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

if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') % Backwards compatability, and also needed for the Affine+DCT co-registration of ASL-T1w
	AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
else
	AffineTransfPath = [];
end

xASL_spm_deformations(x, x.P.Path_SliceGradient, x.P.Pop_Path_SliceGradient, 0, [], AffineTransfPath, x.P.Path_y_ASL); % nearest neighbor

if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') % Backwards compatability, and also needed for the Affine+DCT co-registration of ASL-T1w
	AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
else
	AffineTransfPath = [];
end

xASL_spm_deformations(x, x.P.Path_SliceGradient_extrapolated, x.P.Pop_Path_SliceGradient_extrapolated,1,[], AffineTransfPath, x.P.Path_y_ASL);
% Since smoothing is desirable for this slice gradient, and spline edge artifacts are not, we use trilinear here as interpolation method

%% Add some more extrapolation to be sure
SGim = xASL_io_Nifti2Im(x.P.Pop_Path_SliceGradient_extrapolated);
SGim(SGim==0) = NaN;
% Count the number of NaNs
countNanLast = numel(SGim)+1;
while sum(isnan(SGim(:)))<countNanLast
	% And make sure that the number of NaNs is decreasing with the filtering, if not 
	countNanLast = sum(isnan(SGim(:)));
	SGim = xASL_im_ndnanfilter(SGim,'rect',[8 8 8],2);
end

xASL_io_SaveNifti(x.P.Pop_Path_SliceGradient_extrapolated,x.P.Pop_Path_SliceGradient_extrapolated,SGim,[],0);
PathFig = 'placeholder';
end
