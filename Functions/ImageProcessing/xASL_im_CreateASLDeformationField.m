function xASL_im_CreateASLDeformationField(x, bOverwrite, EstimatedResolution, PathLowResNIfTI)
%xASL_im_CreateASLDeformationField Adapt deformation field for lower resolution
%
% FORMAT:   xASL_im_CreateASLDeformationField(x, bOverwrite, EstimatedResolution)
%
% INPUT:
%   x                   - structure containing fields with all information required to run this submodule (REQUIRED)
%   bOverwrite          - true for overwriting any existing ASL deformation field (OPTIONAL, DEFAULT = false)
%   EstimatedResolution - X Y Z scalar with estimated effective spatial
%                         resolution of images that flowfield will be applied to (OPTIONAL,
%                         DEFAULT = obtain this automaticall)
%   PathLowResNIfTI     - path to the NIfTI file that deformations will be
%                         applied to (OPTIONAL, default=ASL4D)
% OUTPUT:                 n/a
% --------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function smooths a transformation flow field to a lower
%              resolution. Usually, we use a high resolution anatomical
%              image (e.g. {{3D T1w}}) to obtain the flowfields from native
%              space to standard space, and apply these to the lower
%              resolution ASL images. Because of the resolution
%              differences, the flowfields need to be downsampled/smoothed,
%              to avoid deformation effects that are crispier than the
%              functional image that is investigated. This function
%              performs the following steps:
%
%              1. Obtain resolutions
%              2. Fill NaNs at edges y_T1.nii flowfield to prevent interpolation artifact
%              3. Smooth flowfield
%              4. Fill NaNs at edges y_ASL.nii
% 
%              Note that if the resolution of ASL is not significantly (i.e. >0.5 mm in
%              any dimension) lower than T1w, the y_T1.nii is copied to y_ASL.nii
% --------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_im_CreateASLDeformationField(x);
%               xASL_im_CreateASLDeformationField(x, 1, [3 3 7]);
% __________________________________
% Copyright (C) 2015-2021 ExploreASL


%% Admin
if nargin<4 || isempty(PathLowResNIfTI)
    PathLowResNIfTI = x.P.Path_ASL4D;
end
if nargin<2 || isempty(bOverwrite)
    bOverwrite = false; % By default, this function will be skipped if the ASL deformation field already exists
end

if xASL_exist(x.P.Path_y_ASL,'file') && ~bOverwrite
    return; % skip this function
elseif ~xASL_exist(x.P.Path_y_T1,'file')
    warning([x.P.Path_y_T1 ' didnt exist, skipping...']);
    return;
end

%% 1) Obtain resolutions
if nargin<3 || isempty(EstimatedResolution)
	EstimatedResolution = xASL_init_DefaultEffectiveResolution(PathLowResNIfTI, x);
end

TemplateResolution = [1.5 1.5 1.5];
sKernel = (EstimatedResolution.^2 - TemplateResolution.^2).^0.5; % assuming the flow fields are in 1.5x1.5x1.5 mm
sKernel(EstimatedResolution<TemplateResolution) = 0;

%% 2) Fill NaNs at edges y_T1.nii flowfield to prevent interpolation artifact
xASL_im_FillNaNs(x.P.Path_y_T1, 3, x.settings.Quality); % First fill NaNs, to prevent interpolation artifacts

%% 3) Smooth flowfield
% Note that we provide the T1w resolution, not the resolution of y_T1,
% which is the standard space resolution (as it pulls not pushes).
niiT1w = xASL_io_ReadNifti(x.P.Path_T1);
resSrc = niiT1w.hdr.pixdim(2:4);

xASL_im_PreSmooth(PathLowResNIfTI, x.P.Path_y_T1, x.P.Path_y_ASL, [], resSrc); % we need to add the effective resolution here still!
% sKernel, as calculated above, can be used for this. But the major
% rotations need to be taken into account, between the effective
% resolution as specified, and the one in the different NIfTIs
% (e.g. the ASL & T1w are usually acquired transversal & sagittally,
% respectively

%% 4) Fill NaNs at edges y_ASL.nii
xASL_im_FillNaNs(x.P.Path_y_ASL, 3); % Again fill NaNs, if needed


end
