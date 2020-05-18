function [ x ] = xASL_init_PopulationSettings( x )
%xASL_init_PopulationSettings Prepares memory for population module
% ExploreASL 2018


if ~isfield(x,'Quality')
    x.Quality=1;
end

% xASL_adm_CreateDir( x.D.T1_ASLREGDIR);
% xASL_adm_CreateDir(x.D.DICOMparameterDir);

x.STRUCT_TEMPLATE_IM = fullfile(x.D.PopDir, 'Templates', 'pGM_bs-mean.nii');
pGM_exist   = true;
if ~xASL_exist(x.STRUCT_TEMPLATE_IM,'file')
    pGM_exist   = false;
elseif sum(sum(sum(xASL_io_Nifti2Im(x.STRUCT_TEMPLATE_IM))))==0
    pGM_exist   = false;
end
if pGM_exist
    fprintf('%s\n','Using population-specific pGM template');
else
    x.STRUCT_TEMPLATE_IM = fullfile(x.D.TemplateDir,'rc1T1.nii');
end

% x.GradualSkull = xASL_io_Nifti2Im(fullfile(x.D.MapsSPMmodifiedDir,
% 'rbrainmask.nii')); % not useed

% Load T1 background image
load_T1         = fullfile(x.D.PopDir, 'Templates', [x.P.STRUCT '_bs-mean.nii']);
T1_exist        = 1;
if     ~xASL_exist(load_T1,'file')
        T1_exist    = 0;
elseif  sum(sum(sum(xASL_io_Nifti2Im(load_T1))))==0
        T1_exist    = 0;
end
if  T1_exist
    fprintf('%s\n','Using population-specific structural template');
else
    load_T1 = fullfile(x.D.MapsSPMmodifiedDir, ['r' x.P.STRUCT '.nii']);
end


background_mask                     = xASL_io_Nifti2Im(load_T1);
MaxN                                = max(background_mask(:));
% Apply gradual masking
IM1MASK                             = background_mask .* x.GradualSkull + (1-x.GradualSkull).*MaxN; % with white background
% Put black line around mask
InvertMask                          = IM1MASK~=MaxN;
Erosion1                            = logical(InvertMask-xASL_im_DilateErodeFull(InvertMask,'erode',xASL_im_DilateErodeSphere(1)) );
Erosion2                            = logical(InvertMask-xASL_im_DilateErodeFull(InvertMask,'erode',xASL_im_DilateErodeSphere(2)) );
Erosion3                            = logical(InvertMask-xASL_im_DilateErodeFull(InvertMask,'erode',xASL_im_DilateErodeSphere(3)) );
Erosion4                            = logical(InvertMask-xASL_im_DilateErodeFull(InvertMask,'erode',xASL_im_DilateErodeSphere(4)) );
Erosion2(Erosion1)                  = 0;
Erosion3(Erosion1 | Erosion2)       = 0;
Erosion4(Erosion1 | Erosion2 | Erosion3)       = 0;
IM1MASK(Erosion1)                   = 0;
IM1MASK(Erosion2)                   = IM1MASK(Erosion2).*0.25;
IM1MASK(Erosion3)                   = IM1MASK(Erosion3).*0.5;
IM1MASK(Erosion4)                   = IM1MASK(Erosion4).*0.75;
% Erosion5 = 1
InvertMask                          = IM1MASK~=MaxN;
Erosion2                            = xASL_im_DilateErodeFull(InvertMask,'erode',xASL_im_DilateErodeSphere(2));
IM1MASK(~Erosion2)                  = MaxN;




% IM1MASK             = background_mask .* x.GradualSkull; % with black background

%         IM1_outside         = background_mask.*~x.skull;
%
%         sort_IM1            = sort(nonzeros(IM1_outside));
%         Value9              = sort_IM1(round(0.999*length(sort_IM1)));
%         IM1_outside(IM1_outside>Value9)     = Value9;
%         IM1_outside         = IM1_outside ./ max(nonzeros(IM1_outside)) .* 1500;
%
%         sort_IM1            = sort(nonzeros(IM1MASK));
%         Value9              = sort_IM1(round(0.999*length(sort_IM1)));
%         IM1MASK(IM1MASK>Value9)     = Value9;

%% Figure used in Sleep Paper
% x.S.CorSlices   = [];
% x.S.SagSlices   = [];
% x.S.TraSlices   = 34:4:34+15*4; % 34:4:34+17*4 was original, now removed 2 slices

background_mask                     = (IM1MASK ./ max(nonzeros(IM1MASK(:))) ) .* 1750;
background_mask                     = background_mask ./ max( background_mask(:) );

DATA_OUT                            = xASL_im_TransformData2View( background_mask, x );
x.background_view_clr               = double(DATA_OUT ./ (max(DATA_OUT(:)) ));
x.background_view_clr               = repmat(x.background_view_clr,[1 1 3]);
% not used





end
