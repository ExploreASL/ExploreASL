function [EffectiveResolution] = xASL_init_DefaultEffectiveResolution(PathASL, x)
%xASL_init_DefaultEffectiveResolution Educated-guess ASL effective resolution
%
% FORMAT: [EffectiveResolution] = xASL_init_DefaultEffectiveResolution(PathASL, x)
%
% INPUT:
%   PathASL       - path to NIfTI containing raw ASL timeseries (REQUIRED)
%   x             - x structure containing all input parameters (REQUIRED)
%   x.Q.Sequence  - type of ASL readout used (2D_EPI, 3D GRASE or 3D spiral)
%
% OUTPUT:
%   EffectiveResolution     - mm FWHM of estimated Gaussian kernel
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL module provides an educated guess on
% the effective spatial resolution of ASL. This may depend on the
% combination of acquisition PSF, reconstruction filter, head motion.
% Note that the derived/processed images may have a lower effective
% resolution because of smoothing effects from interpolation. Note that
% this remains an educated-guess, the actual FWHM may still differ,
% especially for 3D GRASE sequences, where e.g. the choice of number of
% segments can affect the smoothness.
%
% This function conducts the following steps:
%
% 1. Educated-guess FWHM
% 2. Attempt accounting for in-plane interpolation in reconstruction
% 3. Calculate and report effective spatial resolution
%
% EXAMPLE: EffectiveResolution = xASL_init_DefaultEffectiveResolution('/MyStudy/sub-001/ASL_1/ASL4D.nii.gz, x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: Petr, 2018 MAGMA; Vidorreta 2013 Neuroimage
%
% Copyright 2015-2020 ExploreASL

%% ----------------------------------------------------------------------------------------
%% Admin
tIM = xASL_io_ReadNifti(PathASL);
NativeResolution = tIM.hdr.pixdim(2:4);
if ~isfield(x.Q,'Sequence')
    warning(['Setting x.Q.Sequence missing, skipping']);
    return;
end


%% ----------------------------------------------------------------------------------------
%% 1) Educated-guess FWHM
switch lower(x.Q.Sequence)
    case '2d_epi'
        % for 2D EPI the PSF is effectively negligible
        Estimated_FWHM = [1 1 1];
    case '3d_spiral_wip'
        % educated guess: 3D spiral acquision has large through-plane PSF
        Estimated_FWHM = [1.1158 1.1289 1.9525]; % average between 3D_spiral & 3D_GRASE
    case '3d_spiral'
        % estimation by Petr, MAGMA
        % Paper Jan Petr, 4.9 & 5.1 where X Y, which should be similar
        % difference with 3D_spiral WIP is that this product sequence applies a
        % Fermi filter in reconstruction
        EstimatedEffectiveResolution = [4.3 4.4 10.1];
        BasedOnResolution = [3.8 3.8 4];
        Estimated_FWHM = EstimatedEffectiveResolution./BasedOnResolution;
    case '3d_grase'
        Estimated_FWHM = [1.1 1.1 1.38]; % KISS, Vidoretta, NeuroImage 2013
        % The in-plane 2 translates the reconstruction resolution into
        % acquisition resolution. A lot still depends on the actual readout
        % used, this is perhaps the most variable sequence in terms of PSF
        % FWHM. Roughly speaking, the 3d_grase PSF should be somewhere between the 2d_epi
        % and 3d_spiral.
    otherwise
        warning(['Unknown setting x.Q.Sequence=' xASL_num2str(x.Q.Sequence)]);
end

%% ----------------------------------------------------------------------------------------
%% 2) Attempt accounting for in-plane interpolation in reconstruction
if ~isempty(regexpi(x.Q.Sequence,'3D_spiral')) && NativeResolution(1)<4
        NativeResolution(1:2) = max(NativeResolution(1:2),[3.75 3.75]);
        % to account for in-plane interpolation
elseif ~isempty(regexpi(x.Q.Sequence,'3D_GRASE')) && NativeResolution(1)<3
        NativeResolution(1:2) = max(NativeResolution(1:2),[3.8 3.8]);
end

%% ----------------------------------------------------------------------------------------
%% 3) Calculate and report effective spatial resolution
EffectiveResolution = NativeResolution.*Estimated_FWHM;
fprintf('%s\n',[x.Q.Sequence ' NIfTI has native resolution ' num2str(NativeResolution(1)) ' ' num2str(NativeResolution(2)) ' ' num2str(NativeResolution(3)) ', assuming PSF ' num2str(Estimated_FWHM(1)) ' ' num2str(Estimated_FWHM(2)) ' ' num2str(Estimated_FWHM(3)) ' this gives estimated effective resolution ' num2str(EffectiveResolution(1)) ' ' num2str(EffectiveResolution(1)) ' ' num2str(EffectiveResolution(3))])



end
