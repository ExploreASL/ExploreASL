function [EffectiveResolution] = xASL_init_DefaultEffectiveResolution(PathASL, x)
%xASL_init_DefaultEffectiveResolution Educated-guess ASL effective resolution
%
% FORMAT: [EffectiveResolution] = xASL_init_DefaultEffectiveResolution(PathASL, x)
%
% INPUT:
%   PathASL               - path to NIfTI containing raw ASL timeseries (REQUIRED)
%   x                     - x structure containing all input parameters (REQUIRED)
%   x.Q.PulseSequenceType - type of ASL readout used (EPI, GRASE or spiral)
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
% 1. Basic cases when resolution is not calculated
% 2. Educated-guess FWHM
% 3. Attempt accounting for in-plane interpolation in reconstruction
% 4. Calculate and report effective spatial resolution
%
% EXAMPLE: EffectiveResolution = xASL_init_DefaultEffectiveResolution('/MyStudy/sub-001/ASL_1/ASL4D.nii.gz, x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: Petr, 2018 MAGMA; Vidorreta 2013 Neuroimage
%
% Copyright 2015-2024 ExploreASL

%% ----------------------------------------------------------------------------------------
%% Admin
[tIM, ~, tJson] = xASL_io_ReadNifti(PathASL);
NativeResolution = tIM.hdr.pixdim(2:4);

% Check for missing JSON
if isempty(tJson) || ~isfield(tJson, 'Q') 
	error(['Missing JSON sidecar for ' PathASL]);
end

% Obtain the x.Q.PulseSequenceType and x.Q.MRAcquisition fields from the json
x.Q.PulseSequenceType = tJson.Q.PulseSequenceType;
x.Q.MRAcquisitionType = tJson.Q.MRAcquisitionType;

% Also make sure we read the Vendor B0-strength, number of spirals, and AcquisitionVoxelSize if defined
x.Q.Vendor = tJson.Q.Vendor;
x.MagneticFieldStrength = tJson.MagneticFieldStrength;
if isfield(tJson, 'NumberOfArms')
	x.Q.NumberOfArms = tJson.NumberOfArms;
else
	x.Q.NumberOfArms = [];
end

if isfield(tJson, 'AcquisitionVoxelSize')
	x.Q.AcquisitionVoxelSize = tJson.AcquisitionVoxelSize;
else
	x.Q.AcquisitionVoxelSize = [];
end

%% Obtain the effective resolution
%% 1. Basic cases when resolution is not calculated
if ~isempty(x.Q.AcquisitionVoxelSize) && numel(x.Q.AcquisitionVoxelSize) == 3
	% If estimated effective resolution is given in JSON-sidecar, then we use it
	EffectiveResolution = x.Q.AcquisitionVoxelSize;
	fprintf('%s\n', ['Assume effective resolution ' num2str(EffectiveResolution(1)) ' ' num2str(EffectiveResolution(1)) ' ' num2str(EffectiveResolution(3)) ' as provided in the JSON-sidecar'])
elseif x.MagneticFieldStrength > 3
	% For more than 3T scanners, we don't have this estimated as the resolutions are much higher
	EffectiveResolution = NativeResolution;
	fprintf('%s\n', ['Assume effective resolution equal to voxel size ' num2str(EffectiveResolution(1)) ' ' num2str(EffectiveResolution(1)) ' ' num2str(EffectiveResolution(3))]);
elseif ~isfield(x.Q, 'PulseSequenceType') || ~isfield(x.Q, 'MRAcquisitionType') || ~isfield(x.Q, 'Vendor')
	% If sequence or vendor are still missing we skip this function
	error('Settings of x.Q.PulseSequenceType or x.Q.MRAcquisitionType or x.Q.Vendor are missing');
else	
	%% ----------------------------------------------------------------------------------------
	%% Use precalculated values
	%% 1) Educated-guess FWHM
	if strcmpi(x.Q.PulseSequenceType, 'EPI') && strcmpi(x.Q.MRAcquisitionType, '2D')
		% for 2D EPI the PSF is effectively negligible
		Estimated_FWHM = [1 1 1];
	elseif strcmpi(x.Q.PulseSequenceType, 'spiral_wip') && strcmpi(x.Q.MRAcquisitionType, '3D')
		% educated guess: 3D spiral acquision has large through-plane PSF
		Estimated_FWHM = [1.1158 1.1289 1.9525]; % average between 3D_spiral & 3D_GRASE
	elseif strcmpi(x.Q.PulseSequenceType, 'spiral') && strcmpi(x.Q.MRAcquisitionType, '3D')
		% estimation by Petr, MAGMA
		% Paper Jan Petr, 4.9 & 5.1 where X Y, which should be similar
		% difference with 3D_spiral WIP is that this product sequence applies a
		% Fermi filter in reconstruction
		EstimatedEffectiveResolution = [4.3 4.4 10.1];
		BasedOnResolution = [3.8 3.8 4];
		Estimated_FWHM = EstimatedEffectiveResolution./BasedOnResolution;
	elseif strcmpi(x.Q.PulseSequenceType, 'GRASE') && strcmpi(x.Q.MRAcquisitionType, '3D')
		Estimated_FWHM = [1.1 1.1 1.38]; % KISS, Vidoretta, NeuroImage 2013
		% The in-plane 2 translates the reconstruction resolution into
		% acquisition resolution. A lot still depends on the actual readout
		% used, this is perhaps the most variable sequence in terms of PSF
		% FWHM. Roughly speaking, the 3d_grase PSF should be somewhere between the 2d_epi
		% and 3d_spiral.
	else
		warning(['Unknown setting x.Q.PulseSequenceType=' xASL_num2str(x.Q.PulseSequenceType) ', x.Q.MRAcuqisitionType=' xASL_num2str(x.Q.MRAcquisitionType)]);
	end

	%% ----------------------------------------------------------------------------------------
	%% 2) Attempt accounting for in-plane interpolation in reconstruction
	if strcmpi(x.Q.PulseSequenceType, 'spiral') && strcmpi(x.Q.MRAcquisitionType, '3D')
		% GE tends to upsample their spiral acquisitions to 1.6-1.9mm voxels
		% For non-GE, or in-plane resolution higher than 2mm, we can't assume that the reconstruction was upsampled and we leave the native resolution intact
		if regexpi(x.Q.Vendor, 'GE') && NativeResolution(1) < 2
			% The individual resolution will also depend on the FOV, so we cannot set a fixed resolution but rather modify it based on the native resolution
			if isempty(x.Q.NumberOfArms)
				warning('Number of spirals not defined for a GE sequence');
			elseif x.Q.NumberOfArms == 8
				% For the standard number of spirals, we use 2 times the native resolution
				NativeResolution(1:2) = 2 * NativeResolution(1:2);
			elseif x.Q.NumberOfArms == 4
				NativeResolution(1:2) = 3 * NativeResolution(1:2);
			else
				warning('Number of spirals for a GE sequence is not in the common range of 4 or 8')
			end
		end
	elseif strcmpi(x.Q.PulseSequenceType, 'GRASE') && strcmpi(x.Q.MRAcquisitionType, '3D') && NativeResolution(1)<3
		NativeResolution(1:2) = max(NativeResolution(1:2),[3.8 3.8]);
	end

	%% ----------------------------------------------------------------------------------------
	%% 3) Calculate and report effective spatial resolution
	EffectiveResolution = NativeResolution.*Estimated_FWHM;
	fprintf('%s\n',[x.Q.PulseSequenceType ' NIfTI has native resolution ' num2str(NativeResolution(1)) ' ' num2str(NativeResolution(2)) ' ' num2str(NativeResolution(3)) ', assuming PSF ' num2str(Estimated_FWHM(1)) ' ' num2str(Estimated_FWHM(2)) ' ' num2str(Estimated_FWHM(3)) ' this gives estimated effective resolution ' num2str(EffectiveResolution(1)) ' ' num2str(EffectiveResolution(1)) ' ' num2str(EffectiveResolution(3))])
end

end
