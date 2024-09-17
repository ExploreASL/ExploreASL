function scaleFactor = xASL_adm_GetPhilipsScaling(parms, header)
% Checks the nifti header and the additional parameters in a side-car for a nifti file and extract correction
% scaling factors for Philips.
%
% FORMAT: scaleFactor = xASL_adm_GetPhilipsScaling(parms, header)
% 
% INPUT:
%   parms    - parms read from parms.mat and JSON (REQUIRED)
%   header   - the nifti file header (REQUIRED)
%
% OUTPUT:
% scaleFactor       - a multiplicate factor - multiply the NiFTI image with this number to get the correct scaling
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script provides the correct scaling factors for a NIfTI file. It checks the header of the NIfTI
%              that normally has the same scaling as RescaleSlope in DICOM, it checks if dcm2nii (by the info in JSON)
%              has already converted the scale slopes to floating point. And if not, the derive the correct
%              scaling factor to be applied.
%              The function works with Philips-specific scale-slopes:
%              RWVSlope = Real world value slope (0040,9225) 
%              MRScaleSlope (2005,100e), (2005,110e), (2005,120e)
%              RescaleSlopeOriginal (2005,0x140a), (2005,110a)
%              These are different from the standard slopes in DICOM that are used as display slopes in Philips:
%              RescaleSlope (0028,1053)
%              With new version of dcm2niix-20220720 and ExploreASL-1.10.0, dcm2niix does RWVSlope scalings correctly and we do not have to fix that.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: scaleFactor = xASL_adm_GetPhilipsScaling('ASL4D_parms.mat','ASL.json');
% __________________________________
% Copyright 2015-2022 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

if nargin < 2 || isempty(parms) || isempty(header)
	error('Need 2 input parameters');
	
end
% Read the NIFTI header to check the scaling there
rescaleSlopeNifti = header.dat.scl_slope;
% bApplyScaling should be set to true if we need to apply the scalings, otherwise to false
bApplyScaling = true;
if isfield(parms,'RWVSlope')
	% Obtains dcm2niix version used for conversion and ExploreASL version
	% If ExploreASL version >= 1.10.0 or dcm2niix version >= 20210317
	if isfield(parms, 'ConversionSoftware') && (strcmp(parms.ConversionSoftware, 'dcm2nii') || strcmp(parms.ConversionSoftware, 'dcm2niix')) && isfield(parms, 'ConversionSoftwareVersion')
		dcm2niixVersion = regexp(parms.ConversionSoftwareVersion,'v\d*.\d*.(\d{8}).*','tokens');
		dcm2niixVersion = str2num(dcm2niixVersion{1}{1});
		%dcm2niixVersion = str2num(dcm2niixVersion(end-7:end));
	else
		dcm2niixVersion = 0;
	end
	
	% The correct dcm2niix version has to be used that is able to handle this parameter, and dcm2niix also needs to report the scalings were fixed
	if dcm2niixVersion >= 20210317 && isfield(parms, 'UsePhilipsFloatNotDisplayScaling') && (parms.UsePhilipsFloatNotDisplayScaling == 1)
		scaleFactor = 1;
		bApplyScaling = false;
	else
		% If the RealWorldValue is present, then dcm2nii scales to them and ignores everything else
		fprintf('=========================================== Scaling ==========================================\n');
		if numel(parms.RWVSlope)>1
			[~,idx] = min(abs(parms.RWVSlope - rescaleSlopeNifti));
			parms.RWVSlope = parms.RWVSlope(idx);
			if isfield(parms, 'RescaleSlopeOriginal') && numel(parms.RescaleSlopeOriginal)>1
				parms.RescaleSlopeOriginal = parms.RescaleSlopeOriginal(idx);
			end
			
			if isfield(parms, 'RescaleSlope') && numel(parms.RescaleSlope)>1
				parms.RescaleSlope = parms.RescaleSlope(idx);
			end
			
			if isfield(parms, 'MRScaleSlope') && numel(parms.MRScaleSlope)>1
				parms.MRScaleSlope = parms.MRScaleSlope(idx);
			end
		end
		
		if ~isempty(parms.RWVSlope)
			if ~xASL_stat_IsEqualTol(parms.RWVSlope, rescaleSlopeNifti, parms.RWVSlope/100) && (rescaleSlopeNifti ~= 1)
				fprintf('%s\n', ['RWVSlope (' xASL_num2str(parms.RWVSlope) ') and NIfTI slope (' xASL_num2str(rescaleSlopeNifti) ') differ, using RWVSlope']);
			end
			if isfield(parms,'RescaleSlopeOriginal') && ~xASL_stat_IsEqualTol(parms.RWVSlope, parms.RescaleSlopeOriginal, parms.RWVSlope/100) && (rescaleSlopeNifti ~= 1)
				fprintf('%s\n', ['RWVSlope (' xASL_num2str(parms.RWVSlope) ') and RescaleSlopeOriginal (' xASL_num2str(parms.RescaleSlopeOriginal) ') differ, using RWVSlope']);
			end
			if isfield(parms,'RescaleSlope') && ~xASL_stat_IsEqualTol(parms.RWVSlope, parms.RescaleSlope, parms.RWVSlope/100) && (rescaleSlopeNifti ~= 1)
				fprintf('%s\n', ['RWVSlope (' xASL_num2str(parms.RWVSlope) ') and RescaleSlope (' xASL_num2str(parms.RescaleSlope) ') differ, using RWVSlope']);
			end
		else
			fprintf('Warning: RWVSlope value is empty...\n');
		end
		
		if parms.RWVSlope == 1
			warning('RWVSlope was 1, could be a scale slope issue');
		end
		
		% Set the correction scaling factor
		scaleFactor = parms.RWVSlope;
	end
			
elseif isfield(parms, 'UsePhilipsFloatNotDisplayScaling') && (parms.UsePhilipsFloatNotDisplayScaling == 1)
	% Without the RWVSlope set and with the UsePhilipsFloatNotDisplayScaling set to 1, we know that dcm2nii did the scaling correctly and we don't have to further scale anything
	scaleFactor = 1;
	bApplyScaling = false;
elseif (~isfield(parms, 'RescaleSlope')) && (~isfield(parms, 'MRScaleSlope')) && (~isfield(parms, 'RescaleSlopeOriginal')) && (~isfield(parms, 'UsePhilipsFloatNotDisplayScaling'))
	% All fields are missing - we can ignore the quantification as all the fields have been removed and scaling done properly (we assume)
	% Note that RWVSlope was ruled out previously, and if UsePhilipsFloatNotDisplayScaling exists, it must be 0
	scaleFactor = 1;
	bApplyScaling = false;
else
	% Standard scaling using RescaleSlope/RescaleSlopeOriginal or the Nifti slope - which should all be either non-existing or 1 or all equal
	% Set the scaling to remove
	scaleFactor = 1;
	
	if (rescaleSlopeNifti ~= 1)
		scaleFactor = rescaleSlopeNifti;
	end
	
	if isfield(parms, 'RescaleSlopeOriginal') && (parms.RescaleSlopeOriginal ~= 1)
		if (scaleFactor ~= 1) && ~xASL_stat_IsEqualTol(scaleFactor, parms.RescaleSlopeOriginal, scaleFactor/100)
			warning('%s\n%s\n',...
				['Discrepancy in RescaleSlopeOriginal (' xASL_num2str(parms.RescaleSlopeOriginal) ') and NIFTI Slopes (' xASL_num2str(rescaleSlopeNifti) ')'],...
				'Using RescaleSlopeOriginal');
		end
		scaleFactor = parms.RescaleSlopeOriginal;
	end
	
	if isfield(parms, 'RescaleSlope') && (parms.RescaleSlope ~= 1)
		if (scaleFactor ~= 1) && ~xASL_stat_IsEqualTol(scaleFactor, parms.RescaleSlope, scaleFactor/100)
			if (rescaleSlopeNifti ~= 1) && ~xASL_stat_IsEqualTol(rescaleSlopeNifti, parms.RescaleSlope, scaleFactor/100)
				warning('%s\n%s\n',...
				['Discrepancy in RescaleSlope (' xASL_num2str(parms.RescaleSlope) ') and NIFTI Slopes (' xASL_num2str(rescaleSlopeNifti) ')'],...
				'Using RescaleSlope');
			else
				warning('%s\n%s\n',...
				['Discrepancy in RescaleSlope (' xASL_num2str(parms.RescaleSlope) ') and RescaleSlopeOriginal (' xASL_num2str(parms.RescaleSlopeOriginal) ')'],...
				'Using RescaleSlope');
			end
		end
		scaleFactor = parms.RescaleSlope;
	end
		
	if scaleFactor == 1
		warning('Scale slope was 1, could be a scale slope issue');
	end
end
% In case the scaleFactor was read as zero from DICOM or JSON, then don't apply it
if scaleFactor == 0
	bApplyScaling = false;
	scaleFactor = 1;
end
% Apply the correct scaling
if bApplyScaling
	if ~isfield(parms,'MRScaleSlope')
        warning('MRScaleSlope missing, potential quantification error, skipping');
        return;
	end
	% Make a local copy to potentially modify
    MRScaleSlope = parms.MRScaleSlope;
	
	if numel(MRScaleSlope) > 1
		ssInd = find(MRScaleSlope ~= 1);
		if isempty(ssInd)
			MRScaleSlope = 1;
		else
			if numel(ssInd) > 1
				warning('Philips MR ScaleSlope has more than a single value, could be a scale slope issue');
			end
			MRScaleSlope = MRScaleSlope(ssInd);
		end
	end
	if MRScaleSlope == 1
		warning('Philips MR ScaleSlope was 1, could be a scale slope issue.');
	end
	
	fprintf('%s\n',['Using DICOM (re)scale slopes ' xASL_num2str(scaleFactor) ' * ' xASL_num2str(MRScaleSlope)]);
	scaleFactor = 1./(scaleFactor .* MRScaleSlope);
end
% Set scaleFactor to zero if xASL_adm_GetPhilipsScaling failed?
if isempty(scaleFactor)
    warning('xASL_adm_GetPhilipsScaling failed, scaleFactor will be set to 1...');
    scaleFactor = 1;
end
% Linebreak after printing of user feedback
fprintf('\n');
end
