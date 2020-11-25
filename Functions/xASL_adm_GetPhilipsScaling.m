function scaleFactor = xASL_adm_GetPhilipsScaling(parms,header)
% Checks the nifti header and the additional parameters in a side-car for a nifti file and extract correction
% scaling factors for Philips.
%
% FORMAT: scaleFactor = xASL_adm_GetPhilipsScaling(pathParmsMat,pathNifti)
% 
% INPUT:
%   parms    - parms read from parms.mat and JSON (REQUIRED)
%   header   - the nifti file header(REQUIRED)
%
% OUTPUT:
% scaleFactor       - a multiplicate factor - multiply the NiFTI image with this number to get the correct scaling
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script provides the correct scaling factors for a NIfTI file. It checks the header of the NIfTI
%              that normally has the same scaling as RescaleSlope in DICOM, it checks if dcm2nii (by the info in JSON)
%              has already converted the scale slopes to floating point. And if not, the derive the correct
%              scaling factor to be applied.
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: scaleFactor = xASL_adm_GetPhilipsScaling('ASL4D_parms.mat','ASL.json');
% __________________________________
% Copyright 2015-2020 ExploreASL

if nargin < 2 || isempty(parms) || isempty(header)
	error('Need 2 input parameters');
	
end

% Read the NIFTI header to check the scaling there
rescaleSlopeNifti = header.dat.scl_slope;

if isfield(parms,'RWVSlope')
	% If the RealWorldValue is present, then dcm2nii scales to them and ignores everything else
	if length(parms.RWVSlope)>1
		[~,idx] = min(abs(parms.RWVSlope - rescaleSlopeNifti));
		parms.RWVSlope = parms.RWVSlope(idx);
		if isfield(parms,'RescaleSlopeOriginal') && length(parms.RescaleSlopeOriginal)>1
			parms.RescaleSlopeOriginal = parms.RescaleSlopeOriginal(idx);
		end
		
		if isfield(parms,'RescaleSlope') && length(parms.RescaleSlope)>1
			parms.RescaleSlope = parms.RescaleSlope(idx);
		end
		
		if isfield(parms,'MRScaleSlope') && length(parms.MRScaleSlope)>1
			parms.MRScaleSlope = parms.MRScaleSlope(idx);
		end
	end
	
	if ~isnear(parms.RWVSlope,rescaleSlopeNifti,parms.RWVSlope/100) && (rescaleSlopeNifti ~= 1)
		fprintf('%s\n', ['RWVSlope (' xASL_num2str(parms.RWVSlope) ') and NIfTI slope (' xASL_num2str(rescaleSlopeNifti) ') differ, using RWVSlope']);
	end
	if isfield(parms,'RescaleSlopeOriginal') && ~isnear(parms.RWVSlope,parms.RescaleSlopeOriginal,parms.RWVSlope/100) && (rescaleSlopeNifti ~= 1)
		fprintf('%s\n', ['RWVSlope (' xASL_num2str(parms.RWVSlope) ') and RescaleSlopeOriginal (' xASL_num2str(parms.RescaleSlopeOriginal) ') differ, using RWVSlope']);
	end
	if isfield(parms,'RescaleSlope') && ~isnear(parms.RWVSlope,parms.RescaleSlope,parms.RWVSlope/100) && (rescaleSlopeNifti ~= 1)
		fprintf('%s\n', ['RWVSlope (' xASL_num2str(parms.RWVSlope) ') and RescaleSlope (' xASL_num2str(parms.RescaleSlope) ') differ, using RWVSlope']);
	end
			
	if parms.RWVSlope == 1
		warning('RWVSlope was 1, could be a scale slope issue');
	end
			
	% Set the scaling to remove
	scaleFactor = parms.RWVSlope;
			
elseif isfield(parms,'UsePhilipsFloatNotDisplayScaling') && (parms.UsePhilipsFloatNotDisplayScaling == 1)
	% Without the RWVSlope set and with the UsePhilipsFloatNotDisplayScaling set to 1, we know that dcm2nii did the scaling correctly and we don't have to further scale anything
	scaleFactor = 0;
elseif (~isfield(parms,'RescaleSlope')) && (~isfield(parms,'MRScaleSlope')) && (~isfield(parms,'RescaleSlopeOriginal')) && (~isfield(parms,'UsePhilipsFloatNotDisplayScaling'))
	% All fields are missing - we can ignore the quantification as all the fields have been removed and scaling done properly (we assume)
	% Note that RWVSlope was ruled out previously, and if UsePhilipsFloatNotDisplayScaling exists, it must be 0
	scaleFactor = 0;
else
	% Standard scaling using RescaleSlope/RescaleSlopeOriginal or the Nifti slope - which should all be either non-existing or 1 or all equal
	% Set the scaling to remove
	scaleFactor = 1;
	
	if (rescaleSlopeNifti ~= 1)
		scaleFactor = rescaleSlopeNifti;
	end
	
	if isfield(parms,'RescaleSlopeOriginal') && (parms.RescaleSlopeOriginal ~= 1)
		if (scaleFactor ~= 1) && ~isnear(scaleFactor,parms.RescaleSlopeOriginal,scaleFactor/100)
			warning('%s\n%s\n',...
				['Discrepancy in RescaleSlopeOriginal (' xASL_num2str(parms.RescaleSlopeOriginal) ') and NIFTI Slopes (' xASL_num2str(rescaleSlopeNifti) ')'],...
				'Using RescaleSlopeOriginal');
		end
		scaleFactor = parms.RescaleSlopeOriginal;
	end
	
	if isfield(parms,'RescaleSlope') && (parms.RescaleSlope ~= 1)
		if (scaleFactor ~= 1) && ~isnear(scaleFactor,parms.RescaleSlope,scaleFactor/100)
			if (rescaleSlopeNifti ~= 1) && ~isnear(rescaleSlopeNifti,parms.RescaleSlope,scaleFactor/100)
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
		
% Apply the correct scaling
if scaleFactor
    if ~isfield(parms,'MRScaleSlope')
        warning('MRScaleSlope missing, potential quantification error, skipping');
        return;
	end
	% Make a local copy to potentially modify
    MRScaleSlope = parms.MRScaleSlope;
	
	if length(MRScaleSlope) > 1
		ssInd = find(MRScaleSlope ~= 1);
		if isempty(ssInd)
			MRScaleSlope = 1;
		else
			if length(ssInd) > 1
				warning('Philips MR ScaleSlope has more than a single value, could be a scale slope issue');
			end
			MRScaleSlope = MRScaleSlope(ssInd);
		end
	end
	if MRScaleSlope == 1
		warning('Philips MR ScaleSlope was 1, could be a scale slope issue.');
	end
	
	fprintf('%s\n',['Using DICOM (re)scale slopes ' num2str(scaleFactor) ' * ' num2str(MRScaleSlope)]);
	scaleFactor = 1./(scaleFactor .* MRScaleSlope);
end


end