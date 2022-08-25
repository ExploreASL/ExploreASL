function xASL_bids_parseM0(pathASLNifti)
%xASL_bids_parseM0 Check the ASL file in Legacy format, but with BIDS sidecars for its M0 possibilities and 
% finish the conversion of the ASL NIFTI to the ExploreASL legacy format
%
% FORMAT:      xASL_bids_parseM0(pathASLNifti)
%
% INPUT:       pathASLNifti - Path to a ASL NIFTI file in legacy filename and BIDS sidecars (CHAR ARRAY, REQUIRED)
% 
% OUTPUT:      n/a
%
% DESCRIPTION: Check the .JSON and aslContext.tsv sidecards of an ASL file in BIDS format and find the 
%              specified M0 possibilities. Then it converts the ASL file to ExploreASL legacy format including 
%              splitting of ASL and M0 NIFTIes if needed. Note that the sidecars are in BIDS, but the file-structure
%              is already expected to be in Legacy format
%
% EXAMPLE:     xASL_bids_parseM0('/test/ASL4D.nii')
%
% __________________________________
% Copyright 2015-2021 ExploreASL


%% Check input
if nargin < 1 || isempty(pathASLNifti)
    warning('Empty ASL path...');
    return
end

% Verify that pathASLNifti leads to ASL Nifti
[Fpath, Ffile] = xASL_fileparts(pathASLNifti);
PathJSON = fullfile(Fpath, [Ffile '.json']);

%% Parse & process M0 options

PathM0 = fullfile(Fpath, 'M0.nii');

JSON = spm_jsonread(PathJSON);
if isfield(JSON, 'M0Type')
    switch JSON.M0Type
        %% Option 1. Included
		% The M0 file is included in the ASL timeseries as specified in the aslContext
        case 'Included'
            PathContext = fullfile(Fpath, [Ffile 'context.tsv']);
            if exist(PathContext, 'file')
               TSV = xASL_tsvRead(PathContext);
			   % Checks for the m0scan in ASLContext
               M0Index = find(strcmp(TSV, 'm0scan'))-1;
               if ~isempty(M0Index)
				   % If specified to remove the Dummy ASL scans, then remove them while splitting ASL to ASL and M0
				   fprintf('%s\n', '"M0Type":"Included", splitting ASL and M0 NIfTIs');
                   if isfield(JSON,'DummyScanPositionInASL4D') && ~isempty(JSON.DummyScanPositionInASL4D)
					   xASL_io_SplitASL(pathASLNifti, M0Index,JSON.DummyScanPositionInASL4D);
				   else
					   xASL_io_SplitASL(pathASLNifti, M0Index);
				   end
               else
                   warning(['M0Index missing in ' PathContext]);
               end
            else
                warning([PathContext ' missing']);
            end
        %% Option 2. Separate
		% M0 image is already provided as a separate image - checks for its existence
        case 'Separate'
			if ~xASL_exist(fullfile(PathM0))
                warning(['Missing: ' PathM0]);
			else
				JSON.M0 = 'separate_scan';
				JSON = rmfield(JSON,'M0Type');
				spm_jsonwrite(PathJSON, JSON);
			end
        %% Option 3. Estimate
		% If a single M0 values is provided then keep it, just rename the field
        case 'Estimate'
            if isfield(JSON, 'M0Estimate')
                JSON.M0 = JSON.M0Estimate;
				JSON = rmfield(JSON,'M0Type');
                spm_jsonwrite(PathJSON, JSON);
            else
                warning(['Field M0_value missing in ' PathJSON]);
            end
        %% Option 4. Absent
        % (in this case we use the control image as pseudo-M0, 
        % which happens in xASL_wrp_Resample
        case {'use_control_as_m0', 'UseControlAsM0'}
			JSON.M0 = 'UseControlAsM0';
			JSON = rmfield(JSON,'M0Type');
			spm_jsonwrite(PathJSON, JSON);
			if ~isfield(JSON, 'BackgroundSuppression')
				warning(['M0 is defined as UserControlAsM0, but missing field BackgroundSuppression in ' PathJSON]);
			elseif JSON.BackgroundSuppression == true && ~isfield(JSON,'BackgroundSuppressionPulseTime')
				warning(['M0 is defined as UserControlAsM0, BackgroundSuppression is ON, but BSup timings are undefined in ' PathJSON]);
			end
		case 'Absent'
			JSON.M0 = 'Absent';
			JSON = rmfield(JSON,'M0Type');
				
			if xASL_exist(fullfile(PathM0))
                    warning('"M0Type":"Absent" but separate M0 NIfTI detected, consider setting "M0Type":"Separate"');
			elseif isfield(JSON, 'BackgroundSuppression') && (~JSON.BackgroundSuppression || isfield(JSON,'BackgroundSuppressionPulseTime'))
				NIfTI_ASL = xASL_io_ReadNifti(pathASLNifti);
				if size(NIfTI_ASL.dat,4)>1
					% background suppression was off, so we check if we have
					% multiple ASL volumes for using mean control as pseudo-M0
					warning('Multiple ASL volumes detected, setting M0 to "UseControlAsM0" instead of "Absent"');
					JSON.M0 = 'UseControlAsM0';
				end
			end
			
			spm_jsonwrite(PathJSON, JSON);
        otherwise
            error(['Invalid M0Type in ' PathJSON]);
        end
else
    warning(['Field M0Type missing in ' PathJSON]);
end
        
end
