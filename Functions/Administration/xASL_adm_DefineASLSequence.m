function [x] = xASL_adm_DefineASLSequence(x, bVerbose)
%xASL_adm_DefineASLSequence Obtain ASL sequence type for sequence-specific
%image processing
%
% FORMAT: [x] = xASL_adm_DefineASLSequence(x)
%
% INPUT:
%   x                - x structure containing all input parameters (REQUIRED)
%   x.Q.readoutDim   - dimensionality of readout (2D or 3D) (OPTIONAL)
%   x.Q.Vendor       - Either 'GE', 'Philips', 'Siemens' (OPTIONAL)
% OUTPUT:
%   x               - x structure containing all output parameters
%   x.Q.Sequence    - sequence type (readout)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL function tries to check what ASL sequence is
% being processed, if this was not already defined in x.Q.Sequence.
% It does so by checking known combinations of readout dimensionality
% (x.Q.readoutDim) and Vendor, knowing the product sequences of the Vendors.
%
% EXAMPLE: x = xASL_adm_DefineASLSequence(x);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2023 ExploreASL


%% Check quantification fields readoutDim & x.Q.Sequence
if ~isfield(x.Q, 'readoutDim') || isempty(x.Q.readoutDim)
    warning('x.Q.readoutDim parameter missing');
end

% Check for illegal sequence definitions
if isfield(x.Q,'Sequence')
    if isempty(regexpi(x.Q.Sequence, '(2d_epi|3d_grase|3d_spiral)'))
        warning(['Unknown ASL sequence: ' xASL_num2str(x.Q.Sequence)]);
        fprintf('%s\n', 'Trying to fix this');
        x.Q = rmfield(x.Q, 'Sequence');

    end
end

%% Check vendor field
% Assume vendor=GE for 3D spiral sequences 
% (though this is tricky for special cases, Siemens and Philips are both
% working on a 3D spiral)
if isfield(x.Q,'Sequence') && (~isfield(x.Q, 'Vendor') || isempty(x.Q.Vendor))
    if strcmpi(x.Q.Sequence,'3d_spiral')
        warning('x.Q.Vendor missing but 3D spiral sequence detected, assuming vendor GE');
        x.Q.Vendor = 'GE';
    end
end

if ~isfield(x.Q, 'Vendor') || isempty(x.Q.Vendor)
    warning('x.Q.Vendor missing, skipping determining ASL sequence');
elseif  isempty(regexpi(x.Q.Vendor, 'Gold Standard Phantoms|GE|Philips|Siemens'))
    warning('Unknown Vendor specified in x.Q.Vendor');
elseif ~isempty(regexpi(x.Q.Vendor, 'Gold Standard Phantoms'))
    fprintf('%s\n', 'Digital Reference Object ASL-DRO detected');
end


%% Try to work out which ASL sequence we have
% First assume that 2D is 2D EPI, irrespective of Vendor
if ~isfield(x.Q, 'Sequence') && isfield(x.Q, 'readoutDim') 
	if strcmpi(x.Q.readoutDim,'2D')
		x.Q.Sequence = '2D_EPI';
		if bVerbose
            if ~isempty(regexpi(x.Q.Vendor,'Gold Standard Phantoms'))
			    fprintf('%s\n', 'Processing as if this is a 2D EPI sequence');
			    fprintf('%s\n', 'Though the acquisition is not simulated, this will assume acquisition of multi-slice 2D acquisitions');
			    fprintf('%s\n', 'and heavy geometric distortion and minimal smoothness');
		    else
			    fprintf('%s\n', '2D readout detected, assuming 2D EPI');
            end
        end
	elseif isfield(x.Q, 'Vendor') && strcmpi(x.Q.readoutDim,'3D')
		if  ~isempty(regexpi(x.Q.Vendor, 'Philips')) || ~isempty(regexpi(x.Q.Vendor, 'Siemens'))
			x.Q.Sequence = '3D_GRASE'; % assume that 3D Philips or Siemens is 3D GRASE
			if bVerbose; fprintf('%s\n', '3D readout detected with vendor Philips or Siemens, assuming 3D GRASE'); end
		elseif ~isempty(regexpi(x.Q.Vendor, 'GE'))
			x.Q.Sequence = '3D_spiral'; % assume that 3D GE is 3D spiral
			if bVerbose; fprintf('%s\n', '3D readout detected with vendor GE, assuming 3D spiral'); end
		elseif ~isempty(regexpi(x.Q.Vendor, 'Gold Standard Phantoms'))
			x.Q.Sequence = '3D_GRASE'; % assume that this is simulated 3D GRASE by the DRO
			if bVerbose
                fprintf('%s\n', 'Processing as if this is a 3D GRASE sequence');
			    fprintf('%s\n', 'Though the acquisition is not simulated, this will assume acquisition of a single 3D volume');
			    fprintf('%s\n', 'and intermediate amount of geometric distortion and smoothness');
            end
		end
	end
end

%% Warn if we couldn't detect a sequence, default to 3D spiral
% We default to 3D spiral because this sequence has whole-brain coverage
% and negligible geometric distortion, so we don't do any special
% readout-specific image processing for this sequence
if ~isfield(x.Q,'Sequence') || isempty(x.Q.Sequence)
    warning('No x.Q.Sequence defined');
    fprintf('If there are multiple sequence types, this needs to be implemented yet here\n');
    fprintf('Otherwise, please define x.Q.Sequence\n');
    fprintf('Setting x.Q.Sequence=3D spiral to fool the pipeline here\n');
    fprintf('As this sequence doesnt have any susceptibility masks\n');
    fprintf('Note that this disables any masking of susceptibility signal dropout areas\n');
    x.Q.Sequence = '3D_spiral';
    x.Q.readoutDim = '3D';
else
    if bVerbose; fprintf('%s\n', [x.Q.Sequence ' sequence detected']); end
end


end