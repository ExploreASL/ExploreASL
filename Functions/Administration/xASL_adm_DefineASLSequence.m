function [xQ] = xASL_adm_DefineASLSequence(xQ, bVerbose)
%xASL_adm_DefineASLSequence Obtain ASL sequence type for sequence-specific
%image processing
%
% FORMAT: [xQ] = xASL_adm_DefineASLSequence(xQ, bVerbose)
%
% INPUT:
%   xQ                   - x.Q structure containing all input parameters (REQUIRED)
%   xQ.MRAcquisitionType - dimensionality of readout (2D or 3D) (OPTIONAL)
%   xQ.Vendor            - Either 'GE', 'Philips', 'Siemens' (OPTIONAL)
%   bVerbose             - verbose output (OPTIONAL, DEFAULT=true)
%
% OUTPUT:
%   xQ                   - x structure containing all output parameters
%   xQ.PulseSequenceType - sequence type (readout)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL function tries to check what ASL sequence is
% being processed, if this was not already defined in xQ.PulseSequenceType.
% It does so by checking known combinations of readout dimensionality
% (xQ.MRAcquisitionType) and Vendor, knowing the product sequences of the Vendors.
%
% EXAMPLE: xQ = xASL_adm_DefineASLSequence(xQ);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2024 ExploreASL


%% Check quantification fields MRAcquisitionType & xQ.PulseSequenceType
if nargin<2 || isempty(bVerbose)
    bVerbose = true;
end

if ~isfield(xQ, 'MRAcquisitionType') || isempty(xQ.MRAcquisitionType)
    warning('xQ.MRAcquisitionType parameter missing');
end

% Check for illegal sequence definitions
if isfield(xQ, 'PulseSequenceType')
    if isempty(regexpi(xQ.PulseSequenceType, '(epi|grase|spiral)', 'once'))
        warning(['Unknown ASL sequence: ' xASL_num2str(xQ.PulseSequenceType)]);
        fprintf('%s\n', 'Trying to fix this');
        xQ = rmfield(xQ, 'PulseSequenceType');

    end
end

%% Check vendor field
% Assume vendor=GE for 3D spiral sequences 
% (though this is tricky for special cases, Siemens and Philips are both working on a 3D spiral)
if isfield(xQ, 'PulseSequenceType') && (~isfield(xQ, 'Vendor') || isempty(xQ.Vendor))
    if strcmpi(xQ.PulseSequenceType,'spiral')
        warning('xQ.Vendor missing but 3D spiral sequence detected, assuming vendor GE');
        xQ.Vendor = 'GE';
    end
end

if ~isfield(xQ, 'Vendor') || isempty(xQ.Vendor)
    warning('xQ.Vendor missing, skipping determining ASL sequence');
elseif  isempty(regexpi(xQ.Vendor, 'Gold Standard Phantoms|GE|Philips|Siemens', 'once'))
    warning('Unknown Vendor specified in xQ.Vendor');
elseif ~isempty(regexpi(xQ.Vendor, 'Gold Standard Phantoms', 'once'))
    fprintf('%s\n', 'Digital Reference Object ASL-DRO detected');
end


%% Try to work out which ASL sequence we have
% First assume that 2D is 2D EPI, irrespective of Vendor
if ~isfield(xQ, 'PulseSequenceType') && isfield(xQ, 'MRAcquisitionType') 
	if strcmpi(xQ.MRAcquisitionType, '2D')
		xQ.PulseSequenceType = 'EPI';
        if bVerbose
            if ~isempty(regexpi(xQ.Vendor, 'Gold Standard Phantoms', 'once'))
			    fprintf('%s\n', 'Processing as if this is a 2D EPI sequence');
			    fprintf('%s\n', 'Though the acquisition is not simulated, this will assume acquisition of multi-slice 2D acquisitions');
			    fprintf('%s\n', 'and heavy geometric distortion and minimal smoothness');
		    else
			    fprintf('%s\n', '2D readout detected, assuming 2D EPI');
            end
        end
	elseif isfield(xQ, 'Vendor') && strcmpi(xQ.MRAcquisitionType, '3D')
		if  ~isempty(regexpi(xQ.Vendor, 'Philips', 'once')) || ~isempty(regexpi(xQ.Vendor, 'Siemens', 'once'))
			xQ.PulseSequenceType = 'GRASE'; % assume that 3D Philips or Siemens is 3D GRASE
			if bVerbose; fprintf('%s\n', '3D readout detected with vendor Philips or Siemens, assuming 3D GRASE'); end
		elseif ~isempty(regexpi(xQ.Vendor, 'GE', 'once'))
			xQ.PulseSequenceType = 'spiral'; % assume that 3D GE is 3D spiral
			if bVerbose; fprintf('%s\n', '3D readout detected with vendor GE, assuming 3D spiral'); end
		elseif ~isempty(regexpi(xQ.Vendor, 'Gold Standard Phantoms', 'once'))
			xQ.PulseSequenceType = 'GRASE'; % assume that this is simulated 3D GRASE by the DRO
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
if ~isfield(xQ, 'PulseSequenceType') || isempty(xQ.PulseSequenceType)
    warning('No xQ.PulseSequenceType defined');
    fprintf('If there are multiple sequence types, this needs to be implemented yet here\n');
    fprintf('Otherwise, please define xQ.PulseSequenceType\n');
    fprintf('Setting xQ.PulseSequenceType=3D spiral to fool the pipeline here\n');
    fprintf('As this sequence doesnt have any susceptibility masks\n');
    fprintf('Note that this disables any masking of susceptibility signal dropout areas\n');
    xQ.PulseSequenceType = 'spiral';
    xQ.MRAcquisitionType = '3D';
else
    if bVerbose; fprintf('%s\n', [xQ.PulseSequenceType ' sequence detected']); end
end


end
