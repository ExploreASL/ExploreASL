function [xQ] = xASL_adm_DefineASLReadout(xQ, bVerbose)
%xASL_adm_DefineASLReadout Obtain ASL readout type for readout-specific image processing
%
% FORMAT: [xQ] = xASL_adm_DefineASLReadout(xQ, bVerbose)
%
% INPUT:
%   xQ                   - x.Q structure containing all input parameters (REQUIRED)
%   xQ.MRAcquisitionType - dimensionality of readout (2D or 3D) (OPTIONAL)
%   xQ.Vendor            - Either 'GE', 'Philips', 'Siemens' (OPTIONAL)
%   bVerbose             - verbose output (OPTIONAL, DEFAULT=true)
%
% OUTPUT:
%   xQ                   - x structure containing all output parameters
%   xQ.PulseSequenceType - pulse sequence readout type
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This ExploreASL function tries to check what ASL readout is
% being processed, if this was not already defined in xQ.PulseSequenceType.
% It does so by checking known combinations of readout dimensionality
% (xQ.MRAcquisitionType) and Vendor, knowing the product sequences of the Vendors.
%
% EXAMPLE: xQ = xASL_adm_DefineASLReadout(xQ);
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________



%% Check quantification fields MRAcquisitionType & xQ.PulseSequenceType
if nargin<2 || isempty(bVerbose)
    bVerbose = true;
end

if ~isfield(xQ, 'MRAcquisitionType') || isempty(xQ.MRAcquisitionType)
    warning('xQ.MRAcquisitionType parameter missing');
end

% Check and optionally fix illegal pulse sequence definitions
% Tier 1: explicit pattern table for vendor-specific names
% Tier 2: fuzzy match (Damerau-Levenshtein) for typo tolerance
if isfield(xQ, 'PulseSequenceType') && ~isempty(xQ.PulseSequenceType)
    pstRaw   = xQ.PulseSequenceType;
    pstLower = lower(pstRaw);
    pstCorrected = '';

    %% Tier 1: explicit pattern table for vendor-specific names
    if ~isempty(regexpi(pstLower, 'epi|ep2d|epfid|epse|se\-?epi|ffe\-?epi', 'once'))
        pstCorrected = 'EPI';
    elseif ~isempty(regexpi(pstLower, 'grase|tgse', 'once'))
        pstCorrected = 'GRASE';
    elseif ~isempty(regexpi(pstLower, 'spiral', 'once'))
        pstCorrected = 'spiral';
    end

    %% Tier 2: fuzzy match against normalized keywords (typo fallback)
    if isempty(pstCorrected)
        pstStripped = regexprep(pstLower, '^(\d+d[_\-]?)', ''); % strip "3D_" etc.
        if isempty(pstStripped)
            pstStripped = pstLower;
        end
        keywords = {'epi', 'grase', 'spiral'};
        bestDist = inf;
        bestKey  = '';
        for k = 1:length(keywords)
            d = xASL_stat_EditDistance(pstStripped, keywords{k});
            if d < bestDist
                bestDist = d;
                bestKey  = keywords{k};
            end
        end
        % Keyword length-aware threshold: 1 for short (<=4), 2 for longer
        if bestDist <= 1 || (length(bestKey) > 4 && bestDist <= 2)
            pstCorrected = bestKey;
        end
    end

    %% Apply result
    if ~isempty(pstCorrected)
        if ~strcmpi(pstRaw, pstCorrected)
            warning(['PulseSequenceType "' pstRaw '" corrected to "' pstCorrected '"']);
        end
        xQ.PulseSequenceType = pstCorrected;
    else
        warning(['Unknown ASL readout PulseSequenceType: ' xASL_num2str(pstRaw) '; falling back to inference']);
        xQ = rmfield(xQ, 'PulseSequenceType');
    end
end

%% Check vendor field
% Use BIDS Manufacturer field if present to set Vendor (with fuzzy fallback for typos).
% Only fall back to the spiral->GE assumption when Vendor and Manufacturer are both absent.
knownVendors = {'GE', 'Philips', 'Siemens', 'Gold Standard Phantoms'};

if ~isfield(xQ, 'Vendor') || isempty(xQ.Vendor)
    if isfield(xQ, 'Manufacturer') && ~isempty(xQ.Manufacturer)
        % Manufacturer is BIDS: use it to set Vendor
        mfrLower = xQ.Manufacturer;
        xQ.Vendor = '';
        % Substring match first (handles "Siemens Healthineers" etc.)
        for v = 1:length(knownVendors)
            if ~isempty(regexpi(mfrLower, knownVendors{v}, 'once'))
                xQ.Vendor = knownVendors{v};
                break
            end
        end
        % Fuzzy fallback (handles typos like "Seimens", "PHILPS")
        if isempty(xQ.Vendor)
            bestDist = inf;
            bestV    = '';
            for v = 1:length(knownVendors)
                d = xASL_stat_EditDistance(lower(mfrLower), lower(knownVendors{v}));
                if d < bestDist
                    bestDist = d;
                    bestV    = knownVendors{v};
                end
            end
            % Length-aware threshold: short vendor strings (GE) need tighter threshold
            vLen = length(bestV);
            if (vLen <= 4 && bestDist <= 1) || (vLen > 4 && bestDist <= 2)
                warning(['Manufacturer "' xQ.Manufacturer '" fuzzy-corrected to vendor "' ...
                         bestV '" (edit distance ' xASL_num2str(bestDist) ')']);
                xQ.Vendor = bestV;
            else
                warning(['Unknown Manufacturer "' xQ.Manufacturer '" cannot set Vendor']);
            end
        end
    elseif isfield(xQ, 'PulseSequenceType') && strcmpi(xQ.PulseSequenceType, 'spiral')
        % Last-resort heuristic: spiral implies GE on production systems
        warning('Vendor/Manufacturer missing but spiral readout detected, assuming vendor GE');
        xQ.Vendor = 'GE';
    else
        warning('xQ.Vendor missing and no Manufacturer to infer from');
    end
end

% Warn if Vendor is set but unrecognized
if isfield(xQ, 'Vendor') && ~isempty(xQ.Vendor)
    if isempty(regexpi(xQ.Vendor, 'Gold Standard Phantoms|GE|Philips|Siemens', 'once'))
        warning('Unknown Vendor specified in xQ.Vendor');
    elseif ~isempty(regexpi(xQ.Vendor, 'Gold Standard Phantoms', 'once'))
        fprintf('%s\n', 'Digital Reference Object ASL-DRO detected');
    end
end


%% Try to work out which ASL readout we have
% First assume that 2D is 2D EPI, irrespective of Vendor
if ~isfield(xQ, 'PulseSequenceType') && isfield(xQ, 'MRAcquisitionType') 
	if strcmpi(xQ.MRAcquisitionType, '2D')
		xQ.PulseSequenceType = 'EPI';
        if bVerbose
            if ~isempty(regexpi(xQ.Vendor, 'Gold Standard Phantoms', 'once'))
			    fprintf('%s\n', 'Processing as if this is a 2D EPI readout');
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
                fprintf('%s\n', 'Processing as if this is a 3D GRASE readout');
			    fprintf('%s\n', 'Though the acquisition is not simulated, this will assume acquisition of a single 3D volume');
			    fprintf('%s\n', 'and intermediate amount of geometric distortion and smoothness');
            end
		end
	end
end

%% Consistency check between PulseSequenceType and MRAcquisitionType
if isfield(xQ, 'PulseSequenceType') && ~isempty(xQ.PulseSequenceType) ...
        && isfield(xQ, 'MRAcquisitionType') && ~isempty(xQ.MRAcquisitionType)
    pst = lower(xQ.PulseSequenceType);
    mra = upper(xQ.MRAcquisitionType);
    bDRO = isfield(xQ, 'Vendor') && ~isempty(regexpi(xQ.Vendor, 'Gold Standard Phantoms', 'once'));
    if strcmpi(mra, '2D') && strcmpi(pst, 'spiral')
        warning('Implausible combination: 2D + spiral (no production 2D spiral exists). Check sidecar.');
    elseif strcmpi(mra, '2D') && strcmpi(pst, 'grase')
        warning('Implausible combination: 2D + GRASE (GRASE is a 3D readout by design). Check sidecar.');
    elseif strcmpi(mra, '3D') && strcmpi(pst, 'epi') && ~bDRO
        warning('Implausible combination: 3D + EPI (production 3D EPI is rare). Check sidecar.');
    end
end

%% Warn if we couldn't detect a readout
if ~isfield(xQ, 'PulseSequenceType') || isempty(xQ.PulseSequenceType) || ~isfield(xQ, 'MRAcquisitionType') || isempty(xQ.MRAcquisitionType)
    error('We cannot detect the PulseSequenceType and/or MRAcquisitionType')
else
    if bVerbose; fprintf('%s\n', [xQ.PulseSequenceType ' readout detected']); end
end

end
