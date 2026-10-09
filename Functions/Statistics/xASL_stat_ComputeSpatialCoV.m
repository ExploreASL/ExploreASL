function [sCoVPVC0, sCoVPVC2] = xASL_stat_ComputeSpatialCoV(imCBF, imMask, nMinSize, bOutput, imPVprimary, imPVsecondary)
%xASL_stat_ComputeSpatialCoV calculates spatial coefficient of variation (sCoV) in the image with optional partial volume correction.
%
% FORMAT: [sCoVPVC0, sCoVPVC2] = xASL_stat_ComputeSpatialCoV(imCBF[, imMask, nMinSize, bOutput, imPVprimary, imPVsecondary])
%
% INPUT:
%   imCBF       - input CBF volume (REQUIRED)
%   imMask      - mask for the calculation (OPTIONAL, DEFAULT finite part of imCBF)
%   nMinSize    - minimal size of the ROI in voxels, if not big enough, then return NaN (OPTIONAL, DEFAULT 0)
%   bOutput     - vector of length 1 or 2 that specifies which output is requested
%                 [sCoVPVC0, sCoVPVC2] (OPTIONAL, DEFAULT [1 0 ])
%   imPVprimary        - Primary partial volume map with the same size as imCBF (OPTIONAL, but REQUIRED for bPVC==2)
%   imPVsecondary        - Secondary partial volume map with the same size as imCBF (OPTIONAL, but REQUIRED for bPVC==2)
%
% OUTPUT:
%   sCoVPVC0   - calculated spatial coefficient of variation
%   sCoVPVC2   - calculated spatial coefficient of variation with partial volume correction
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It calculates the spatial CoV value on finite part of imCBF. Optionally a mask IMMASK is provide, 
%              ROIs of size < NMINSIZE are ignored, and PVC is done for for sCoVPVC2 using imPVprimary and imPVsecondary masks and constructing
%              pseudoCoV from pseudoCBF image. The values are calculated over IMMASK. imPVprimary and imPVsecondary are only used for PVC2
%
% 1. Admin
% 2. Create masks
% 3. sCoV computation
%
% EXAMPLE: sCoVPVC0 = xASL_stat_ComputeSpatialCoV(imCBF)
%          sCoVPVC0 = xASL_stat_ComputeSpatialCoV(imCBF, imMask, [])
%          sCoVPVC0 = xASL_stat_ComputeSpatialCoV(imCBF, [], [], [])
%          sCoVPVC0 = xASL_stat_ComputeSpatialCoV(imCBF, imMask, 290, [1 0])
%          [sCoVPVC0, sCoVPVC2] = xASL_stat_ComputeSpatialCoV(imCBF, imMask, [], [1 1], imPVprimary, imPVsecondary)
%          [~, sCoVPVC2] = xASL_stat_ComputeSpatialCoV(imCBF, [], [], [0 1], imPVprimary, imPVsecondary)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: Mutsaerts HJ, Petr J, Vaclavu L, van Dalen JW, Robertson AD, Caan MW, Masellis M, Nederveen AJ, Richard E, MacIntosh BJ. The spatial 
%             coefficient of variation in arterial spin labeling cerebral blood flow images. Journal of Cerebral Blood Flow & Metabolism. 
%             2017 Sep;37(9):3184-92.
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________



%% 1. Admin
if nargin < 1
	error('imCBF is a required parameter');
end

if nargin < 2 || isempty(imMask)
	imMask = ones(size(imCBF));
end

if nargin < 3 || isempty(nMinSize)
	nMinSize = 0;
end

if nargin < 4 || isempty(bOutput)
	bOutput = [1 0];
else
	bOutput((end+1):2) = 0;
end

if nargin < 5
	imPVprimary = [];
end
if nargin < 6
	imPVsecondary = [];
end

% Initialize the output
if nargout<1
	bOutput(1) = 0;
end

if nargout>=2 && bOutput(2)
	sCoVPVC2 = NaN;
	% If running PVC, then need imPVprimary and imPVsecondary of the same size as imCBF
	if ~isequal(size(imCBF),size(imPVprimary)) || ~isequal(size(imCBF),size(imPVsecondary))
		warning('When running PVC, need imPVprimary and imPVsecondary of the same size as imCBF');
	end
else
	bOutput(2) = 0;
end

if ~any(imCBF>0, 'all')
    warning('CBF image is empty, skipping');
    return;
end


%% 2. Create masks
% Only compute in real data
imMask = (imMask>0) & isfinite(imCBF) & (imCBF~=0);

% Constrain calculation to the mask and to finite values
imCBF = imCBF(imMask);

if ~isempty(imPVprimary)
	imPVprimary = imPVprimary(imMask); 
end
if ~isempty(imPVsecondary)
	imPVsecondary = imPVsecondary(imMask); 
end
    
maskSize = sum(imMask, 'all');

if maskSize < nMinSize || maskSize <= 0 
	sCoVPVC0  = NaN;
    sCoVPVC2  = NaN;
    return;
end

%% 3. sCoV computation

sCoVPVC0 = std(imCBF, 'omitnan') / mean(imCBF, 'omitnan');

if bOutput(2)
    % Partial volume correction by normalizing by expected variance because of the structural data
    PseudoCoV = imPVprimary + 0.3.*imPVsecondary;
    PseudoCoV = std(PseudoCoV) / mean(PseudoCoV);
    sCoVPVC2 = sCoVPVC0./PseudoCoV;
end

end

