function [sCoVPVC0, sCoVPVC2] = xASL_stat_ComputeSpatialCoV(imCBF, imMask, nMinSize, bOutput, imGM, imWM)
%xASL_stat_ComputeSpatialCoV calculates spatial coefficient of variation (sCoV) in the image with optional partial volume correction.
%
% FORMAT: [sCoVPVC0, sCoVPVC2] = xASL_stat_ComputeSpatialCoV(imCBF[, imMask, nMinSize, bOutput, imGM, imWM])
%
% INPUT:
%   imCBF       - input CBF volume (REQUIRED)
%   imMask      - mask for the calculation (OPTIONAL, DEFAULT finite part of imCBF)
%   nMinSize    - minimal size of the ROI in voxels, if not big enough, then return NaN (OPTIONAL, DEFAULT 0)
%   bOutput     - vector of length 1 to 2 that specifies if each of the outputs is provided in the order
%                 [sCoVPVC0, sCoVPVC2] (OPTIONAL, DEFAULT [1 0 ])
%   imGM        - GM partial volume map with the same size as imCBF (OPTIONAL, but REQUIRED for bPVC==2)
%   imWM        - WM partial volume map with the same size as imCBF (OPTIONAL, but REQUIRED for bPVC==2)
%
% OUTPUT:
%   sCoVPVC0   - calculated spatial coefficient of variation
%   sCoVPVC2   - calculated spatial coefficient of variation with partial volume correction
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It calculates the spatial CoV value on finite part of imCBF. Optionally a mask IMMASK is provide, 
%              ROIs of size < NMINSIZE are ignored, and PVC is done for for sCoVPVC2 using imGM and imWM masks and constructing
%              pseudoCoV from pseudoCBF image. The values are calculated over IMMASK. imGM and imWM are only used for PVC2
%
% 1. Admin
% 2. Create masks
% 3. sCoV computation
%
% EXAMPLE: sCoVPVC0 = xASL_stat_ComputeSpatialCoV(imCBF)
%          sCoVPVC0 = xASL_stat_ComputeSpatialCoV(imCBF, imMask, [])
%          sCoVPVC0 = xASL_stat_ComputeSpatialCoV(imCBF, [], [], [])
%          sCoVPVC0 = xASL_stat_ComputeSpatialCoV(imCBF, imMask, 290, [1 0])
%          [sCoVPVC0, sCoVPVC2] = xASL_stat_ComputeSpatialCoV(imCBF, imMask, [], [1 1], imGM, imWM)
%          [~, sCoVPVC2] = xASL_stat_ComputeSpatialCoV(imCBF, [], [], [0 1], imGM, imWM)
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
end

if nargin < 5
	imGM = [];
end
if nargin < 6
	imWM = [];
end

% Initialize the output
if nargout>=1 && bOutput(1)
	sCoVPVC0 = NaN;
else
	sCoVPVC0 = [];
end

if nargout>=2 && length(bOutput)>1 && bOutput(2)
	sCoVPVC2 = NaN;
	% If running PVC, then need imGM and imWM of the same size as imCBF
	if ~isequal(size(imCBF),size(imGM)) || ~isequal(size(imCBF),size(imWM))
		warning('When running PVC, need imGM and imWM of the same size as imCBF');
	end
else
	sCoVPVC2 = [];
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

if ~isempty(imGM)
	imGM = imGM(imMask); 
end
if ~isempty(imWM)
	imWM = imWM(imMask); 
end
    
maskSize = sum(imMask, 'all');

if maskSize < nMinSize || maskSize <= 0 
	sCoVPVC0  = NaN;
    sCoVPVC2  = NaN;
    return;
end

%% 3. sCoV computation

sCoVPVC0 = std(imCBF, 'omitnan') / mean(imCBF, 'omitnan');

if ~isempty(sCoVPVC2)
    % Partial volume correction by normalizing by expected variance because of the structural data
    PseudoCoV = imGM + 0.3.*imWM;
    PseudoCoV = std(PseudoCoV) / mean(PseudoCoV);
    sCoVPVC2 = sCoVPVC0./PseudoCoV;
end

end

