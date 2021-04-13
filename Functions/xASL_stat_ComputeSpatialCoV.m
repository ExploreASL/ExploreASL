function sCov = xASL_stat_ComputeSpatialCoV(imCBF, imMask, nMinSize, bPVC, bParametric, imGM, imWM)
%xASL_stat_ComputeSpatialCoV calculates spatial coefficient of variation (sCoV) in the image with optional partial volume correction.
%
% FORMAT: sCov = xASL_stat_ComputeSpatialCoV(imCBF[, imMask, nMinSize, bPVC, bParametric, imGM, imWM])
%
% INPUT:
%   imCBF       - input CBF volume (REQUIRED)
%   imMask      - mask for the calculation (OPTIONAL, DEFAULT finite part of imCBF)
%   nMinSize    - minimal size of the ROI in voxels, if not big enough, then return NaN (OPTIONAL, DEFAULT 0)
%   bPVC        - perform PV-correction (OPTIONAL, DEFAULT 0)
%                 0 - do not do partial volume correction
%                 2 - partial volume correction by using pseudoCov calculated from imGM, imWM
%   bParametric - performs parametric statistics (1 mean) or non-parametric when turned off (0 median) (OPTIONAL, DEFAULT 1) 
%   imGM        - GM partial volume map with the same size as imCBF (OPTIONAL, but REQUIRED for bPVC==2)
%   imWM        - WM partial volume map with the same size as imCBF (OPTIONAL, but REQUIRED for bPVC==2)
%
% OUTPUT:
%   sCov   - calculated spatial coefficient of variation
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It calculates the spatial CoV value on finite part of imCBF. Optionally a mask IMMASK is provide, 
%              ROIs of size < NMINSIZE are ignored, and PVC is done for bPVC==2 using imGM and imWM masks and constructing
%              pseudoCoV from pseudoCBF image. For bPVC~=2, imGM and imWM are ignored
%
% 1. Admin
% 2. Create masks
% 3. sCoV computation
%
% EXAMPLE: sCov = xASL_stat_ComputeSpatialCoV(imCBF)
%          sCov = xASL_stat_ComputeSpatialCoV(imCBF, imMask, [])
%          sCov = xASL_stat_ComputeSpatialCoV(imCBF, [], [], 0)
%          sCov = xASL_stat_ComputeSpatialCoV(imCBF, imMask, 290, 0)
%          sCov = xASL_stat_ComputeSpatialCoV(imCBF, imMask, [], 2, imGM, imWM)
%          sCov = xASL_stat_ComputeSpatialCoV(imCBF, [], [], 2, imGM, imWM)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: Mutsaerts HJ, Petr J, Vaclavu L, van Dalen JW, Robertson AD, Caan MW, Masellis M, Nederveen AJ, Richard E, MacIntosh BJ. The spatial 
%             coefficient of variation in arterial spin labeling cerebral blood flow images. Journal of Cerebral Blood Flow & Metabolism. 
%             2017 Sep;37(9):3184-92.
% __________________________________
% Copyright (C) 2015-2021 ExploreASL


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

if nargin < 4 || isempty(bPVC)
	bPVC = 0;
end

if nargin<5 || isempty(bParametric)
	bParametric = 1;
end

if (bParametric == 0) || (bPVC == 1)
	error('This option is not implemented yet')
end

if nargin < 6
	imGM = [];
end
if nargin < 7
	imWM = [];
end

if bPVC == 2
	% If running PVC, then need imGM and imWM of the same size as imCBF
	if ~isequal(size(imCBF),size(imGM)) || ~isequal(size(imCBF),size(imWM))
		warning('When running PVC, need imGM and imWM of the same size as imCBF');
	end
end

%% 2. Create masks
% Only compute in real data
imMask = (imMask>0) & isfinite(imCBF);

imMask = imMask & (imCBF~=0); % Exclude zero values as well

% Constrain calculation to the mask and to finite values
imCBF = imCBF(imMask);

if ~isempty(imGM)
	imGM = imGM(imMask); 
end
if ~isempty(imWM)
	imWM = imWM(imMask); 
end
    
if sum(imMask(:))<nMinSize 
    sCov  = NaN;
    return;
end

%% 3. sCoV computation

sCov = std(imCBF(:)) / mean(imCBF(:));

if bPVC==2
    % Partial volume correction by normalizing by expected variance because of the structural data
    PseudoCoV = imGM + 0.3.*imWM;
    PseudoCoV = std(PseudoCoV(:)) / mean(PseudoCoV(:));
    sCov = sCov./PseudoCoV;
end
    

end

