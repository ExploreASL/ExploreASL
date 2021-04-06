function sCov = xASL_stat_ComputeSpatialCoV(imCBF,imMask,nMinSize,bPVC,imGM,imWM)
% Calculates spatial coefficient of variation (sCoV) in the image with optional partial volume correction.
%
% FORMAT: sCov = xASL_stat_ComputeSpatialCoV(imCBF)
%         sCov = xASL_stat_ComputeSpatialCoV(imCBF,imMask)
%         sCov = xASL_stat_ComputeSpatialCoV(imCBF,imMask,nMinSize)
%         sCov = xASL_stat_ComputeSpatialCoV(imCBF,imMask,nMinSize,bPVC,imGM,imWM)
%
% INPUT:
%   imCBF  - input CBF volume
%   imMask - mask for the calculation (DEFAULT finite part of imCBF)
%   nMinSize - minimal size of the ROI in voxels (DEFAULT 0), if not big enough, then return NaN
%   bPVC   - perform PV-correction (DEFAULT 0)
%            0 - don't do partial volume correction
%            2 - partial volume correction by using pseudoCov calculated from imGM, imWM
%   imGM   - GM partial volume map with the same size as imCBF, mandatory for bPVC==2
%   imWM   - WM partial volume map with the same size as imCBF, mandatory for bPVC==2
% OUTPUT:
%   sCov   - calculated spatial coefficient of variation
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It calculates the spatial CoV value on finite part of imCBF. Optionally a mask IMMASK is provide, 
%              ROIs of size < NMINSIZE are ignored, and PVC is done for bPVC==2 using imGM and imWM masks and constructing
%              pseudoCoV from pseudoCBF image. For bPVC~=2, imGM and imWM are ignored
%
% EXAMPLE: sCov = xASL_stat_ComputeSpatialCoV(imCBF)
%          sCov = xASL_stat_ComputeSpatialCoV(imCBF,imMask,[])
%          sCov = xASL_stat_ComputeSpatialCoV(imCBF,[],[],0)
%          sCov = xASL_stat_ComputeSpatialCoV(imCBF,imMask,290,0)
%          sCov = xASL_stat_ComputeSpatialCoV(imCBF,imMask,[],2,imGM,imWM)
%          sCov = xASL_stat_ComputeSpatialCoV(imCBF,[],[],2,imGM,imWM)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: Mutsaerts HJ, Petr J, V??clav?? L, van Dalen JW, Robertson AD, Caan MW, Masellis M, Nederveen AJ, Richard E, MacIntosh BJ. The spatial 
%             coefficient of variation in arterial spin labeling cerebral blood flow images. Journal of Cerebral Blood Flow & Metabolism. 
%             2017 Sep;37(9):3184-92.
% __________________________________
% Copyright ?? 2015-2019 ExploreASL
%
% 2017-00-00 HJ


%% Admin
if nargin < 2 || isempty(imMask)
	imMask = ones(size(imCBF));
end
if nargin < 3 || isempty(nMinSize)
	nMinSize = 0;
end
if nargin < 4 || isempty(bPVC)
	bPVC = 0;
end
if nargin < 5
	imGM = [];
end
if nargin < 6
	imWM = [];
end

if bPVC == 2
	% If running PVC, then need imGM and imWM of the same size as imCBF
	if ~isequal(size(imCBF),size(imGM)) || ~isequal(size(imCBF),size(imWM))
		warning('xASL_stat_ComputeSpatialCoV: When running PVC, need imGM and imWM of the same size as imCBF');
	end
end

%% Only compute in real data
imMask = (imMask>0) & isfinite(imCBF);

imMask = imMask & (imCBF~=0); % Exclude zero values as well

%% Constrain calculation to the mask and to finite values
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

%% Computation

sCov = std(imCBF(:)) / mean(imCBF(:));

if bPVC==2
    % Partial volume correction by normalizing by expected variance because of the structural data
    PseudoCoV = imGM + 0.3.*imWM;
    PseudoCoV = std(PseudoCoV(:)) / mean(PseudoCoV(:));
    sCov = sCov./PseudoCoV;
end
    

end

