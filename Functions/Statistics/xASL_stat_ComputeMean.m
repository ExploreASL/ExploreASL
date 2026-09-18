function [resMeanPVC0, resMeanPVC1, resMeanPVC2GM, resMeanPVC2WM, resMedian] = xASL_stat_ComputeMean(imCBF, imMask, nMinSize, bOutput, imGM, imWM)
%xASL_stat_ComputeMean calculates mean or median of CBF in the image across a mask with an optional partial volume correction.
%
% FORMAT:  [resMeanPVC0, resMeanPVC1, resMeanPVC2GM, resMeanPVC2WM, resMedian] = xASL_stat_ComputeMean(imCBF[, imMask, nMinSize, imGM, imWM])
%
% INPUT:
%   imCBF  - input CBF volume (REQUIRED)
%   imMask - mask for the calculation (OPTIONAL, DEFAULT = finite part of imCBF)
%   nMinSize - minimal size of the ROI in voxels, if not big enough, then return NaN
%            - ignore when 0 (OPTIONAL, default = 0)
%   bOutput - vector of length 1 to 5 that specifies if each of the outputs is provided in the order
%             [resMeanPVC0, resMeanPVC1, resMeanPVC2GM, resMeanPVC2WM, resMedian] (OPTIONAL, DEFAULT [1 0 0 0 0])
%   imGM   - GM partial volume map with the same size as imCBF
%            (OPTIONAL, REQUIRED for bPVC==2 and bPVC==1)
%   imWM   - WM partial volume map with the same size as imCBF
%            (OPTIONAL, REQUIRED for bPVC==2)
% OUTPUT:
%   resMeanPVC0 - mean value with PVC0
%   resMeanPVC1 - mean value with PVC1
%   resMeanPVC2GM - mean value with PVC2 in GM
%   resMeanPVC2WM - mean value with PVC2 in WM
%   resMedian   - median value
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It calculates mean or median of CBF over the mask imMask if the mask volume exceeds nMinSize. It calculates either
%              a mean, a median, or a mean after PVC - depending on what outputs. For the PVC options, it needs also imGM and imWM and returns the
%              separate PV-corrected values calculated over the entire ROI. It is a single function call that assigns all value. The reason is that we can share 
%              joint calculations that take most of the calculation time. PVC options are: 0 - don't do partial volume correction, just calculate a mean or median on imMask
%              1 - simple partial volume correction by normalizaton by the GM volume - see Petr et al. 2018
%              2 - partial volume correction using linear regression and imGM, imWM maps according to Asllani et al. 2008
%
% 1. Admin
% 2. Mask calculations
% 3. Calculate the ROI statistics
% 3a. No PVC and simple mean
% 3b. No PVC and median
% 3c. Simple PVC
% 3d. Full PVC on a region
%
% EXAMPLE: resMeanPVC0 = xASL_stat_ComputeMean(imCBF)
%          resMeanPVC0 = xASL_stat_ComputeMean(imCBF,imMask,[])
%          [resMeanPVC0,~,~,~,resMedian] = xASL_stat_ComputeMean(imCBF,[],[],[1,0,0,0,1])
%          [~,~,~,~,resMedian] = xASL_stat_ComputeMean(imCBF,imMask,290,[0 0 0 0 1])
%          [resMeanPVC0, resMeanPVC1] = xASL_stat_ComputeMean(imCBF,imMask,[],[1 1 0 0 0],imGM)
%          [~,~,resMeanPVC2GM,~] = xASL_stat_ComputeMean(imCBF,imMask,[],[0 0 1 0 0],imGM,imWM)
%          [resMeanPVC0, resMeanPVC1, resMeanPVC2GM, resMeanPVC2WM, resMedian] = xASL_stat_ComputeMean(imCBF,[],[],[1 1 1 1 1],imGM,imWM)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: Asllani I, Borogovac A, Brown TR. Regression algorithm correcting for partial volume effects in arterial spin labeling MRI. Magnetic 
%             Resonance in Medicine: An Official Journal of the International Society for Magnetic Resonance in Medicine. 2008 Dec;60(6):1362-71.
% 
%             Petr J, Mutsaerts HJ, De Vita E, Steketee RM, Smits M, Nederveen AJ, Hofheinz F, van den Hoff J, Asllani I. Effects of systematic partial 
%             volume errors on the estimation of gray matter cerebral blood flow with arterial spin labeling MRI. Magnetic Resonance Materials in 
%             Physics, Biology and Medicine. 2018 Dec 1;31(6):725-34.
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% ExploreASL; see permissions and limitations at https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% __________________________________


%% 1. Admin
if nargin<1
	error('imCBF is a required input parameter');
end

if nargin<2 || isempty(imMask)
	imMask = ones(size(imCBF));
end

if nargin<3 || isempty(nMinSize)
	nMinSize = 0;
end

if nargin<4 || isempty(bOutput)
	bOutput = 1;
end

if nargin<5
	imGM = [];
end

if nargin<6 
	imWM = [];
end

% Initialize the output
if nargout>=1 && bOutput(1)
	resMeanPVC0 = NaN;
else
	resMeanPVC0 = [];
end

if nargout>=2 && length(bOutput)>1 && bOutput(2)
	resMeanPVC1 = NaN;
	if isempty(imGM)
		error('Cannot calculate resMeanPVC1 when imGM is not provided');
	end
else
	resMeanPVC1 = [];
end

if nargout>=3 && length(bOutput)>2 && bOutput(3)
	resMeanPVC2GM = NaN;
	if isempty(imWM) || isempty(imGM)
		error('Cannot calculate resMeanPVC2GM when imGM and imWM are not provided');
	end
else
	resMeanPVC2GM = [];
end

if nargout>=4 && length(bOutput)>3 && bOutput(4)
	resMeanPVC2WM = NaN;
	if isempty(imWM) || isempty(imGM)
		error('Cannot calculate resMeanPVC2WM when imGM and imWM are not provided');
	end
else
	resMeanPVC2WM = [];
end

if nargout>=5 && length(bOutput)>4 && bOutput(5)
	resMedian = NaN;
else
	resMedian = [];
end

if ~any(imCBF>0, 'all')
    warning('CBF image is empty, skipping');
    return;
end

if ~isempty(resMeanPVC2GM) || ~isempty(resMeanPVC2WM)
	% If running PVC, then need imGM and imWM of the same size as imCBF
	if ~isequal(size(imCBF),size(imGM)) || ~isequal(size(imCBF),size(imWM))
		warning('When running PVC, the GM and WM maps should have the same size as the CBF image');
	elseif size(imGM, 4)>2
        warning('Invalid GM map size');
    elseif size(imWM, 4)>2
        warning('Invalid WM map size');
	end
end

%% 2. Mask calculations

% Only compute in real data
imMask = imMask>0 & isfinite(imCBF) & (imCBF~=0);

% Constrain calculation to the mask and to finite values
imCBF = imCBF(imMask);

% Limit imGM and imWM to imMask, if provided
if ~isempty(imGM)
	imGM = imGM(imMask); 
end

if ~isempty(imWM)
	imWM = imWM(imMask); 
end

maskSize = sum(imMask, 'all');

if maskSize < nMinSize || maskSize <= 0 
    return;
end

%% 3. Calculate the ROI statistics
sumCBF = sum(imCBF, 1, 'omitnan');

%% 3a. No PVC and simple mean
if ~isempty(resMeanPVC0)
	resMeanPVC0 = sumCBF/maskSize;
end

%% 3b. No PVC and median
if ~isempty(resMedian)
	resMedian = median(imCBF, 1, 'omitnan'); % this is non-parametric
end

%% 3c. Simple PVC
if ~isempty(resMeanPVC1)
	if isempty(imGM)
		error('imGM needs to be provided for bPVC == 1');
	end
	resMeanPVC1 = sumCBF/sum(imGM, 1, 'omitnan');
end
	
%% 3d. Full PVC on a region
if ~isempty(resMeanPVC2GM) || ~isempty(resMeanPVC2WM)
	% although assuming that CBF in CSF = 0, that maps are optimally resampled (cave
	% smoothing of c1T1 & c2T1 to ASL smoothness!) and that TotalVolume-GM-WM = CSF
	% The current absence of modulation here will not change a lot according to Jan Petr

	% Real original Partial Volume Error Correction (PVEC)
	% Normal matrix inverse solves a system of linear equations.
	% If the matrix is not square = more equations than unknowns, then the pseudo-inverse gives solution in the least-square sense - meaning the sum of squares
	% of the error (CBF-CBF*inv(PV)) is minimized.
	%
	% you can see how close you get:
	% gwcbf*gwpv'
	% gwcbf*gwpv' - cbf'

	gwpv                       = imGM;
	gwpv(:,2)                  = imWM;
	gwcbf                      = (imCBF')*pinv(gwpv');
	if ~isempty(resMeanPVC2GM)
		resMeanPVC2GM = gwcbf(1);
	end
	if ~isempty(resMeanPVC2WM)
		resMeanPVC2WM = gwcbf(2);
	end
end


    % % Print histograms to check validity
    % IMPLEMENT THIS LATER ON GROUP LEVEL FOR EACH ROI                            
    % pGMmask         = logical((GMmask) .*TempCurrentMask);
    % pWMmask         = logical((WMmask) .*TempCurrentMask);
    % pCSFmask        = logical((CSFmask).*TempCurrentMask);
    % [Xpgm  Npgm]   = hist(temp(pGMmask));
    % [Xpwm  Npwm]   = hist(temp(pWMmask));
    % [Xcsf  Ncsf]   = hist(temp(pCSFmask));
    % [Xfull Nfull]  = hist(temp(TempCurrentMask));
    % 
    % fig     = figure('Visible','off');
    % plot(Npgm, Xpgm,'r'); %NaNs will be plotted as zeros
    % hold on
    % plot(Npwm, Xpwm,'b');
    % hold on
    % plot(Ncsf, Xcsf,'g');
    % hold on
    % plot(Nfull, Xfull,'k');
    % 
    % xlabel('CBF (mL/100g/min)');
    % ylabel('Norm frequency');
    % title('ROI histogram. Black = full ROI, red = pGM>0.7, blue = pWM>0.7, green = pCSF>0.7, black = full ROI');
    % OutputFile      = fullfile(OutputDir,[x.S.Measurements{iMeas} '_' x.SUBJECTS{iSubject} '_' x.SESSIONS{iSession} '.jpg']);
    % print(gcf,'-djpeg','-r200', OutputFile);
    % close    
end

