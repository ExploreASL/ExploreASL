function [CBF_GM, CBF_WM] = xASL_stat_ComputeMean(imCBF, imMask, nMinSize, bPVC, bParametric, imGM, imWM)
%xASL_stat_ComputeMean calculates mean or median of CBF in the image across a mask with an optional partial volume correction.
%
% FORMAT:  [CBF_GM CBF_WM] = xASL_stat_ComputeMean(imCBF[, imMask, nMinSize, bPVC, bParametric, imGM, imWM])
%
% INPUT:
%   imCBF  - input CBF volume (REQUIRED)
%   imMask - mask for the calculation (OPTIONAL, DEFAULT = finite part of imCBF)
%   nMinSize - minimal size of the ROI in voxels, if not big enough, then return NaN
%            - ignore when 0 (OPTIONAL, default = 0)
%   bPVC   - perform PV-correction (OPTIONAL, DEFAULT = 0)
%            0 - don't do partial volume correction, just calculate a mean or median on imMask
%            1 - simple partial volume correction by normalizaton by the GM volume - see Petr et al. 2018
%            2 - partial volume correction using linear regression and imGM, imWM maps according to Asllani et al. 2008
%   bParametric - performs parametric statistics (1 mean) or non-parametric when turned off (0 median) (OPTIONAL, DEFAULT 1) 
%   imGM   - GM partial volume map with the same size as imCBF
%            (OPTIONAL, REQUIRED for bPVC==2 and bPVC==1)
%   imWM   - WM partial volume map with the same size as imCBF
%            (OPTIONAL, REQUIRED for bPVC==2)
% OUTPUT:
%   CBF_GM - calculated CBF for options bPVC==0:2
%   CBF_WM - calculated WM CBF for bPVC==2
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It calculates mean or median of CBF over the mask imMask if the mask volume exceeds nMinSize. It calculates either
%              a mean, a median, or a mean after PVC, depending on the settings of bPVC. For the PVC options, it needs also imGM and imWM and returns the
%              separate PV-corrected values calculated over the entire ROI.
%
% 1. Admin
% 2. Mask calculations
% 3. Calculate the ROI statistics
% 3a. No PVC and simple mean
% 3b. No PVC and median
% 3c. Simple PVC
% 3d. Full PVC on a region
%
% EXAMPLE: CBF_GM = xASL_stat_ComputeMean(imCBF)
%          CBF_GM = xASL_stat_ComputeMean(imCBF,imMask,[])
%          CBF_GM = xASL_stat_ComputeMean(imCBF,[],[],0,0)
%          CBF_GM = xASL_stat_ComputeMean(imCBF,imMask,290,0,1)
%          CBF_GM = xASL_stat_ComputeMean(imCBF,imMask,[],1,1,imGM)
%          CBF_GM = xASL_stat_ComputeMean(imCBF,imMask,[],2,1,imGM,imWM)
% [CBF_GM CBF_WM] = xASL_stat_ComputeMean(imCBF,[],[],2,1,imGM,imWM)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: Asllani I, Borogovac A, Brown TR. Regression algorithm correcting for partial volume effects in arterial spin labeling MRI. Magnetic 
%             Resonance in Medicine: An Official Journal of the International Society for Magnetic Resonance in Medicine. 2008 Dec;60(6):1362-71.
% 
%             Petr J, Mutsaerts HJ, De Vita E, Steketee RM, Smits M, Nederveen AJ, Hofheinz F, van den Hoff J, Asllani I. Effects of systematic partial 
%             volume errors on the estimation of gray matter cerebral blood flow with arterial spin labeling MRI. Magnetic Resonance Materials in 
%             Physics, Biology and Medicine. 2018 Dec 1;31(6):725-34.
% __________________________________
% Copyright (C) 2015-2021 ExploreASL

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

if nargin<4 || isempty(bPVC)
	bPVC = 0;
end

if nargin<5 || isempty(bParametric)
	bParametric = 1;
end

if bPVC>0 && bParametric == 0
	error('Non-parametric PVC is not implemented yet')
end

if nargin<6 
	imGM = [];
end

if nargin<7 
	imWM = [];
end

% Initialize the output
CBF_GM = NaN;
CBF_WM = NaN;

if nargout>1 && bPVC~=2
	warning('Calculates CBF_WM only for PVC==2 option');
end

if (nargout>1) && isempty(imWM)
	warning('Cannot calculate CBF_WM when imWM is not provided');
end

if sum(isfinite(imCBF(:)))==0
    warning('CBF image is empty, skipping');
    return;
end

if bPVC == 2
	% If running PVC, then need imGM and imWM of the same size as imCBF
	if ~isequal(size(imCBF),size(imGM)) || ~isequal(size(imCBF),size(imWM))
		warning('When running PVC, need imGM and imWM of the same size as imCBF');
	end
end

%% 2. Mask calculations

% Only compute in real data
imMask = imMask>0 & isfinite(imCBF);

imMask = imMask & (imCBF~=0); % Exclude zero values as well

% Constrain calculation to the mask and to finite values
imCBF = imCBF(imMask);

% Limit imGM and imWM to imMask, if provided
if ~isempty(imGM)
	imGM = imGM(imMask); 
end

if ~isempty(imWM)
	imWM = imWM(imMask); 
end

if sum(imMask(:))<nMinSize 
    CBF_GM = NaN;
    CBF_WM = NaN;
    return;
end

%% 3. Calculate the ROI statistics
switch (bPVC)
	case 0
		if bParametric
			%% 3a. No PVC and simple mean
			CBF_GM = xASL_stat_MeanNan(imCBF);
		
		else
			%% 3b. No PVC and median
			CBF_GM = xASL_stat_MedianNan(imCBF); % this is non-parametric
		end

	case 1
		%% 3c. Simple PVC
		if isempty(imGM)
			error('imGM needs to be provided for bPVC == 1');
		end
		CBF_GM = xASL_stat_SumNan(imCBF)/xASL_stat_SumNan(imGM);
	
	case 2
		%% 3d. Full PVC on a region
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
		CBF_GM                     = gwcbf(1);
		CBF_WM                     = gwcbf(2);
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

