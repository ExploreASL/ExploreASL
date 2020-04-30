function [CBF_GM, CBF_WM] = xASL_stat_ComputeMean(imCBF,imMask,nMinSize,bPVC,imGM,imWM)
%Calculates mean CBF in the image with optional partial volume correction.
%
% FORMAT:  [CBF_GM CBF_WM] = xASL_stat_ComputeMean(imCBF [,imMask,nMinSize,bPVC,imGM,imWM])
%
% INPUT:
%   imCBF  - input CBF volume (REQUIRED)
%   imMask - mask for the calculation (OPTIONAL, default = finite part of imCBF)
%   nMinSize - minimal size of the ROI in voxels, if not big enough, then return NaN
%            - ignore when 0 (OPTIONAL, default = 0)
%   bPVC   - perform PV-correction (OPTIONAL, default = 0)
%            0 - don't do partial volume correction, just calculate a median on imMask (+imGM and imWM)
%            1 - simple partial volume correction by dividing by the GM mask
%            2 - partial volume correction using linear regression and imGM, imWM masks
%   imGM   - GM partial volume map with the same size as imCBF
%            (OPTIONAL, mandatory for bPVC==2 and bPVC==1)
%   imWM   - WM partial volume map with the same size as imCBF
%            (OPTIONAL, mandatory for bPVC==2, if provided for bPVC=0 gives also CBF_WM)
% OUTPUT:
%   CBF_GM - calculated CBF for options bPVC==0:2
%   CBF_WM - calculated WM CBF for bPVC==2, and for ==0 if imWM provided
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It behaves in a similar way as VAR.
%
% EXAMPLE: CBF_GM = xASL_stat_ComputeMean(imCBF)
%          CBF_GM = xASL_stat_ComputeMean(imCBF,imMask,[])
%          CBF_GM = xASL_stat_ComputeMean(imCBF,[],[],0)
%          CBF_GM = xASL_stat_ComputeMean(imCBF,imMask,290,0)
%          CBF_GM = xASL_stat_ComputeMean(imCBF,imMask,290,0,imGM)
% [CBF_GM CBF_WM] = xASL_stat_ComputeMean(imCBF,imMask,290,0,imGM,imWM)
%          CBF_GM = xASL_stat_ComputeMean(imCBF,imMask,[],1,imGM)
%          CBF_GM = xASL_stat_ComputeMean(imCBF,imMask,[],2,imGM,imWM)
% [CBF_GM CBF_WM] = xASL_stat_ComputeMean(imCBF,[],[],2,imGM,imWM)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% REFERENCES: Asllani I, Borogovac A, Brown TR. Regression algorithm correcting for partial volume effects in arterial spin labeling MRI. Magnetic 
%             Resonance in Medicine: An Official Journal of the International Society for Magnetic Resonance in Medicine. 2008 Dec;60(6):1362-71.
% 
%             Petr J, Mutsaerts HJ, De Vita E, Steketee RM, Smits M, Nederveen AJ, Hofheinz F, van den Hoff J, Asllani I. Effects of systematic partial 
%             volume errors on the estimation of gray matter cerebral blood flow with arterial spin labeling MRI. Magnetic Resonance Materials in 
%             Physics, Biology and Medicine. 2018 Dec 1;31(6):725-34.
% __________________________________
% Copyright (C) 2015-2019 ExploreASL
%
% 2017-00-00 HJ+JP

% ------------------------------------------------------------------------------
% Admin
% ------------------------------------------------------------------------------

if nargin<2 || isempty(imMask)
	imMask = ones(size(imCBF));
end
if nargin<3 || isempty(nMinSize)
	nMinSize = 0;
end
if nargin<4 || isempty(bPVC)
	bPVC = 0;
end

if nargin<5 
	imGM = [];
end

if nargin<6 
	imWM = [];
end

if (nargout>1) && isempty(imWM)
	warning('Cannot calculate CBF_WM when imWM is not provided');
end

% ------------------------------------------------------------------------------
% Mask calculations
% ------------------------------------------------------------------------------

% only compute in real data
imMask = imMask>0 & isfinite(imCBF);

% Constrain calculation to the mask and to finite values
imCBF = imCBF(imMask);

% Limits imGM and imWM to imMask, if provided
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

switch (bPVC)
    case -1
    % No PVC and mean
	if isempty(imGM)
		CBF_GM = xASL_stat_MeanNan(imCBF); 
	else
		CBF_GM = xASL_stat_MeanNan(imCBF(imGM>0.7)); 
	end
	if ~isempty(imWM)
		CBF_WM = xASL_stat_MeanNan(imCBF(imWM>0.7));
	end
    case 0
    % No PVC
	if isempty(imGM)
		CBF_GM = xASL_stat_MedianNan(imCBF); % this is non-parametric
	else
		CBF_GM = xASL_stat_MedianNan(imCBF(imGM>0.7)); % this is non-parametric
	end
	if ~isempty(imWM)
		CBF_WM = xASL_stat_MedianNan(imCBF(imWM>0.7)); % this is non-parametric
	end
 case 1
    % Simple PVC
	if isempty(imGM)
		error('imGM needs to be provided for bPVC == 1');
	end
    CBF_GM = xASL_stat_SumNan(imCBF)/xASL_stat_SumNan(imGM);
	
	if nargout > 1
		error('Only CBF_GM can be provided for bPVC == 1');
	end
case 2
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

    % When we assume a GM ROI, we want to include some more voxels with a bit of WM content
	if isempty(imGM) || isempty(imWM)
		error('imGM and imWM need to be provided for bPVC == 2');
	end
	
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

