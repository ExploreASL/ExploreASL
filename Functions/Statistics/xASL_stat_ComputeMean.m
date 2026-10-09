function [meanPVC0, meanPVC1, meanPVC2primary, meanPVC2secondary, medianPVC0] = xASL_stat_ComputeMean(imCBF, imMask, nMinSize, bOutput, imPVprimary, imPVsecondary)
%xASL_stat_ComputeMean calculates mean or median of CBF in the image across a mask with an optional partial volume correction.
%
% FORMAT:  [meanPVC0, meanPVC1, meanPVC2primary, meanPVC2secondary, medianPVC0] = xASL_stat_ComputeMean(imCBF[, imMask, nMinSize, bOutput, imPVprimary, imPVsecondary])
%
% INPUT:
%   imCBF  - input CBF volume (REQUIRED)
%   imMask - mask for the calculation (OPTIONAL, DEFAULT = finite part of imCBF)
%   nMinSize - minimal size of the ROI in voxels, if not big enough, then return NaN
%            - ignore when 0 (OPTIONAL, default = 0)
%   bOutput - vector of length between 1 and 5 that specifies which output is requested
%             [meanPVC0, meanPVC1, meanPVC2primary, meanPVC2secondary, medianPVC0] (OPTIONAL, DEFAULT [1 0 0 0 0])
%   imPVprimary   - Primary partial volume map with the same size as imCBF
%            (OPTIONAL, REQUIRED for bPVC==2 and bPVC==1)
%   imPVsecondary   - Secondary partial volume map with the same size as imCBF
%            (OPTIONAL, REQUIRED for bPVC==2)
% OUTPUT:
%   meanPVC0 - mean value with PVC0
%   meanPVC1 - mean value with PVC1
%   meanPVC2primary - mean value with PVC2 in the primary PV
%   meanPVC2secondary - mean value with PVC2 in the secondary PV
%   medianPVC0   - median value
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: It calculates mean or median of CBF over the mask imMask if the mask volume exceeds nMinSize. It calculates either
%              a mean, a median, or a mean after PVC - depending on what outputs. For the PVC options, it needs also imPVprimary and imPVsecondary and returns the
%              separate PV-corrected values calculated over the entire ROI. It is a single function call that assigns all value. The reason is that we can share 
%              joint calculations that take most of the calculation time. PVC options are: 0 - don't do partial volume correction, just calculate a mean or median on imMask
%              1 - simple partial volume correction by normalizaton by the imPVprimary volume - see Petr et al. 2018
%              2 - partial volume correction using linear regression and imPVprimary, imPVsecondary maps according to Asllani et al. 2008
%
% 1. Admin
% 2. Mask calculations
% 3. Calculate the ROI statistics
% 3a. No PVC and simple mean
% 3b. No PVC and median
% 3c. Simple PVC
% 3d. Full PVC on a region
%
% EXAMPLE: meanPVC0 = xASL_stat_ComputeMean(imCBF)
%          meanPVC0 = xASL_stat_ComputeMean(imCBF,imMask,[])
%          [meanPVC0,~,~,~,medianPVC0] = xASL_stat_ComputeMean(imCBF,[],[],[1,0,0,0,1])
%          [~,~,~,~,medianPVC0] = xASL_stat_ComputeMean(imCBF,imMask,290,[0 0 0 0 1])
%          [meanPVC0, meanPVC1] = xASL_stat_ComputeMean(imCBF,imMask,[],[1 1 0 0 0],imPVprimary)
%          [~,~,meanPVC2primary,~] = xASL_stat_ComputeMean(imCBF,imMask,[],[0 0 1 0 0],imPVprimary,imPVsecondary)
%          [meanPVC0, meanPVC1, meanPVC2primary, meanPVC2secondary, medianPVC0] = xASL_stat_ComputeMean(imCBF,[],[],[1 1 1 1 1],imPVprimary,imPVsecondary)
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
	bOutput = [1 0 0 0 0];
else
	bOutput(end+1:5) = 0;
end

if nargin<5
	imPVprimary = [];
end

if nargin<6 
	imPVsecondary = [];
end

% Initialize the output
meanPVC0 = NaN;
meanPVC1 = NaN;
meanPVC2primary = NaN;
meanPVC2secondary = NaN;
medianPVC0 = NaN;

if nargout>=1 && bOutput(1)
	
else
	bOutput(1) = 0; % Do not calculate that output when it will not be assigned to output parameters
end

if nargout>=2 && bOutput(2)
	if isempty(imPVprimary)
		error('Cannot calculate meanPVC1 when imPVprimary is not provided');
	end
else
	bOutput(2) = 0;
end

if nargout>=3 && bOutput(3)
	if isempty(imPVsecondary) || isempty(imPVprimary)
		error('Cannot calculate meanPVC2primary when imPVprimary and imPVsecondary are not provided');
	end
else
	bOutput(3) = 0;
end

if nargout>=4 && bOutput(4)
	if isempty(imPVsecondary) || isempty(imPVprimary)
		error('Cannot calculate meanPVC2secondary when imPVprimary and imPVsecondary are not provided');
	end
else
	bOutput(4) = 0;
end

if nargout<5
	bOutput(5) = 0;
end

if ~any(imCBF>0, 'all')
    warning('CBF image is empty, skipping');
    return;
end

if bOutput(3) || bOutput(4)
	% If running PVC, then need imPVprimary and imPVsecondary of the same size as imCBF
	if ~isequal(size(imCBF),size(imPVprimary)) || ~isequal(size(imCBF),size(imPVsecondary))
		warning('When running PVC, the primary and secondary PV maps should have the same size as the CBF image');
	elseif size(imPVprimary, 4)>2
        warning('Invalid imPVprimary map size');
    elseif size(imPVsecondary, 4)>2
        warning('Invalid imPVsecondary map size');
	end
end

%% 2. Mask calculations

% Only compute in real data
imMask = imMask>0 & isfinite(imCBF) & (imCBF~=0);

% Constrain calculation to the mask and to finite values
imCBF = imCBF(imMask);

% Limit imPVprimary and imPVsecondary to imMask, if provided
if ~isempty(imPVprimary)
	imPVprimary = imPVprimary(imMask); 
end

if ~isempty(imPVsecondary)
	imPVsecondary = imPVsecondary(imMask); 
end

maskSize = sum(imMask, 'all');

if maskSize < nMinSize || maskSize <= 0 
    return;
end

%% 3. Calculate the ROI statistics
sumCBF = sum(imCBF, 1, 'omitnan');

%% 3a. No PVC and simple mean
if bOutput(1)
	meanPVC0 = sumCBF/maskSize;
end

%% 3b. No PVC and median
if bOutput(5)
	medianPVC0 = median(imCBF, 1, 'omitnan'); % this is non-parametric
end

%% 3c. Simple PVC
if bOutput(2)
	if isempty(imPVprimary)
		error('imPVprimary needs to be provided for bPVC == 1');
	end
	meanPVC1 = sumCBF/sum(imPVprimary, 1, 'omitnan');
end
	
%% 3d. Full PVC on a region
if bOutput(3) || bOutput(4)
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

	gwpv                       = imPVprimary;
	gwpv(:,2)                  = imPVsecondary;
	gwcbf                      = (imCBF')*pinv(gwpv');
	if bOutput(3)
		meanPVC2primary = gwcbf(1);
	end
	if bOutput(4)
		meanPVC2secondary = gwcbf(2);
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

