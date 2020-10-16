function [M0IM] = xASL_quant_M0(M0IM, x)
%xASL_quant_M0 Quantification of M0
%
% FORMAT: [M0IM] = xASL_quant_M0(M0IM, x)
%
% INPUT:
%   M0IM        - M0 image matrix or path
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
% OUTPUT:
%   x           - struct containing pipeline environment parameters, useful when only initializing ExploreASL/debugging
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function quantifies the M0, except for the difference in voxel size
% between the M0 and ASL source data (which is scaled in
% xASL_wrp_ProcessM0.m). This function runs the following steps:
%
% 1. Correct scale slopes, if Philips
% 2. Skip M0 quantification if ~x.ApplyQuantification(4)
% 3. Set TR specifically for GE
% 4. Check for correct TR values
% 5. Quantify the M0, either for single 3D volume or slice-wise
% 6. Apply custom scalefactor if requested (x.M0_GMScaleFactor)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: M0IM = xASL_quant_M0('MyStudy/sub-001/ASL_1/rM0.nii', x);
% __________________________________
% Copyright 2015-2020 ExploreASL


%% ------------------------------------------------------------------------------------------------------
% Admin
if nargin<2 || isempty(x)
    warning('x input parameter missing, skipping');
elseif isempty(M0IM) || sum(isfinite(M0IM(:)))==0
    warning('Invalid M0 image loaded, skipping');
end

if strcmpi(x.M0,'no_background_suppression')
    warning('This option will be obsolete in the future, please use M0=UseControlAsM0 instead');
    x.M0 = 'UseControlAsM0'; % backward compatibility
end

if strcmpi(x.M0, 'UseControlAsM0')
    fprintf('%s\n','x.M0_usesASLtiming didnt exist, set to 1, because of absence background suppression, same timing in mean control (used as M0) as in ASL');
    x.M0_usesASLtiming = 1;
    M0ParmsMat  = x.P.Path_ASL4D_parms_mat;
elseif ~isfield(x,'M0_usesASLtiming')
    x.M0_usesASLtiming = 0; % default value
    M0ParmsMat  = x.P.Path_M0_parms_mat;
end

% Allow inputting path instead of image
M0IM = xASL_io_Nifti2Im(M0IM);

%% ------------------------------------------------------------------------------------------------------
% 1. Correct scale slopes, if Philips
if ~x.ApplyQuantification(2)
    fprintf('%s\n','M0 ScaleSlopes skipped');
else
    if ~isempty(regexpi(x.Vendor,'Philips'))

		scaleFactor = xASL_adm_GetPhilipsScaling(xASL_adm_LoadParms(M0ParmsMat, x),xASL_io_ReadNifti(x.P.Path_M0));

		if scaleFactor
			M0IM = M0IM .* scaleFactor;
		end
    end
end

%% ------------------------------------------------------------------------------------------------------
% 2. Skip M0 quantification if ~x.ApplyQuantification(4)
if ~x.ApplyQuantification(4)
    fprintf('%s\n','M0 quantification for incomplete T1 relaxation skipped');
    % M0 quantification here is only for the incomplete T1 recovery
    % One potential reason of skipping this, is when a pseudo-M0 is used
    % that has background suppression. In this case, there is no incomplete
    % inversion recovery, the background suppression pulses have saturated
    % the signal instead.
else

    % Correction for incomplete T1 relaxation of M0 because of short TR.
    % Control image/M0 reference value is smaller than should be, paired subtraction is not
    % (i.e. the labeled blood doesn't have this incomplete T1 relaxation
    % issue from short TR).
    % Therefore, correct this (~2% decrease)
    % Using WM value here, since we used the M0 biasfield correction method
    % in which the M0 is roughly eroded to the WM
    if ~isfield(x.Q,'TissueT1')
        fprintf('%s\n','x.Q.TissueT1 did not exist, default=1240 used');
        x.Q.TissueT1 = 800; % 800 ms for WM, average of frontal & occipital WM from Lu et al., JMRI 2005 & 3 studies they refer to
        % for GM, this would be 1300 ms
        % Here we use the WM T1, since we smooth a lot
    end

    %% ------------------------------------------------------------------------------------------------------
    % 3. Set TR specifically for GE
    if ~isempty(regexpi(x.Vendor,'GE')) && isempty(regexpi(x.Vendor,'Siemens')) &&  isempty(regexpi(x.Vendor,'Philips'))
        TR = 2000; % GE does an inversion recovery, which takes 2 s and hence signal has decayed 2 s
        fprintf('%s\n','GE M0 scan, so using 2 s as TR (GE inversion recovery M0)');
    else
        M0_parms = xASL_adm_LoadParms(M0ParmsMat, x);
        TR = M0_parms.RepetitionTime; % This will be either the separate M0-scan value or the ASL scan value (if UseControlAsM0)
    end

    %% ------------------------------------------------------------------------------------------------------
    % 4. Check for correct TR values
    if TR<1000
        error(['Unusually small TR for ASL M0: ' num2str(TR)]);
    end

    if length(TR)>1
        warning(['Multiple M0 TRs found: ' num2str(TR(:)') ', ' num2str(min(TR)) ' used']);
        TR  = min(TR);
    end

    %% ------------------------------------------------------------------------------------------------------
    % 5. Quantify the M0, either for single 3D volume or slice-wise
    if strcmpi(x.readout_dim,'3D') % for 3D readout, we assume the M0 readout is at the end of the TR
            NetTR = TR;
            fprintf('%s\n','Single 3D M0 readout assumed');
    elseif  strcmpi(x.readout_dim,'2D') % for 2D readouts, there are slice timing differences

            % Calculate SliceReadoutTime for shortestTR
            [x] = xASL_quant_SliceReadoutTime_ShortestTR(x);

            SliceIM = zeros(size(M0IM));
            for iZ=1:size(M0IM,3)
                SliceIM(:,:,iZ) = iZ;
            end

            if  x.M0_usesASLtiming % in this case, the M0 readout has the exact same timing as the ASL readout
                                   % this is the case e.g. for Philips 3D GRASE
                NetTR = x.Q.LabelingDuration+x.Q.Initial_PLD+x.Q.SliceReadoutTime.*(SliceIM-1);
                fprintf('%s\n','2D sliceWise M0 readout assumed, same timing as ASL slices readout used');
            elseif ~x.M0_usesASLtiming
                NetTR = TR - ((max(SliceIM(:))-SliceIM).*x.Q.SliceReadoutTime);
                % Here we assume the M0 slices readout were at the end of the
                % TR, with the same time between slices as for the ASL readout
                fprintf('%s\n','2D SliceWise M0 readout, assumed M0 slices readout at end of TR');
            else
                error('Unknown x.M0_usesASLtiming specified');
            end

    else
        error('Unknown x.readout_dim specified');
    end

    corr_T1 = 1 ./ (1-exp(-NetTR/x.Q.TissueT1));
    M0IM = M0IM .* corr_T1;

    min_TR = min(NetTR( NetTR~=0 & isfinite(NetTR)));
    max_TR = max(NetTR( NetTR~=0 & isfinite(NetTR)));
    minCorr_T1 = min(corr_T1( corr_T1~=0 & isfinite(corr_T1)));
    maxCorr_T1 = max(corr_T1( corr_T1~=0 & isfinite(corr_T1)));

    fprintf('%s\n', ['M0 correction of incomplete T1 relaxation: TR range ' num2str(min_TR) '-' num2str(max_TR) ' ms, T1 WM tissue ' num2str(x.Q.TissueT1) ' ms, gives factor ' num2str(minCorr_T1) '-' num2str(maxCorr_T1)]);
end

%% ------------------------------------------------------------------------------------------------------
% 6. Apply custom scalefactor if requested (x.M0_GMScaleFactor)
if ~isfield(x,'M0_GMScaleFactor') || isempty(x.M0_GMScaleFactor)
    x.M0_GMScaleFactor = 1; % no scaling
else
    fprintf('%s\n',['M0 scaling corrected by GMScaleFactor ' xASL_num2str(x.M0_GMScaleFactor)]);
    M0IM = M0IM.*x.M0_GMScaleFactor;
end

% x.M0_GMScaleFactor - add additional scale factor to multiply the M0 image by (OPTIONAL, default = 1)
% This can be useful when you have background suppression but no control/M0
% image without background suppression. If you then know the M0 scalefactor
% for the GM, you can use the control image as M0 and use this parameter to
% scale back what was suppressed by background suppression.
% Note that there is no option for separate tissue scaling (e.g. WM & GM),
% because ExploreASL pragmatically smooths the M0 a lot, assuming that
% head motion and registration between M0 & ASL4D will differ between
% patients and controls.



% For the SliceWise M0 TR correction, we have to define an equation that
% would serve most examples without a too large offset.
% There are sequences that include M0 as a control-label pair within the
% ASL sequence without labeling or background suppression. But also M0
% scans exist which have decreased the TR to e.g. 2 ms, as a separate scan.
% In the former example, calculating the slice timing is the same as in
% ASL, in the latter it is not. In the latter, we assume that the slices
% are acquired near the end of the TR.

% Here we compare the latter equation with the former, and the error is
% small. We will use the latter in general here, and the former is supplied
% as option: x.M0_usesASLtiming. E.g. this is the case for Philips 3D
% GRASE.



% Here, we compare the difference between both equations.
% E.g. suppose a Philips 3D GRASE M0 with TR=4035, x.Q.LabelingDuration =
% 1800, x.Q.Initial_PLD = 1525, x.Q.SliceReadoutTime = 30,
% x.Q.TissueT1 = 1240, nSlices=14, SliceN = [1:nSlices].

% Default Philips equation on the scanner:
% NetTR = x.Q.LabelingDuration+x.Q.Initial_PLD+x.Q.SliceReadoutTime.*(SliceN-1)
% gives [3600,3640,3680,3720,3760,3800,3840,3880,3920,3960,4000,4040,4080,4120;]

% giving correction factor with equationa:
% where corr_T1             = 1 ./ (1-exp(-NetTR/x.Q.TissueT1))
% gives
% [1.05802865189666,1.05608332293690,1.05420654956043,1.05239569649126,1.05064824413736,1.04896178258070,1.04733400593878,1.04576270707089,1.04424577260449,1.04278117825921,1.04136698444764,1.04000133213376,1.03868243893151,1.03740859542711;]

% The more general equation:
% NetTR     = TR - ((nSlices-SliceN).*qnt_PLDslicereadout)
% gives
% [3645,3675,3705,3735,3765,3795,3825,3855,3885,3915,3945,3975,4005,4035;]
% where corr_T1             = 1 ./ (1-exp(-NetTR/x.Q.TissueT1));
% gives
% [1.05584503232182,1.05443748830939,1.05306720374486,1.05173310385921,1.05043414927497,1.04916933462559,1.04793768723894,1.04673826588166,1.04557015956089,1.04443248638046,1.04332439244853,1.04224505083405,1.04119366056948,1.04016944569733;]
% which is a overestimation of
% [0.997936143250067,0.998441567448549,0.998919238534379,0.999370395912622,0.999796225936145,1.00019786425810,1.00057639807047,1.00093286823500,1.00126827131231,1.00158356149466,1.00187965244733,1.00215741906377,1.00241769913868,1.00266129496361;]


end
