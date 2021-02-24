function [ScaleImage, CBF] = xASL_quant_SinglePLD(PWI, M0_im, SliceGradient, x)
%xASL_quant_SinglePLD % Perform a multi-step quantification
% FORMAT: [ScaleImage[, CBF]] = xASL_quant_SinglePLD(PWI, M0_im, SliceGradient, x)
%
% INPUT:
%   PWI             - image matrix of perfusion-weighted image (REQUIRED)
%   M0_im           - M0 image (can be a single number or image matrix) (REQUIRED)
%   SliceGradient   - image matrix showing slice number in current ASL space (REQUIRED for 2D multi-slice)
%   x               - struct containing pipeline environment parameters (REQUIRED)
%
% OUTPUT:
% ScaleImage        - image matrix containing net/effective quantification scale factor
% CBF               - Quantified CBF image
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script performs a multi-step quantification, by
%              initializing a ScaleImage that travels through this script & gets changed by the following quantification
%              factors:
%
%              1.    {{PLD scalefactor}} (gradient if 2D multi-slice) (if x.ApplyQuantification(3))
%              2.    {{Label decay scale factor}} for single (blood T1) - or dual-compartment (blood+tissue T1) model, CASL or PASL
%                    Single-compartment model: Alsop MRM 2014
%                    Dual-compartment model: Wang MRM 2002: Gevers JMRI 2012 (if x.ApplyQuantification(3))
%              3.    {{Scaling to physiological units}} [ml/gr/ms =>ml/100gr/min =>(60,000 ms=>min)(1 gr=>100gr)]
%                    (if x.ApplyQuantification(3))
%              4.    {{Vendor-specific scalefactor}} (if x.ApplyQuantification(1) -> future move to dcm2niiX stage)
%              Finally, we:
%              5.    Divide PWI/M0 (if x.ApplyQuantification(5))
%              6.    Print parameters used
%
%              Note that the output always goes to the CBF image (in the
%              future this could go to different stages, e.g. dcm2niiX or
%              PWI stage)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [ScaleImage, CBF] = xASL_quant_SinglePLD(PWI, M0_im, SliceGradient, x);
% __________________________________
% Copyright 2015-2019 ExploreASL


%% Admin
fprintf('%s\n','Quantification CBF single PLD:');

if  xASL_stat_SumNan(M0_im(:))==0
    error('Empty M0 image, something went wrong in M0 processing');
end

x.Q.Initial_PLD = unique(x.Q.Initial_PLD);
if numel(x.Q.Initial_PLD)>1
    warning('Multiple PLDs detected, selecting maximal value');
    fprintf('%s\n', 'A multi-PLD quantification may provide more accurate results');
    x.Q.Initial_PLD = max(x.Q.Initial_PLD);
end

ScaleImage = 1; % initializing (double data format by default in Matlab)

% Convert to double precision to increase quantification precision
% This is especially useful with the large factors that we can multiply and
% divide with in ASL quantification
PWI = double(PWI);
M0_im = double(M0_im);

if ~x.ApplyQuantification(3)
    fprintf('%s\n','We skip the scaling of a.u. to label intensity');
else
	% Calculate the vector of SliceTime
	imASL = xASL_io_ReadNifti(x.P.Path_ASL4D);
	nSlices = size(imASL.dat,3);
	SliceTime = xASL_quant_SliceTimeVector(x,nSlices);
	
    %% 1    PLD scalefactor (gradient if 2D multi-slice)
    switch lower(x.readout_dim)
        case '3d'
            fprintf('%s\n','3D sequence, not accounting for SliceReadoutTime (homogeneous PLD for complete volume)');
            ScaleImage = ScaleImage.*x.Q.Initial_PLD;

        case '2d' % Load slice gradient
            fprintf('%s\n','2D sequence, accounting for SliceReadoutTime');

            SliceGradient = double(SliceGradient);
            % Correct NaNs
            % Fix reslicing errors:
            if sum(isnan(SliceGradient(:)))>0
                SliceGradient = xASL_im_ndnanfilter(SliceGradient,'gauss',[16 16 16],2); % code doesn't smooth, only extrapolates
            end

            SliceGradient(~isfinite(SliceGradient)) = 1;
            SliceGradient(SliceGradient<1) = 1;
			
			imGradientRounded = round(SliceGradient);
			imGradientRounded(imGradientRounded<1) = 1;
			imGradientRounded(imGradientRounded>length(SliceTime)) = length(SliceTime);
			ScaleImage = ScaleImage.*(x.Q.Initial_PLD + SliceTime(imGradientRounded)); % effective/net PLD
        otherwise
            error('Wrong x.readout_dim value!');
    end

    if xASL_stat_SumNan(ScaleImage(:))==0
        error('Wrong PLD definition!');
    end



    %% 2    Label decay scale factor for single (blood T1) - or dual-compartment (blood+tissue T1) model, CASL or PASL
    switch x.Q.nCompartments
        case 1 % single-compartment model
            switch lower(x.Q.LabelingType)
                case 'pasl'
                    DivisionFactor = x.Q.LabelingDuration;
                    fprintf('%s\n','Using single-compartment PASL');
                case 'casl'
                    DivisionFactor = x.Q.BloodT1 .* (1 - exp(-x.Q.LabelingDuration./x.Q.BloodT1));
                    fprintf('%s\n','Using single-compartment CASL');
            end

            ScaleImage = exp((ScaleImage./x.Q.BloodT1)) ./ (2.*x.Q.LabelingEfficiency.* DivisionFactor);

        case 2 % dual-compartment model
            switch lower(x.Q.LabelingType)
                case 'pasl'
                    DivisionFactor = x.Q.LabelingDuration;
                    ScaleImage = exp((x.Q.ATT./x.Q.BloodT1)).*exp(((ScaleImage-x.Q.ATT)./x.Q.TissueT1))./ (2.*x.Q.LabelingEfficiency.*DivisionFactor);
                    fprintf('%s\n','Using dual compartment PASL');
                case 'casl'
                    DivisionFactor = x.Q.TissueT1 .* (exp((min(x.Q.ATT-ScaleImage,0))./x.Q.TissueT1) - exp((min(x.Q.ATT-x.Q.LabelingDuration-ScaleImage,0))./x.Q.TissueT1));
                    ScaleImage = exp((x.Q.ATT./x.Q.BloodT1))./ (2.*x.Q.LabelingEfficiency.* DivisionFactor);
                    fprintf('%s\n','Using dual compartment CASL');
            end

        otherwise
            error('Unknown x.Q.nCompartments');
    end


    %% 3    Scaling to physiological units
    ScaleImage = ScaleImage.*60000.*100.*x.Q.Lambda;
    % (For some reason, GE sometimes doesn't need the 1 gr->100 gr conversion)
    % & old Siemens sequence also didn't need the 1 gr->100 gr conversion
end


%% 4    Vendor-specific scalefactor
if ~x.ApplyQuantification(1)
    fprintf('%s\n','We skip the vendor-specific scalefactors');
else
    % Load the stored parameters
	ASL_parms = xASL_adm_LoadParms(x.P.Path_ASL4D_parms_mat, x);

	% Throw warning if no Philips scans, but some of the scale slopes are not 1:
	if isempty(regexpi(x.Vendor,'Philips'))
		if isfield(ASL_parms,'RescaleSlopeOriginal') && ASL_parms.RescaleSlopeOriginal~=1
			warning('We detected a RescaleSlopeOriginal~=1, verify that this is not a Philips scan!!!');
		end
		if isfield(ASL_parms,'MRScaleSlope') && ASL_parms.MRScaleSlope~=1
			warning('We detected a ScaleSlope~=1, verify that this is not a Philips scan!!!');
		end
		if isfield(ASL_parms,'RWVSlope') && ASL_parms.RWVSlope~=1
			warning('We detected a RWVSlope~=1, verify that this is not a Philips scan!!!');
		end
	end

	% Set GE specific scalings
	if ~isempty(regexpi(x.Vendor,'GE'))
		if ~isfield(x.Q,'NumberOfAverages')
			% GE accumulates signal instead of averaging by NEX, therefore division by NEX is required
			error('GE-data expected, "NumberOfAverages" should be a dicom-field, but was not found!!!')
		else
			x.Q.NumberOfAverages = max(x.Q.NumberOfAverages); % fix for combination of M0 & PWI in same nifti, for GE quantification
		end

		switch lower(x.Vendor)
			% For some reason the older GE Alsop Work in Progress (WIP) version
			% has a different scale factor than the current GE product sequence

			case {'ge_product','ge'} % GE new version
				%                 qnt_R1gain = 1/32;
				%                 qnt_C1 = 6000; % GE constant multiplier

				%                 qnt_GEscaleFactor = (qnt_C1*qnt_R1gain)/(x.Q.NumberOfAverages); % OLD incorrect
				qnt_R1gain = 32; %  PWI is scaled up by 32 (default GE scalefactor)
				qnt_GEscaleFactor = qnt_R1gain*x.Q.NumberOfAverages;
				% division by x.Q.NumberOfAverages as GE sums difference image instead of averaging

			case 'ge_wip' % GE old version
				qnt_RGcorr = 45.24; % Correction for receiver gain in PDref (but not used apparently?)
				% or should this be 6000/45.24?
				qnt_GEscaleFactor = qnt_RGcorr*x.Q.NumberOfAverages;
			otherwise
				error('Please set x.Vendor to GE_product or GE_WIP');
		end

		ScaleImage = ScaleImage./qnt_GEscaleFactor;
		fprintf('%s\n',['Quantification corrected for GE scale factor ' num2str(qnt_GEscaleFactor) ' for NSA=' num2str(x.Q.NumberOfAverages)]);

		% Set Philips specific scaling
	elseif ~isempty(regexpi(x.Vendor,'Philips'))
		% Philips has specific scale & rescale slopes
		% If these are not corrected for, only relative CBF quantification can be performed,
		% i.e. scaled to wholebrain, the wholebrain perfusion cannot be calculated.

		scaleFactor = xASL_adm_GetPhilipsScaling(xASL_adm_LoadParms(x.P.Path_ASL4D_parms_mat,x),xASL_io_ReadNifti(x.P.Path_ASL4D));

		if scaleFactor
			ScaleImage = ScaleImage .* scaleFactor;
		end

		% Siemens specific scalings
	elseif strcmpi(x.Vendor,'Siemens')
		if ~strcmpi(x.Vendor,'Siemens_JJ_Wang') && strcmpi(x.M0,'separate_scan')
			% Some Siemens readouts divide M0 by 10, others don't
			ScaleImage = ScaleImage./10;
			fprintf('%s\n','M0 corrected for Siemens scale factor 10')
		end
	end
end


%% 5    Divide PWI/M0
% Match sizes
MatchSizeM0 = round([size(PWI,1)./size(M0_im,1) size(PWI,2)./size(M0_im,2) size(PWI,3)./size(M0_im,3) size(PWI,4)./size(M0_im,4) size(PWI,5)./size(M0_im,5) size(PWI,6)./size(M0_im,6) size(PWI,7)./size(M0_im,7)]);
MatchSizeSI = round([size(PWI,1)./size(ScaleImage,1) size(PWI,2)./size(ScaleImage,2) size(PWI,3)./size(ScaleImage,3) size(PWI,4)./size(ScaleImage,4) size(PWI,5)./size(ScaleImage,5) size(PWI,6)./size(ScaleImage,6) size(PWI,7)./size(ScaleImage,7)]);

if sum(MatchSizeM0==0) || sum(MatchSizeSI==0)
    error('PWI dimensions too small compared to M0 and/or ScaleImage dimensions');
end

M0_im = repmat(M0_im,MatchSizeM0);
ScaleImage = repmat(ScaleImage,MatchSizeSI);

if ~x.ApplyQuantification(5)
    fprintf('%s\n','We skip the PWI/M0 division');
else
    ScaleImage = ScaleImage./M0_im;
end

%% 6    Apply quantification
CBF = PWI.*ScaleImage;


%% 6    Print parameters used
fprintf('%s\n',' model with parameters:');

if x.ApplyQuantification(3)
    switch lower(x.Q.LabelingType)
        case 'pasl'
            fprintf('%s',['TI1 = ' num2str(x.Q.LabelingDuration) ' ms, ']);
            fprintf('%s',['TI (ms) = ' num2str(x.Q.Initial_PLD)]);
        case 'casl'
            fprintf('%s',['LabelingDuration = ' num2str(x.Q.LabelingDuration) ' ms, ']);
            fprintf('%s',['PLD (ms) = ' num2str(x.Q.Initial_PLD)]);
	end

	if max(SliceTime)>0 && strcmpi(x.readout_dim,'2D')
		fprintf('%s',[' + ' num2str(SliceTime(2)-SliceTime(1)) ' ms*(slice-1)']);
	end

    fprintf('\n%s',['labeling efficiency (neck*Bsup) = ' num2str(x.Q.LabEff_Orig) ' * ' num2str(x.Q.LabEff_Bsup) ', ']);
    fprintf('\n%s','assuming ');
    fprintf('%s',['labda = ' num2str(x.Q.Lambda) ', ']);
    fprintf('%s\n',['T1 arterial blood = ' num2str(x.Q.BloodT1) ' ms']);

    if x.Q.nCompartments==2
        fprintf('%s',['ATT = ' num2str(x.Q.ATT) ' ms, ']);
        fprintf('%s\n',['T1tissue = ' num2str(x.Q.TissueT1) ' ms']);
    end
end


end
