function [ScaleImage, CBF, ATT, ABV, Tex] = xASL_quant_ASL(PWI4D, M0_im, imSliceNumber, x, bUseBasilQuantification)
%xASL_quant_ASL Perform a multi-step quantification of single or multi-PLD with or without BASIL
% FORMAT: [ScaleImage[, CBF, ATT, ABV, Tex]] = xASL_quant_ASL(PWI4D, M0_im, imSliceNumber, x[, bUseBasilQuantification])
%
% INPUT:
%   PWI4D           - 4D timeseries of (control-label subtracted) perfusion-weighted images (REQUIRED)
%   M0_im           - M0 image (can be a single number or image matrix) (REQUIRED)
%   imSliceNumber   - image matrix showing slice number in current ASL space (REQUIRED for 2D multi-slice)
%   x               - struct containing pipeline environment parameters (REQUIRED)
%   bUseBasilQuantification - boolean, true for using FSL BASIL for
%                             quantification, false for using ExploreASL's
%                             own quantification (OPTIONAL, DEFAULT = false for singlePLD, true for multiPLD)
%
% OUTPUT:
% ScaleImage        - image matrix containing net/effective quantification scale factor
% CBF               - Quantified CBF image
% ATT               - Estimated ATT map (if multi-PLD, otherwise empty)
% ABV               - Estimated arterial blood volume map (if multi-PLD, otherwise empty)
% Tex               - Estimated map of time of exchange across BBB (if multi-TE is available, otherwise empty)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script performs a multi-step quantification, by
%              initializing a ScaleImage that travels through this script & gets changed by the following quantification
%              factors:
%
%              1.    {{PLD scalefactor}} (gradient if 2D multi-slice) (if x.modules.asl.ApplyQuantification(3))
%              2.    {{Label decay scale factor}} for single (blood T1) - or dual-compartment (blood+tissue T1) model, CASL or PASL
%                    Single-compartment model: Alsop MRM 2014
%                    Dual-compartment model: Wang MRM 2002: Gevers JMRI 2012 (if x.modules.asl.ApplyQuantification(3))
%              3.    {{Scaling to physiological units}} [ml/gr/ms =>ml/100gr/min =>(60,000 ms=>min)(1 gr=>100gr)]
%                    (if x.modules.asl.ApplyQuantification(3))
%              4.    {{Manufacturer-specific scalefactor}} (if x.modules.asl.ApplyQuantification(1) -> future move to dcm2niiX stage)
%              Finally, we:
%              5.    Divide PWI/M0 (if x.modules.asl.ApplyQuantification(5)) & apply all scaling (if x.modules.asl.ApplyQuantification(6))
%              6.    Print parameters used
%
%              Note that the output always goes to the CBF image (in the
%              future this could go to different stages, e.g. dcm2niiX or
%              PWI stage)
%
%              Note that BASIL is also implemented, but it doesn't allow a
%              standard space quantification yet (it would need to use
%              imSliceNumber)
%
%              Note that in Matlab we only quantify single-PLD here,
%              multi-PLD or multi-TE quantifications are performed with BASIL or FABBER here, respectively.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [ScaleImage, CBF, ATT, ABV, Tex] = xASL_quant_ASL(PWI4D, M0_im, imSliceNumber, x, bUseBasilQuantification);
% __________________________________
% Copyright (c) 2015-2023 ExploreASL


%% Admin
fprintf('%s\n', 'Performing quantification');

ATT = [];
Tex = [];
ABV = [];

if nargin<5 || isempty(bUseBasilQuantification)
	if x.modules.asl.bMultiPLD
		bUseBasilQuantification = true;
	else
		bUseBasilQuantification = false;
	end
end

if nargin < 4
	error('Four input parameters required');
elseif isempty(x) || ~isstruct(x)
    error('Illegal x structure');
end

if  xASL_stat_SumNan(M0_im(:))==0 
	if x.modules.asl.ApplyQuantification(5)
		error('Empty M0 image, something went wrong in M0 processing');
	else
		warning('Empty M0 image detected');
	end
end

if x.modules.asl.bMultiPLD && ~bUseBasilQuantification
	error('Multi-PLD quantification currently works only with BASIL');
end

ScaleImage = 1; % initializing (double data format by default in Matlab)
ScaleImageABV = 1;

% Convert to double precision to increase quantification precision
% This is especially useful with the large factors that we can multiply and
% divide with in ASL quantification
PWI4D = double(PWI4D);
M0_im = double(M0_im);

if ~x.modules.asl.ApplyQuantification(3)
    fprintf('%s\n','We skip the scaling of a.u. to label intensity');
	% In case we don't run quantification, we still need to average the PWI image
	PWI = xASL_im_ASLSubtractionAveraging(x, [], PWI4D);
else
    
    % If we have multiple PLDs but requested a single PLD quantification, we use the highest PLD
    if ~x.modules.asl.bMultiPLD && numel(x.Q.uniqueInitial_PLD)>1
		warning('Multiple PLDs detected, taking the highest PLD only');
		fprintf('%s\n', 'A multi-PLD quantification may provide more accurate results');
		x.Q.uniqueInitial_PLD = max(x.Q.uniqueInitial_PLD);
	end
    
    % Calculate the vector of SliceReadoutTime
    SliceReadoutTime = xASL_quant_SliceTiming(x,x.P.Path_ASL4D);
    
    
    %% 1    PLD scalefactor (gradient if 2D multi-slice)
    % For BASIL the x.Q.SliceReadoutTime is used internally, otherwise
    % x.Q.SliceReadoutTime is added to ScaleImage
    
    switch lower(x.Q.readoutDim)
        case '3d'
            fprintf('%s\n','3D sequence, not accounting for SliceReadoutTime (homogeneous PLD for complete volume)');
            x.Q.SliceReadoutTime = 0;
            if bUseBasilQuantification
                x.Q.BasilSliceReadoutTime = 0;
            else
                ScaleImage = ScaleImage.*x.Q.uniqueInitial_PLD;
            end

        case '2d' % Load slice gradient
            fprintf('%s\n','2D sequence, accounting for SliceReadoutTime');

            imSliceNumber = double(imSliceNumber);
            % Correct NaNs
            % Fix reslicing errors:
            if sum(isnan(imSliceNumber(:)))>0
                imSliceNumber = xASL_im_ndnanfilter(imSliceNumber,'gauss',[16 16 16],2); % code doesn't smooth, only extrapolates
            end

            imSliceNumber(~isfinite(imSliceNumber)) = 1;
            imSliceNumber(imSliceNumber<1) = 1;
			
			imSliceNumber = round(imSliceNumber);
			imSliceNumber(imSliceNumber<1) = 1;
			imSliceNumber(imSliceNumber>length(SliceReadoutTime)) = length(SliceReadoutTime);
            
            % BASIL doesn't use a vector but a difference between slices
			if bUseBasilQuantification
				if max(SliceReadoutTime)>0 && length(SliceReadoutTime) > 1
					x.Q.BasilSliceReadoutTime = SliceReadoutTime(2)-SliceReadoutTime(1);
				else
					x.Q.BasilSliceReadoutTime = 0;
				end
			else
				ScaleImage = ScaleImage.*(x.Q.uniqueInitial_PLD + SliceReadoutTime(imSliceNumber)); % effective/net PLD            
			end
            
        otherwise
            error('Wrong x.Q.readoutDim value!');
    end

    if xASL_stat_SumNan(ScaleImage(:))==0
        error('Wrong PLD definition!');
    end


    %% 2. Run BASIL quantification
    if bUseBasilQuantification
        % Here we perform FSL quantification

        % #1543 We may want to move the NaN extrapolation here

		[PWI, ATT, ABV, Tex] = xASL_quant_FSL(x.P.Path_PWI4D, x); % #1543 Quick and dirty fix for the NaN extrapolation at the start of xASL_quant_FSL
		
		% If resultFSL is not 0, something went wrong
        % This will issue a warning inside xASL_quant_FSL
	else
        % This part should only run when we don't use FSL BASIL/FABBER
        % or as a fallback when FSL BASIL/FABBER crashed
        
        % First, we average PWI4D into PWI, for single-PLD only
        % (later, when we have multi-PLD, multi-echo, or multi-labeling quantification here as well,
        % then we could use PWI3D here)

		if isfield(x.Q, 'SaveCBF4D') && x.Q.SaveCBF4D
            PWI = PWI4D;
            fprintf('%s\n', 'Quantifying CBF in 4D');
            % In this case, we want to save a quantified version of PWI4D
		else
            PWI = xASL_im_ASLSubtractionAveraging(x, [], PWI4D);
		end

        %% 3    Label decay scale factor for single (blood T1) - or dual-compartment (blood+tissue T1) model, CASL or PASL
        if isfield(x.Q,'LabelingType') && isfield(x.Q,'uniqueLabelingDuration')
            ScaleImage = xASL_sub_ApplyLabelDecayScaleFactor(x, ScaleImage);
        else
            warning('Please define both LabelingType and LabelingDuration of this dataset...');
        end

        %% 4    Scaling to physiological units
        % (For some reason, GE sometimes doesn't need the 1 gr->100 gr conversion)
        % (& old Siemens sequence also didn't need the 1 gr->100 gr conversion)
        ScaleImage = ScaleImage.*60000.*100.*x.Q.Lambda;
    end
end


%% 4    Manufacturer-specific scalefactor
if ~x.modules.asl.ApplyQuantification(1)
    fprintf('%s\n','We skip the Manufacturer-specific scalefactors');
else
    % Load the stored parameters
	ASL_parms = xASL_adm_LoadParms(x.P.Path_ASL4D_parms_mat, x);

	% Throw warning if no Philips scans, but some of the scale slopes are not 1:
	if isempty(regexpi(x.Q.Vendor, 'Philips', 'once'))
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
	if ~isempty(regexpi(x.Q.Vendor, 'GE', 'once'))
		if ~isfield(x.Q,'NumberOfAverages')
			% GE accumulates signal instead of averaging by NEX, therefore division by NEX is required
			error('GE-data assumed, "NumberOfAverages" should be a BIDS field, but was not found')
		else
			x.Q.NumberOfAverages = max(x.Q.NumberOfAverages); % fix for combination of M0 & PWI in same nifti, for GE quantification
		end

		% For some reason the older GE Alsop Work in Progress (WIP) version
		% has a different scale factor than the current GE product sequence
		if isfield(x.Q, 'SoftwareVersions')
			softwareVersions = xASL_str2num(x.Q.SoftwareVersions(1:2));
		else
			softwareVersions = 0;
		end
		
		if ~softwareVersions || softwareVersions > 15
			qnt_ReceiverGain = 32; % PWI is scaled up by 32 (default GE scalefactor)
		else
			% GE older WIP version
			qnt_ReceiverGain = 45.24;
		end
        
        qnt_GEscaleFactor = qnt_ReceiverGain*x.Q.NumberOfAverages;
        % division by x.Q.NumberOfAverages as GE sums difference image instead of averaging

		ScaleImage = ScaleImage./qnt_GEscaleFactor;
		ScaleImageABV = ScaleImageABV./qnt_GEscaleFactor;
		fprintf('%s\n',['Quantification corrected for GE scale factor ' xASL_num2str(qnt_GEscaleFactor) ' for NSA=' xASL_num2str(x.Q.NumberOfAverages)]);

		% Set Philips specific scaling
	elseif ~isempty(regexpi(x.Q.Vendor,'Philips', 'once'))
		% Philips has specific scale & rescale slopes
		% If these are not corrected for, only relative CBF quantification can be performed,
		% i.e. scaled to wholebrain, the wholebrain perfusion cannot be calculated.

		scaleFactor = xASL_adm_GetPhilipsScaling(xASL_adm_LoadParms(x.P.Path_ASL4D_parms_mat,x),xASL_io_ReadNifti(x.P.Path_ASL4D));

		if scaleFactor
			ScaleImage = ScaleImage .* scaleFactor;
			ScaleImageABV = ScaleImageABV .* scaleFactor;
		end

		% Siemens specific scalings
	elseif strcmpi(x.Q.Vendor,'Siemens')
		if ~strcmpi(x.Q.Vendor,'Siemens_JJ_Wang') && strcmpi(x.Q.M0,'separate_scan')
			% Some Siemens readouts divide M0 by 10, others don't
			ScaleImage = ScaleImage./10;
			ScaleImageABV = ScaleImageABV./10;
			fprintf('%s\n','M0 corrected for Siemens scale factor 10')
		end
	end
end


%% 5    Divide PWI/M0
% Match sizes
MatchSizeM0 = round([size(PWI,1)./size(M0_im,1) size(PWI,2)./size(M0_im,2) size(PWI,3)./size(M0_im,3) size(PWI,4)./size(M0_im,4) size(PWI,5)./size(M0_im,5) size(PWI,6)./size(M0_im,6) size(PWI,7)./size(M0_im,7)]);
MatchSizeSI = round([size(PWI,1)./size(ScaleImage,1) size(PWI,2)./size(ScaleImage,2) size(PWI,3)./size(ScaleImage,3) size(PWI,4)./size(ScaleImage,4) size(PWI,5)./size(ScaleImage,5) size(PWI,6)./size(ScaleImage,6) size(PWI,7)./size(ScaleImage,7)]);
MatchSizeSIABV = round([size(PWI,1)./size(ScaleImageABV,1) size(PWI,2)./size(ScaleImageABV,2) size(PWI,3)./size(ScaleImageABV,3) size(PWI,4)./size(ScaleImageABV,4) size(PWI,5)./size(ScaleImageABV,5) size(PWI,6)./size(ScaleImageABV,6) size(PWI,7)./size(ScaleImageABV,7)]);

if sum(MatchSizeM0==0) || sum(MatchSizeSI==0)
    error('PWI dimensions too small compared to M0 and/or ScaleImage dimensions');
end

if ~x.modules.asl.ApplyQuantification(3)
	% Convert to double precision if not done previously
	M0_im = double(M0_im);
	PWI = double(PWI);
    ABV = double(ABV);
end

M0_im = repmat(M0_im, MatchSizeM0);
ScaleImage = repmat(ScaleImage, MatchSizeSI);
ScaleImageABV = repmat(ScaleImageABV, MatchSizeSIABV);

if ~x.modules.asl.ApplyQuantification(5)
    fprintf('%s\n','We skip the PWI/M0 division');
else
    ScaleImage = ScaleImage./M0_im;
	ScaleImageABV = ScaleImageABV./M0_im;
end


%% 6    Apply scaling/quantification
if x.modules.asl.ApplyQuantification(6)
    CBF = PWI.*ScaleImage;
	if numel(ABV) > 1
		ABV = ABV.*ScaleImageABV;
	end
else
    CBF = PWI;
end

%% 6    Print parameters used
fprintf('%s\n','with parameters:');

if x.modules.asl.ApplyQuantification(3)
	% Print the prepared parameters
	switch lower(x.Q.LabelingType)
		case 'pasl'
			fprintf('%s\n',['TI1 = ' xASL_num2str(x.Q.uniqueLabelingDuration) ' ms,']);
			fprintf('%s\n',['TI  = ' xASL_num2str(x.Q.uniqueInitial_PLD) ' ms.']);
		case 'casl'
			fprintf('%s\n',['LabelingDuration = ' xASL_num2str(x.Q.uniqueLabelingDuration) ' ms,']);
			fprintf('%s\n',['PLD              = ' xASL_num2str(x.Q.uniqueInitial_PLD) ' ms.']);
		otherwise
			fprintf('Labeling duration and/or PLD undefined...');
	end

	if max(SliceReadoutTime)>0 && strcmpi(x.Q.readoutDim,'2D')
		fprintf('%s\n',[' + ' xASL_num2str(SliceReadoutTime(2)-SliceReadoutTime(1)) ' ms*(slice-1).']);
    else
        fprintf('.\n');
	end

    fprintf('\n%s',['labeling efficiency (neck*Bsup) = ' xASL_num2str(x.Q.LabEff_Orig) ' * ' xASL_num2str(x.Q.LabEff_Bsup) ', ']);
    fprintf('\n%s','assuming ');
    fprintf('%s',['labda = ' xASL_num2str(x.Q.Lambda) ', ']);
    fprintf('%s\n',['T1 arterial blood = ' xASL_num2str(x.Q.BloodT1) ' ms']);

    if x.Q.nCompartments==2
        fprintf('%s',['ATT = ' xASL_num2str(x.Q.ATT) ' ms, ']);
        fprintf('%s\n',['T1tissue = ' xASL_num2str(x.Q.TissueT1) ' ms']);
    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determine Label Decay Scale Factor
function ScaleImage = xASL_sub_ApplyLabelDecayScaleFactor(x, ScaleImage)

    if numel(x.Q.uniqueLabelingDuration)>1
        warning('Quantification with multiple Labeling Durations not implemented yet, taking the first');
        x.Q.uniqueLabelingDuration = x.Q.uniqueLabelingDuration(1);
    end

    switch x.Q.nCompartments
        case 1 % single-compartment model
            switch lower(x.Q.LabelingType)
                case 'pasl'
                    DivisionFactor = x.Q.uniqueLabelingDuration;
                    fprintf('%s\n','Using a single-compartment PASL model');
                case 'casl'
                    DivisionFactor = x.Q.BloodT1 .* (1 - exp(-x.Q.uniqueLabelingDuration./x.Q.BloodT1));
                    fprintf('%s\n','Using a single-compartment CASL model');
            end

            ScaleImage = exp((ScaleImage./x.Q.BloodT1)) ./ (2.*x.Q.LabelingEfficiency.* DivisionFactor);

        case 2 % dual-compartment model
            switch lower(x.Q.LabelingType)
                case 'pasl'
                    DivisionFactor = x.Q.uniqueLabelingDuration;
                    ScaleImage = exp((x.Q.ATT./x.Q.BloodT1)).*exp(((ScaleImage-x.Q.ATT)./x.Q.TissueT1))./ (2.*x.Q.LabelingEfficiency.*DivisionFactor);
                    fprintf('%s\n','Using a dual-compartment PASL model');
                case 'casl'
                    DivisionFactor = x.Q.TissueT1 .* (exp((min(x.Q.ATT-ScaleImage,0))./x.Q.TissueT1) - exp((min(x.Q.ATT-x.Q.uniqueLabelingDuration-ScaleImage,0))./x.Q.TissueT1));
                    ScaleImage = exp((x.Q.ATT./x.Q.BloodT1))./ (2.*x.Q.LabelingEfficiency.* DivisionFactor);
                    fprintf('%s\n','Using a dual-compartment CASL model');
            end

        otherwise
            error('Unknown x.Q.nCompartments');
    end

end
