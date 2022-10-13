function [ScaleImage, CBF, ATT, TExch] = xASL_quant_MultiPLD(PWI, M0_im, imSliceNumber, x, bUseBasilQuantification)
%xASL_quant_MultiPLD % Perform a multi-step quantification using BASIL
% FORMAT: [ScaleImage[, CBF, ATT, TExch]] = xASL_quant_MultiPLD(PWI, M0_im, imSliceNumber, x[, bUseBasilQuantification])
%
% INPUT:
%   PWI             - 4D (4th dimension per PLD) image matrix of perfusion-weighted image (REQUIRED)
%   M0_im           - M0 image (can be a single number or image matrix) (REQUIRED)
%   imSliceNumber   - image matrix showing slice number in current ASL space (REQUIRED for 2D multi-slice)
%   x               - struct containing pipeline environment parameters (REQUIRED)
%   bUseBasilQuantification - boolean, true for using FSL BASIL for
%                             quantification, false for using ExploreASL's
%                             own quantification (OPTIONAL, DEFAULT = true)
%
% OUTPUT:
% ScaleImage        - image matrix containing net/effective quantification scale factor
% CBF               - Quantified CBF image
% ATT               - Estimated ATT map
% TExch             - Estimated map of time of exchange across BBB (if multi-TE is available)
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
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [ScaleImage, CBF, ATT] = xASL_quant_MultiPLD(PWI, M0_im, imSliceNumber, x, bUseBasilQuantification);
% __________________________________
% Copyright (c) 2015-2022 ExploreASL


%% Admin
fprintf('Quantification CBF multi-PLD:\n');

if nargin < 4
	error('Four input parameters required.');
end

if  xASL_stat_SumNan(M0_im(:))==0 && x.modules.asl.ApplyQuantification(5)
    error('Empty M0 image, something went wrong in M0 processing');
elseif xASL_stat_SumNan(M0_im(:))==0
    warning('Weird M0 image detected');    
end

if nargin<5 || isempty(bUseBasilQuantification)
    bUseBasilQuantification = true;
end

if ~bUseBasilQuantification
	error('Multi-PLD quantification currently works only with BASIL');
end

ScaleImage = 1; % initializing (double data format by default in Matlab)

% Convert to double precision to increase quantification precision
% This is especially useful with the large factors that we can multiply and
% divide with in ASL quantification
PWI = double(PWI);
M0_im = double(M0_im);

if ~x.modules.asl.ApplyQuantification(3)
    fprintf('%s\n','We skip the scaling of a.u. to label intensity');
else
    % Calculate the vector of SliceReadoutTime
    SliceReadoutTime = xASL_quant_SliceTiming(x,x.P.Path_ASL4D);
    
 %% 1    PLD scalefactor (gradient if 2D multi-slice)
    % For BASIL the x.Q.SliceReadoutTime is used internally, otherwise
    % x.Q.SliceReadoutTime is added to ScaleImage
    
    switch lower(x.Q.readoutDim)
        case '3d'
            fprintf('%s\n','3D sequence, not accounting for SliceReadoutTime (homogeneous PLD for complete volume)');
            x.Q.SliceReadoutTime = 0;
			x.Q.BasilSliceReadoutTime = 0;

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
			if max(SliceReadoutTime)>0 && length(SliceReadoutTime) > 1
				x.Q.BasilSliceReadoutTime = SliceReadoutTime(2)-SliceReadoutTime(1);
			else
				x.Q.BasilSliceReadoutTime = 0;
			end
            
        otherwise
            error('Wrong x.Q.readoutDim value!');
    end

    if xASL_stat_SumNan(ScaleImage(:))==0
        error('Wrong PLD definition!');
    end
    
    % BASIL/FABBER quantification
    [PWI, ATT, TExch] = xASL_quant_Basil(PWI, x);

end


%% 4    Manufacturer-specific scalefactor
if ~x.modules.asl.ApplyQuantification(1)
    fprintf('%s\n','We skip the Manufacturer-specific scalefactors');
else
    % Load the stored parameters
	ASL_parms = xASL_adm_LoadParms(x.P.Path_ASL4D_parms_mat, x);

	% Throw warning if no Philips scans, but some of the scale slopes are not 1:
	if isempty(regexpi(x.Q.Vendor,'Philips'))
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
	if ~isempty(regexpi(x.Q.Vendor,'GE'))
		if ~isfield(x.Q,'NumberOfAverages')
			% GE accumulates signal instead of averaging by NEX, therefore division by NEX is required
			error('GE-data expected, "NumberOfAverages" should be a dicom-field, but was not found!!!')
		else
			x.Q.NumberOfAverages = max(x.Q.NumberOfAverages); % fix for combination of M0 & PWI in same nifti, for GE quantification
		end

		% For some reason the older GE Alsop Work in Progress (WIP) version
		% has a different scale factor than the current GE product sequence
		if isfield(x.Q,'SoftwareVersions')
			softwareVersions = xASL_str2num(x.Q.SoftwareVersions(1:2));
		else
			softwareVersions = 0;
		end
		
		if ~softwareVersions || softwareVersions > 15
			% GE new version
			% division by x.Q.NumberOfAverages as GE sums difference image instead of averaging
			qnt_R1gain = 32; %  PWI is scaled up by 32 (default GE scalefactor)
			qnt_GEscaleFactor = qnt_R1gain*x.Q.NumberOfAverages;
		else
			% GE older WIP version
			qnt_RGcorr = 45.24; % Correction for receiver gain in PDref (but not used apparently?)
			qnt_GEscaleFactor = qnt_RGcorr*x.Q.NumberOfAverages;
		end

		ScaleImage = ScaleImage./qnt_GEscaleFactor;
		fprintf('%s\n',['Quantification corrected for GE scale factor ' xASL_num2str(qnt_GEscaleFactor) ' for NSA=' xASL_num2str(x.Q.NumberOfAverages)]);

		% Set Philips specific scaling
	elseif ~isempty(regexpi(x.Q.Vendor,'Philips'))
		% Philips has specific scale & rescale slopes
		% If these are not corrected for, only relative CBF quantification can be performed,
		% i.e. scaled to wholebrain, the wholebrain perfusion cannot be calculated.

		scaleFactor = xASL_adm_GetPhilipsScaling(xASL_adm_LoadParms(x.P.Path_ASL4D_parms_mat,x),xASL_io_ReadNifti(x.P.Path_ASL4D));

		if scaleFactor
			ScaleImage = ScaleImage .* scaleFactor;
		end

		% Siemens specific scalings
	elseif strcmpi(x.Q.Vendor,'Siemens')
		if ~strcmpi(x.Q.Vendor,'Siemens_JJ_Wang') && strcmpi(x.Q.M0,'separate_scan')
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

if ~x.modules.asl.ApplyQuantification(3)
	% Convert to double precision if not done previously
	M0_im = double(M0_im);
	PWI = double(PWI);
end

M0_im = repmat(M0_im,MatchSizeM0);
ScaleImage = repmat(ScaleImage,MatchSizeSI);

if ~x.modules.asl.ApplyQuantification(5)
    fprintf('%s\n','We skip the PWI/M0 division');
else
    ScaleImage = ScaleImage./M0_im;
end

%% 6    Apply quantification
if x.modules.asl.ApplyQuantification(6)
    CBF = PWI.*ScaleImage;
else
    CBF = PWI;
end

CBF = xASL_stat_MeanNan(CBF, 4);

%% 6    Print parameters used
fprintf('%s\n',' model with parameters:');

if x.modules.asl.ApplyQuantification(3)
	% Obtain the unique PLDs and LDs
	[PLDs, ind] = unique(x.Q.Initial_PLD, 'stable');
	if length(x.Q.LabelingDuration)>1
		LDs = (x.Q.LabelingDuration(ind))';
	else
		LDs = ones(size(PLDs))*x.Q.LabelingDuration;
	end
	% Skip the first PLD of the block for Time-Encoded
	if x.modules.asl.bTimeEncoded
		numberBlocks = numel(PLDs)/x.Q.TimeEncodedMatrixSize;
		ind = (ones(numberBlocks,1)*(2:x.Q.TimeEncodedMatrixSize) + (0:(numberBlocks-1))' * x.Q.TimeEncodedMatrixSize * ones(1,x.Q.TimeEncodedMatrixSize-1))';
		PLDs = PLDs(ind(:)');
		LDs = LDs(ind(:)');
	end
	
    switch lower(x.Q.LabelingType)
        case 'pasl'
            fprintf('%s\n',['TI1 = ' xASL_num2str(LDs) ' ms,']);
            fprintf('%s\n',['TI  = ' xASL_num2str(PLDs) ' ms.']);
        case 'casl'
            fprintf('%s\n',['LabelingDuration = ' xASL_num2str(LDs) ' ms,']);
            fprintf('%s\n',['PLD              = ' xASL_num2str(PLDs) ' ms.']);
    end

	if max(SliceReadoutTime)>0 && strcmpi(x.Q.readoutDim,'2D')
		fprintf('%s',[' + ' xASL_num2str(SliceReadoutTime(2)-SliceReadoutTime(1)) ' ms*(slice-1)']);
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