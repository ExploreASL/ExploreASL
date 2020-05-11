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
%              1)    PLD scalefactor (gradient if 2D multi-slice) (if x.ApplyQuantification(3))
%              2)    Label decay scale factor for single (blood T1) - or dual-compartment (blood+tissue T1) model, CASL or PASL
%                    Single-compartment model: Alsop MRM 2014
%                    Dual-compartment model: Wang MRM 2002: Gevers JMRI 2012 (if x.ApplyQuantification(3))
%              3)    Scaling to physiological units [ml/gr/ms =>ml/100gr/min =>(60,000 ms=>min)(1 gr=>100gr)]
%                    (if x.ApplyQuantification(3))
%              4)    Vendor-specific scalefactor (if x.ApplyQuantification(1) -> future move to dcm2niiX stage)
%              Finally, we:
%              5)    Divide PWI/M0 (if x.ApplyQuantification(5))
%              6)    Print parameters used
%              Note that the output always goes to the CBF image (in the
%              future this could go to different stages, e.g. dcm2niiX or
%              PWI stage)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [ScaleImage, CBF] = xASL_quant_SinglePLD(PWI, M0_im, SliceGradient, x);
% __________________________________
% Copyright 2015-2019 ExploreASL


%% Admin
fprintf('%s\n','Quantification CBF single PLD:');

if  xASL_stat_SumNan(M0_im(:))==0
    error('Empty M0 image, something went wrong in M0 processing');
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
    %% 1    PLD scalefactor (gradient if 2D multi-slice)
    switch x.readout_dim
        case '3D'
            fprintf('%s\n','3D sequence, not accounting for SliceReadoutTime (homogeneous PLD for complete volume)');
            ScaleImage = ScaleImage.*x.Q.Initial_PLD;

        case '2D' % Load slice gradient
            fprintf('%s\n','2D sequence, accounting for SliceReadoutTime');

            if ~isnumeric(x.Q.SliceReadoutTime)
                if strcmp(x.Q.SliceReadoutTime,'shortestTR')
                    ASL_parms = xASL_adm_LoadParms(x.P.Path_ASL4D_parms_mat, x);
                    if  isfield(ASL_parms,'RepetitionTime')
                        %  Load original file to get nSlices
                        nSlices = size(xASL_io_Nifti2Im(x.P.Path_ASL4D),3);
                        x.Q.SliceReadoutTime = (ASL_parms.RepetitionTime-x.Q.LabelingDuration-x.Q.Initial_PLD)/nSlices;
                    else
                        error('ASL_parms.RepetitionTime was expected but did not exist!');
                    end
                end
            end

            SliceGradient = double(SliceGradient);
            % Correct NaNs
            % Fix reslicing errors:
            if sum(isnan(SliceGradient(:)))>0
                SliceGradient = xASL_im_ndnanfilter(SliceGradient,'gauss',[16 16 16],2); % code doesn't smooth, only extrapolates
            end

            SliceGradient(~isfinite(SliceGradient)) = 1;
            SliceGradient(SliceGradient<1) = 1;

            ScaleImage = ScaleImage.*(x.Q.Initial_PLD + ((SliceGradient-1) .* x.Q.SliceReadoutTime)); % effective/net PLD
        otherwise
            error('Wrong x.readout_dim value!');
    end

    if xASL_stat_SumNan(ScaleImage(:))==0
        error('Wrong PLD definition!');
    end



    %% 2    Label decay scale factor for single (blood T1) - or dual-compartment (blood+tissue T1) model, CASL or PASL
    switch x.Q.nCompartments
        case 1 % single-compartment model
            switch x.Q.LabelingType
                case 'PASL'
                    DivisionFactor = x.Q.LabelingDuration;
                    fprintf('%s\n','Using single-compartment PASL');
                case 'CASL'
                    DivisionFactor = x.Q.BloodT1 .* (1 - exp(-x.Q.LabelingDuration./x.Q.BloodT1));
                    fprintf('%s\n','Using single-compartment CASL');
            end

            ScaleImage = exp((ScaleImage./x.Q.BloodT1)) ./ (2.*x.Q.LabelingEfficiency.* DivisionFactor);

        case 2 % dual-compartment model
            switch x.Q.LabelingType
                case 'PASL'
                    DivisionFactor = x.Q.LabelingDuration;
                    ScaleImage = exp((x.Q.ATT./x.Q.BloodT1)).*exp(((ScaleImage-x.Q.ATT)./x.Q.TissueT1))./ (2.*x.Q.LabelingEfficiency.*DivisionFactor);
                    fprintf('%s\n','Using dual compartment PASL');
                case 'CASL'
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
	if isempty(regexp(x.Vendor,'Philips'))
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
	if ~isempty(regexp(x.Vendor,'GE'))
		if ~isfield(x.Q,'NumberOfAverages')
			% GE accumulates signal instead of averaging by NEX, therefore division by NEX is required
			error('GE-data expected, "NumberOfAverages" should be a dicom-field, but was not found!!!')
		else
			x.Q.NumberOfAverages = max(x.Q.NumberOfAverages); % fix for combination of M0 & PWI in same nifti, for GE quantification
		end
		
		switch x.Vendor
			% For some reason the older GE Alsop Work in Progress (WIP) version
			% has a different scale factor than the current GE product sequence
			
			case 'GE_product' % GE new version
				%                 qnt_R1gain = 1/32;
				%                 qnt_C1 = 6000; % GE constant multiplier
				
				%                 qnt_GEscaleFactor = (qnt_C1*qnt_R1gain)/(x.Q.NumberOfAverages); % OLD incorrect
				qnt_R1gain = 32; %  PWI is scaled up by 32 (default GE scalefactor)
				qnt_GEscaleFactor = qnt_R1gain*x.Q.NumberOfAverages;
				% division by x.Q.NumberOfAverages as GE sums difference image instead of averaging
				
			case 'GE_WIP' % GE old version
				qnt_RGcorr = 45.24; % Correction for receiver gain in PDref (but not used apparently?)
				% or should this be 6000/45.24?
				qnt_GEscaleFactor = qnt_RGcorr*x.Q.NumberOfAverages;
			otherwise
				error('Please set x.Vendor to GE_product or GE_WIP');
		end
		
		ScaleImage = ScaleImage./qnt_GEscaleFactor;
		fprintf('%s\n',['Quantification corrected for GE scale factor ' num2str(qnt_GEscaleFactor) ' for NSA=' num2str(x.Q.NumberOfAverages)]);
		
		% Set Philips specific scaling
	elseif ~isempty(regexp(x.Vendor,'Philips'))
		% Philips has specific scale & rescale slopes
		% If these are not corrected for, only relative CBF quantification can be performed,
		% i.e. scaled to wholebrain, the wholebrain perfusion cannot be calculated.
		
		% Read the NIFTI header to check the scaling there
		HeaderTemp = xASL_io_ReadNifti(x.P.Path_ASL4D); % original NIfTI
		Rescale_NIfTI = HeaderTemp.dat.scl_slope;
		
		if isfield(ASL_parms,'RWVSlope')
			% If the RealWorldValue is present, then dcm2nii scales to them and ignores everything else
			if ~isnear(ASL_parms.RWVSlope,Rescale_NIfTI,ASL_parms.RWVSlope/100) && (Rescale_NIfTI ~= 1)
				fprintf('%s\n', ['RWVSlope (' xASL_num2str(ASL_parms.RWVSlope) ') and NIfTI slope (' xASL_num2str(Rescale_NIfTI) ') differ, using RWVSlope']);
			end
			if isfield(ASL_parms,'RescaleSlopeOriginal') && ~isnear(ASL_parms.RWVSlope,ASL_parms.RescaleSlopeOriginal,ASL_parms.RWVSlope/100) && (Rescale_NIfTI ~= 1)
				fprintf('%s\n', ['RWVSlope (' xASL_num2str(ASL_parms.RWVSlope) ') and RescaleSlopeOriginal (' xASL_num2str(ASL_parms.RescaleSlopeOriginal) ') differ, using RWVSlope']);
			end
			if isfield(ASL_parms,'RescaleSlope') && ~isnear(ASL_parms.RWVSlope,ASL_parms.RescaleSlope,ASL_parms.RWVSlope/100) && (Rescale_NIfTI ~= 1)
				fprintf('%s\n', ['RWVSlope (' xASL_num2str(ASL_parms.RWVSlope) ') and RescaleSlope (' xASL_num2str(ASL_parms.RescaleSlope) ') differ, using RWVSlope']);
			end
			
			if ASL_parms.RWVSlope == 1
				warning('RWVSlope was 1, could be a scale slope issue');
			end
			
			% Set the scaling to remove
			bDoPhilipsScaling = 1;
			sPhilipsScaleSlope = ASL_parms.RWVSlope;
		elseif isfield(ASL_parms,'UsePhilipsFloatNotDisplayScaling') && (ASL_parms.UsePhilipsFloatNotDisplayScaling == 1)
			% Without the RWVSlope set and with the UsePhilipsFloatNotDisplayScaling set to 1, we know that dcm2nii did the scaling correctly and we don't have to further scale anything
			bDoPhilipsScaling = 0;
			sPhilipsScaleSlope = 0;
		elseif (~isfield(ASL_parms,'RescaleSlope')) && (~isfield(ASL_parms,'MRScaleSlope')) && (~isfield(ASL_parms,'RescaleSlopeOriginal')) && (~isfield(ASL_parms,'UsePhilipsFloatNotDisplayScaling'))
			% All fields are missing - we can ignore the quantification as all the fields have been removed and scaling done properly (we assume)
			% Note that RWVSlope was ruled out previously, and if UsePhilipsFloatNotDisplayScaling exists, it must be 0
			bDoPhilipsScaling = 0;
			sPhilipsScaleSlope = 0;
		else 
			% Standard scaling using RescaleSlope/RescaleSlopeOriginal or the Nifti slope - which should all be either non-existing or 1 or all equal
			% Set the scaling to remove
			bDoPhilipsScaling = 1;
			sPhilipsScaleSlope = 1;
			
			if isfield(ASL_parms,'RescaleSlopeOriginal') && (ASL_parms.RescaleSlopeOriginal ~= 1)
				if (sPhilipsScaleSlope ~= 1) && ~isnear(sPhilipsScaleSlope,ASL_parms.RescaleSlopeOriginal,sPhilipsScaleSlope/100)
					warning('%s\n', ['Discrepancy in RescaleSlopeOriginal (' xASL_num2str(ASL_parms.RescaleSlopeOriginal) ')']);
				end
				sPhilipsScaleSlope = ASL_parms.RescaleSlopeOriginal;
			end
			
			if isfield(ASL_parms,'RescaleSlope') && (ASL_parms.RescaleSlope ~= 1)
				if (sPhilipsScaleSlope ~= 1) && ~isnear(sPhilipsScaleSlope,ASL_parms.RescaleSlope,sPhilipsScaleSlope/100)
					warning('%s\n', ['Discrepancy in RescaleSlope (' xASL_num2str(ASL_parms.RescaleSlope) ')']);
				end
				sPhilipsScaleSlope = ASL_parms.RescaleSlope;
			end
			
			if (Rescale_NIfTI ~= 1)
				if (sPhilipsScaleSlope ~= 1) && ~isnear(sPhilipsScaleSlope,Rescale_NIfTI,sPhilipsScaleSlope/100)
					warning('%s\n', ['Discrepancy in NIFTI Slopes (' xASL_num2str(Rescale_NIfTI) ')']);
				end
				sPhilipsScaleSlope = Rescale_NIfTI;
			end

			if sPhilipsScaleSlope == 1
				warning('Scale slope was 1, could be a scale slope issue');
			end
		end
		
		% Apply the correct scaling
		if bDoPhilipsScaling
			if ASL_parms.MRScaleSlope == 1
				warning('Philips MR ScaleSlope was 1, could be a scale slope issue.');
			end
			ScaleImage = ScaleImage./(sPhilipsScaleSlope.*ASL_parms.MRScaleSlope);
			fprintf('%s',['Using DICOM (re)scale slopes ' num2str(sPhilipsScaleSlope) ' * ' num2str(ASL_parms.MRScaleSlope)]);
		end

		% Siemens specific scalings
	elseif strcmp(x.Vendor,'Siemens')
		if ~strcmp(x.Vendor,'Siemens_JJ_Wang') && strcmp(x.M0,'separate_scan')
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
    switch x.Q.LabelingType
        case 'PASL'
            fprintf('%s',['TI1 = ' num2str(x.Q.LabelingDuration) ' ms, ']);
            fprintf('%s',['TI (ms) = ' num2str(x.Q.Initial_PLD)]);
        case 'CASL'
            fprintf('%s',['LabelingDuration = ' num2str(x.Q.LabelingDuration) ' ms, ']);
            fprintf('%s',['PLD (ms) = ' num2str(x.Q.Initial_PLD)]);
    end

    if isfield(x.Q,'SliceReadoutTime')
        if x.Q.SliceReadoutTime>0 && strcmp(x.readout_dim,'2D')
            fprintf('%s',[' + ' num2str(x.Q.SliceReadoutTime) ' ms*(slice-1)']);
        end
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