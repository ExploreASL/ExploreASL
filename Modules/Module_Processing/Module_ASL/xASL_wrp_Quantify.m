function xASL_wrp_Quantify(x, PWI4D_Path, pathOutputCBF, M0Path, SliceGradientPath, bSaveCBF4D)
%xASL_wrp_Quantify Submodule of ExploreASL ASL Module, that performs quantfication
%
% FORMAT: xASL_wrp_Quantify(x [, PWI4D_Path, pathOutputCBF, M0Path, SliceGradientPath, bSaveCBF4D])
%
% INPUT:
%   x                   - structure containing fields with all information required to run this submodule (REQUIRED)
%   PWI4D_Path          - path to NIfTI with perfusion-weighted image (PWI4D) (OPTIONAL, default = x.P.Pop_Path_PWI4D)
%   pathOutputCBF       - path to NifTI to create, with the quantified CBF map (OPTIONAL, DEFAULT = x.P.Pop_Path_qCBF)
%   M0Path              - path to NifTI containing M0 image (OPTIONAL, default = x.Pop_Path_M0)
%   SliceGradientPath   - path to Slice gradient NIfTI (OPTIONAL, default = x.P.Pop_Path_SliceGradient_extrapolated)
%   bSaveCBF4D          - Boolean to save CBF quantified in 4D (OPTIONAL, default = false)
%
% OUTPUT: n/a
% OUTPUT FILES: NIfTI containing quantified CBF map in native or standard space (depending on input NIfTI),
% or other derivatives that need a quantification, e.g. FEAST
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This submodule converts PWIs to quantified CBF maps (or
% related derivatives). Note that we don't delete x.P.Path_PWI4D here, as
% this NIfTI file may be needed by xASL_wrp_VisualQC_ASL.m
%
%           0. Admin
%           1. Load PWI
%           2. Prepare M0
%           3. Hematocrit & blood T1 correction
%           4. ASL & M0 parameters comparisons
%           5. Load SliceGradient
%           6. Initialize & define quantification parameters
%           7. Define labeling efficiency
%           8. Perform quantification
%           9. Save files
%          10. Perform FEAST quantification (if exist)
%          11. Create standard space masked image to visualize masking effect
%
% EXAMPLE: xASL_wrp_Quantify(x);
%
% REFERENCES: 
%     Ivanov D, Gardumi A, Haast RAM, Pfeuffer J, Poser BA, Uludag K. 
%     Comparison of 3 T and 7 T ASL techniques for concurrent functional perfusion and BOLD studies
%     Neuroimage. 2017; 156:363-376.
% __________________________________
% Copyright (C) 2015-2024 ExploreASL

%% ------------------------------------------------------------------------------------------------
%% 0.   Administration

% by default, we use standard space NIfTIs
if nargin<2 || isempty(PWI4D_Path)
    PWI4D_Path = x.P.Pop_Path_PWI4D;
end

if nargin<3 || isempty(pathOutputCBF)
    pathOutputCBF = x.P.Pop_Path_qCBF;
end
% Define output path for ATT map based on the CBF output path
% Replace CBF with ATT in the output path
iStringCBF = regexpi(pathOutputCBF, 'CBF');
iStringCBF = iStringCBF(end);
pathOutputATT = [pathOutputCBF(1:(iStringCBF-1)) 'ATT' pathOutputCBF((iStringCBF+3):end)];
pathOutputTex = [pathOutputCBF(1:(iStringCBF-1)) 'Tex' pathOutputCBF((iStringCBF+3):end)];
pathOutputABV = [pathOutputCBF(1:(iStringCBF-1)) 'ABV' pathOutputCBF((iStringCBF+3):end)];

if nargin<4 || isempty(M0Path)
    M0Path = x.P.Pop_Path_M0;
end
if nargin<5 || isempty(SliceGradientPath)
    SliceGradientPath = x.P.Pop_Path_SliceGradient_extrapolated;
end

if ~isfield(x.modules.asl,'bUseBasilQuantification') || isempty(x.modules.asl.bUseBasilQuantification)
   x.modules.asl.bUseBasilQuantification = false;
end

if nargin<6 || isempty(bSaveCBF4D)
	bSaveCBF4D = false;
end

if x.modules.asl.bUseBasilQuantification && bSaveCBF4D
	warning('Cannot save CBF4D when using BASIL');
end

if ~isfield(x,'Q')
    x.Q = struct;
end

% Define quantification parameters
x = xASL_quant_DefineQuantificationParameters(x);

% Check if PWI4D exists
if ~xASL_exist(PWI4D_Path, 'file')
    warning('Skipped xASL_wrp_Quantify: files missing, please rerun step 4: xASL_wrp_ResampleASL');

    fprintf('%s\n', ['Missing: ' PWI4D_Path]);
    return;
end

% 0b. Remove output paths if they exist, to avoid confusion between different ExploreASL repetitions

xASL_delete(pathOutputCBF);
xASL_delete(pathOutputATT);
xASL_delete(pathOutputTex);
xASL_delete(pathOutputABV);

% For BASIL, only native data are processed and standard space data are not directly quantified, but only transformed
% So that's why we need to delete both the native and standard space data at once for BASIL
if x.modules.asl.bUseBasilQuantification
    xASL_delete(x.P.Pop_Path_qCBF);
    xASL_delete(x.P.Pop_Path_ATT);
    xASL_delete(x.P.Pop_Path_Tex);
    xASL_delete(x.P.Pop_Path_ABV);
end

%% ------------------------------------------------------------------------------------------------
%% 1.   Load PWI4D
fprintf('%s\n','Loading PWI4D & M0 images');

% Load ASL single PWI
[~, jsonASL] = xASL_io_Nifti2Im(PWI4D_Path); % Load CBF nifti
jsonASL = xASL_bids_parms2BIDS([], jsonASL, 0); % BIDS to Legacy conversion

if isempty(jsonASL)
    % If jsonASL is missing, we have to throw an error. Because all x.Q.PLD are for the raw unsubtracted data.
	% We could have recalculated the correct PLD vectors, but that should not really be done here but in the 
	% previous functions.
	error(['JSON file missing for: ' PWI4D_Path]);
end

% Assign the shortest minimal positive TE for ASL
ASLshortestTE = min(jsonASL.Q.EchoTime(jsonASL.Q.EchoTime > 0));

%% ------------------------------------------------------------------------------------------------
%% 2.   Prepare M0
if x.modules.asl.ApplyQuantification(5)==0
    % M0 division disabled, so we use a dummy M0 value only
    M0_im = NaN;

elseif isnumeric(x.Q.M0)
        % Single value per scanner
        % In this case we assume that this nifti value has been properly acquired,
        % does not need any corrections, and whole ASL PWI will be divided by this single value.

        M0_im = x.Q.M0;
        fprintf('%s\n',['Single M0 value ' num2str(M0_im) ' used']);

        if x.modules.asl.ApplyQuantification(4)
            % in case of separate M0, or M0 because of no background suppression,
            % T2* effect is similar in both images and hence removed by division
			if ~isempty(ASLshortestTE)
				T2_star_factor = exp(ASLshortestTE/x.Q.T2star);
				M0_im = M0_im./T2_star_factor;
				fprintf('%s\n',['M0 image corrected for T2* decay during TE in PWI, TE was ' xASL_num2str(ASLshortestTE) ' ms, using T2* ' xASL_num2str(x.Q.T2star) ' ms, this resulting in factor ' xASL_num2str(T2_star_factor)]);
				% If obtained by e.g. CSF inversion recovery, make sure that this is corrected for blood-water partition coefficient (0.76) and density of brain tissue (1.05 g/mL)
			else
				error('EchoTime unknown for ASL, but it is needed for the quantification to correct for ASL vs M0 signal differences.');
			end
        end

else
    % In this case we have a proper M0 image that should have voxel-wise
    % orientation agreement with the ASL PWI. We may need to correct it"s values.

    % Load M0 data

    fprintf('%s\n','M0 scan used');
    M0_im = xASL_io_Nifti2Im(M0Path);

    if xASL_stat_SumNan(M0_im(:))==0
        warning(['Empty M0:' M0Path]);
    end

    % NB: M0 quantification has largely been done in previous script,
    % load M0-parms only to check that ASL & M0-parms are identical
    % Scale slopes & incomplete T1 relaxation were already corrected in M0 module
end

%% ------------------------------------------------------------------------------------------------
%% 3.   Hematocrit & blood T1 correction
% Here, we check if the user has provided a hematocrit value for this
% subject_session_run. Only then, we create a x.Q.BloodT1.
% Below, at the quantification section, this is only taken into account
% when x.Q.BloodT1 exists, otherwise default Blood T1 values are used based
% on MagneticFieldStrength.
if isfield(x, 'hematocrit') && isfield(x, 'Hematocrit') % eventually this should be x.Q.Hematocrit, or x.S.SetsID from participants.tsv
    warning('Two hematocrit fields, not sure which one to use, please provide one only!');
elseif isfield(x, 'hematocrit')
    x.Hematocrit = x.hematocrit;
	x = rmfield(x, 'hematocrit');
end
if isfield(x,'Hematocrit')
    x.Q.BloodT1 = xASL_quant_Hct2BloodT1(x.Hematocrit, [], x.MagneticFieldStrength);
end

%% ------------------------------------------------------------------------------------------------
%% 4)   ASL & M0 parameters comparisons (e.g. TE, these should be the same with a separate M0 scan, for similar T2 & T2*-related quantification effects, and for similar geometric distortion)
if strcmpi(x.Q.M0,'separate_scan')
	[~, jsonM0] = xASL_io_Nifti2Im(x.P.Path_M0);
	jsonM0 = xASL_bids_parms2BIDS([], jsonM0, 0);
    %M0_parms = xASL_adm_LoadParms(x.P.Path_M0_parms_mat, x);
	
	% Assigns the shortest minimal positive TE for M0
	M0shortestTE = [];
	if isfield(jsonM0.Q, 'EchoTime')
		M0shortestTE = min(jsonM0.Q.EchoTime(jsonM0.Q.EchoTime > 0));
    else
        warning('Missing M0 Echo Time');
	end

    % Check echo times
    if  ~isempty(ASLshortestTE) && ~isempty(M0shortestTE)
        
		% Check equality of TE, but allow them to be 1% different, % Throw error if TE of ASL and M0 are not exactly the same!
		if max(~isnear(ASLshortestTE, M0shortestTE, 0.05*abs(ASLshortestTE)))
			% Here we allow for a 5% difference in TE, before giving the warning, which equals to 0.75 ms on 14 ms
			warning('TE of ASL and M0 are unequal. Check geometric distortion...');
		end
		
        % Correction factor and name for 3D spiral sequences
        if strcmpi(x.Q.Sequence,'3D_spiral')
			CorrFactor = x.Q.T2;
			CorrName = 'T2';
        else % assume T2* signal decay 2D_EPI or 3D GRASE
			CorrFactor = x.Q.T2star;
			CorrName = 'T2star';
        end

        % Correct M0 for any EchoTime differences between ASL & M0
        if x.modules.asl.ApplyQuantification(4)
			% If the shortest TEs are unequal, then we have to compensate for this
			if isempty(ASLshortestTE)
				error('EchoTime unknown for ASL, but it is needed for the quantification to correct for ASL vs M0 signal differences.');
			end
			if isempty(M0shortestTE)
				error('EchoTime unknown for M0, but it is needed for the quantification to correct for ASL vs M0 signal differences.');
			end
			if ASLshortestTE ~= M0shortestTE
				ScalingASL = exp(ASLshortestTE/CorrFactor);
				ScalingM0 = exp(M0shortestTE/CorrFactor);

				M0_im = M0_im.*ScalingM0./ScalingASL;
				fprintf('Delta TE between ASL %s ms & M0 %s ms, for %s, assuming %s decay of arterial blood, factor applied to M0: %s\n', ...
					num2str(ASLshortestTE),num2str(M0shortestTE),...
					x.Q.Sequence, CorrName, num2str(ScalingM0/ScalingASL));
			end
        end
        
    else
        warning('Could not compare TEs from ASL & M0, JSON fields missing...');
    end
end



if ~x.modules.asl.ApplyQuantification(3) % if conversion PWI for label units is not requested
    SliceGradient = [];
else

    %% ------------------------------------------------------------------------------------------------
    %% 5    Load SliceGradient
    if  strcmpi(x.Q.readoutDim,'2D')
        SliceGradient = xASL_io_Nifti2Im(SliceGradientPath);
    else
        SliceGradient = [];
    end


    %% ------------------------------------------------------------------------------------------------
    %% 6.   Initialize quantification parameters
	if ~isfield(x.Q,'nCompartments') || isempty(x.Q.nCompartments)
        x.Q.nCompartments = 1; % by default, we use a single-compartment model, as proposed by the Alsop et al. MRM 2014 concensus paper
    elseif x.Q.nCompartments>2 || x.Q.nCompartments<1
		warning(['Unknown x.Q.nCompartments: ' xASL_numstr2(x.Q.nCompartments)]);
        fprintf('%s\n', 'Now x.Q.nCompartments set to 1 (single compartment model)');
        x.Q.nCompartments = 1;
	end

	if ~isfield(x.Q,'ATT')
        x.Q.ATT = 1800; % ms as default micro-vascular ATT
	end



    % Check correct order of magnitude blood T1 (this value should be around 1700, or ~ 1000-3000)
    if x.Q.BloodT1<10
        x.Q.BloodT1 = x.Q.BloodT1.*1000;
    end

	if ~isfield(x.Q,'LabelingType')
           error('Unknown LabelingType, needed for quantification');
    elseif isempty(regexpi(x.Q.LabelingType, '^(PC|P|C)ASL$'))
           error('x.Q.LabelingType was invalid, should be PASL|CASL|PCASL');
    elseif strcmpi(x.Q.LabelingType,'PCASL')
           x.Q.LabelingType = 'CASL';
	end

	if isfield(x.Q,'BackgroundSuppression') && ~x.Q.BackgroundSuppression
		% In case BSup is defined and turned off, then we set the number of pulses to 0
		x.Q.BackgroundSuppressionNumberPulses = 0;
		% otherwise, Bsup is either not defined or BSup == true
	elseif ~isfield(x.Q,'BackgroundSuppressionNumberPulses') || isempty(x.Q.BackgroundSuppressionNumberPulses)
		% number of pulses or presence of BSup is not specified, then assign defaults
        if strcmpi(x.Q.readoutDim, '3d') && strcmpi(x.Q.Vendor, 'ge')
            warning('Unknown number of background suppression pulses, assuming 5 pulses for this GE 3D sequence (including the pre-pulse)');
			x.Q.BackgroundSuppression = true; % We can set to TRUE since we know that this is either empty or TRUE already
            x.Q.BackgroundSuppressionNumberPulses = 5;
        elseif strcmpi(x.Q.readoutDim, '3d') && strcmpi(x.Q.Vendor, 'philips')
            warning('Unknown number of background suppression pulses, assuming 4 pulses for this Philips 3D sequence');
			x.Q.BackgroundSuppression = true; % We can set to TRUE since we know that this is either empty or TRUE already
            x.Q.BackgroundSuppressionNumberPulses = 4;
        elseif strcmpi(x.Q.readoutDim, '3d') && strcmpi(x.Q.Vendor, 'siemens')
            warning('Unknown number of background suppression pulses, assuming 4 pulses for this Siemens 3D sequence');
			x.Q.BackgroundSuppression = true; % We can set to TRUE since we know that this is either empty or TRUE already
            x.Q.BackgroundSuppressionNumberPulses = 4;
			% For 3D cases, the defaults are rather consistent for product sequences, for 2D below, we are not so sure
		else
			if isfield(x.Q,'BackgroundSuppression')
				% So if BSup was defined, we set the default to 2
				warning('Unknown number of background suppression pulses for a 2D, assuming 2 pulses');
				x.Q.BackgroundSuppressionNumberPulses = 2;
			else
				% otherwise we assume noBsup
				warning('For 2D sequences, the Bsup is turned off by default, assuming 0 pulses');
				x.Q.BackgroundSuppression = false;
				x.Q.BackgroundSuppressionNumberPulses = 0;
			end
        end
	end


    %% 7.   Labeling efficiency
    if ~isfield(x.Q,'LabelingEfficiency') || isempty(x.Q.LabelingEfficiency)
		switch lower(x.Q.LabelingType)
			case 'pasl'
				if x.MagneticFieldStrength == 7
					x.Q.LabelingEfficiency = 0.95; % FAIR at 7T, Ivanov, Neuroimage, 2017
				else
					x.Q.LabelingEfficiency = 0.98; % (concensus paper, Wong, MRM 1998)
				end
			case 'casl'
				x.Q.LabelingEfficiency = 0.85; % (concensus paper, Dai, MRM 2008)
		end
    end
    x.Q.LabEff_Bsup = 1; % default for no background suppression
    % Apply the effect of background suppression pulses on labeling efficiency
    switch x.Q.BackgroundSuppressionNumberPulses
        case 0 % when you have an M0, but no background suppression used for ASL
            % Then the labeling efficiency doesn't change by background suppression
        case 2 % e.g. Philips 2D EPI or Siemens 3D GRASE
            x.Q.LabEff_Bsup = 0.83; % 0.83 = 2 background suppression pulses (Garcia et al., MRM 2005)
        case 4 % e.g. Philips 3D GRASE
            x.Q.LabEff_Bsup = 0.81; % 0.81 = as implemented by Philips
        case 5 % e.g. GE 3D spiral
            x.Q.LabEff_Bsup = 0.75; % 0.75 = 5 background suppression pulses (GE FSE) (Garcia et al., MRM 2005)
    end

    x.Q.LabEff_Orig = x.Q.LabelingEfficiency;
    x.Q.LabelingEfficiency = x.Q.LabelingEfficiency*x.Q.LabEff_Bsup;
end

%% ------------------------------------------------------------------------------------------------
%% 8.   Perform Quantification
if ~x.modules.asl.bQuantifyMultiPLD || x.modules.asl.bUseBasilQuantification % multi-PLD with BASIL or single-PLD
    [~, CBF, ATT, ABV, Tex] = xASL_quant_ASL(PWI4D_Path, M0_im, SliceGradient, x, x.modules.asl.bUseBasilQuantification, bSaveCBF4D); % also runs BASIL, but only in native space!
else
    % multi-PLD quantification without BASIL
    error('Multi PLD quantification without BASIL is not yet implemented.');
end

if x.modules.asl.ApplyQuantification(5)==0
    MeanCBF = xASL_stat_MeanNan(CBF(:));
    if MeanCBF>666 % this is the average including air
        CBF = CBF .* (10./MeanCBF); % protection against extremely high values
        warning('M0 division was disabled & CBF image had too high values');
        fprintf('%s\n',['mean whole image CBF normalized from ' xASL_num2str(MeanCBF) ' to 10 mL/100g/min']);
    end
end

%% ------------------------------------------------------------------------------------------------
%% 9.	Save files
% Both ExploreASL and BASIL-quantified maps will be saved similarly here
fprintf('%s\n','Saving PWI & CBF niftis');

xASL_io_SaveNifti(PWI4D_Path, pathOutputCBF, CBF, 32, 0);

if numel(ATT) > 1
	% Save the ATT file
	xASL_io_SaveNifti(PWI4D_Path, pathOutputATT, ATT, 32, 0);
end

if numel(ABV) > 1 
	% Save the ABV file
	xASL_io_SaveNifti(PWI4D_Path, pathOutputABV, ABV, 32, 0);
end

if numel(Tex) > 1
	% Save the Tex file
	xASL_io_SaveNifti(PWI4D_Path, pathOutputTex, Tex, 32, 0);
end

%% 9.b Save files in standard space for BASIL native space output
% Transform BASIL CBF to standard space as BASIL only quantifies in native space
if x.modules.asl.bUseBasilQuantification && strcmp(x.P.Path_CBF, pathOutputCBF)
    if exist(x.P.Path_mean_PWI_Clipped_sn_mat, 'file') % Backwards compatability, and also needed for the Affine+DCT co-registration of ASL-T1w
        AffineTransfPath = x.P.Path_mean_PWI_Clipped_sn_mat;
    else
        AffineTransfPath = [];
    end
    
    xASL_spm_deformations(x, {x.P.Path_CBF}, {x.P.Pop_Path_qCBF}, [], [], AffineTransfPath, x.P.Path_y_ASL);
	if xASL_exist(x.P.Path_ATT,'file')
		xASL_spm_deformations(x, {x.P.Path_ATT}, {x.P.Pop_Path_ATT}, [], [], AffineTransfPath, x.P.Path_y_ASL);
	end
    
	if xASL_exist(x.P.Path_Tex,'file')
		xASL_spm_deformations(x, {x.P.Path_Tex}, {x.P.Pop_Path_Tex}, [], [], AffineTransfPath, x.P.Path_y_ASL);
	end

	if xASL_exist(x.P.Path_ABV,'file')
		xASL_spm_deformations(x, {x.P.Path_ABV}, {x.P.Pop_Path_ABV}, [], [], AffineTransfPath, x.P.Path_y_ASL);
	end
end


%% ------------------------------------------------------------------------------------------------
%% 10.   FEAST quantification
% run FEAST quantification if crushed & non-crushed ASL sessions exist
if (x.dataset.nSessions>1 && isfield(x,'session') && isfield(x.session,'options') && strcmp(x.session.options{1},'crushed') && strcmp(x.session.options{2},'non-crushed'))
    % run FEASTS quantification if current session =session 2
    if strcmp(x.dir.SESSIONDIR(length(x.dir.SUBJECTDIR)+2:end),'ASL_2')
        xASL_quant_FEAST(x);
    end
end


%% ------------------------------------------------------------------------------------------------
%% 11.  Create standard space masked image to visualize masking effect
if xASL_exist(x.P.Pop_Path_qCBF, 'file') && (strcmp(pathOutputCBF, x.P.Pop_Path_qCBF) || x.modules.asl.bUseBasilQuantification) && ~bSaveCBF4D
    % Load CBF image
    MaskedCBF = xASL_io_Nifti2Im(x.P.Pop_Path_qCBF);
    % Mask vascular voxels (i.e. set them to NaN)
    MaskVascularMNI = xASL_io_Nifti2Im(x.P.Pop_Path_MaskVascular);
    MaskedCBF(~MaskVascularMNI) = NaN;
    
    % Mask susceptibility voxels (i.e. set them to NaN)
    if ~xASL_exist(x.P.Pop_Path_MaskSusceptibility,'file')
        warning([x.P.Pop_Path_MaskSusceptibility ' missing, cannot create ' x.P.Pop_Path_qCBF_masked]);
        fprintf('Please rerun xASL_wrp_CreateIndividualMask (7th step of ASL module)...\n');
        return;
    else
        xASL_io_SaveNifti(x.P.Pop_Path_qCBF, x.P.Pop_Path_qCBF_masked, MaskedCBF, [], false);
    end
end

end

