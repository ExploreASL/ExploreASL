function [ScaleImage, CBF] = xASL_quant_Basil(PWI, M0_im, SliceGradient, x)
%xASL_quant_Basil Perform quantification using FSL BASIL
% FORMAT: [ScaleImage[, CBF]] = xASL_quant_Basil(PWI, M0_im, SliceGradient, x)
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
% DESCRIPTION: This script performs quantification of the PWI using the FSL Basil pipeline. Final calibration to
%              physiological units is performed by dividing the quantified PWI by the M0 image/value
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE: [ScaleImage, CBF] = xASL_quant_Basil(PWI, M0_im, SliceGradient, x);
% __________________________________
% Copyright 2015-2020 ExploreASL

fprintf('%s\n','Quantification CBF using FSL Basil:');

%
% For now we think ExploreASL is only handling single PLD data
%
SinglePLD = 1;

if  xASL_stat_SumNan(M0_im(:))==0
    error('Empty M0 image, something went wrong in M0 processing');
end

% Convert to double precision to increase quantification precision
% This is especially useful with the large factors that we can multiply and
% divide with in ASL quantification
PWI = double(PWI);
M0_im = double(M0_im);

if ~x.ApplyQuantification(3)
    fprintf('%s\n','We skip the scaling of a.u. to label intensity');
else
    %
    % For a 2D readout, determine the slice time gradient
    %
    switch x.readout_dim
        case '3D'
            fprintf('%s\n','3D sequence, not accounting for SliceReadoutTime (homogeneous PLD for complete volume)');
            x.Q.SliceReadoutTime = 0.0

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
        otherwise
            error('Wrong x.readout_dim value!');
    end

    %
    % Write the PWI out as  Nifti file for Basil to read as input
    % FIXME would be good to have a brain mask too at this point
    %
    xASL_io_CreateNifti('basil_input_data.nii', PWI)

    %
    % option_file contains options which are passed to Fabber
    % basil_options is a character array containing CL args for the Basil command
    %
    option_file = fopen('model_options.txt', 'w');
    basil_options = '';

    % Basic acquisition and tissue parameters
    fprintf(option_file, '# Basil options written by ExploreASL\n');
    fprintf(option_file, '--repeats=%i\n', size(PWI, 4));
    fprintf(option_file, '--ti1=%f\n', (x.Q.LabelingDuration + x.Q.Initial_PLD)/1000);
    fprintf(option_file, '--t1b=%f\n', x.Q.BloodT1/1000);
    fprintf(option_file, '--tau=%f\n', x.Q.LabelingDuration/1000); % FIXME tau could be list
    fprintf(option_file, '--slicedt=%f\n', x.Q.SliceReadoutTime/1000);
    fprintf(option_file, '--save-model-fit\n');

    %
    % FIXME Aquisition options we might be able to use in the future
    %
    %fprintf(option_file, '--FA=%f\n', fa);
    %fprintf('Basil: Flip angle for look-locker readout: %f\n', fa);
    %fprintf(option_file, '--sliceband=%i\n', sliceband);
    %fprintf('Basil: Multi-band setup with number of slices per band: %i\n', slicedband);

    % This helps avoid failure on the structural-space image
    fprintf(option_file, '--allow-bad-voxels\n');

    % FIXME is a user-specified T1 map possible in ExploreASL?
    %fprintf(option_file, '--t1im=%s\n', t1im)
    %fprintf('Basil: Using supplied T1 (tissue) image in BASIL: %s\n', $t1im)

    % Labelling type - PASL or pCASL
    switch x.Q.LabelingType
        case 'PASL'
            fprintf('Basil: pASL model\n');
        case 'CASL'
            fprintf(option_file, '--casl\n');
            fprintf('Basil: cASL/pcASL model\n');
    end

    %
    % Model option - 1 or 2 compartment
    % The 1-compartment model is what we call 'White paper mode' in oxford-asl.
    % This means zero ATT (all bolus delivered by imaging time) and blood T1 only
    % The 2-compartment model is the 'standard' buxton model which takes into account
    % the ATT and tissue T1 value.
    %
    switch x.Q.nCompartments
        case 1
            fprintf(option_file, '--bat=0\n');
            fprintf(option_file, '--t1=%f\n', x.Q.BloodT1/1000);
            fprintf('Basil: Single-compartment (white paper mode)\n');
        case 2
            fprintf(option_file, '--bat=%f\n', x.Q.ATT/1000);
            fprintf(option_file, '--t1=%f\n', x.Q.TissueT1/1000);
            fprintf('Basil: 2-compartment - ATT=%fs\n', x.Q.ATT/1000);
    end

    %
    %  Bolus duration typically fixed but can be inferred
    %
    if isfield(x.Q,'BasilInferTau') & x.Q.BasilInferTau
         fprintf(option_file, '--infertau\n');
         fprintf('Basil: Infer bolus duration component\n')
    else
        fprintf('Basil: Fixed bolus duration component\n')
    end

    %
    % ATT and arterial component inference only possible with multi-PLD
    %
    if SinglePLD
        fprintf('Basil: Single-delay data - cannot infer ATT or arterial component\n');
        x.Q.BasilInferATT = 0;
        x.Q.BasilInferArt = 0;
    end

    %
    % Infer arterial transit time
    %
    if isfield(x.Q,'BasilInferATT') & x.Q.BasilInferATT
        if ~isfield(x.Q, 'BasilATTSD')
            x.Q.BasilATTSD = 1.0;
        end
        fprintf(option_file, '--inferbat\n');
        fprintf(option_file, '--batsd=%f\n', x.Q.BasilATTSD);
	    fprintf('Basil: Setting std dev of the (tissue) BAT prior std.dev. to %f\n', x.Q.BasilATTSD);
    else
        basil_options = [basil_options ' --fixbat'];
        fprintf('Basil: Fixed arterial arrival time\n');
    end

    %
    % Infer arterial component
    %
    if isfield(x.Q,'BasilInferATT') & x.Q.BasilInferATT
        fprintf(option_file, '--inferart\n');
        fprintf('Basil: Infer arterial component');
        fprintf('Basil: Variable arterial component arrival time');
    end

    %
    % Noise specification. For small numbers of time points we need informative
    % noise prior. User can specify assumed SNR for this, or give noise std.dev
    % directly.
    %
    if ~isfield(x.Q,'BasilSNR') | ~x.Q.BasilSNR
        x.Q.BasilSNR = 10;
    end

    if size(PWI, 4) < 5
        x.Q.BasilNoisePrior = 1;
        fprintf('Basil: Small number of volumes (%i < 5): informative noise prior will be used\n', size(PWI, 4));
    end

    if isfield(x.Q,'BasilNoisePrior') & x.Q.BasilNoisePrior
        % Use an informative noise prior
        if ~isfield(x.Q,'BasilNoiseSD') | ~x.Q.BasilNoiseSD
            fprintf('Basil: Using SNR of %f to set noise std dev\n', x.Q.BasilSNR);
            % Estimate signal magntiude FIXME brain mask assume half of voxels
            mag_max = max(PWI, [], 4);
            brain_mag = 2*mean(mag_max, 'all');
            fprintf('Basil: Mean maximum signal across brain: %f\n', brain_mag);
            % This will correspond to whole brain CBF (roughly) - about 0.5 of GM
            x.Q.BasilNoiseSD = sqrt(brain_mag * 2 / x.Q.BasilSNR);
        end
        fprintf('Basil: Using a prior noise std.dev. of: %f\n', x.Q.BasilNoiseSD);
        fprintf(option_file, '--prior-noise-stddev=%f\n', x.Q.BasilNoiseSD);
    end

    %
    % Various optional features
    %

    if isfield(x,'BasilSpatial') & x.Q.BasilSpatial
        fprintf('Basil: Instructing BASIL to use automated spatial smoothing\n');
        basil_options = [basil_options ' --spatial'];
    end

    if isfield(x.Q,'BasilInferT1') & x.Q.BasilInferT1
        fprintf(option_file, '--infert1\n');
        fprintf('Basil: Instructing BASIL to infer variable T1 values\n');
    end

    if isfield(x.Q,'BasilExch')
        fprintf('Basil: Using exchange model: %s\n', x.Q.BasilExch);
        fprintf(option_file, '--exch=%s\n', x.Q.BasilExch);
    end

    if isfield(x.Q,'BasilDisp')
        fprintf('Basil: Using dispersion model: %s\n', x.Q.BasilDisp);
        fprintf(option_file, '--disp=%s\n', x.Q.BasilDisp);
    end

    if isfield(x.Q,'BasilDebug') & x.Q.BasilDebug
        basil_options = [basil_options ' --devel'];
    end
    fclose(option_file);

    %
    % Run Basil and retrieve CBF output
    % FIXME might not be step1 if using spatial mode
    %
    args.bAutomaticallyDetectFSL=1;
    xASL_fsl_RunFSL(['basil -i basil_input_data -@ model_options.txt -o basil_out' basil_options], args);

    ftiss = xASL_io_ReadNifti('basil_out/step1/mean_ftiss.nii');
    CBF_nocalib = ftiss.dat(:, :, :);

    %
    % Scaling to physiological units
    % Note different to xASL_quant_SinglePLD since Fabber has T1 in seconds
    % and does not take into account labelling efficiency
    %
    CBF_nocalib = CBF_nocalib .* 6000 .* x.Q.Lambda ./ x.Q.LabelingEfficiency;
    % (For some reason, GE sometimes doesn't need the 1 gr->100 gr conversion)
    % & old Siemens sequence also didn't need the 1 gr->100 gr conversion
end


%% 4    Vendor-specific scalefactor
if ~x.ApplyQuantification(1)
    fprintf('%s\n','We skip the vendor-specific scalefactors');
else
    % Throw warning if no Philips scans, but scale slopes are not 1:
    ASL_parms = xASL_adm_LoadParms(x.P.Path_ASL4D_parms_mat, x);
    if isempty(regexp(x.Vendor,'Philips')) && isfield(ASL_parms,'RescaleSlopeOriginal')
        if ASL_parms.RescaleSlopeOriginal~=1
            warning('We detected a RescaleSlope~=1, verify that this is not a Philips scan!!!');
        end
    elseif isempty(regexp(x.Vendor,'Philips')) && isfield(ASL_parms,'MRScaleSlope')
        if ASL_parms.MRScaleSlope~=1
            warning('We detected a ScaleSlope~=1, verify that this is not a Philips scan!!!');
        end
    end


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

        CBF_nocalib = CBF_nocalib./qnt_GEscaleFactor;
        fprintf('%s\n',['Quantification corrected for GE scale factor ' num2str(qnt_GEscaleFactor) ' for NSA=' num2str(x.Q.NumberOfAverages)]);

    elseif ~isempty(regexp(x.Vendor,'Philips'))
        % Philips has specific scale & rescale slopes
        % If these are not corrected for, only relative CBF quantification can be performed,
        % i.e. scaled to wholebrain, the wholebrain perfusion cannot be calculated.

        % QUICK & DIRTY FIX, NEEDS TO BE CHECKED
        % We assume that the original NIfTI rescale slope is correct
        % since the newer dcm2niiX version sometimes scales extra to use
        % the full 16 bit width.

        HeaderTemp = xASL_io_ReadNifti(x.P.Path_ASL4D); % original NIfTI
        if isfield(ASL_parms,'RescaleSlopeOriginal')
            % here we check whether these rescale slopes differ
            % significantly (i.e. more than 5%)
            Rescale_NIfTI = HeaderTemp.dat.scl_slope;
            Rescale_JSON = ASL_parms.RescaleSlopeOriginal;
            if Rescale_JSON<0.95*Rescale_NIfTI || Rescale_JSON>1.05*Rescale_NIfTI
                warning('Philips rescaleslopes differed >5% between JSON & NIfTI!');
                fprintf('%s\n', ['JSON RescaleSlope is ' xASL_num2str(Rescale_JSON)]);
                fprintf('%s\n', ['NIfTI RescaleSlope is ' xASL_num2str(Rescale_NIfTI)]);
                fprintf('%s\n', 'If NIfTI rescaleslope is 0.5<RescaleSlope<20 & non-one, we use this');
                if Rescale_NIfTI>0.5 && Rescale_NIfTI<20 && Rescale_NIfTI~=1
                    ASL_parms.RescaleSlopeOriginal = HeaderTemp.dat.scl_slope;
                    fprintf('%s\n', 'Using the NIfTI RescaleSlope');
                else
                    fprintf('%s\n', 'Using the JSON RescaleSlope');
                end
            end
        end

        if ASL_parms.RescaleSlopeOriginal==1
            warning('ASL_parms.RescaleSlopeOriginal was 1, could be a scale slope issue');
        end

        % Correct scale slopes
        CBF_nocalib = CBF_nocalib./(ASL_parms.RescaleSlopeOriginal.*ASL_parms.MRScaleSlope);
        % NifTI scale slope has already been corrected for by SPM's nifti function
        fprintf('%s',['Using DICOM (re)scale slopes ' num2str(ASL_parms.RescaleSlopeOriginal) ' * ' num2str(ASL_parms.MRScaleSlope)]);

    elseif strcmp(x.Vendor,'Siemens') && ~strcmp(x.Vendor,'Siemens_JJ_Wang') && strcmp(x.M0,'separate_scan')
          % Some Siemens readouts divide M0 by 10, others don't
          CBF_nocalib = CBF_nocalib./10;
          fprintf('%s\n','M0 corrected for Siemens scale factor 10')
    end
end

%% 5    Divide PWI/M0
% Match sizes
MatchSizeM0 = round([size(PWI,1)./size(M0_im,1) size(PWI,2)./size(M0_im,2) size(PWI,3)./size(M0_im,3) size(PWI,4)./size(M0_im,4) size(PWI,5)./size(M0_im,5) size(PWI,6)./size(M0_im,6) size(PWI,7)./size(M0_im,7)]);
MatchSizeSI = round([size(PWI,1)./size(CBF_nocalib,1) size(PWI,2)./size(CBF_nocalib,2) size(PWI,3)./size(CBF_nocalib,3) size(PWI,4)./size(CBF_nocalib,4) size(PWI,5)./size(CBF_nocalib,5) size(PWI,6)./size(CBF_nocalib,6) size(PWI,7)./size(CBF_nocalib,7)]);

if sum(MatchSizeM0==0) || sum(MatchSizeSI==0)
    error('PWI dimensions too small compared to M0 and/or CBF_nocalib dimensions');
end

M0_im = repmat(M0_im,MatchSizeM0);
CBF_nocalib = repmat(CBF_nocalib,MatchSizeSI);
    
if ~x.ApplyQuantification(5)
    fprintf('%s\n','We skip the PWI/M0 division');
else    
    CBF = CBF_nocalib./M0_im;
end

%
% Unsure if ScaleImage is needed elsewhere but calculate it anyway
%
ScaleImage = CBF ./ PWI;

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