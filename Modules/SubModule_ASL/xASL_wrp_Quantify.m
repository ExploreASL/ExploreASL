function xASL_wrp_Quantify(x, PWI_Path, OutputPath, M0Path, SliceGradientPath, bUseBasilQuantification)
%xASL_wrp_Quantify Submodule of ExploreASL ASL Module, that performs quantfication
%
% FORMAT: xASL_wrp_Quantify(x)
%
% INPUT:
%   x                   - structure containing fields with all information required to run this submodule (REQUIRED)
%   PWI_Path            - path to NIfTI with perfusion-weighted image (PWI) (OPTIONAL, default = x.P.Pop_Path_PWI)
%   OutputPath          - path to NifTI to create, with the quantified CBF map
%   M0Path              - path to NifTI containing M0 image (OPTIONAL, default = x.Pop_Path_M0)
%   SliceGradientPath   - path to Slice gradient NIfTI (OPTIONAL, default = x.P.Pop_Path_SliceGradient_extrapolated)
%   bUseBasilQuantification - boolean, true for using FSL BASIL for
%                             quantification, false for using ExploreASL's
%                             own quantification (OPTIONAL, DEFAULT = false)
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
% __________________________________
% Copyright (C) 2015-2020 ExploreASL


%% ------------------------------------------------------------------------------------------------
%% 0.   Administration

% Use either original or motion estimated ASL4D
% Use despiked ASL only if spikes were detected and new file has been created
% Otherwise, despiked_raw_asl = same as original file
if ~xASL_exist(x.P.Path_despiked_ASL4D,'file')
    x.P.Path_despiked_ASL4D = x.P.Path_ASL4D;
end

% by default, we use standard space NIfTIs
if nargin<2 || isempty(PWI_Path)
    PWI_Path = x.P.Pop_Path_PWI;
end
if nargin<3 || isempty(OutputPath)
    OutputPath = x.P.Pop_Path_qCBF;
end
if nargin<4 || isempty(M0Path)
    M0Path = x.P.Pop_Path_M0;
end
if nargin<5 || isempty(SliceGradientPath)
    SliceGradientPath = x.P.Pop_Path_SliceGradient_extrapolated;
end
if nargin<6 || isempty(bUseBasilQuantification)
    bUseBasilQuantification = false;
end
if ~isfield(x,'Q')
    x.Q = struct;
end
if ~isfield(x.Q,'T2star') || isempty(x.Q.T2star)
    x.Q.T2star = 48; % default for 3T
end
if ~isfield(x.Q,'T2') || isempty(x.Q.T2)
    x.Q.T2 = 180; % default for 3T (ref Jean Chen, MRM 2009)
end

if ~xASL_exist(PWI_Path, 'file')
    warning('Skipped xASL_wrp_Quantify: files missing, please rerun step 4: xASL_wrp_ResampleASL');
    fprintf('%s\n', ['Missing: ' PWI_Path]);
    return;
end

%% ------------------------------------------------------------------------------------------------
%% 1.   Load PWI
fprintf('%s\n','Loading PWI & M0 images');

% Load ASL PWI
PWI = xASL_io_Nifti2Im(PWI_Path); % Load CBF nifti
ASL_parms = xASL_adm_LoadParms(x.P.Path_ASL4D_parms_mat, x);

if xASL_stat_SumNan(PWI(:))==0
    warning(['Empty PWI image:' PWI_Path]);
end

% % % Philips dcm2niiX scaling fix:
% % % in the new dcm2niiX with the default Philips scaling,
% % % there can be a discrepancy in the rescale slope that was stored as
% % % DICOM rescale slope, and the NIfTI rescale slope applied by dcm2niiX
% % % to use the full 16 bit spectrum. This discrepancy is not normally
% % % taken into account by SPM's NIfTI function.
% % HeaderTemp = xASL_io_ReadNifti(x.P.Path_ASL4D);
% %
% % if HeaderTemp.dat.scl_slope~=1 && ASL_parms.RescaleSlope~=1
% %     ReRescaleFactor = ASL_parms.RescaleSlope./HeaderTemp.dat.scl_slope;
% %     PWI = PWI.*ReRescaleFactor;
% % end


%% ------------------------------------------------------------------------------------------------
%% 2.   Prepare M0
if isnumeric(x.M0)
        % Single value per scanner
        % In this case we assume that this nifti value has been properly acquired,
        % does not need any corrections, and whole ASL PWI will be divided by this single value.

        M0_im = x.M0;
        fprintf('%s\n',['Single M0 value ' num2str(M0_im) ' used']);

        if x.ApplyQuantification(4)
            % in case of separate M0, or M0 because of no background suppression,
            % T2* effect is similar in both images and hence removed by division
            T2_star_factor = exp(ASL_parms.EchoTime/x.Q.T2star);
            M0_im = M0_im./T2_star_factor;
            fprintf('%s\n',['M0 image corrected for T2* decay during TE in PWI, TE was ' num2str(ASL_parms.EchoTime) ' ms, using T2* ' num2str(x.Q.T2star) ' ms, this resulting in factor ' num2str(T2_star_factor)]);
            % If obtained by e.g. CSF inversion recovery, make sure that this is corrected for blood-water partition coefficient (0.76) and density of brain tissue (1.05 g/mL)
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

    M0_parms = xASL_adm_LoadParms(x.P.Path_M0_parms_mat, x);

% %     % Philips dcm2niiX scaling fix:
% %     % in the new dcm2niiX with the default Philips scaling,
% %     % there can be a discrepancy in the rescale slope that was stored as
% %     % DICOM rescale slope, and the NIfTI rescale slope applied by dcm2niiX
% %     % to use the full 16 bit spectrum. This discrepancy is not normally
% %     % taken into account by SPM's NIfTI function.
% %     HeaderTemp = xASL_io_ReadNifti(x.P.Path_M0);
% %
% %     if HeaderTemp.dat.scl_slope~=1 && M0_parms.RescaleSlope~=1
% %         ReRescaleFactor = M0_parms.RescaleSlope./HeaderTemp.dat.scl_slope;
% %         M0_im = M0_im.*ReRescaleFactor;
% %     end



    % NB: M0 quantification has largely been done in previous script,
    % load M0-parms only to check that ASL & M0-parms are identical
    % Scale slopes & incomplete T1 relaxation were already corrected in M0 module
end


%% ------------------------------------------------------------------------------------------------
%% 3.   Hematocrit & blood T1 correction
if isfield(x,'Hematocrit')
        x.Q.BloodT1 = xASL_quant_Hct2BloodT1(x.Hematocrit);
elseif isfield(x,'hematocrit')
        x.Q.BloodT1 = xASL_quant_Hct2BloodT1(x.hematocrit);
end




%% ------------------------------------------------------------------------------------------------
%% 4)   ASL & M0 parameters comparisons (e.g. TE, these should be the same with a separate M0 scan, for similar T2 & T2*-related quantification effects, and for similar geometric distortion)
if strcmpi(x.M0,'separate_scan')
    if  isfield(ASL_parms,'EchoTime') && isfield(M0_parms,'EchoTime')
        % Check equality of TE, but allow them to be 1% different, % Throw error if TE of ASL and M0 are not exactly the same!
        if  ASL_parms.EchoTime<(M0_parms.EchoTime*0.95) || ASL_parms.EchoTime>(M0_parms.EchoTime*1.05)
            % Here we allow for a 5% difference in TE, before giving the warning, which equals to 0.75 ms on 14 ms
            warning('TE of ASL and M0 were unequal! Check geometric distortion');
        end

        if strcmpi(x.Sequence,'3D_spiral')
           CorrFactor = x.Q.T2;
           CorrName = 'T2';
        else % assume T2* signal decay 2D_EPI or 3D GRASE
            CorrFactor = x.Q.T2star;
            CorrName = 'T2star';
        end

        % Correct M0 for any EchoTime differences between ASL & M0
        if x.ApplyQuantification(4)
            ScalingASL = exp(ASL_parms.EchoTime/CorrFactor);
            ScalingM0 = exp(M0_parms.EchoTime/CorrFactor);
            M0_im = M0_im.*ScalingM0./ScalingASL;
            fprintf('%s\n', ['Delta TE between ASL ' num2str(ASL_parms.EchoTime) 'ms & M0 ' num2str(M0_parms.EchoTime) 'ms, for ' x.Sequence ', assuming ' CorrName ' decay of arterial blood, factor applied to M0: ' num2str(ScalingM0/ScalingASL)]);
        end
    else
        warning('Could not compare TEs from ASL & M0, JSON fields missing!');
    end
end



if ~x.ApplyQuantification(3) % if conversion PWI for label units is not requested
    SliceGradient = [];
else

    %% ------------------------------------------------------------------------------------------------
    %% 5    Load SliceGradient
    if  strcmpi(x.readout_dim,'2D')
        SliceGradient = xASL_io_Nifti2Im(SliceGradientPath);
    else
        SliceGradient = [];
    end


    %% ------------------------------------------------------------------------------------------------
    %% 6.   Initialize quantification parameters
    if ~isfield(x.Q,'nCompartments') || isempty(x.Q.nCompartments) || x.Q.nCompartments>2
        x.Q.nCompartments = 1; % by default, we use a single-compartment model, as proposed by the Alsop et al. MRM 2014 concensus paper
    end

    if ~isfield(x.Q,'ATT')
        x.Q.ATT = 1800; % ms as default micro-vascular ATT
    end
    if ~isfield(x.Q,'TissueT1')
        x.Q.TissueT1 = 1240; % T1 GM tissue @ 3T
    end
    if ~isfield(x.Q,'BloodT1')
        x.Q.BloodT1  = 1650; % T1 relaxation time of arterial blood DEFAULT=1650 @ 3T
    end
    if ~isfield(x.Q,'T2art')
        x.Q.T2art = 50; % T2* of arterial blood, only used when no M0 image
    end
    if ~isfield(x.Q,'Lambda')
        x.Q.Lambda = 0.9; % Brain/blood water coefficient (mL 1H/ mL blood)
    end
    % Check correct order of magnitude blood T1 (this value should be around 1700, or ~ 1000-3000)
    if x.Q.BloodT1<10
        x.Q.BloodT1 = x.Q.BloodT1.*1000;
    end

    if ~isfield(x.Q,'LabelingType')
           error('Unknown LabelingType, needed for quantification');
    elseif isempty(regexpi(x.Q.LabelingType,'^(PC|P|C)ASL$'))
           error('x.Q.LabelingType was invalid, should be PASL|CASL|PCASL');
    elseif strcmpi(x.Q.LabelingType,'PCASL')
           x.Q.LabelingType = 'CASL';
    end

    if ~isfield(x.Q,'BackgroundSuppressionNumberPulses')
        warning('No background suppression pulses known, assuming no background suppression');
        x.Q.BackgroundSuppressionNumberPulses = 0;
    end


    %% 7.   Labeling efficiency
    if ~isfield(x.Q,'LabelingEfficiency') || isempty(x.Q.LabelingEfficiency)
        if ~isfield(x.Q,'LabelingEfficiency')
            switch lower(x.Q.LabelingType)
                case 'pasl'
                    x.Q.LabelingEfficiency = 0.98; % (concensus paper, Wong, MRM 1998)
                case 'casl'
                    x.Q.LabelingEfficiency = 0.85; % (concensus paper, Dai, MRM 2008)
            end
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
[~, CBF] = xASL_quant_SinglePLD(PWI, M0_im, SliceGradient, x, bUseBasilQuantification); % also runs BASIL, but only in native space!

if x.ApplyQuantification(5)==0
    MeanCBF = xASL_stat_MeanNan(CBF(:));
    if MeanCBF>666 % this is the average including air
        CBF = CBF .* (10./MeanCBF);
        warning('M0 division was disabled & CBF image had too high values');
        fprintf('%s\n',['mean whole image CBF normalized from ' xASL_num2str(MeanCBF) ' to 10 mL/100g/min']);
    end
end




%% ------------------------------------------------------------------------------------------------
%% 9.	Save files

fprintf('%s\n','Saving PWI & CBF niftis');

% Unmasked CBF
if bUseBasilQuantification
    [Fpath, Ffile, Fext] = xASL_fileparts(OutputPath);
    OutputPath = fullfile(Fpath, [Ffile '_Basil' Fext]);
end

xASL_io_SaveNifti(PWI_Path, OutputPath, CBF, 32, 0);


%% ------------------------------------------------------------------------------------------------
%% 10.   FEAST quantification
% run FEAST quantification if crushed & non-crushed ASL sessions exist
xASL_quant_FEAST(x);


%% ------------------------------------------------------------------------------------------------
%% 11.  Create standard space masked image to visualize masking effect
if strcmp(OutputPath, x.P.Pop_Path_qCBF)
    % Load CBF image
    MaskedCBF = xASL_io_Nifti2Im(x.P.Pop_Path_qCBF);
    % Mask vascular voxels (i.e. set them to NaN)
    MaskVascularMNI = xASL_io_Nifti2Im(x.P.Pop_Path_MaskVascular);
    MaskedCBF(~MaskVascularMNI) = NaN;
    
    % Mask susceptibility voxels (i.e. set them to NaN)
    if strcmpi(x.Sequence, '2d_epi') || strcmpi(x.Sequence, '3d_grase')
        if ~xASL_exist(x.P.Path_Pop_MaskSusceptibility)
            warning([x.P.Path_Pop_MaskSusceptibility ' missing, cannot create ' x.P.Pop_Path_qCBF_masked]);
            return;
        else
            MaskSusceptibility = xASL_io_Nifti2Im(x.P.Path_Pop_MaskSusceptibility);
            MaskedCBF(~MaskSusceptibility) = NaN;
        end
    end
    xASL_io_SaveNifti(x.P.Pop_Path_qCBF, x.P.Pop_Path_qCBF_masked, MaskedCBF, [], false);
end


end