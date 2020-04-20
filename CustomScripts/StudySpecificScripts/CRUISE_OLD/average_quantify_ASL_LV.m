function average_quantify_ASL( x, InputPath  ) % InputPath = rtemp_ASL4D.nii
%average_quantify_ASL % Quantification is performed here according to ASL consensus paper (Alsop, MRM 2016)
%
%
%
%
% 1    Prepare M0 image
% 2    Prepare CBF image
% 3    Load slice gradient if 2D
% 4    CBF quantification equation
% 5    Outlier rejection
% 6    Division by M0 & scale slopes
% 7    Remove non-perfusion values
%
%       BACKGROUND INFORMATION


%% Administration

% Cave: file can either be equal to x.P.ASL4D, or to ['despiked_' x.P.ASL4D]

[path file ext]                 = fileparts(x.despiked_raw_asl);
[path_dummy x.P.SessionID dummy] 	= fileparts(path);
[path_dummy x.P.SubjectID dummy] 	= fileparts(path_dummy);
clear dummy path_dummy path ext

CBF_nii                 = InputPath; % rtemp_ASL4D.nii LENA. Otherwise its fullfile(x.D.PopDir, ['PWI_' x.P.SubjectID '_' x.P.SessionID '.nii']);
%CBF_nii                 = fullfile(x.D.ROOT, x.subject_regexp, {'ASL_1'}, 'rtemp_PWI4D.nii'); % Added by Lena for 4D quantification
%CBF_nii                 = fullfile(x.D.PopDir, ['PWI_' x.P.SubjectID '_' x.P.SessionID '.nii']);
slice_gradient_file     = ['slice_gradient_' x.P.SubjectID '_' x.P.SessionID '.nii'];
slice_gradient_load     = fullfile(x.D.PopDir, slice_gradient_file);

qCBF_nii        = fullfile(x.D.PopDir,['q' x.P.CBF '_' x.P.SubjectID '_' x.P.SessionID  '.nii']);
UnTreat_nii     = fullfile(x.D.PopDir,['qCBF_untreated_'  x.P.SubjectID '_' x.P.SessionID '.nii']);
Diff_nii        = fullfile(x.D.PopDir,['diff_vasc_treat_'  x.P.SubjectID '_' x.P.SessionID '.nii']);

%% Get identification
% Identify current subject
for iS=1:x.nSubjects
    if  strcmp(x.SUBJECTS{iS},x.P.SubjectID)
        iSubj=iS;
    end
end

% Identify current session
for iS=1:x.nSessions
    if  strcmp(x.SESSIONS{iS},x.P.SessionID)
        iSess=iS;
    end
end

iSubjSess   = (iSubj-1)*x.nSessions + iSess;

%% Ability to skip quantification, if images were already quantified

if ~isfield(x,'DisableQuantification')
    % When DisableQuantification is activated,
    % ASL_im & PWI will be the unchanged input image but
    % ScaleImage will still apply the scaleslopes (for Philips).
    % Other vendor-specific scaling will not be applied, as this mostly is
    % only used in the raw ASL data
    x.DisableQuantification   = 0;
end

if ~isfield(x.Q,'nCompartments')
    x.Q.nCompartments   = 1; % by default, we use a single-compartment model, as proposed by the Alsop et al. MRM 2014 concensus paper
end

if ~isfield(x.Q,'ATT')
    x.Q.ATT             = 1800; % ms as default micro-vascular ATT
end

if ~isfield(x.Q,'TissueT1')
    x.Q.TissueT1             = 1240; % T1 GM tissue @ 3T
end

% Load ASL PWI
tnii            = nifti(CBF_nii); % Load CBF nifti
ASL_im          = single(tnii.dat(:,:,:,:)); % Added by Lena for 4D quantification

if  xASL_stat_SumNan(ASL_im(:))==0
    error('Empty M0 image, something went wrong in M0 processing');
end

% ADDED BY LENA for rtemp_ASL4D.nii input from 4D
if size(ASL_im,4)>1
    [control_im label_im]   = Check_control_label_order( ASL_im );
    ASL_im                  = control_im - label_im;
    PWI                     = xASL_stat_MeanNan(ASL_im,4);
end

%% Hematocrit correction x.T1a

for iS=1:length(x.S.SetsName)
    if  strcmp(x.S.SetsName{iS},'Hematocrit')
        x.Q.BloodT1 = calc_blood_t1_from_Hct( x.S.SetsID(iSubjSess,iS) );
    end
end

%% Blood T1 correction
for iS=1:length(x.S.SetsName)
    if  strcmp(x.S.SetsName{iS},'BloodT1')
        x.Q.BloodT1     = x.S.SetsID(iSubjSess,iS);
    end
end

%% Correct order of magnitude blood T1 (this value should be around 1700, or ~ 1000-3000)
if  x.Q.BloodT1<10
    x.Q.BloodT1          = x.Q.BloodT1.*1000;
end

%% 1    Prepare M0 image
fprintf('%s\n','Preparing M0 image');
if      isnumeric(x.M0)
        % Single value per scanner
        % In this case we assume that this nifti value has been properly acquired,
        % does not need any corrections, and whole ASL PWI will be divided by this single value.

        % Check if M0 values are defined per SubjectSession:
        for iS=1:length(x.S.SetsName)
            if  strcmp(x.S.SetsName{iS},'M0')
                x.M0     = x.S.SetsID(iSubjSess,iS);
            end
        end


        fprintf('%s\n',['Single M0 value ' num2str(x.M0) ' used']);
        % If obtained by e.g. CSF inversion recovery, make sure that this is corrected for blood-water partition coefficient (0.76) and density of brain tissue (1.05 g/mL)

        M0_im   = x.M0;
else
    % In this case we have a proper M0 image that should have voxel-wise
    % orientation agreement with the ASL PWI. We may need to correct it"s values.

    fprintf('%s\n','M0 scan used');
    M0_nii  = fullfile(x.D.PopDir, [x.P.M0 '_' x.P.SubjectID '_' x.P.SessionID '.nii']);
    M0_im   = xASL_io_ReadNifti(M0_nii);

    % Load M0 parameter file
    % NB: M0 quantification has largely been done in previous script,
    % load M0-parms only to check that ASL & M0-parms are identical
    M0_parms_file= fullfile(x.SESSIONDIR, x.M0PARMSFILE);
    if exist(M0_parms_file,'file')
        M0_P = load(M0_parms_file,'-mat');
        if isfield(M0_P,'parms')
            M0_parms = M0_P.parms;
        else
            fprintf(2,'ERROR in module_ASL: could not load M0 parameters from:\n%s\n',M0_parms_file);
            return % exit with an error
        end
        clear M0_P M0_parms_file
    end

    if  isfield(M0_parms,'RescaleSlope') && isfield(M0_parms,'RescaleSlopeOriginal')
        if  M0_parms.RescaleSlopeOriginal~=M0_parms.RescaleSlope
            warning('M0 RescaleSlope & RescaleSlopeOriginal are not identical, check if scaling was applied already on Philips scanner');
        end
    end

    % Confirm that niftis have same matrix size
    if tnii.dat.dim(1)~=M0_im.dat.dim(1) || tnii.dat.dim(2)~=M0_im.dat.dim(2) || tnii.dat.dim(3)~=M0_im.dat.dim(3)
        error('M0 has different matrix size than ASL-scan, should be the same!');
    end
    % Confirm that nifti is properly averaged, having only 3 dimensions
    if length(M0_im.dat.dim)>3
        error('M0 image has more than 3 dimensions! Inproperly averaged?');
    end

    M0_im           = single(M0_im.dat(:,:,:)); % Get M0 image
    M0_im           = M0_im ./ x.Q.Lambda;
    fprintf('%s\n',['M0 image corrected for Labda: ' num2str(x.Q.Lambda)]);
    % qnt_labda (0.9) = brain-blood partition coefficient.
    % Used for M0 image
    % Scale slopes & incomplete T1 relaxation were already corrected in M0 module
end

if  xASL_stat_SumNan(M0_im(:))==0
    error('Empty M0 image, something went wrong in M0 processing');
end

%% 2    Prepare PWI image
fprintf('%s\n','Preparing PWI image');

% Load ASL parameter file
ASL_parms_file= fullfile(x.SESSIONDIR, x.ASLPARMSFILE);
if exist(ASL_parms_file,'file')
    ASL_P = load(ASL_parms_file,'-mat');
    if isfield(ASL_P,'parms')
        ASL_parms = ASL_P.parms;
    else
        fprintf(2,'ERROR in module_ASL: could not load ASL parameters from:\n%s\n',M0_parms_file);
        return % exit with an error
    end
    clear ASL_P ASL_parms_file
end

if  isfield(ASL_parms,'RescaleSlope') && isfield(ASL_parms,'RescaleSlopeOriginal')
    if  ASL_parms.RescaleSlopeOriginal~=ASL_parms.RescaleSlope
        warning('ASL RescaleSlope & RescaleSlopeOriginal are not identical, check if scaling was applied already on Philips scanner');
    end
end

switch x.M0
    case 'separate_scan'
        % Check equality of TE, but allow them to be 1% different, % Throw error if TE of ASL and M0 are not exactly the same!
        if  ASL_parms.EchoTime<(M0_parms.EchoTime*0.99) || ASL_parms.EchoTime>(M0_parms.EchoTime*1.01)
            fprintf('%s\n','With a separate M0, echotimes of ASL and M0 should be identical');
            fprintf('%s\n','But they were not!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        end
end

%% Possibility for individual initial PLD (different scans have different PLDs)
qnt_init_PLD    = x.Q.Initial_PLD;

for iSet=1:length(x.S.SetsName)
    if  strcmp(x.S.SetsName{iSet},'Init_PLD')
        qnt_init_PLD    = x.S.SetsID(iSubjSess,iSet);
    end
end

% Create factor, single compartment model, including T2* correction (if M0 = number)
% Take care of PLD
switch x.readout_dim
    case '3D'
        fprintf('%s\n','3D sequence, not accounting for SliceReadoutTime (homogeneous PLD for complete volume)');
        qnt_slice_gradient         = qnt_init_PLD;

    case '2D' % Load slice gradient
%% 3    Load slice gradient if 2D
        fprintf('%s\n','2D sequence, accounting for SliceReadoutTime (inhomogeneous/slice-specific PLD)');

        if ~isfield(x,'qnt_PLDslicereadout')
            error('x.Q.SliceReadoutTime was not defined!');

        else
            if  isnumeric(x.Q.SliceReadoutTime)

            else
                switch x.Q.SliceReadoutTime
                    case 'shortestTR'
                        if  isfield(ASL_parms,'RepetitionTime')
                            %  Load original file to get nSlices
                            tORI        = xASL_io_ReadNifti(x.P.ASL4D);
                            nSlices     = size(tORI.dat,3);
                            clear tORI

                            x.Q.SliceReadoutTime     = (ASL_parms.RepetitionTime-x.Q.LabelingDuration-qnt_init_PLD)/nSlices;
                        else error('ASL_parms.RepetitionTime expected but did not exist!');
                        end
                    otherwise;   error('Invalid x.Q.SliceReadoutTime!');
                end
            end

            %% Possibility individual SliceReadoutTime
            for iSet=1:length(x.S.SetsName)
                if  strcmp(x.S.SetsName{iSet},'SliceReadoutTime')
                    x.Q.SliceReadoutTime    = x.S.SetsID(iSubjSess,iSet);
                end
            end


            fprintf('%s\n',['Using SliceReadoutTime ' num2str(x.Q.SliceReadoutTime) ' ms']);
        end

        qnt_slice_gradient      = xASL_io_Nifti2Im( slice_gradient_load );

        % Fix reslicing errors:
        for ii=1:12 % remove NaNs to be sure, this performs an extrapolation 12 times, with 3.5 mm FWHM
            if  sum(isnan(qnt_slice_gradient(:)))>0
                qnt_slice_gradient      = xASL_im_ndnanfilter(qnt_slice_gradient,'gauss',[1.885 1.885 1.885],2);
            end
        end

        qnt_slice_gradient(~isfinite(qnt_slice_gradient))   = 0;
        qnt_slice_gradient                                  = max(qnt_slice_gradient,1);

        if  max(qnt_slice_gradient(:))>50 % not expected in 2D sequences
            error('Erroneous values in slice_gradient_images!');
        end

        qnt_slice_gradient      = qnt_init_PLD + ((qnt_slice_gradient-1) .* x.Q.SliceReadoutTime); % effective PLD

        qnt_slice_gradient      = repmat(qnt_slice_gradient,[1 1 1 size(ASL_im,4)]);

        % Check for erroneous values in slice_gradient_image for 2D
        % sequence

        % Now slice gradient has effective PLD numbers
    otherwise; error('Wrong x.readout_dim value!');
end

if  xASL_stat_SumNan(qnt_slice_gradient(:))==0
    error('Wrong PLD definition!');
end

%% 4    CBF quantification: labeling efficiency
switch x.Q.LabelingType
    case 'PASL'
        x.Q.LabelingEfficiency                     = 0.98; % (concensus paper, Wong, MRM 1998)
    case 'CASL'
        x.Q.LabelingEfficiency                     = 0.85; % (concensus paper, Dai, MRM 2008)
    otherwise error(['Unknown ' x.Q.LabelingType ', should be PASL or CASL (latter also for PCASL)'])
end

% If the labeling efficiency is specified per subject, use this
for iS=1:length(x.S.SetsName)
    if  strcmp(x.S.SetsName{iS},'qnt_lab_eff')
        x.Q.LabelingEfficiency     = x.S.SetsID(iSubjSess,iS);
    end
end

% Apply the effect of background suppression pulses on labeling efficiency
switch x.M0
    case 'no_background_suppression'
        % don't change it
    otherwise
        switch x.Q.BackGrSupprPulses
            case 0 % when you have an M0, but no background suppression used for ASL
                % Then the labeling efficiency doesn't change by background suppression
            case 2 % e.g. Philips 2D EPI or Siemens 3D GRASE
                x.Q.LabelingEfficiency             = x.Q.LabelingEfficiency*0.83; % 0.83 = 2 background suppression pulses (Garcia et al., MRM 2005)
            case 4 % e.g. Philips 3D GRASE
				x.Q.LabelingEfficiency             = x.Q.LabelingEfficiency*0.81; % 0.81 = as implemented by Philips
            case 5 % e.g. GE 3D spiral
                x.Q.LabelingEfficiency             = x.Q.LabelingEfficiency*0.75; % 0.75 = 5 background suppression pulses (GE FSE) (Garcia et al., MRM 2005)
        end
end


%% 4    CBF quantification: CASL & PASL
fprintf('%s\n','Getting CBF quantification factor');

switch x.Q.nCompartments

    case 1 % single-compartment model (Alsop et al MRM 2014 concensus paper)
        switch x.Q.LabelingType
            case 'PASL'
                DivisionFactor     = x.Q.LabelingDuration;

                fprintf('%s\n','Using single-compartment PASL model');
            case 'CASL'
                DivisionFactor     = x.Q.BloodT1 .* (1 - exp(-x.Q.LabelingDuration./x.Q.BloodT1));

                fprintf('%s\n','Using single-compartment CASL model');
        end;

        qnt_slice_gradient = exp((qnt_slice_gradient./x.Q.BloodT1)) ./ (2.*x.Q.LabelingEfficiency.* DivisionFactor );



    case 2 % (Wang 2002) dual-compartment model. For equation, see Gevers et al 2012
        switch x.Q.LabelingType
            case 'PASL'
                DivisionFactor = x.Q.LabelingDuration;
                qnt_slice_gradient = exp((x.Q.ATT./x.Q.BloodT1))*exp(((qnt_slice_gradient-x.Q.ATT)./x.Q.TissueT1))./ (2.*x.Q.LabelingEfficiency.*DivisionFactor );

                fprintf('%s\n','Using Dual compartment PASL model');
            case 'CASL'
                DivisionFactor     = x.Q.TissueT1 .* (exp((min(x.Q.ATT-qnt_slice_gradient,0))./x.Q.TissueT1) - exp((min(x.Q.ATT-x.Q.LabelingDuration-qnt_slice_gradient,0))./x.Q.TissueT1));
                qnt_slice_gradient = exp((x.Q.ATT./x.Q.BloodT1))./ (2.*x.Q.LabelingEfficiency.* DivisionFactor );

                fprintf('%s\n','Using Dual compartment CASL model');
        end;

    otherwise
        error('Wrong x.Q.nCompartments setting');
end

%                                 correcting T1-decay from PLD to label_duration
fprintf('%s\n','ASL quantification:');
fprintf('%s\n',['labda = ' num2str(x.Q.Lambda)]);
fprintf('%s\n',['labeling efficiency = ' num2str(x.Q.LabelingEfficiency)]);
fprintf('%s\n',['T1 arterial blood = ' num2str(x.Q.BloodT1) ' ms']);

switch x.Q.LabelingType
    case 'PASL'
        fprintf('%s\n',['TI1 = ' num2str(x.Q.LabelingDuration)   ' ms']);
        fprintf('%s\n',['initial TI  = ' num2str(qnt_init_PLD) ' ms']);
    case 'CASL'
        fprintf('%s\n',['labeling duration = ' num2str(x.Q.LabelingDuration) ' ms']);
        fprintf('%s\n',['initial post-label delay = ' num2str(qnt_init_PLD) ' ms']);
end

% Now slice gradient has quantification factor numbers
% if  strcmp(x.Vendor,'Siemens') && strcmp(x.M0, 'no_background_suppression') && strcmp(x.readout_dim,'3D')
% % Siemens sometimes seems to multiply by 100, sometimes not, see also below where there is another Siemens-specific scaleslope
%     qnt_slice_gradient              = qnt_slice_gradient*60000;
% else
qnt_slice_gradient              = qnt_slice_gradient*60000*100;
% end

% Now slice gradient has quantification factor numbers scaled to physiological units ( ml/gr/ms =>ml/100gr/min =>(60,000 ms=>min)(1 gr=>100gr) )
% (For some reason, GE sometimes doesn't need the 1 gr->100 gr conversion)
% & old Siemens sequence also didn't need the 1 gr->100 gr conversion

if  x.DisableQuantification % In case the input was a CBF image
    qnt_slice_gradient          = ones(size(ASL_im));
end

PWI                             = ASL_im .* qnt_slice_gradient;

clear qnt_slice_gradient
% Now ASL_image has been scaled with the quantification factors, acknowledging multi-slice acquisitions







%% 5    Division by M0 & scale slopes

ScaleImage          = ones(size(PWI));

% Correct scale slopes
ScaleImage          = ScaleImage./(ASL_parms.RescaleSlopeOriginal.*ASL_parms.MRScaleSlope);
fprintf('%s\n',['ASL image corrected for dicom scale slopes ' num2str(ASL_parms.RescaleSlopeOriginal) ' and ' num2str(ASL_parms.MRScaleSlope)]);
% nifti scale slope has already been corrected for by SPM nifti

% Division by M0
if  isnumeric(x.M0) & ~x.DisableQuantification % do T2* correction of arterial blood
    % in case of separate M0, or M0 because of no background suppression,
    % T2* effect is similar in both images and hence removed by division
        T2_star_factor            = exp(ASL_parms.EchoTime/x.Q.T2art);
        ScaleImage                = ScaleImage .* T2_star_factor;

        fprintf('%s\n',['ASL image corrected for T2* decay during TE, TE was ' num2str(ASL_parms.EchoTime) ' ms, using T2* ' num2str(x.Q.T2art) ' ms, this resulting in factor ' num2str(T2_star_factor)]);
end

if ~x.DisableQuantification
    %ScaleImage          = ScaleImage./M0_im;% uncommented by Lena because of dimensions must agree error
    ScaleImage          = bsxfun(@rdivide,ScaleImage,M0_im);
    fprintf('%s\n','ASL image divided by M0 image');
end

clear M0_im


%% 6    Vendor-specific quantification correction

if ~x.DisableQuantification
    if  ~isempty(findstr(x.Vendor,'GE'))

            if ~isfield(ASL_parms,'NumberOfAverages')
                % GE accumulates signal instead of averaging by NEX, therefore division by NEX is required
                error('GE-data expected, "NumberOfAverages" should be a dicom-field, but was not found!!!')
            end

            ASL_parms.NumberOfAverages  = max(ASL_parms.NumberOfAverages); % fix for combination of M0 & PWI in same nifti, for GE quantification

            switch x.Vendor
                % For some reason the older GE Alsop Work in Progress (WIP) version
                % has a different scale factor than the current GE product sequence

                case 'GE_product' % GE new version

                    qnt_R1gain          = 1/32;          % R1 analogue gain/ simple multiplier
                    qnt_C1              = 6000;          % GE constant multiplier

                    qnt_GEscaleFactor   = (qnt_C1*qnt_R1gain)/(ASL_parms.NumberOfAverages);

                case 'GE_WIP' %     GE old version

                    qnt_RGcorr          = 45.24;         % Correction for receiver gain in PDref (but not used apparently?)
                    qnt_GEscaleFactor   = qnt_RGcorr*ASL_parms.NumberOfAverages;
            end

            ScaleImage            = ScaleImage./qnt_GEscaleFactor;
            fprintf('%s\n',['Quantification corrected for GE scale factor ' num2str(qnt_GEscaleFactor)]);
%               fprintf('%s\n','No GE scale factor applied');

    elseif    strcmp(x.Vendor,'Siemens') && ~strcmp(x.Vendor,'Siemens_JJ_Wang')
                % The Siemens readout divides M0 by 10
                % But JJ Wang's sequence doesn't
              ScaleImage            = ScaleImage./10;
              fprintf('%s\n','M0 corrected for Siemens 3D scale factor 10')
    end
end

%  8    Outlier rejection was removed; this does not work voxel-wise in
%  high-res. It is possible when smoothing (i.e. removing spatially
%  correlated temporal variability) but this remains to be validated


%% 7    Repair negative vascular artifacts by switching sign
% Rationale: vascular contamination, CBF values can turn out far too
% large. However, by virtue of perfusion fluctuation, values may turn out
% negatively. One solution is to create an absolute perfusion map (i.e.
% turn negative values into positive). However, this will also introduce
% noise.

% One property that we use is that vascular artifacts will be
% regions with relatively large negative signal whereas noise will be regions with relatively small
% negative signal. Negative signal by wrong background suppression timing
% might be fixed with this as well.

% We first correct the negative artifacts, and then we fix the large
% macro-vascular signal (which could include negative perfusion artifacts
% that have been post-processed/switched sign

% Process whole image (not to have registration issues)
% but get negative region distribution from GM CBF only

% See \\ExploreASL\Literature\Macrovascular signal processing.docx

GMmask                  = fullfile( x.D.MapsDir, 'rgrey.nii');
GMmask                  = xASL_io_ReadNifti(GMmask);
GMmask                  = GMmask.dat(:,:,:)>0.5;

% NegativeMask            = PWI<0 & GMmask; Uncommented to allow 4D
NegativeMask            = bsxfun(@lt,PWI,0);
NegativeMask            = bsxfun(@and,NegativeMask,GMmask);

if  sum(NegativeMask(:))>0 % if there are subzero voxels

    % Label (i.e. assign a number) to each subzero clusters
    labeled                 = spm_bwlabel(double(NegativeMask),26);

    % Get the mean value of each subzero clusters
    for iP=1:max(labeled(:))
        MeanValue(iP)       = mean(PWI(labeled==iP));
    end

    % Define the threshold below which we will switch sign
    medianValue             = median(MeanValue);
    madValue                = xASL_stat_MadNan(MeanValue,1);
    ClipThr                 = medianValue - (4*madValue);

else
    ClipThr                 = -max(PWI(:)); % don't clip
end

TreatedPWI                      = PWI;
TreatedPWI(TreatedPWI<ClipThr)  = abs(TreatedPWI(TreatedPWI<ClipThr));

% Clip at zero
TreatedPWI(TreatedPWI<0)= 0;
clear medianValue madValue ClipThr GMmask



%% 8   Compress high macro-vascular signal
% We search for voxels that are higher than median + 4 MAD in the PWI image
% & in the PWI/M0 image. Vascular artifacts will have a high intensity in
% both images, whereas errors by division by M0 will only have a high
% intensity on the M0 image, and high values due to a biasfield will only
% be noticeable on the PWI image

TreatedPWI              = VascularCompression( TreatedPWI, ScaleImage );

qCBF                    = TreatedPWI;
TreatedPWI              = TreatedPWI./ScaleImage;

%% 9	Save files

fprintf('%s\n','Saving PWI & CBF niftis');

% Untreated PWI, this can be used for spatial CoV or for WM perfusion
xASL_io_SaveNifti(CBF_nii, UnTreat_nii ,PWI.*ScaleImage,32);

% CBF image that has been treated, and is best fit for statistical analysis
xASL_io_SaveNifti(CBF_nii, qCBF_nii ,qCBF,32); % single precision for better precision

% % Save difference image to be able to check what happened in vascular processing
% % Brainmasked & clipped
%
% DiffPWI                 = x.skull .* (qCBF-(PWI.*ScaleImage));
% DiffPWI                 = ClipVesselImage( DiffPWI,0.95,1,0);
% xASL_io_SaveNifti(CBF_nii, Diff_nii ,DiffPWI,32); % single precision for better precision



%% 10   Computation FEAST-based transit times (uses images that were not vasculary treated)
% if there is a crushed & non-crushed scan, then transit times will be
% computed by division of these scans, provided sessions are exactly named
% as defined below.
% Computation is performed if CurrentSession=session 2

if  x.nSessions>1
    if  strcmp(x.session.options{1},'crushed') && strcmp(x.session.options{2},'non-crushed')
        if  strcmp(x.SESSIONDIR(length(x.SUBJECTDIR)+2:end),'ASL_2')

            fprintf('%s\n','Saving TT nifti');

            xASL_adm_CreateDir(x.D.TTCheckDir);

            %% load data

            MASKnii     = x.skull; % this is not masked for WM or GM, we could do that

            CBF_1Name   = fullfile( x.D.PopDir, ['qCBF_untreated_' x.P.SubjectID '_' x.SESSIONS{1} '.nii'] );
            CBF_2Name   = fullfile( x.D.PopDir, ['qCBF_untreated_' x.P.SubjectID '_' x.SESSIONS{2} '.nii'] );
            gradient_1  = fullfile( x.D.PopDir, ['slice_gradient_' x.P.SubjectID '_' x.SESSIONS{1} '.nii'] );
            gradient_2  = fullfile( x.D.PopDir, ['slice_gradient_' x.P.SubjectID '_' x.SESSIONS{2} '.nii'] );

            CBF_1       = xASL_io_Nifti2Im(CBF_1Name);
            CBF_2       = xASL_io_Nifti2Im(CBF_2Name);
            gradient_1  = xASL_io_Nifti2Im(gradient_1);
            gradient_2  = xASL_io_Nifti2Im(gradient_2);

            % Clip CBF images below zero


            %% Correct different PLD scales & combine them
            PLD_1       = qnt_init_PLD + ((gradient_1-1) .* x.Q.SliceReadoutTime);
            PLD_2       = qnt_init_PLD + ((gradient_2-1) .* x.Q.SliceReadoutTime);

            gradient_1  = exp(PLD_1./x.Q.BloodT1) / (2.*x.Q.LabelingEfficiency.*x.Q.BloodT1 .* (1- exp(-x.Q.LabelingDuration./x.Q.BloodT1)) );
            gradient_2  = exp(PLD_2./x.Q.BloodT1) / (2.*x.Q.LabelingEfficiency.*x.Q.BloodT1 .* (1- exp(-x.Q.LabelingDuration./x.Q.BloodT1)) );
            CBF_1       = CBF_1./gradient_1;
            CBF_2       = CBF_2./gradient_2;

            PLD_combined= (PLD_1+PLD_2)./2;

            N           = 8;



%             % Masking
%             CBF_1(MASKnii==0)         = NaN;
%             CBF_2(MASKnii==0)         = NaN;

            CBF_1                   = CBF_1.*MASKnii;
            CBF_2                   = CBF_2.*MASKnii;


            % smooth maps, ignoring NaNs
            % CAVE: NaNs are interpolated with data,
            % hence the maps should be masked later!
            % gausswin with even kernel (e.g. SD=4) performs best

            CBF_1                   = xASL_im_ndnanfilter(CBF_1,'gauss',[N N N]./1.06,0);
            CBF_2                   = xASL_im_ndnanfilter(CBF_2,'gauss',[N N N]./1.06,0);

            FEAST_ratio             = CBF_1./CBF_2; % crushed/non-crushed
            FEAST_ratio(FEAST_ratio>1)  = 1; % clip @ 1 (== ATT = PLD)
            FEAST_ratio(FEAST_ratio<0)  = 0; % clip @ 0 (== infinite ATT)

            % Calculate
            qnt_PLDdecay            = exp(-PLD_combined/x.Q.BloodT1);
            qnt_combidecay          = exp( (-x.Q.LabelingDuration - PLD_combined) / x.Q.BloodT1);
            TT                      = -x.Q.BloodT1 .* reallog( FEAST_ratio .* (qnt_PLDdecay  - qnt_combidecay ) + qnt_combidecay );

            xASL_io_SaveNifti(CBF_1Name,fullfile(x.D.PopDir,['TT_'  x.P.SubjectID '.nii']),TT,32);
        end
    end
end


%% Housekeeping
if  x.DELETETEMP
    if ~isempty(strfind(x.despiked_raw_asl,'despiked'))
        delete(x.despiked_raw_asl);
        if  xASL_exist( fullfile( x.SESSIONDIR, 'wmask_ICV.nii'))
            xASL_delete( fullfile(x.SESSIONDIR, 'wmask_ICV.nii') );
        end
        if  xASL_exist( fullfile( x.SESSIONDIR, 'rwmask_ICV.nii'))
            xASL_delete( fullfile(x.SESSIONDIR, 'rwmask_ICV.nii') );
        end
    end
end

end
