function x = DataParTemplate(x)
% DataParTemplate Template for study-specific parameter table
%
% FORMAT: x = DataParTemplate(x)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: 
% This template provides an overview of the table of parameters used by ExploreASL, 
% which is very specific to the study to be processed.
% Most of it is ASL quantification-related, some of it are image processing
% parameters also applicable to the structural module.
% It allows only for minimal pipeline modifications (e.g. x.Quality, x.DELETETEMP)
% as most of the ExploreASL environment parameters are loaded through ExploreASL_Initialize.m
%
% Here we list potential data parameters. Most are optional, and are ignored if not provided.
% For an example, see the DataPar*.json file(s) in the TestDataSet.
% For legacy reasons & ease of use, these can be set up in a DataPar*.m,
% which will be converted to the DataPar*.json per BIDS.
% The compiled version of ExploreASL only allows for JSON input.
%
% For the DataPar.json:
% Make sure that booleans are inputted as numbers (e.g. 0 or 1) and not as words (e.g. true or false)
% Scalars can be inputted as scalars, but vectors need to be inputted as string, e.g. '[1 2]'
% This is to allow for valid JSONs. The conversion is carried out internally.
%
% EXAMPLE: ExploreASL_Master('//StudyFolder/analysis/DataPar_HiQ.json');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2021 ExploreASL

% -----------------------------------------------------------------------------------------------------------------------------------------------------
% ENVIRONMENT PARAMETERS
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% x.bAutomaticallyDetectFSL - Boolean to automatically detect the FSL version
%                             if disabled, this function will try to use the system-initialized FSL 
%                             and throw an error if FSL is not initialized
%                             (OPTIONAL, DEFAULT = disabled)
% x.MakeNIfTI4DICOM - Boolean to output CBF native space maps resampled and/or registered to the original T1w/ASL, and contrast adapted and in 12 bit
% 					  range allowing to convert the NIfTI to a DICOM file, e.g. for implementation in PACS or other DICOM archives
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% STUDY PARAMETERS
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% x.name   - string for the name of the study (OPTIONAL)
%          - example: 'AD-study'
% x.D.ROOT - path to analysis root folder where data is stored (OPTIONAL, DEFAULT = pwd) 
%          - example: '/home/hjmutsaerts/TestDataSet'
% x.subject_regexp - string with regular expression for ExploreASL to find subjects by foldername (REQUIRED)
%                  - example: '^\d{3}$' for three digits
% x.exclusion - cell with list of subjects to exclude (OPTIONAL, DEFAULT = empty)
%             - example: {'005' '018'}
% x.SESSIONS  - use this to define sessions (OPTIONAL, DEFAULT = {'ASL_1'})
%             - example: {'ASL_1' 'ASL_2'}
%             - Specific options: for FEAST: 1=crushed, 2=not crushed. This used to be other way around, 
%               but the crushed image registers better with the pGM image
% x.session.options - this is how the sessions will be called (OPTIONAL)
%                   - example: {'baseline' 'drug'}
%                   - For FEAST, this should be {'non-crushed' 'crushed'}
% x.ForceInclusionList - Use this field if you want to use a selection of subjects rather than taking all available subjects
%                        from directories (OPTIONAL, DEFAULT = use all subjects)
%                      - example: load(fullfile(x.D.ROOT,'LongitudinalList.mat'))

% -----------------------------------------------------------------------------------------------------------------------------------------------------
% M0 PARAMETERS and OPTIONS
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% x.M0_conventionalProcessing - boolean - use the conventional M0 processing (per consensus paper) (OPTIONAL, DEFAULT = 0)
%                             - options: 
%                               - 1 = standard processing
%                               - 0 = new image processing (improved masking & smoothing)
% x.M0 - choose which M0 option to use (REQUIRED)
%      - options:
%        - 'separate_scan' = for a separate M0 NIfTI (needs to be in the same folder called M0.nii)
%        - 3.7394*10^6 = single M0 value to use 
%        - 'UseControlAsM0' = will copy the mean control image as M0.nii and process
%                             as if it was a separately acquired M0 image (taking TR etc from the
%                             ASL4D.nii). Make sure that no background suppression was used, otherwise
%                             this option is invalid
% x.M0_GMScaleFactor - add additional scale factor to multiply the M0 image by (OPTIONAL, default = 1)
%                      This can be useful when you have background suppression but no control/M0
%                      image without background suppression. If you then know the M0 scalefactor
%                      for the GM, you can use the control image as M0 and use this parameter to
%                      scale back what was suppressed by background suppression.
%                      Note that there is no option for separate tissue scaling (e.g. WM & GM),
%                      because ExploreASL pragmatically smooths the M0 a lot, assuming that
%                      head motion and registration between M0 & ASL4D will differ between
%                      patients and controls. 
% x.M0PositionInASL4D - indicates the position of M0 in TimeSeries, if it is integrated by the vendor in the 
%                       DICOM export. Will move this from ASL4D.nii to M0.nii(OPTIONAL, DEFAULT = [] (no M0 in timeseries))
%                     - Note that the x.M0PositionInASL4D parameter is independent from the x.M0 parameter choice.
%                     - example for Philips 3D GRASE = '[1 2]' % (first control-label pair)
%                     - example for Siemens 3D GRASE = 1 % first image
%                     - example for GE 3D spiral = 2 % where first image is PWI & last = M0

% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SEQUENCE PARAMETERS
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% x.Q.BackgroundSuppressionNumberPulses - used to estimate decrease of labeling efficiency (REQUIRED)
%                       - options:
%                         - 0 = (no background suppression)
%                         - 2 = labeling efficiency factor 0.83 (e.g. Philips 2D EPI & Siemens 3D GRASE)
%                         - 4 = labeling efficiency factor 0.81 (e.g. Philips 3D GRASE)
%                         - 5 = labeling efficiency factor 0.75 (e.g. GE 3D spiral)
% x.Q.BackgroundSuppressionPulseTime - Vector containing timing, in seconds, 
%                                      of the background suppression pulses
%                                      before the start of the readout (per
%                                      BIDS) (REQUIRED when
%                                      x.Q.UseControlAsM0 &
%                                      x.Q.BackgroundSuppressionNumberPulses>0)
% x.Q.PresaturationTime - time in ms before the start of the readout, scalar, when the slice has been saturated (90 degree flip)
%                    this has to come before all the bSup pulses, but doesn't need to be always specified 
%                    (OPTIONAL, defaults to PLD (PASL) or PLD+LabDur ((P)CASL)
% x.readout_dim - string specifying the readout type (REQUIRED)
%               - options:
%                 - '2D' for slice-wise readout
%                 - '3D' for volumetric readout
% x.Vendor - string containing the vendor used. This parameter is used to apply the vendor-specific 
%            scale factors (REQUIRED for ASL)
%          - options:
%            - 'GE_product' 
%            - 'GE_WIP' 
%            - 'Philips' 
%            - 'Siemens'. 
% x.Sequence - string containing the sequence used (REQUIRED for ASL)
%            - options:
%              - '3D_spiral'
%              - '3D_GRASE'
%              - '2D_EPI'
% x.Q.LabelingType - string containing the labeling strategy used (REQUIRED for ASL)
%                  - options:
%                    - 'PASL' (pulsed Q2-TIPS)
%                    - 'CASL' (CASL/PCASL)
%                  - note: pulsed without Q2TIPS cannot be reliably quantified because the bolus width 
%                    cannot be identified CASL & PCASL are both continuous ASL methods, identical quantification
% x.Q.Initial_PLD - value of PLD (ms), for 3D this is fixed for whole brain, for 2D this is the PLD of first 
%                   acquired slice (REQUIRED for ASL)
%                 - example: 1800
% x.Q.LabelingDuration - value of labeling duration (ms) (REQUIRED for ASL)
%                      - example: 1800
% x.Q.SliceReadoutTime - value (ms) of time added to the PLD after reading out each slice (REQUIRED for 2D ASL sequences)
%                      - example: 31
%                      - Other option = 'shortestTR'; % shortest TR enabled gives each sequence the minimal TR. 
%                        This enables calculating slice delay per subject

% -----------------------------------------------------------------------------------------------------------------------------------------------------
% QUANTIFICATION PARAMETERS
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% x.Q.Lambda - Brain/blood water coefficient (mL 1H/ mL blood), (OPTIONAL, DEFAULT = 0.9) 
%            - example: 0.32 (for GSP phantom)
% x.Q.T2art - T2* of arterial blood at 3T, only used when no M0 image (ms) (OPTIONAL, DEFAULT = 50)
% x.Q.BloodT1 - T1 relaxation time of arterial blood (ms) 
%             - (OPTIONAL, DEFAULT = 1650 @ 3T, 1350 @ 1.5 T, (Alsop MRM 2014), 1800 for GSP phantom)
% x.Q.TissueT1 - T1 relaxation time of GM tissue (ms) (OPTIONAL, DEFAULT=1240 @ 3T, 920 @ 1.5 T (Alsop MRM 2014))
% x.Q.nCompartments - number of modeled compartments for quantification (OPTIONAL, DEFAULT = 1)
%                   - options:
%                     - 1 = a single-compartment quantification model (default by concensus paper)
%                     - 2 = a dual-compartment quantification model 
% x.ApplyQuantification - a vector of 1x5 logical values specifying which types on quantified images should be
%                         calculated and saved (OPTIONAL, DEFAULT = '[1 1 1 1 1]' = all enabled)
%                       - fields:
%                         - 1) Apply ScaleSlopes ASL4D (xASL_wrp_Quantify, future at dcm2niiX stage)
%                         - 2) Apply ScaleSlopes M0 (xASL_quant_M0, future at dcm2niiX stage)
%                         - 3) Convert PWI a.u. to label (xASL_wrp_Quantify, future at xASL_wrp_Reslice?)
%                         - 4) Quantify M0 a.u. (xASL_quant_M0, corrects for incomplete T1 relaxation)
%                         - 5) Perform division by M0
%                       - examples: 
%                         - ASL4D is an already quantified CBF image, disable all quantification '[0 0 0 0 0]'
%                         - To compare label but not CBF (e.g. label in vessels or sinus vs tissue): '[1 1 1 1 0]''
%                       - Note that the output always goes to CBF.nii
% x.Q.SaveCBF4D - boolean, true to also save 4D CBF timeseries, if ASL4D had timeseries (OPTIONAL, DEFAULT=false)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% GENERAL PROCESSING PARAMETERS
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% x.Quality - boolean specifying on which quality the pipeline should be run (OPTIONAL, DEFAULT = 1)
%           - options:
%             - 1 = normal quality 
%             - 0 = lower quality, fewer iterations and lower resolution of processing for a fast try-out
% x.DELETETEMP - boolean for removing the temporary files (OPTIONAL, DEFAULT = 1)
%              - options:
%                - 0 = keeping all files
%                - 1 = delete temporary files created by the pipeline
% x.SkipIfNoFlair - boolean to skip processing of subjects that do not have a FLAIR image (OPTIONAL, DEFAULT = 0)
%                   These parameters can be useful when some data is still complete, but one
%                   would like to start image processing already.
%                 - options:
%                   - 1 = skip processing of a subject that does not have a FLAIR image
%                   - 0 = do not skip anything
% x.SkipIfNoASL   - boolean to skip processing of subjects that do not have a ASL image (OPTIONAL, DEFAULT = 0)
%                 - options:
%                   - 1 = skip processing of a subject that does not have a ASL image
%                   - 0 = do not skip anything
% x.SkipIfNoM0    - boolean to skip processing of subjects that do not have a M0 image (OPTIONAL, DEFAULT = 0)
%                 - options:
%                   - 1 = skip processing of a subject that does not have a M0 image
%                   - 0 = do not skip anything
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% STRUCTURAL PROCESSING PARAMETERS
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% x.bRunModule_LongReg - run longitudinal registration (OPTIONAL, DEFAULT = 0)
% x.bRunModule_DARTEL - run between-subject registration/create templates (OPTIONAL, DEFAULT = 0)
% x.SegmentSPM12 - boolean to specify if SPM12 segmentation is run instead of CAT12 (OPTIONAL, DEFAULT = 0)
%                - options:
%                  - 1 = run SPM12
%                  - 0 = run CAT12
% x.bHammersCAT12 - boolean specifying if CAT12 should provide Hammers volumetric ROI results (OPTIONAL, DEFAULT = 0)
% x.bFixResolution - resample to a resolution that CAT12 accepts (OPTIONAL, DEFAULT=false)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% ASL PROCESSING PARAMETERS
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% x.motion_correction - boolean to perform motion correction in case of timeseries (OPTIONAL, DEFAULT = 1)
%                     - options:
%                       - 1 = on
%                       - 0 = off
% x.SpikeRemovalThreshold - minimal t-stat improval needed to remove motion spikes (OPTIONAL, DEFAULT = 0.01) 
%                         - examples:
%                           - 1 = effectively disabling spike removal
% x.bRegistrationContrast - specifies the image contrast used for
%                                 registration (OPTIONAL, DEFAULT = 2):
%                           - 0 = Control->T1w
%                           - 1 = CBF->pseudoCBF from template/pGM+pWM
%                                 (skip if sCoV>0.667)
%                           - 2 = automatic (mix of both)
%                           - 3 = option 2 & force CBF->pseudoCBF irrespective of sCoV
% x.bAffineRegistration - specifies if the ASL-T1w rigid-body
%                         registration is followed up by an affine
%                         registration (OPTIONAL, DEFAULT = 2)
%                  - 0 = affine registration disabled
%                  - 1 = affine registration enabled
%                  - 2 = affine registration automatically chosen based on
%                        spatial CoV of PWI
% x.bDCTRegistration -  Specifies if to include the DCT registration on top of Affine, all other 
%                            requirements for affine are thus also taken into account (OPTIONAL, DEFAULT = 0)
%                            the x.bAffineRegistration must be >0 for DCT to run
%                          - 0 = DCT registration disabled
%                          - 1 = DCT registration enabled if affine enabled and conditions for affine passed
%                          - 2 = DCT enabled as above, but use PVC on top of it to get the local intensity scaling right
% x.bRegisterM02ASL - boolean specifying whether M0 is registered to
%                     mean_control image (or T1w if no control image exists)
%                     It can be useful to disable M0 registration if the
%                     ASL registration is done based on the M0, and little
%                     motion is expected between the M0 and ASL
%                     acquisition.
%                     If no separate M0 image is available, this parameter
%                     will have no effect. This option is disabled
%                     automatically for 3D spiral
%                     (OPTIONAL, DEFAULT = 0)
%                     - 0 = M0 registration disabled
%                     - 1 = M0 registration enabled (DEFAULT)
% x.bUseMNIasDummyStructural - When structural (e.g. T1w) data is missing, copy population-average
%  							   MNI templates as dummy structural templates. With this option, the
% 							   ASL module copies the structural templates to fool the pipeline,
% 							   resulting in ASL registration to these templates. While the rigid-body
% 							   parameters might still be found somewhat correctly, with this option
% 							   it is advised to enable affine registration for ASL as well, since ASL
% 							   and these dummy structural images will differ geometrically.
% 							   When disabled, an error will be issued instead when the structural images
% 							   are missing.
% 							   (OPTIONAL, DEFAULT = 0).
% 							   - 1 = enabled
% 							   - 0 = disabled
%
%% Masking parameters
%   x.S.bMasking        - vector specifying if we should mask a ROI with a subject-specific mask
%                       (1 = yes, 0 = no)
%                       [1 0 0 0] = susceptibility mask (either population-or subject-wise)
%                       [0 1 0 0] = vascular mask (only subject-wise)
%                       [0 0 1 0] = subject-specific tissue-masking (e.g. pGM>0.5)
%                       [0 0 0 1] = WholeBrain masking (used as memory compression)
%                       [0 0 0 0] = no masking at all
%                       [1 1 1 1] = apply all masks
%                       Can also be used as boolean, where 
%                       1 = [1 1 1 1]
%                       0 = [0 0 0 0]
%                       Can be useful for e.g. loading lesion masks outside the GM
%                       (OPTIONAL, DEFAULT=1)

x.name = ExampleDataSet;
x.subject_regexp = '^Sub-\d{3}$';
x.M0 = 'separate_scan';
x.M0PositionInASL4D = '[1 2]';
x.readout_dim = '2D';
x.Quality = 0;
x.DELETETEMP = 1;
x.Vendor = 'Philips';
x.Q.BackgroundSuppressionNumberPulses = 2;
x.LabelingType = 'CASL';
x.Initial_PLD = 1525;
x.LabelingDuration = 1650;
x.SliceReadoutTime = 43.7647;

end
