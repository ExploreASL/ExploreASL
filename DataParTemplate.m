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
% EXAMPLE: ExploreASL_Master('//StudyFolder/analysis/DataPar_HiQ.json');
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% Copyright 2015-2019 ExploreASL

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
%                       DICOM export. Will move this from ASL4D.nii to M0.nii(OPTIONAL, DEFAULT = no M0 in timeseries) 
%                     - example for Philips 3D GRASE = [1 2] % (first control-label pair)
%                     - example for Siemens 3D GRASE = 1 % first image
%                     - example for GE 3D spiral = 2 % where first image is PWI & last = M0

% -----------------------------------------------------------------------------------------------------------------------------------------------------
% SEQUENCE PARAMETERS
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% x.Q.BackGrSupprPulses - used to estimate decrease of labeling efficiency (REQUIRED)
%                       - options:
%                         - 0 = (no background suppression)
%                         - 2 = labeling efficiency factor 0.83 (e.g. Philips 2D EPI & Siemens 3D GRASE)
%                         - 4 = labeling efficiency factor 0.81 (e.g. Philips 3D GRASE)
%                         - 5 = labeling efficiency factor 0.75 (e.g. GE 3D spiral)
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
%                         calculated and saved (OPTIONAL, DEFAULT = [1 1 1 1 1] = all enabled)
%                       - fields:
%                         - 1) Apply ScaleSlopes ASL4D (xASL_wrp_Quantify, future at dcm2niiX stage)
%                         - 2) Apply ScaleSlopes M0 (xASL_quant_M0, future at dcm2niiX stage)
%                         - 3) Convert PWI a.u. to label (xASL_wrp_Quantify, future at xASL_wrp_Reslice?)
%                         - 4) Quantify M0 a.u. (xASL_quant_M0)
%                         - 5) Perform division by M0
%                       - examples: 
%                         - ASL4D is an already quantified CBF image, disable all quantification [0 0 0 0 0]
%                         - To compare label but not CBF (e.g. label in vessels or sinus vs tissue): [1 1 1 1 0]
%                       - Note that the output always goes to CBF.nii
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% PROCESSING PARAMETERS
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% x.SpikeRemovalThreshold - minimal t-stat improval needed to remove motion spikes (OPTIONAL, DEFAULT = 0.01) 
%                         - examples:
%                           - 1 = effectively disabling spike removal
% x.motion_correction - boolean to perform motion correction in case of timeseries (OPTIONAL, DEFAULT = 1)
%                     - options:
%                       - 1 = on
%                       - 0 = off
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
% x.SegmentSPM12 - boolean to specify if SPM12 is run instead of CAT12 (OPTIONAL, DEFAULT = 0)
%                - options:
%                  - 1 = run SPM12
%                  - 0 = run CAT12
% x.bFixResolution - resample to a resolution that CAT12 accepts (OPTIONAL, DEFAULT=false)
% x.bRegistrationContrast - specifies the image contrast used for
%                           registration (OPTIONAL, DEFAULT = 2):
%                           - 0 = Control-T1w
%                           - 1 = CBF - pseudoCBF from template/pGM+pWM
%                           - 2 = automatic (mix of both)
% x.bAffineRegistration - specifies if the ASL-T1w rigid-body
%                         registration is followed up by an affine
%                         registration (OPTIONAL, DEFAULT = 2)
%                  - 0 = affine registration disabled
%                  - 1 = affine registration enabled
%                  - 2 = affine registration automatically chosen based on
%                        spatial CoV of PWI

x.name = ExampleDataSet;
x.subject_regexp = '^Sub-\d{3}$';
x.M0 = 'separate_scan';
x.readout_dim = '2D';
x.Quality = 0;
x.DELETETEMP = 1;
x.Vendor = 'Philips';
x.Q.BackGrSupprPulses = 2;
x.LabelingType = 'CASL';
x.Initial_PLD = 1525;
x.LabelingDuration = 1650;
x.SliceReadoutTime = 43.7647;

end
