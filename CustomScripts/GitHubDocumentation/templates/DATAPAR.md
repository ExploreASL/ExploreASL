
----
## Import Parameters

Input parameters of ExploreASL can be specified in the `DataPar.json` file of the corresponding study. Required and optional parameters are displayed in the table below and in the ExploreASL `DataParTemplate.json`.


----
### Environment Parameters

----
#### FSL
Parameter    | x.bAutomaticallyDetectFSL
----------   | ----------------------------
Description  | Boolean to automatically detect the FSL version if disabled, this function will try to use the system-initialized FSL and throw an error if FSL is not initialized.
Example      | n/a
Required     | optional
Default      | disabled

----
#### Make NIfTI 4D DICOM
Parameter    | x.MakeNIfTI4DICOM
----------   | ----------------------------
Description  | Boolean to output CBF native space maps resampled and/or registered to the original T1w/ASL, and contrast adapted and in 12 bit  range allowing to convert the NIfTI to a DICOM file, e.g. for implementation in PACS or other DICOM archives.
Example      | n/a
Required     | n/a
Default      | n/a

----
### Study Parameters

----
#### Name
Parameter    | x.name
----------   | ----------------------------
Description  | String for the name of the study.
Example      | `AD-study`
Required     | optional
Default      | n/a

----
#### Root directory
Parameter    | x.D.ROOT
----------   | ----------------------------
Description  | Path to analysis root folder where data is stored.
Example      | `/home/hjmutsaerts/TestDataSet`
Required     | optional
Default      | `pwd`

----
#### Regular Expression
Parameter    | x.subject_regexp
----------   | ----------------------------
Description  | String with regular expression for ExploreASL to find subjects by foldername.
Example      | `^\d{3}$` for three digits
Required     | required
Default      | n/a

----
#### Exclusion
Parameter    | x.exclusion
----------   | ----------------------------
Description  | Cell with list of subjects to exclude.
Example      | `{'005' '018'}`
Required     | optional
Default      | empty

----
#### Sessions
Parameter    | x.SESSIONS
----------   | ----------------------------
Description  | Use this to define sessions.
Example      | `{'ASL_1' 'ASL_2'}`
Required     | optional
Default      | `{'ASL_1'}`

----
#### Session Options
Parameter    | x.session.options
----------   | ----------------------------
Description  | This is how the sessions will be called. For FEAST, this should be `{'non-crushed' 'crushed'}`.
Example      | `{'baseline' 'drug'}`
Required     | optional
Default      | n/a

----
#### Inclusion List
Parameter    | x.ForceInclusionList
----------   | ----------------------------
Description  | Use this field if you want to use a selection of subjects rather than taking all available subjects from directories.
Example      | `load(fullfile(x.D.ROOT,'LongitudinalList.mat'))`
Required     | optional
Default      | use all subjects

----
### M0 Parameters and Options

----
#### Conventional Processing
Parameter    | x.M0_conventionalProcessing
----------   | ----------------------------
Description  | Boolean, use the conventional M0 processing (per consensus paper). 1 = standard processing, 0 = new image processing (improved masking & smoothing).
Example      | n/a
Required     | optional
Default      | 0

----
#### M0
Parameter    | x.M0
----------   | ----------------------------
Description  | Choose which M0 option to use. Options: `separate_scan` = for a separate M0 NIfTI (needs to be in the same folder called M0.nii), `3.7394*10^6` = single M0 value to use, `UseControlAsM0` = will copy the mean control image as M0.nii and process as if it was a separately acquired M0 image (taking TR etc from the ASL4D.nii). Make sure that no background suppression was used, otherwise this option is invalid.
Example      | n/a
Required     | required
Default      | n/a

----
#### GM Scale Factor
Parameter    | x.M0_GMScaleFactor
----------   | ----------------------------
Description  | Add additional scale factor to multiply the M0 image by. This can be useful when you have background suppression but no control/M0 image without background suppression. If you then know the M0 scalefactor for the GM, you can use the control image as M0 and use this parameter to scale back what was suppressed by background suppression. Note that there is no option for separate tissue scaling (e.g. WM & GM), because ExploreASL pragmatically smooths the M0 a lot, assuming that head motion and registration between M0 & ASL4D will differ between patients and controls. 
Example      | n/a
Required     | optional
Default      | 1

----
#### M0 Position in ASL Scan
Parameter    | x.M0PositionInASL4D
----------   | ----------------------------
Description  | Indicates the position of M0 in TimeSeries, if it is integrated by the vendor in the DICOM export. Will move this from ASL4D.nii to M0.nii.
Example      | Philips 3D GRASE = `[1 2]` (first control-label pair), Siemens 3D GRASE = `1` (first image), GE 3D spiral = `2` (where first image is PWI & last = M0)
Required     | optional
Default      | no M0 in timeseries

----
### Sequence Parameters

----
#### Background Suppression Pulses
Parameter    | x.Q.BackGrSupprPulses
----------   | ----------------------------
Description  | Used to estimate decrease of labeling efficiency. Options: 0 = (no background suppression), 2 = labeling efficiency factor 0.83 (e.g. Philips 2D EPI & Siemens 3D GRASE), 4 = labeling efficiency factor 0.81 (e.g. Philips 3D GRASE), 5 = labeling efficiency factor 0.75 (e.g. GE 3D spiral).
Example      | n/a
Required     | required
Default      | n/a

----
#### Readout Dimension
Parameter    | x.readout_dim
----------   | ----------------------------
Description  | String specifying the readout type. Options: `2D` for slice-wise readout,  `3D` for volumetric readout.
Example      | n/a
Required     | required
Default      | n/a

----
#### Vendor
Parameter    | x.Vendor
----------   | ----------------------------
Description  | String containing the vendor used. This parameter is used to apply the vendor-specific scale factors. Options: `GE_product`, `GE_WIP`, `Philips`, `Siemens`. 
Example      | n/a
Required     | required
Default      | n/a

----
#### Sequence
Parameter    | x.Sequence
----------   | ----------------------------
Description  | String containing the sequence used. Options: `3D_spiral`, `3D_GRASE`, `2D_EPI`.
Example      | n/a
Required     | required
Default      | n/a

----
#### Labeling Type
Parameter    | x.Q.LabelingType
----------   | ----------------------------
Description  | String containing the labeling strategy used. Options: 'PASL' (pulsed Q2-TIPS), 'CASL' (CASL/PCASL), Note: pulsed without Q2TIPS cannot be reliably quantified because the bolus width cannot be identified CASL & PCASL are both continuous ASL methods, identical quantification.
Example      | n/a
Required     | required
Default      | n/a

----
#### Initial Post Labeling Delay
Parameter    | x.Q.Initial_PLD
----------   | ----------------------------
Description  | Value of PLD (ms), for 3D this is fixed for whole brain, for 2D this is the PLD of first acquired slice.
Example      | 1800
Required     | required
Default      | n/a

----
#### Labeling Duration
Parameter    | x.Q.LabelingDuration
----------   | ----------------------------
Description  | Value of labeling duration (ms).
Example      | 1800
Required     | required
Default      | n/a

----
#### Slice Readout Time
Parameter    | x.Q.SliceReadoutTime
----------   | ----------------------------
Description  | Value (ms) of time added to the PLD after reading out each slice. Other option = `shortestTR`; shortest TR enabled gives each sequence the minimal TR. This enables calculating slice delay per subject.
Example      | 31
Required     | required for 2D ASL sequences
Default      | n/a

----
### Quantification Parameters

----
#### Lambda
Parameter    | x.Q.Lambda
----------   | ----------------------------
Description  | Brain/blood water coefficient (mL 1H/ mL blood).
Example      | 0.32 (for GSP phantom)
Required     | optional
Default      | 0.9

----
#### T2* of arterial blood
Parameter    | x.Q.T2art
----------   | ----------------------------
Description  | T2* of arterial blood at 3T, only used when no M0 image (ms).
Example      | n/a
Required     | optional
Default      | 50

----
#### Blood T1
Parameter    | x.Q.BloodT1
----------   | ----------------------------
Description  | T1 relaxation time of arterial blood (ms).
Example      | n/a
Required     | optional
Default      | 1650 @ 3T, 1350 @ 1.5 T, (Alsop MRM 2014), 1800 for GSP phantom

----
#### Tissue T1
Parameter    | x.Q.TissueT1
----------   | ----------------------------
Description  | T1 relaxation time of GM tissue (ms).
Example      | n/a
Required     | optional
Default      | 1240 @ 3T, 920 @ 1.5 T (Alsop MRM 2014)

----
#### Compartments
Parameter    | x.Q.nCompartments
----------   | ----------------------------
Description  | Number of modeled compartments for quantification. Options: 1 = a single-compartment quantification model (default by concensus paper), 2 = a dual-compartment quantification model.
Example      | n/a
Required     | optional
Default      | 1

----
#### Apply Quantification
Parameter    | x.ApplyQuantification
----------   | ----------------------------
Description  | A vector of 1x5 logical values specifying which types on quantified images should be calculated and saved. Fields: 1) Apply ScaleSlopes ASL4D (xASL_wrp_Quantify, future at dcm2niiX stage), 2) Apply ScaleSlopes M0 (xASL_quant_M0, future at dcm2niiX stage), 3) Convert PWI a.u. to label (xASL_wrp_Quantify, future at xASL_wrp_Reslice?), 4) Quantify M0 a.u. (xASL_quant_M0, corrects for incomplete T1 relaxation), 5) Perform division by M0.
Example      | ASL4D is an already quantified CBF image, disable all quantification '[0 0 0 0 0]'. To compare label but not CBF (e.g. label in vessels or sinus vs tissue): `[1 1 1 1 0]`. Note that the output always goes to CBF.nii.
Required     | optional
Default      | `[1 1 1 1 1]` = all enabled

----
### Processing Parameters

----
#### Spike Removal Threshold
Parameter    | x.SpikeRemovalThreshold
----------   | ----------------------------
Description  | Minimal t-stat improval needed to remove motion spikes.
Example      | 1 = effectively disabling spike removal
Required     | optional
Default      | 0.01

----
#### Motion Correction
Parameter    | x.motion_correction 
----------   | ----------------------------
Description  | Boolean to perform motion correction in case of timeseries. Options: 1 = on, 0 = off.
Example      | n/a
Required     | optional
Default      | 1

----
#### Quality
Parameter    | x.Quality
----------   | ----------------------------
Description  | Boolean specifying on which quality the pipeline should be run. Options: 1 = normal quality, 0 = lower quality, fewer iterations and lower resolution of processing for a fast try-out.
Example      | n/a
Required     | optional
Default      | 1

----
#### Delete Temporary Files
Parameter    | x.DELETETEMP
----------   | ----------------------------
Description  | Boolean for removing the temporary files. Options: 0 = keeping all files, 1 = delete temporary files created by the pipeline.
Example      | n/a
Required     | optional
Default      | 1

----
#### Skip if no FLAIR
Parameter    | x.SkipIfNoFlair
----------   | ----------------------------
Description  | Boolean to skip processing of subjects that do not have a FLAIR image. These parameters can be useful when some data is still complete, but one would like to start image processing already. Options: 1 = skip processing of a subject that does not have a FLAIR image, 0 = do not skip anything.
Example      | n/a
Required     | optional
Default      | 0

----
#### Skip if no ASL
Parameter    | x.SkipIfNoASL
----------   | ----------------------------
Description  | Boolean to skip processing of subjects that do not have a ASL image. Options: 1 = skip processing of a subject that does not have a ASL image, 0 = do not skip anything.
Example      | n/a
Required     | optional
Default      | 0

----
#### Skip if no M0
Parameter    | x.SkipIfNoM0
----------   | ----------------------------
Description  | Boolean to skip processing of subjects that do not have a M0 image. Options: 1 = skip processing of a subject that does not have a M0 image, 0 = do not skip anything.
Example      | n/a
Required     | optional
Default      | 0

----
#### Segment SPM12
Parameter    | x.SegmentSPM12
----------   | ----------------------------
Description  | Boolean to specify if SPM12 segmentation is run instead of CAT12 (OPTIONAL, DEFAULT = 0). Options: 1 = run SPM12, 0 = run CAT12.
Example      | n/a
Required     | optional
Default      | 0

----
#### Hammers CAT 12
Parameter    | x.bHammersCAT12
----------   | ----------------------------
Description  | Boolean specifying if CAT12 should provide Hammers volumetric ROI results.
Example      | n/a
Required     | optional
Default      | 0

----
#### Fix Resolution
Parameter    | x.bFixResolution
----------   | ----------------------------
Description  | Resample to a resolution that CAT12 accepts.
Example      | n/a
Required     | optional
Default      | false

----
#### Registration Contrast
Parameter    | x.bRegistrationContrast
----------   | ----------------------------
Description  | specifies the image contrast used for registration: 0 = Control->T1w, 1 = CBF->pseudoCBF from template/pGM+pWM (skip if sCoV>0.667), 2 = automatic (mix of both), 3 = option 2 & force CBF->pseudoCBF irrespective of sCoV.
Example      | n/a
Required     | optional
Default      | 2

----
#### Affine Registration
Parameter    | x.bAffineRegistration
----------   | ----------------------------
Description  | Specifies if the ASL-T1w rigid-body registration is followed up by an affine registration: 0 = affine registration disabled, 1 = affine registration enabled, 2 = affine registration automatically chosen based on spatial CoV of PWI.
Example      | n/a
Required     | optional
Default      | 2

----
#### DCT Registration
Parameter    | x.bDCTRegistration
----------   | ----------------------------
Description  | Specifies if to include the DCT registration on top of Affine, all other requirements for affine are thus also taken into account the x.bAffineRegistration must be >0 for DCT to run: 0 = DCT registration disabled, 1 = DCT registration enabled if affine enabled and conditions for affine passed, 2 = DCT enabled as above, but use PVC on top of it to get the local intensity scaling right.
Example      | n/a
Required     | optional
Default      | 0

----
#### Register M0 to ASL
Parameter    | x.bRegisterM02ASL
----------   | ----------------------------
Description  | Boolean specifying whether M0 is registered to mean_control image (or T1w if no control image exists). It can be useful to disable M0 registration if the ASL registration is done based on the M0, and little motion is expected between the M0 and ASL acquisition. If no separate M0 image is available, this parameter will have no effect. This option is disabled automatically for 3D spiral: 0 = M0 registration disabled, 1 = M0 registration enabled (DEFAULT).
Example      | n/a
Required     | optional
Default      | 0

----
#### Use MNI as dummy Structural
Parameter    | x.bUseMNIasDummyStructural
----------   | ----------------------------
Description  | When structural (e.g. T1w) data is missing, copy population-average MNI templates as dummy structural templates. With this option, the ASL module copies the structural templates to fool the pipeline, resulting in ASL registration to these templates. While the rigid-body parameters might still be found somewhat correctly, with this option it is advised to enable affine registration for ASL as well, since ASL and these dummy structural images will differ geometrically. When disabled, an error will be issued instead when the structural images are missing. 1 = enabled, 0 = disabled.
Example      | n/a
Required     | optional
Default      | 0

----
#### Masking
Parameter    | x.S.bMasking
----------   | ----------------------------
Description  | Vector specifying if we should mask a ROI with a subject-specific mask (1 = yes, 0 = no), `[1 0 0 0]` = susceptibility mask (either population-or subject-wise), `[0 1 0 0]` = vascular mask (only subject-wise),  `[0 0 1 0]` = subject-specific tissue-masking (e.g. pGM>0.5), `[0 0 0 1]` = WholeBrain masking (used as memory compression), `[0 0 0 0]` = no masking at all, `[1 1 1 1]` = apply all masks. Can also be used as boolean, where: `1 = [1 1 1 1]`, `0 = [0 0 0 0]`. Can be useful for e.g. loading lesion masks outside the GM.
Example      | n/a
Required     | optional
Default      | 1


