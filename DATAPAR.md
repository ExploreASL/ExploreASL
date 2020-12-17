
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
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Readout Dimension
Parameter    | x.readout_dim
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Vendor
Parameter    | x.Vendor
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Sequence
Parameter    | x.Sequence
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Labeling Type
Parameter    | x.Q.LabelingType
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Initial Post Labeling Delay
Parameter    | x.Q.Initial_PLD
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Labeling Duration
Parameter    | x.Q.LabelingDuration
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Slice Readout Time
Parameter    | x.Q.SliceReadoutTime
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
### Quantification Parameters

----
#### Lambda
Parameter    | x.Q.Lambda
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### T2 art
Parameter    | x.Q.T2art
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Blood T1
Parameter    | x.Q.BloodT1
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Tissue T1
Parameter    | x.Q.TissueT1
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Compartments
Parameter    | x.Q.nCompartments
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Apply Quantification
Parameter    | x.ApplyQuantification
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
### Processing Parameters

----
#### Spike Removal Threshold
Parameter    | x.SpikeRemovalThreshold
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Motion Correction
Parameter    | x.motion_correction 
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Quality
Parameter    | x.Quality
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Delete Temporary Files
Parameter    | x.DELETETEMP
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Skip if no FLAIR
Parameter    | x.SkipIfNoFlair
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Skip if no ASL
Parameter    | x.SkipIfNoASL
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Skip if no M0
Parameter    | x.SkipIfNoM0
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Segment SPM12
Parameter    | x.SegmentSPM12
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Fix Resolution
Parameter    | x.bFixResolution
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Registration Contrast
Parameter    | x.bRegistrationContrast
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Affine Registration
Parameter    | x.bAffineRegistration
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### DCT Registration
Parameter    | x.bDCTRegistration
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Register M0 to ASL
Parameter    | x.bRegisterM02ASL
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Use MNI as dummy Structural
Parameter    | x.bUseMNIasDummyStructural
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a

----
#### Masking
Parameter    | x.S.bMasking
----------   | ----------------------------
Description  | 
Example      | n/a
Required     | n/a
Default      | n/a


