
# DataParTemplate

Template for study-specific parameter tables

This template provides an overview of the table of parameters used by ExploreASL, which is very specific to the study to be processed.
Most of it is ASL quantification-related, some of it are image processing parameters also applicable to the structural module.
It allows only for minimal pipeline modifications (e.g. `x.settings.Quality`, `x.settings.DELETETEMP`) as most of the ExploreASL environment parameters are loaded through `ExploreASL_Initialize.m`.

Here we list potential data parameters. Most are optional, and are ignored if not provided.
For an example, see the `DataPar*.json` file(s) in the TestDataSet.
For legacy reasons & ease of use, these can be set up in a `DataPar*.m`,
which will be converted to the `DataPar*.json` per BIDS.
The compiled version of ExploreASL only allows for JSON input.

For the DataPar.json:
Make sure that booleans are inputted as numbers (e.g. `1` or `0`) and not as words (e.g. `true` or `false`).
Scalars can be inputted as scalars. For vectors/arrays of character arrays, we recommend to insert 
vectors of strings, e.g. `["optionA" "optionB"]` instead of `['optionA', 'optionB']`.
This is to allow for valid JSONs. The conversion is carried out internally.


## ENVIRONMENT PARAMETERS

|                                      | Description                                   | Defaults           |
| ------------------------------------ |:---------------------------------------------:|:------------------:|
| x.external.bAutomaticallyDetectFSL   | Boolean to automatically detect the FSL version if disabled, this function will try to use the system-initialized FSL and throw an error if FSL is not initialized. | OPTIONAL, DEFAULT = disabled |



## STUDY PARAMETERS

|                                      | Description                                   | Defaults           |
| ------------------------------------ |:---------------------------------------------:|:------------------:|
| x.dataset.name                       | String for the name of the study, example: `AD-study` | OPTIONAL |
| x.D.ROOT                             | Path to analysis root folder where data is stored, example: `/home/hjmutsaerts/TestDataSet` | OPTIONAL, DEFAULT = pwd |
| x.dataset.subjectRegexp              | String with regular expression for ExploreASL to find subjects by foldername, example: `^\d{3}$` for three digits | REQUIRED |
| x.dataset.exclusion                  | Cell with list of subjects to exclude, example: `{'005' '018'}`| OPTIONAL, DEFAULT = empty |
| x.SESSIONS                           | Use this to define sessions. Specific options: for FEAST: `1=crushed`, `2=not crushed`. This used to be other way around, but the crushed image registers better with the pGM image. Example (`'.json' file): ["ASL_1","ASL_2"]` | OPTIONAL, DEFAULT = `{'ASL_1'}` |
| x.session.options                    | This is how the sessions will be called, example: `{'baseline' 'drug'}`. For FEAST, this should be `{'non-crushed' 'crushed'}`. | OPTIONAL |
| x.dataset.ForceInclusionList | Use this field if you want to use a selection of subjects rather than taking all available subjects from directories. Example: `load(fullfile(x.D.ROOT,'LongitudinalList.mat')`). | OPTIONAL, DEFAULT = use all subjects |



## M0 PARAMETERS and OPTIONS

|                                      | Description                                   | Defaults           |
| ------------------------------------ |:---------------------------------------------:|:------------------:|
| x.settings.M0_conventionalProcessing | Boolean - use the conventional M0 processing (per consensus paper), options: 1 = standard processing, 0 = new image processing (improved masking & smoothing). | OPTIONAL, DEFAULT = 0 |
| x.M0                                 | Choose which M0 option to use: `'separate_scan'` = for a separate M0 NIfTI (needs to be in the same folder called `M0.nii`), `3.7394*10^6` = single M0 value to use, `'UseControlAsM0'` = will copy the mean control image as M0.nii and process as if it was a separately acquired M0 image (taking TR etc from the `ASL4D.nii`). Make sure that no background suppression was used, otherwise this option is invalid. | REQUIRED |
| x.M0_GMScaleFactor                   | Add additional scale factor to multiply the M0 image by This can be useful when you have background suppression but no control/M0 image without background suppression. If you then know the M0 scalefactor for the GM, you can use the control image as M0 and use this parameter to scale back what was suppressed by background suppression. Note that there is no option for separate tissue scaling (e.g. WM & GM), because ExploreASL pragmatically smooths the M0 a lot, assuming that head motion and registration between M0 & ASL4D will differ between patients and controls. | OPTIONAL, default = 1 |
| x.M0PositionInASL4D                  | A vector of integers that indicates the position of M0 in TimeSeries, if it is integrated by the Vendor in the DICOM export. Will move this from ASL4D.nii to M0.nii Note that the x.M0PositionInASL4D parameter is independent from the x.M0 parameter choice. Example for Philips 3D GRASE = '[1 2]' (first control-label pair). Example for Siemens 3D GRASE = 1 first image. Example for GE 3D spiral = 2 where first image is PWI & last = M0. Empty vector should be given (= [] or = null (in JSON)) if no action is to be taken and nothing is removed. | OPTIONAL, DEFAULT = `[] (no M0 in timeseries)` |
| x.DummyScanPositionInASL4D           | A vector of integers that indicates the position of Dummy scans in TimeSeries if they are integrated by the Vendor in the DICOM export. This allows to remove the dummy scans or noise scans that are part of the Timeseries. A new ASL4D.nii is saved with dummy scans removed and the original is backed-up. Works in a similar way as M0PositionInASL4D, both can be entered at the same time and both indicate the original position in the Timeseries independend of each other. Example for Siemens 2D EPI = `[79 80]` Skip the control-label pair used for noise measurements. Example for certain Siemens 3D GRASE = 2 Skip the first dummy control image. Empty vector should be given (= [] or = null (in JSON)) if no action is to be taken and nothing is removed. | OPTIONAL, DEFAULT = `[] (no M0 in timeseries)` |



## SEQUENCE PARAMETERS

|                                       | Description                                   | Defaults           |
| ------------------------------------- |:---------------------------------------------:|:------------------:|
| x.Q.BackgroundSuppressionNumberPulses | Used to estimate decrease of labeling efficiency. Options: 0 = (no background suppression), 2 = labeling efficiency factor `0.83` (e.g. Philips 2D EPI & Siemens 3D GRASE), 4 = labeling efficiency factor `0.81` (e.g. Philips 3D GRASE), 5 = labeling efficiency factor `0.75` (e.g. GE 3D spiral). | REQUIRED |
| x.Q.BackgroundSuppressionPulseTime    | Vector containing timing, in ms, of the background suppression pulses before the start of the readout (per BIDS). | REQUIRED when x.Q.UseControlAsM0 & x.Q.BackgroundSuppressionNumberPulses>0 |
| x.Q.PresaturationTime                 | Time in ms before the start of the readout, scalar, when the slice has been saturated (90 degree flip) this has to come before all the bSup pulses, but doesn't need to be always specified. | OPTIONAL, defaults to PLD (PASL) or PLD+LabDur ((P)CASL) |
| x.Q.readoutDim                        | String specifying the readout type. Options: `'2D'` for slice-wise readout, `'3D'` for volumetric readout. | REQUIRED |
| x.Q.Vendor                            | String containing the Vendor used. This parameter is used to apply the Vendor-specific scale factors, options: 'GE_product', 'GE_WIP', 'Philips', 'Siemens'. | REQUIRED for ASL |
| x.Q.Sequence                          | String containing the sequence used. Options: `'3D_spiral', '3D_GRASE', '2D_EPI'`. | REQUIRED for ASL |
| x.Q.LabelingType                      | String containing the labeling strategy used. Options: `'PASL'` (pulsed Q2-TIPS), `'CASL'` (CASL/PCASL). Note: pulsed without Q2TIPS cannot be reliably quantified because the bolus width cannot be identified CASL & PCASL are both continuous ASL methods, identical quantification. | REQUIRED for ASL |
| x.Q.Initial_PLD                       | Value of PLD (ms), for 3D this is fixed for whole brain, for 2D this is the PLD of first acquired slice, example: 1800. | REQUIRED for ASL |
| x.Q.LabelingDuration                  | Value of labeling duration (ms), example: 1800. | REQUIRED for ASL |
| x.Q.SliceReadoutTime                  | Value (ms) of time added to the PLD after reading out each slice, example: 31. Other option = `'shortestTR'`; shortest TR enabled gives each sequence the minimal TR. This enables calculating slice delay per subject. | REQUIRED for 2D ASL sequences |



## QUANTIFICATION PARAMETERS

|                                       | Description                                   | Defaults           |
| ------------------------------------- |:---------------------------------------------:|:------------------:|
| x.Q.bUseBasilQuantification           | True for using BASIL quantification in addition to ExploreASL's quantification. |  |
| x.Q.Lambda                            | Brain/blood water coefficient (mL 1H/ mL blood). Example: `0.32` (for GSP phantom). | OPTIONAL, DEFAULT = 0.9 |
| x.Q.T2art                             | `T2*` of arterial blood at 3T, only used when no M0 image (ms). | OPTIONAL, DEFAULT = 50 |
| x.Q.BloodT1                           | T1 relaxation time of arterial blood (ms). Defaults (Alsop MRM 2014), 1800 for GSP phantom. | OPTIONAL, DEFAULT = 1650 @ 3T, 1350 @ 1.5 T |
| x.Q.TissueT1                          | T1 relaxation time of GM tissue (ms). Defaults (Alsop MRM 2014). | OPTIONAL, DEFAULT=1240 @ 3T, 920 @ 1.5 T |
| x.Q.nCompartments                     | Number of modeled compartments for quantification. Options: 1 = a single-compartment quantification model (default by concensus paper), 2 = a dual-compartment quantification model. | OPTIONAL, DEFAULT = 1) |
| x.Q.ApplyQuantification               | A vector of 1x5 logical values specifying which types on quantified images should be calculated and saved. Fields: **1)** Apply ScaleSlopes ASL4D (xASL_wrp_Quantify, future at dcm2niiX stage), **2)** Apply ScaleSlopes M0 (xASL_quant_M0, future at dcm2niiX stage), **3)** Convert PWI a.u. to label (xASL_wrp_Quantify, future at xASL_wrp_Reslice?), **4)** Quantify M0 a.u. (xASL_quant_M0, corrects for incomplete T1 relaxation), **5)** Perform division by M0. Examples: ASL4D is an already quantified CBF image, disable all quantification `'[0 0 0 0 0]'`. To compare label but not CBF (e.g. label in vessels or sinus vs tissue): `[1 1 1 1 0]'`. Note that the output always goes to CBF.nii. | OPTIONAL, DEFAULT = `'[1 1 1 1 1]'` = all enabled |
| x.Q.SaveCBF4D                         | Boolean, true to also save 4D CBF timeseries, if ASL4D had timeseries. | OPTIONAL, DEFAULT=false |


## GENERAL PROCESSING PARAMETERS

|                                       | Description                                   | Defaults           |
| ------------------------------------- |:---------------------------------------------:|:------------------:|
| x.settings.Quality                    | Boolean specifying on which quality the pipeline should be run, options: `1` = normal quality, `0` = lower quality, fewer iterations and lower resolution of processing for a fast try-out. | OPTIONAL, DEFAULT = 1 |
| x.settings.DELETETEMP                 | Boolean for removing the temporary files. Options: `0` = keeping all files, `1` = delete temporary files created by the pipeline. | OPTIONAL, DEFAULT = 1 |
| x.settings.SkipIfNoFlair              | Boolean to skip processing of subjects that do not have a FLAIR image. These parameters can be useful when some data is still complete, but one would like to start image processing already. Options: `1` = skip processing of a subject that does not have a FLAIR image `0` = do not skip anything. | OPTIONAL, DEFAULT = 0 |
| x.settings.SkipIfNoASL                | Boolean to skip processing of subjects that do not have a ASL image. Options: `1` = skip processing of a subject that does not have a ASL image, `0` = do not skip anything. | OPTIONAL, DEFAULT = 0 |
| x.settings.SkipIfNoM0                 | Boolean to skip processing of subjects that do not have a M0 image. Options:  `1` = skip processing of a subject that does not have a M0 image, `0` = do not skip anything. | OPTIONAL, DEFAULT = 0 |




## STRUCTURAL PROCESSING PARAMETERS

|                                       | Description                                   | Defaults           |
| ------------------------------------- |:---------------------------------------------:|:------------------:|
| x.modules.bRunLongReg                 | Run longitudinal registration. | OPTIONAL, DEFAULT = 0 |
| x.modules.bRunDARTEL                  | Run between-subject registration/create templates. | OPTIONAL, DEFAULT = 0 |
| x.modules.structural.bSegmentSPM12    | Boolean to specify if SPM12 segmentation is run instead of CAT12. Options: 1 = run SPM12, 0 = run CAT12. | OPTIONAL, DEFAULT = 0 |
| x.modules.structural.bHammersCAT12    | Boolean specifying if CAT12 should provide Hammers volumetric ROI results. | OPTIONAL, DEFAULT = 0 |
| x.modules.structural.bFixResolution   | Resample to a resolution that CAT12 accepts. | OPTIONAL, DEFAULT=false |



## ASL PROCESSING PARAMETERS

|                                        | Description                                   | Defaults           |
| -------------------------------------- |:---------------------------------------------:|:------------------:|
| x.modules.asl.motionCorrection         | Boolean to perform motion correction in case of timeseries. Options: `1` = on, `0` = off. | OPTIONAL, DEFAULT = 1 |
| x.modules.asl.SpikeRemovalThreshold    | Minimal t-stat improval needed to remove motion spikes. Examples: `1` = effectively disabling spike removal. | OPTIONAL, DEFAULT = 0.01 |
| x.modules.asl.bRegistrationContrast    | Specifies the image contrast used for registration: `0` = Control->T1w, `1` = CBF->pseudoCBF from template/pGM+pWM (skip if sCoV>0.667), `2` = automatic (mix of both), `3` = option 2 & force CBF->pseudoCBF irrespective of sCoV. | OPTIONAL, DEFAULT = 2 |
| x.modules.asl.bAffineRegistration      | Specifies if the ASL-T1w rigid-body registration is followed up by an affine registration: `0` = affine registration disabled, `1` = affine registration enabled, `2` = affine registration automatically chosen based on spatial CoV of PWI. | OPTIONAL, DEFAULT = 0 |
| x.modules.asl.bDCTRegistration         | Specifies if to include the DCT registration on top of Affine, all other requirements for affine are thus also taken into account the x.modules.asl.bAffineRegistration must be >0 for DCT to run: `0` = DCT registration disabled `1` = DCT registration enabled if affine enabled and conditions for affine passed, `2` = DCT enabled as above, but use PVC on top of it to get the local intensity scaling right. | OPTIONAL, DEFAULT = 0 |
| x.modules.asl.bRegisterM02ASL          | Boolean specifying whether M0 is registered to mean_control image (or T1w if no control image exists). It can be useful to disable M0 registration if the ASL registration is done based on the M0, and little motion is expected between the M0 and ASL acquisition. If no separate M0 image is available, this parameter will have no effect. This option is disabled automatically for 3D spiral: `0` = M0 registration disabled, `1` = M0 registration enabled (DEFAULT). | OPTIONAL, DEFAULT = 0 |
| x.modules.asl.bUseMNIasDummyStructural | When structural (e.g. T1w) data is missing, copy population-average MNI templates as dummy structural templates. With this option, the ASL module copies the structural templates to fool the pipeline, resulting in ASL registration to these templates. While the rigid-body parameters might still be found somewhat correctly, with this option it is advised to enable affine registration for ASL as well, since ASL and these dummy structural images will differ geometrically. When disabled, an error will be issued instead when the structural image are missing. `1` = enabled, `0` = disabled. | OPTIONAL, DEFAULT = 0 |
| x.modules.asl.bPVCNativeSpace          | Performs partial volume correction (PVC) in ASL native space using the GM and WM maps obtained from previously segmented T1-weighted images. Skipped with warning when those maps do not exist and are not resampled to the ASL space. PVC can take several minutes for larger scans (e.g. 128x128x30), so it is deactivated by default. `1` = enabled, `0` = disabled. | OPTIONAL, DEFAULT = 0 |
| x.modules.asl.PVCNativeSpaceKernel     | Kernel size for the ASL native space PVC. This is ignored when x.modules.asl.bPVCNativeSpace is set to 0. Equal weighting of all voxels within the kernel is assumed. 3D kernel can be used, but any of the dimension can be also set to 1. Only odd number of voxels can be used in each dimension (e.g. `[3 7 5]` not `[2 3 1]`). | OPTIONAL, DEFAULT = `[5 5 1]` for bPVCGaussianMM==0, `[10 10 4]` for bPVCGaussianMM==1 |
| x.modules.asl.bPVCGaussianMM           | If set to 1, PV-correction with a Gaussian weighting is used instead of the equal weights of all voxels in the kernel ('flat' kernel) as per Asllani's original method. Ignored when x.modules.asl.bPVCNativeSpace is set to 0. Unlike with the flat kernel when the size is defined in voxels, here the FWHM of the Gaussian in mm is defined in each dimension. The advantage is twofold - continuous values can be added and a single value can be entered which is valid for datasets with different voxel-sizes without having a kernel of different effective size.`1` = enabled, use Gaussian kernel with FWHM in mm given in PVCNativeSpaceKernel, `0` = disabled, use 'flat' kernel with voxels given in PVCNativeSpaceKernel. | OPTIONAL, DEFAULT = 0 |
| x.modules.asl.bMakeNIfTI4DICOM         | Boolean to output CBF native space maps resampled and/or registered to the original T1w/ASL, and contrast adapted and in 12 bit range allowing to convert the NIfTI to a DICOM file, e.g. for implementation in PACS or other DICOM archives. If set to true, an additional CBF image will be created with modifications that allow it to be easily implemented back into a DICOM for e.g. PACS: 1. Remove peak & valley signal, remove NaNs, rescale to 12 bit integers, apply original orientation (2 copies saved, with original ASL and T1w orientation). |  |



## MASKING & ATLAS PARAMETERS


|                                       | Description                                   | Defaults           |
| ------------------------------------- |:---------------------------------------------:|:------------------:|
| x.S.bMasking                          | Vector specifying if we should mask a ROI with a subject-specific mask (1 = yes, 0 = no): `[1 0 0 0]` = susceptibility mask (either population-or subject-wise), `[0 1 0 0]` = vascular mask (only subject-wise), `[0 0 1 0]` = subject-specific tissue-masking (e.g. pGM>0.5), `[0 0 0 1]` = WholeBrain masking (used as memory compression) `[0 0 0 0]` = no masking at all, `[1 1 1 1]` = apply all masks, Can also be used as boolean, where 1 = `[1 1 1 1]`, 0 = `[0 0 0 0]`. Can be useful for e.g. loading lesion masks outside the GM. | OPTIONAL, DEFAULT=1 |
| x.S.Atlases                           | Vector specifying the atlases which should be used within the population module. Default definition within the Population Module: `x.S.Atlases = {'TotalGM','DeepWM'}`. Available atlases (please check the atlas NIfTI and accompanying files for more information): **Free atlases**: `TotalGM`: Mask of the entire GM `'./External/SPMmodified/MapsAdded/TotalGM.nii'`, `TotalWM`: Mask of the entire WM `'./External/SPMmodified/MapsAdded/TotalWM.nii'`, `DeepWM`: Mask of the deep WM `'./External/SPMmodified/MapsAdded/DeepWM.nii'`, `WholeBrain`: Mask of the entire brain `'./External/SPMmodified/MapsAdded/WholeBrain.nii'`, `MNI_Structural`: MNI cortical atlas '`./External/SPMmodified/MapsAdded/MNI_Structural.nii'`, `Tatu_ACA_MCA_PCA`: Original vascular territories by Tatu et al. `'./External/SPMmodified/MapsAdded/VascularTerritories/CortVascTerritoriesTatu.nii.nii'`, `Tatu_ICA_PCA`: Tatu - only ICA and PCA `'./External/SPMmodified/MapsAdded/VascularTerritories/TatuICA_PCA.nii'`, `Tatu_ICA_L_ICA_R_PCA`: `'./External/SPMmodified/MapsAdded/VascularTerritories/LabelingTerritories.nii'`, `Tatu_ACA_MCA_PCA_Prox_Med_Dist`: Tatu separated to distal/medial/proximal of ACA/MCA/PCA `'./External/SPMmodified/MapsAdded/VascularTerritories/ATTbasedFlowTerritories.nii.nii'`, `Mindboggle_OASIS_DKT31_CMA`: Mindboggle-101 cortical atlas `'./External/Atlases/Mindboggle_OASIS_DKT31_CMA.nii.gz'`. **Free for non-commercial use only**: `HOcort_CONN`: Harvard-Oxford cortical atlas `'./External/Atlases/HOcort_CONN.nii.gz'`, `HOsub_CONN`: Harvard-Oxford subcortical atlas `'./External/Atlases/HOsub_CONN.nii.gz'`, `Hammers`: Alexander Hammers's brain atlas `'./External/Atlases/Hammers.nii.gz'`, `HammersCAT12`: Hammers atlas adapted to DARTEL template of IXI550 space `'./External/Atlases/HammersCAT12.nii'`, `Thalamus`: Harvad-Oxford thalamus atlas `'./External/Atlases/Thalamus.nii.gz'`. | OPTIONAL, DEFAULT=`{'TotalGM','DeepWM'}` |



### Matlab structure example

To create the **JSON** file of this matlab structure, you can use `spm_jsonwrite()`.

```matlab
x.dataset.name = ExampleDataSet;
x.dataset.subjectRegexp = '^Sub-\d{3}$';
x.Q.readoutDim = '2D';
x.Q.Vendor = 'Philips';
x.Q.BackgroundSuppressionNumberPulses = 2;
x.Q.LabelingType = 'CASL';
x.Q.Initial_PLD = 1525;
x.Q.LabelingDuration = 1650;
x.Q.SliceReadoutTime = 43.7647;
x.settings.Quality = 0;
x.settings.DELETETEMP = 1;
x.M0 = 'separate_scan';
x.M0PositionInASL4D = '[1 2]';
x.DummyScanPositionInASL4D = [];
```


