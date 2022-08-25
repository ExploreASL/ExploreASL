//External/SPMmodified is a modified version of SPM12, CAT12 and LST, which is redistributed with the same license (GPL, see license files inside folders). This was developed parallel with ExploreASL, but is a separate project, with its own license.

In this textfile we list the modifications. This concerns the following versions:
SPM12, r7219 (/External/SPMmodified/Contents.txt)
CAT12, r1615 (/External/SPMmodified/toolbox/cat12/CHANGES.txt)
LST, r2.0.15 (/External/SPMmodified/toolbox/LST/lst-version.txt)

Aside from the below list of code modifications,
we have added new maps (//External/SPMmodified/MapsAdded/README_Maps.txt)
as well as new code (//External/SPMmodified/xASL/README_Functions.txt)

*********************************************************************************************************
MAP CHANGES
DATE+Name: 2021-07-07 HM
DESCRIPTION:
/External/SPMmodified/toolbox/cat12/templates_volumes/cat_wmh_miccai2017.nii
Was cleaned up (lower slices contained some messy low probabilities)

*********************************************************************************************************
ASL-SPECIFIC HACKS

DATE+Name: 2019-11-25 JP
DESCRIPTION:
Merged the spm_realign and the xASL version with Zig-zag. Now the zig-zag is added as an extra third option
FILE:
spm_realign.m
ADDED:
Line 1: added the extra option bZigzag
Line 117,137,158,181,571: added extra routines for the Zigzag

*********************************************************************************************************
COST FUNCTION MASKING FOR LESIONS

DATE+NAME:2020-07-08, JP
DESCRIPTION:
hotfix minor bug in loading NIfTIs containing lesion masks in CAT12 #28
FILE:
cat_run_job.m at 577

*********************************************************************************************************
ENABLING LOW QUALITY MODE (FOR QUICK TESTING, RUN EVERYTHING BUT WITH LOWER ITERATIONS AND/OR SPATIAL RESOLUTION)

DATE+NAME:2018_12_29, HM
DESCRIPTION:
Speeding up low quality option, by reduction of iterations & resolution. All steps are the same as quality=1, but simply lower quality.
Useful for performing multiple reproducibility tests.
FILE:
External/SPMmodified/toolbox/cat12/cat_conf_extopts.m
525     : Added xasl_quality as a valid input for CAT12
725,727,762 : Added xasl_quality as a valid input for CAT12
External/SPMmodified/toolbox/cat12/cat_main.m:
905      : set xasl_quality on default to 1
207      : if job.extopts.xasl_quality is on 0, run Geodesic Shooting on fewer iterations & lower resolution
1446    : if job.extopts.xasl_quality is on 0, run Geodesic Shooting on fewer iterations & lower resolution
External/SPMmodified/toolbox/cat12/cat_run_job.m:
126     : if job.extopts.xasl_quality is on 0, increase segmentation sampling distance from 3 to 6
215     : if job.extopts.xasl_quality is on 0, decrease denoising strenght
302     : if job.extopts.xasl_quality is on 0, increase segmentation voxel size to [1.5 1.5 1.5]
652 : if job.extopts.xasl_quality is on 0, increase segmentation sampling distance from 3 to 9


DATE+NAME:2018_12_29, HM
DESCRIPTION:
LST hacks to speed up low quality mode.
FILE:
External/SPMmodified/toolbox/LST/tbx_cfg_LST.m:
218: Added the xasl_quality parameter specification, description, and default value
232,247: Added xasl_quality as a valid input for LPA and LGA algorithms
External/SPMmodified/toolbox/LST/ps_LST_lga.m:
165: First loads the xasl_quality parameter from the job, or set to default == 1
312: if xasl_quality is set to 0, load ps_LST_lga_preproc_default_LowQ.m instead of ps_LST_lga_preproc_default.m
External/SPMmodified/toolbox/LST/ps_LST_lpa.m
115: First loads the xasl_quality parameter from the job, or set to default == 1
218: if xasl_quality is set to 0, load ps_LST_lpa_preproc_default_LowQ.m instead of ps_LST_lpa_preproc_default.m,
likewise if xasl_quality is set to 2 (when WMH_SEGM pre-exists), load ps_LST_lpa_preproc_default_UltraLowQ.m

*********************************************************************************************************
REDUCE CODE SIZE

DATE+Name: 2020-11-29 HM (#243)
DESCRIPTION: Avoid initializing unused toolboxes
FILE:
External/SPMmodified/config/spm_cfg.m @ several locations

DATE+Name: 2020-11-29 HM (#243)
DESCRIPTION: Move SPM templates used by CAT12 out of Fieldmap folder
FILE:
External/SPMmodified/toolbox/cat12/cat_defaults.m @ 224
(also same changes in:
External/SPMmodified/toolbox/cat12/cat_run1173plus/cat_defaults1173plus.m
External/SPMmodified/toolbox/cat12/cat_run1585/cat_defaults1585.m
External/SPMmodified/toolbox/cat12/cat_vol_groupwise_ls.m)

DATE+Name: 2020-07-04 HM
DESCRIPTION:
Disable previous CAT12 versions for increased stability
FILE:
cat_conf_output.m @ 394
cat_conf_tools.m @ 135, 168
tbx_cfg_cat.m @ 11, 103, 107, 112, 145, 239, 242, 245

DATE+NAME:2019_11_13, HM
DESCRIPTION:
Remove unused SPM stuff to reduce data size ExploreASL/SPM compilation:
//SPM/toolbox/cat12/Folders:
/atlases_surfaces
/atlases_surfaces_32k
/CAT.glnx86
/CAT.maci64
/CAT.w32
/html
/templates_surfaces
/templates_surfaces_32k

/templates_volumes:
Removed atlases, partly
left brainmask.nii, cat.nii,
& all Template_._IXI555_MNI152_(GS|)\.nii
SPM/toolbox/DARTEL/icbm152.nii
cat12/templates_volumes/TPM_Age11.5.nii (but returned by JP later?)

DATE+NAME:2019_11_13, HM
DESCRIPTION:
Remove unused SPM stuff to reduce data size ExploreASL/SPM compilation:
SPM12/batches
SPM12/canonical
SPM12/man
SPM12/tests

DATE+Name: 2019-01-16 JP
DESCRIPTION:
Removed the EEF toolbox from SPM. The EEG toolbox initialization thus had to be disabled.
FILE:
External/SPMmodified/config/spm_cfg.m
ADDED:
Line 41 - removed 'spm_cfg_dcm_meeg'
Line 132 - removed 'EEG'
Line 135 - removed 'spm_cfg_eeg'


*********************************************************************************************************
IMAGE PROCESSING IMPROVEMENT

DATE+NAME:2020_07_01, HM
DESCRIPTION: Remove feedback missing log-file (we don't use catlog_txt in xASL)
FILE:
cat_run.m @ 459
cat_run_job.m @ 79

DATE+NAME:2020_07_01, HM
DESCRIPTION: Add progress tracking
FILE:
cat_vol_imcalc.m @ 207

DATE+NAME:2020_07_01, HM
DESCRIPTION: Add progress tracking
FILE:
cat_vol_imcalc.m @ 207
cat_main_reportfig.m @ 35, 91, 159, 231, 300, 400, 485, 594, 615, 617

DATE+NAME:2020_07_01, HM
DESCRIPTION: Improved feedback on ROI creation
FILE:
cat_main_roi.m @ 50

DATE+NAME:2020_07_01, HM
DESCRIPTION: Improved feedback on PDF creation
FILE:
cat_main_reportfig.m @ 34, 622

DATE+NAME:2020_07_01, HM
DESCRIPTION: Save CSF as well in newer CAT12 version
FILE: cat_main.m @ 674

DATE+NAME:2019_11_11, HM
DESCRIPTION:
Replaced TPM.nii by enhanced TPM.nii, from Lorio 2016 Neuroimage. Have more accurate thalamus.
TPM.nii needed to be in single precision for SPM though.
FILE:
//ExploreASL/External/SPMmodified/tpm/TPM.nii

DATE+NAME:2020_03_14, HM
DESCRIPTION:
Hack to output TLV filename as output argument
FILE:
ps_LST_tlv @ 1, end

DATE+NAME:2019_09_02, HM
DESCRIPTION:
Hack to allow setting MinimalLesionSize
FILE:
ps_LST_tlv @ 46, 67-71, 170

DATE+NAME:2019_05_02, HM
DESCRIPTION:
Some bugfixes:
FILE:
cat_main.m:
set NaNs to zeroes:
915: for nT1.nii
1456: for edges of Yy transformation fields
1476: for edges of trans.warped.y transformation fields

DATE+NAME:2019_03_27, Jan Petr
DESCRIPTION:
Enabled save of CSF as c3T1.nii from CAT12
FILE:
External/SPMmodified/toolbox/cat12/cat_conf_opts.m
297     : activated the SAMP option for non-expert mode

DATE+NAME:2020_05_26, Jan Petr
DESCRIPTION:
Remove the link to the older versions like 1173, 1173plus, ,1445, 1585
FILE:
External/SPMmodified/toolbox/cat12/tbx_cfg_cat.m


DATE+NAME:2019_02_08, Jan Petr
DESCRIPTION:
Fixed the option for using and removing lesions from the CAT12 segmentation. Now the lesions are merged and prepared beforehand and the filename is passed to CAT, where this is loaded and all transformation steps are applied to it directly in CAT. Plus, CAT has some extra routines for working with the lesions, some of them were bypassed, because they were not working efficiently, when the lesion was set to NaN in the T1w image.
FILE:
External/SPMmodified/toolbox/cat12/cat_conf_extopts.m
537     : Added xasl_lesion as a valid input for CAT12 containing the name of the merged lesion file
725,727,762 : Added xasl_lesion as a valid input for CAT12
External/SPMmodified/toolbox/cat12/cat_defaults.m
197: Added the default value '' for the xasl_lesion
External/SPMmodified/toolbox/cat12/cat_main_registration.m
407,638: The resliced lesion - ls - can contain NaNs due to reslicing - these must be removed before the mask is applied.
External/SPMmodified/toolbox/cat12/cat_main_updateSPM.m
88 - We save our original lesion as LesionFull for later use, because the other lesion gets stripped (and this does not work when T1w has the lesion set to NaN by ExploreASL).
External/SPMmodified/toolbox/cat12/cat_main.m
708 - calling the registration with the Full Lesion, that is not stripped (as the stripped one does not work when NaNs were set to T1w - because the WM segmentation is then missing there).
1390 - removing the lesion also when the xasl_lesion parameter was set. Few lines below - need to add apply the same transformations also to the non-stripped lesion. Do not divide it by
       255, as this makes it virtually unusable - no idea why CAT does that.
1434 - Not loading the lesions here anymore, but instead only saving original segmentation with xASL_im_SaveOriginal4CAT.
1448 - calling the segmentation again with the non-stripped lesions.
External/SPMmodified/toolbox/cat12/cat_run_job1070.m
523 - We need to load our provided lesion, not use the CAT one (that's usually empty)

DATE+NAME:2019_01_30, Jan Petr
DESCRIPTION:
Added an option for saving the intermediate steps of the CAT12 registration at the end of the segmentation
(i.e. saving the lower quality intermediate DARTEL and GS transformation fields).
FILE:
External/SPMmodified/toolbox/cat12/cat_conf_extopts.m
537     : Added xasl_savesteps as a valid input for CAT12
725,727,762 : Added xasl_savesteps as a valid input for CAT12
External/SPMmodified/toolbox/cat12/cat_main_registration.m
459: Added a call to function xASL_wrp_DARTELSaveReg that saves the current DARTEL registration flow field to a file - this function is ripped out of this CAT file
794: Added a call to function xASL_wrp_GSSaveReg that saves the current GS registration flow field to a file (only for the full resolution one) - this function is ripped out of this CAT file



*********************************************************************************************************
IMAGE PROCESSING REPRODUCIBILITY

DATE+NAME:2020_07_03, Jan Petr
DESCRIPTION:
Reseeding the random number generator for certain functions to make sure they run reproducibly.
Removing the legacy function for that and using RNG instead. Rather than using a 'default' random number generator that might differ
potentially between Matlab version, we have used 'twister' everywhere.
FILEs:
External/SPMmodified/spm_coreg.m - line 291
External/SPMmodified/spm_maff8.m line 79
External/SPMmodified/spm_preproc.m line 132
External/SPMmodified/spm_preproc8.m line 100
External/SPMmodified/spm_realign.m line 236
External/SPMmodified/toolbox/LST/ps_LST_spm_coreg.m line 287
External/SPMmodified/toolbox/OldNorm/spm_affreg.m line 94
External/SPMmodified/toolbox/cat12/cat_run_job.m lin 733 - switch off turning of the warnings for this reseeding reasons
External/SPMmodified/toolbox/cat12/cat_spm_affreg.m line 94
External/SPMmodified/toolbox/cat12/cat_spm_preproc8.m line 110
External/SPMmodified/toolbox/cat12/cat_vol_sanlm.m line 85

DATE+NAME:2019_03_28, Jan Petr
DESCRIPTION:
Added rounding to inversions in CAT12 to avoid having reproducibility issues between OS and Matlab versions
FILE:
cat_main_registration.m
Added round(~,12) to lines 323, 333, 340, 554, 558, 686, 860 -> replaced by HM by xASL_round for backward compatibility

DATE+Name: 2018-12-11 JP
DESCRIPTION:
Rounded the results of matrix inversion at 10e-12 to avoid differences between OS and Matlab versions
FILE:
spm_realign.m
ADDED:
Line 455 - added xASL_round(xxx,12)

DATE+Name: 2018-10-08 JP
DESCRIPTION:
Use own Euclidean distance transformation instead of bwdist, disable the bwdist alternative as well.
Increase between-os reproducibility because the same option is used for all
FILE:
ps_LST_lpa.m
ADDED:
Lines 569-589 commented as this was the alternative option.
Line 637 - bwdist replaced with xASL_im_DistanceTransform

DATE+NAME: 2018-08-22, HM
FILE:
ps_LST_lga (358), ps_LST_lpa (217)
DESCRIPTION:
disable coregistration in LGA/LPA (already performed by ExploreASL), increases reproducibility

DATE+Name:2018_09_05, JP
DESCRIPTION:
replacing convn to improve between-matlab version reproducibility
FILE:
cat_vol_iscale.m (line 490-500)

DATE+NAME: 2018-08-16 Guillaume Flandin
FILE: spm_powell.m
DESCRIPTION:
117: replaced to solve issues with inter-version repeatability
several other replacements, update to version 2017

DATE+NAME:2018-08-30 JP
DESCRIPTION:
Adding our own function for convolution that would not have repeatability issues between Matlab version in the following files:
FILE:
spm_coreg.m line 193
spm_maff8.m line 193
toolbox/LST/ps_LST_spm_coreg line 178
spm_smoothto8bit.m 57
External/SPMmodified/toolbox/OldSeg/spm_maff.m line 106
toolbox/cat12/cat_vol_correct_slice_scaling line 425

*********************************************************************************************************
BIDS/JSON

DATE+NAME:2020_03_12, JP
DESCRIPTION:
In DICOM header reader, removed the extra processing for CSASeriesHeaderInfo as the Phoenix Siemens protocol is saved under that DICOM tag and need a simple conversion to string. The type is redefined as LT in the External/SPMmodified/spm_dicom_dict.txt. Moreover, new DICOM fields are added at the end. The .mat dictionary is generated using a spm_dicom_text_to_dict('spm_dicom_dict.txt') command.
FILE:
spm_dicom_header.m at line 102
External/SPMmodified/spm_dicom_dict.txt at line 3607 and 3665 till the end


DATE+NAME:2020_08_11, JP
DESCRIPTION:
Error in the NIfTI header. 3D NIfTI files (dim[0] == 3) that have the size of the fourth dimension NT==1 (dim[4] == 1), but have TR defined (pixdim[4] ~= 1)
do not pass through the BIDS validator. The SPM NIfTI writing function needs to be modified to save dim[4]==4 in this case even though the size of the fourth dimension is 1.
FILE:
External/SPMmodified/@nifti/private/write_hdr_raw.m at line 20

DATE+NAME:2020_04_29, HM
DESCRIPTION:
Instead of error for an imaginary number, put it in a string
FILE:
spm_jsonwrite.m at line 198


DATE+NAME:2020_01_17, JP
DESCRIPTION:
Edited the spm_jsonwrite so that it makes a newline after each field
FILE:
spm_jsonwrite.m at line 157

DATE+NAME:2019_10_25, JP
DESCRIPTION:
Edited the jsmn.c so that spm_jsonread.c can read JSONs that have empty string as value
FILE:
jsmn.c at line 269

DATE+NAME:2021_01_12, JP
DESCRIPTION:
Edited the spm_jsonread.c so that spm_jsonread MEX correctly reads files that contain 0-character inside the text. A character with value 0 appears in some - otherwise empty - strings.
Which causes that the JSON is not read entirely. We now read the entire JSON and replace 0 by a space.
FILE:
spm_jsonread.c at line 496

DATE+NAME:2022_08_25, JP
DESCRIPTION:
Edited the spm_jsonread.c so that spm_jsonread MEX reports that an error happened within jsonread and also print out the file name being read to be able to identify the source of error.
FILE:
spm_jsonread.c at line 588 and all lines with mexErrMsgTxt, mexWarnMsgTxt, and mexPrintf.

*********************************************************************************************************
OTHER CODE HACKS

DATE+NAME: 2021-03-10 MS (issue #408)
DESCRIPTION: Missing version_matlab field in compiled xASL.
FILE: cat_main_reportstr, 50-54

DATE+NAME:2021-02-04 HM (issue #302)
DESCRIPTION: Deactivate (comment out) calls of SPM/LST/CAT12 to server
FILE:
Cat12 at line 81
cat_update at line 22
spm_update at line 36
cat_io_send_to_server at line 19
ps_LST_update at line 5

DATE+NAME:2022-08-03 HM (#1140)
DESCRIPTION: Attempt repairing broken CSV/TSV files by removing empty cells at the end of rows
FILE: spm_load.m, 102, 162, 191

DATE+NAME:2021-01-20 HM (issue #276)
DESCRIPTION: Manage trailing \t on header only
FILE: spm_load.m, 131, 165

DATE+NAME:2020-10-26 HM (issue #190)
DESCRIPTION: Add atlas ROI creation comments and add creation of catROI_T1.tsv
FILE: cat_main.m, 127

DATE+NAME:2020-09-02 MS (issue #114)
DESCRIPTION:
Make Matlab version information robust for both deployed and undeployed mode (bugfix)
FILE:
cat_io_report.m, 176
cat_vol_qa.m, 507, 574

DATE+NAME:2020-07-06 HM (issue #1)
DESCRIPTION:
Allow running CAT12 without JVM
FILE:
cat_run at 466, 860
cat_main at 884

DATE+NAME:2020-03-29 HM
DESCRIPTION:
288 cosmetic hack, for xASL_TrackProgress
210 improvement xASL_adm_CreateDir instead of mkdir
FILE:
ps_LST_lpa.m at 210 & 288

DATE+NAME:2020-05-26 JP
DESCRIPTION:
Small fix for the boundary conditions
FILE:
replace round with ceil at cat_vol_qa - Line 627

DATE+NAME:2019_10_20, HM
DESCRIPTION:
Hack to avoid warning of absence DEM toolbox
FILE:
spm @ 4919-421

DATE+NAME:2019_10_01, HM
DESCRIPTION:
Hack to insert comments that explain WMH volumetrics
FILE:
ps_LST_tlv @ 164-174

DATE+Name: 2018-09-14 JP
DESCRIPTION:
overrides the settings of the CAT12 log-files to continue writing to the xASL log-file
FILE:
cat_io_report line 41
cat_main line 2629
cat_run_job line 69
cat_run_job1070 line 66
ADDED:
Commented out the lines with 'diary'

DATE+NAME:2018_spring, HM
DESCRIPTION:
Changes to CAT12
FILE:
cat_main.m
18  : debug on
91  : insert disableDARTEL option
315 :
1443: Cost Function Masking -> if Lesion_*.nii are found, these are masked out
1665: here create flow field + DARTEL
2155: Disable GUI output

DATE+NAME:2019_02_13, JP
DESCRIPTION:
Changes to CAT12 to disable DARTEL
FILE:
External/SPMmodified/toolbox/cat12/cat_conf_extopts.m
568: the xasl_disabledartel option is introduced
cat_main.m
107  : insert xasl_disabledartel option

DATE+NAME:2018_09_28, HM
DESCRIPTION:
Changes to LST
FILE:
ps_LST_lpa.m
498: if WMH_SEGM.nii already exist (e.g. from third-party/manual), load this,
     and perform the CleanUp on this pre-existing WMH_SEGM.nii

DATE+Name: 2019-01-24 Michael Stritt
DESCRIPTION:
There are four files that are created due to spm_make_standalone / ExploreASL_make_standalone.
FOLDER:
External/SPMmodified
FILES:
cfg_mlbatch_appcfg_1.m
cfg_mlbatch_appcfg_2.m
cfg_mlbatch_appcfg_master.m
spm_cfg_static_tools.m

DATE+NAME: 2018-08-22 HM
DESCRIPTION:
and some small fixes
FILE:
ps_LST_lga (249)
ADDED:
cd_tmp  = xASL_adm_ConvertSlash(cd_tmp,1);
if ~isdir(tmpFolder)
    mkdir(tmpFolder);
end
FILE:
ps_LST_lga (290)
ADDED:
if ~isdir(tmpFolder)
    mkdir(tmpFolder);
end
FILE:
ps_LST_lpa (183)
ADDED:
cd_tmp = xASL_adm_ConvertSlash( cd_tmp, 1);

DATE+NAME:
2018_spring HM
FILE:
ps_LST_lpa.m
DESCRIPTION:
184: flipping bugfix
698: comments added, which were here in removed/unused function ps_LST_lpa_2.m. otherwise, these functions were identical
754: bugfix slashes (per ConvertSlash.m)



*********************************************************************************************************
MATLAB BACKWARDS COMPATIBILITY

DATE+NAME:2019_05_06, HM
DESCRIPTION:
Fixing rounding backward compatibility
(replacing round by xASL_round)
FILE:
cat_main_register.m, lines 341, 554, 560, 687, 861, 865, 1168






*********************************************************************************************************
REDUCE GRAPHICAL OUTPUT & NON-SPECIFIC WARNINGS. ADD TRACKING PROGRESS AT COMMAND LINE: xASL_TrackProgress

DATE+NAME:2020_06_22, HM
DESCRIPTION:
No need to mention subfolders CAT12, xASL moves files from here
FILE:
cat_main_reportcmd @ 33

DATE+NAME:2020_06_22, HM
DESCRIPTION:
Skip SPM progress bar visualization (use xASL_TrackProgress for CLI progress instead)
FILE:
spm_progress_bar @ 26

DATE+NAME:2020_06_22, HM
DESCRIPTION:
Print instead of warn when nojvm for storing XML
FILE:
cat_io_xml @ 51

DATE+NAME:2019_10_13, HM
DESCRIPTION:
Hack to reinforce command line output when SPM warning/error messages occur
FILE:
spm @ 1156-1157

DATE+NAME:2019_10_10, HM
DESCRIPTION:
Hack to remove redundant output figures with registration
FILE:
spm_coreg @ 169-171
ps_LST_spm_coreg @ 153-155

DATE+NAME:
2018_spring
FILE:
spm_defaults.m:
DESCRIPTION:
73: increase maxmem
75: set resmem to true
OLD:
defaults.stats.maxmem      = 2^29;
defaults.stats.resmem      = false;
NEW:
defaults.stats.maxmem      = 2^31;
defaults.stats.resmem      = true;


DATE+NAME: 2018_spring, HM
FILE:
spm_chi2_plot.m (22)
spm_atranspa.m (18)
spm_imcalc.m (113)
spm_preproc8.m (602)
DESCRIPTION:
Disable warnings:
OLD:
warning
NEW:
disp
NOTE:
spm_imcalc.m only comments out the warning

DATE+NAME:2018-09-04 HM
DESCRIPTION:
Turned warnings off temporarily
FILE:
cat_run_job.m (64, 67)
cat_run_job1070.m (63, 65)

DATE+NAME:2018_spring, HM
DESCRIPTION:
Disable Figure output SPM:
FILE:spm_DesRep.m (702-703, 877)

DATE+NAME:2018_spring, HJ
FILE:
cfg_util.m
DESCRIPTION:
1641, 1686: shorter screenprint

DATE+NAME: 2018 HM
DESCRIPTION:
xASL_TrackProgress Edits -> count percentage completion without graphical output: function (lines)
And replace progress bars
FILE:
spm_groupwise_ls (143,145)
spm_coreg (155, 159)
spm_preproc8 (167, 308)
spm_realign (294)
spm_deformations (409)
cat_vol_sanlm (345)
cat_main (231, 633, 641, 649, 2654, 2674, 2693)
cat_main_gintnorm (62, 150, 363, 390)
ADDED:
xASL_TrackProgress

DATE+NAME:2018_spring, HM
DESCRIPTION:
Handle screen messages: reduces nr of messages, more ExploreASL-specific, replaced by xASL_TrackProgress
FILE:
cat_run_job.m & cat_run_job1070 (163, 181)
cat_main_gintnorm.m (388)
cat_main_updateSPM.m (303-323)
cat_main.m (315, 332, 488,  1035, 1903, 2652, 2668, 2690)
cat_spm_preproc_write8.m & spm_preproc_write8.m (190, 192, 291, 293, 304, 306, 609, 611, 627, 629, 667)
spm_proc8 (196)

DATE+NAME:2018_09_04, HM
DESCRIPTION:
Handle screen messages: reduces nr of messages, more ExploreASL-specific, replaced by xASL_TrackProgress
FILE:
spm_preproc8.m: line 196
spm_preproc_write8.m & cat_spm_preproc_write8.m: lines 178, 180, 275, 277, 287, 289, 588, 590, 606, 608, 646 (lines shown for spm_preproc_write8.m, but are same scripts)
cat_run_job.m (lines 355
cat_spm_preproc_run (lines 140)
cat_run.m (lines 393)
cat_main.m (lines 516, 525, 533, 1763, 1777)
cat_run_job.m & cat_run_job1070.m:  remove warning line 64
