function LST = tbx_cfg_LST
% Configuration file for LST jobs
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience
%
% Paul Schmidt, 2015/06/
%
%#ok<*AGROW>
 
addpath(fileparts(which(mfilename)));
%_______________________________________________________________________


%------------------------------------------------------------------------

data = cfg_files;
data.tag  = 'data';
data.name = 'Images in native space';
data.help = {[...
'Select any image that is in alignment with the lesion map, i.e. T1, FLAIR ...']};
data.filter = 'image';
data.ufilter = '.*';
data.num     = [1 Inf];

data_T1 = cfg_files;
data_T1.tag  = 'data_T1';
data_T1.name = 'T1 images';
data_T1.help = {'Select T1 images.'};
data_T1.filter = 'image';
data_T1.ufilter = '.*';
data_T1.num     = [1 Inf];

data_F2 = cfg_files;
data_F2.tag  = 'data_F2';
data_F2.name = 'FLAIR images';
data_F2.help = {'Select FLAIR images.'};
data_F2.filter = 'image';
data_F2.ufilter = '.*';
data_F2.num     = [1 Inf];

data_coreg         = cfg_files;
data_coreg.tag     = 'data_coreg';
data_coreg.name    = 'Reference images (optional)';
data_coreg.val     = {{''}};
data_coreg.help    = {['Select optional reference images to which the ' ...
    'FLAIR images are coregistered prior to lesion segmentation.']};
data_coreg.filter  = 'image';
data_coreg.ufilter = '.*';
data_coreg.num     = [0 Inf];

data_lm = cfg_files;
data_lm.tag  = 'data_lm';
data_lm.name = 'Lesion maps';
data_lm.help = {'Select lesion probability maps.'};
data_lm.filter = 'image';
data_lm.ufilter = '.*';
data_lm.num     = [1 Inf];

data_plm = cfg_files;
data_plm.tag  = 'data_plm';
data_plm.name = 'Probability lesion maps';
data_plm.help = {[...
'Select lesion probability lesion maps.']};
data_plm.filter = 'image';
data_plm.ufilter = '.*';
data_plm.num     = [1 Inf];

data_ref = cfg_files;
data_ref.tag  = 'data_ref';
data_ref.name = 'Reference Segmentation';
data_ref.help = {[...
'Select binary reference segmentations. These images needed to be in the ' ...
'same folder as the corresponding probability lesion maps. The algorithm ' ...
'automatically searches the folder of each reference image for lesion ' ...
'probability maps obtained by the LGA.']};
data_ref.filter = 'image';
data_ref.ufilter = '.*';
data_ref.num     = [1 Inf];

data_img1 = cfg_files;
data_img1.tag  = 'data_img1';
data_img1.name = 'Image';
data_img1.help = {[...
'Select main image.']};
data_img1.filter = 'image';
data_img1.ufilter = '.*';
data_img1.num     = [1 1];

%{
data_over = cfg_files;
data_over.tag  = 'data_over';
data_over.name = 'Overlay images';
data_over.help = {[...
'Select overlay images.']};
data_over.filter = 'image';
data_over.ufilter = '.*';
data_over.num     = [1 Inf];
%}

data_rep = cfg_files;
data_rep.tag  = 'data_rep';
data_rep.name = 'HTML reports';
data_rep.help = {'Select HTML reports produced by LST.'};
data_rep.filter = 'any';
data_rep.ufilter = '^report';
data_rep.num     = [1 Inf];

LSTmat         = cfg_files;
LSTmat.tag     = 'LSTmat';
LSTmat.name    = 'Select LST.mat';
LSTmat.help    = {'Select the LST.mat file that contains details for the segmentations.'};
LSTmat.filter  = 'mat';
LSTmat.ufilter = '^LST\.mat$';
LSTmat.num     = [1 1];

LSTmat01         = cfg_files;
LSTmat01.tag     = 'LSTmat01';
LSTmat01.name    = 'Select LST.mat for t = 1';
LSTmat01.help    = {'Select the LST.mat file for the first time point.'};
LSTmat01.filter  = 'mat';
LSTmat01.ufilter = '^LST\.mat$';
LSTmat01.num     = [1 1];

LSTmat02         = cfg_files;
LSTmat02.tag     = 'LSTmat02';
LSTmat02.name    = 'Select LST.mat for t = 2';
LSTmat02.help    = {'Select the LST.mat file for the second time point.'};
LSTmat02.filter  = 'mat';
LSTmat02.ufilter = '^LST\.mat$';
LSTmat02.num     = [1 1];

data_long = cfg_files;
data_long.tag  = 'data_long_tmp';
data_long.name = 'Probability lesion maps for a new time point';
data_long.help = {'Select the lesion probability maps for this time point.'};
data_long.filter = 'image';
data_long.ufilter = '.*';
data_long.num     = [1 Inf];



%------------------------------------------------------------------------

initial = cfg_entry;
initial.tag  = 'initial';
initial.name = 'Initial threshold';
initial.strtype = 'e';
initial.num = [1 Inf];
initial.def  = @(val) 0.3;
initial.help = {[...
'Threshold (kappa) for the initial lesion mask. See the Schmidt et al.' ...
'(2012) for details and the manual on how to determine the optimal ' ...
'value for this threshold.']};

mrf = cfg_entry;
mrf.tag  = 'mrf';
mrf.name = 'MRF parameter';
mrf.strtype = 'e';
mrf.num = [1 1];
mrf.def  = @(x) 1;
mrf.help = {'Parameter for the Markov Random Field.'};

maxiter = cfg_entry;
maxiter.tag  = 'maxiter';
maxiter.name = 'Maximum iterations';
maxiter.strtype = 'e';
maxiter.num = [1 1];
maxiter.def  = @(val) 50; 
maxiter.help = {'Number of maximum iterations.'};

opts_lga      = cfg_branch;
opts_lga.tag  = 'opts_lga';
opts_lga.name = 'Options for lesion segmentation';
opts_lga.val  = {initial, mrf, maxiter}; % wmh, atlas, pc, te};
opts_lga.help = {[...
'The main parameter that needs to be set by the user is kappa (Initial ' ...
'threshold). In addition, the user can specify the strength of the ' ...
'Markov random field (MRF parameter) as well as the maximum number of ' ...
'iterations (Maximum iterations) for the LGA. While we do not recommend ' ...
'to change the former parameter the maximum number of iterations should ' ...
'be increased if it reaches its limit during the segmentation process. ' ...
'Finally, the user can choose if the results of the segmentation should ' ...
'be summarized in a HTML report. Although it needs a while to produce ' ...
'this report we recommend it as it makes it easier to check the results.']};

%------------------------------------------------------------------------

bin_thresh = cfg_entry;
bin_thresh.tag = 'bin_thresh';
bin_thresh.name = 'Threshold for binary lesion maps';
bin_thresh.strtype = 'e';
bin_thresh.num = [1 1];
bin_thresh.def  = @(val) .5; 
bin_thresh.help = {[...
'The threshold for computing binary lesion maps. This value should be a ' ...
'number in the interval (0,1].']};

doit = cfg_exbranch;
doit.tag = 'doit';
doit.name = 'Determination of the optimal initial threshold';
doit.val = {data_ref,bin_thresh};%{data_plm,data_ref,bin_thresh};
doit.prog = @ps_LST_doit;
doit.help = {['With this module we offer the opportunity to determine '...
    ' the optimal kappa based on the Dice coefficient, see Schmidt et al. (2012) '...
    'for details. This requires the existence of a reference segmentation '...
    'to be compared with the lesion map. This reference image is a binary '...
    'image in the space of the T1 image where a 1 indicates a lesion.']};

html_report    = cfg_menu;
html_report.tag = 'html_report';
html_report.name = 'Produce HTML report';
html_report.labels = {'no','yes'};
html_report.values = {0 1};
html_report.def = @(val) 1; 
html_report.help    = {[...
    'Should the results of this routine be summarized in a HTML report?']};

xasl_quality    = cfg_menu;
xasl_quality.tag = 'xasl_quality';
xasl_quality.name = 'Run xASL on high quality';
xasl_quality.labels = {'no','yes'};
xasl_quality.values = {0 1};
xasl_quality.def = @(val) 1; 
xasl_quality.help    = {[...
    'Runs xASL on high quality, if 0 then run at faster settings with lower resolution and fewer iterations']};

%------------------------------------------------------------------------

lga        = cfg_exbranch;
lga.tag    = 'lga';
lga.name   = 'LST: Lesion segmentation (LGA)';
lga.val    = {data_T1,data_F2,opts_lga,html_report,xasl_quality};
lga.prog   = @ps_LST_lga;
lga.vout   = @vout;
lga.help   = {['Lesion segmentation by a lesion growth algorithm (LGA). ' ...
'This routine requires a T1 and a FLAIR image. Furthermore, the user has ' ...
'to specify an initial threshold (kappa). See Schmidt et al. (2012) for ' ...
'details. This algorithm produces lesion probability maps (ples...), ' ...
'coregistered bias corrected versions of the FLAIR inputs, a .mat-file ' ...
'with information about the segmentation that is needed for a re-run of ' ...
'the algorithm, and a HTML report along with a folder if this option has ' ...
'been chosen desired.']};

lpa        = cfg_exbranch;
lpa.tag    = 'lpa';
lpa.name   = 'LST: Lesion segmentation (LPA)';
lpa.val    = {data_F2, data_coreg, html_report,xasl_quality};
lpa.prog   = @ps_LST_lpa;
lpa.vout   = @vout;
lpa.help   = {['Lesion segmentation by the LPA requires a FLAIR image ' ...
    'only. However, the user is free to choose an additional image that ' ...
    'serves as a reference image during a coregistration step before ' ...
    'the main lesion segmentation. This may be helpful if the dimension ' ...
    'of the FLAIR image is low or if the goal of the lesion segmentation ' ...
    'is to fill lesions in T1 images. Beside that no additional parameter ' ...
    'needs to be set.']};

generic        = cfg_repeat;
generic.tag    = 'generic';
generic.name   = 'Data';
generic.help   = {'Create as many items as time points.'};
generic.values = {data_long};
generic.num    = [0 Inf];

long        = cfg_exbranch;
long.tag    = 'long';
long.name   = 'LST: Longitudinal lesion segmentation';
%long.val    = {data_plm1, data_plm2,html_report};
long.val    = {generic, html_report};
long.prog   = @ps_LST_long;
long.help   = {[...
'The longitudinal pipeline is able to compare segmented lesion maps ' ...
'for different time points. The segmentations can be derived either by ' ...
'the LGA or the LPA, but not both. Furthermore, it is required that ' ...
'the FLAIR images and .mat-files that are saved during the lesion ' ...
'segmentation process are available in the same folders as the lesions ' ...
'probability maps. The pipeline proceeds by comparing all consecutive ' ...
'time points in an iterative manner. It decides if changes in lesion ' ...
'structure are significant or due to natural variations of the FLAIR ' ...
'signal. Non-significant changes are labeled as lesions in both ' ...
'probability maps, thus, probability lesion maps are corrected within ' ...
'this procedure and may differ from the ones that served as input. As ' ... 
'a final result, lesion change labels are produced for all consecutive ' ...
'time points. In these images the three possible cases decrease, no ' ...
'change and increase are labeled by the numbers 1, 2, and 3, ' ...
'repsectively. In addition, a lesion change plot is constructed. This ' ...
'plot shows the lesion volumes for both time points of all segmented ' ...
'lesions. In this way it is easy to recognize how the lesions structure ' ...
'has been changed, i.e. if the change occured by the appearence of new ' ...
'lesions, the disappearence of old lesions, by the change of already ' ...
'existing lesions, or a combination of these possibilities.']};

merge_reports        = cfg_exbranch;
merge_reports.tag    = 'merge_reports';
merge_reports.name   = 'LST: Merge HTML reports';
merge_reports.val    = {data_rep};
merge_reports.prog   = @ps_LST_merge_reports;
merge_reports.help   = {[...
'This function allows to merge multiple HTML reports that have been ' ...
'obtained by different functions of this toolbox. The resulting report ' ...
'can be moved anywhere on your desk as long as it is provided that the ' ...
'folders for the original reports stay where they are. A report that ' ...
'is readable on different platforms can be obtained by exporting/' ...
'printing the HTML report as a PDF document, a function that is ' ...
'included in all modern browsers.']};

filling = cfg_exbranch;
filling.tag = 'filling';
filling.name = 'Lesion filling';
filling.val = {data, data_plm, html_report};
filling.prog = @ps_LST_lesfill;
filling.vout = @vout_lesfill;
filling.help = {['Lesion filling can be applied to any image that is in ' ...
    'alignment with the lesion probability map. However, it is required ' ...
    'that the .mat-files that are saved during the lesion segmentation ' ...
    'process are available in the same folders as the lesions probability maps.']};

thresholding = cfg_exbranch;
thresholding.tag = 'thresholding';
thresholding.name = 'Thresholding of lesion probabilities';
thresholding.val = {data_plm, bin_thresh};
thresholding.prog = @ps_LST_thresholding;
thresholding.vout = @vout_lesthresh;
thresholding.help = {['This module can be used to threshold probability ' ...
    'lesion maps in order to compute binary lesion maps.']};

%{
gif = cfg_exbranch;
gif.tag = 'gif';
gif.name = 'Create animated GIF and PNGs.';
gif.val = {data_img1, data_over};
gif.prog = @ps_LST_create_gif;
%other.vout = @vout_lesthresh;
gif.help = {['This module can be used to threshold probability ' ...
    'lesion maps to compute binary lesion maps. ... Make sure you have ' ...
    'enough space on your hard disk.']};
%}
    
tlv = cfg_exbranch;
tlv.tag = 'tlv';
tlv.name = 'Extract values of interest';
tlv.val = {data_lm,bin_thresh};
tlv.prog = @ps_LST_tlv;
tlv.help = {['computes the total lesion volume (TLV) and the number of ' ...
    'lesions for a given set of lesion probability maps. The result ' ...
    'of this function is a CSV-file with the name LST_tlv_[date].csv' ...
    'with the following four columns: Path (path of image), ' ...
    'FileName (name of image), LGA (dummy that indicate weather the ' ...
    'lesion map was obtained by LGA (1) or LPA (0)), TLV (total lesion ' ...
    'volume in ml), and N (Number of lesions).']};

%------------------------------------------------------------------------

LST        = cfg_choice;
LST.name   = 'LST';
LST.tag    = 'LST';
LST.values = {lga,doit,lpa,long,merge_reports,filling,thresholding,tlv};

%------------------------------------------------------------------------

function dep = vout(job)

%opts  = job.output_segment;
cdep = cfg_dep;

cdep(end).sname      = 'Lesion probability map';
cdep(end).src_output = substruct('.','lpm','()',{':'});
cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
%{
if opts.ples
    cdep(end).sname      = 'Lesion probability map';
    cdep(end).src_output = substruct('.','lpm','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

if opts.wples
    cdep(end).sname      = 'Normalized lesion probability map';
    cdep(end).src_output = substruct('.','wlpm','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});        
end
%}
dep = cdep;

function dep = vout_lesfill(job)

cdep = cfg_dep;

cdep(end).sname      = 'Filled images in native space';
cdep(end).src_output = substruct('.','Img_filled','()',{':'});
cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
dep = cdep;

function dep = vout_lesthresh(job)

cdep = cfg_dep;

cdep(end).sname      = 'Binary lesion maps';
cdep(end).src_output = substruct('.','bles','()',{':'});
cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
dep = cdep;
