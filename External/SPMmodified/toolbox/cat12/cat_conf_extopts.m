function extopts = cat_conf_extopts(expert,spm)
% Configuration file for extended CAT options
%
% Christian Gaser
% $Id: cat_conf_extopts.m 1354 2018-08-13 10:21:10Z gaser $
%#ok<*AGROW>

if ~exist('expert','var')
  expert = 0; % switch to de/activate further GUI options
end

if ~exist('spm','var')
  spm = 0; % SPM segmentation input
end

%_______________________________________________________________________

% options for output
%-----------------------------------------------------------------------

vox         = cfg_entry;
vox.tag     = 'vox';
vox.name    = 'Voxel size for normalized images';
vox.strtype = 'r';
vox.num     = [1 1];
vox.def     = @(val)cat_get_defaults('extopts.vox', val{:});
vox.help    = {
  'The (isotropic) voxel sizes of any spatially normalised written images. A non-finite value will be replaced by the average voxel size of the tissue probability maps used by the segmentation.'
  ''
};
if expert > 1
  vox.num  = [1 inf];
  vox.help = [vox.help; { 
    'Developer option: '
    '  For multiple values the first value is used for the final output, whereas results for the other values are saved in separate sub directories. '
    ''
    }];
end

%------------------------------------------------------------------------
% SPM, Dartel, Shooting Template Maps
% e.g. for other species
%------------------------------------------------------------------------

darteltpm         = cfg_files;
darteltpm.tag     = 'darteltpm';
darteltpm.name    = 'Dartel Template';
darteltpm.def     = @(val)cat_get_defaults('extopts.darteltpm', val{:});
darteltpm.num     = [1 1];
darteltpm.filter  = 'image';
darteltpm.ufilter = 'Template_1'; 
darteltpm.help    = {
  'Select the first of six images (iterations) of a Dartel template.  The Dartel template must be in multi-volume (5D) nifti format and should contain GM and WM segmentations. '
  ''
  'Please note that the use of an own Dartel template will result in deviations and unreliable results for any ROI-based estimations because the atlases will differ and any ROI processing will be therefore deselected.'
  ''
};

%---------------------------------------------------------------------

shootingtpm         = cfg_files;
shootingtpm.tag     = 'shootingtpm';
shootingtpm.name    = 'Shooting Template';
shootingtpm.def     = @(val)cat_get_defaults('extopts.shootingtpm', val{:});
shootingtpm.num     = [1 1];
shootingtpm.filter  = 'image';
shootingtpm.ufilter = 'Template_0'; 
shootingtpm.help    = {
  'Select the first of five images (iterations) of a Shooting template.  The Shooting template must be in multi-volume (5D) nifti format and should contain GM, WM, and background segmentations and have to be saved with at least 16 bit. '
  ''
  'Please note that the use of an own Shooting template will result in deviations and unreliable results for any ROI-based estimations because the atlases will differ and any ROI processing will be therefore deselected.'
  ''
};

%------------------------------------------------------------------------

cat12atlas         = cfg_files;
cat12atlas.tag     = 'cat12atlas';
cat12atlas.name    = 'CAT12 ROI atlas';
cat12atlas.filter  = 'image';
cat12atlas.ufilter = 'cat';
cat12atlas.def     = @(val)cat_get_defaults('extopts.cat12atlas', val{:});
cat12atlas.num     = [1 1];
cat12atlas.help    = {
  'CAT12 atlas file to handle major regions.'
};

%------------------------------------------------------------------------

brainmask         = cfg_files;
brainmask.tag     = 'brainmask';
brainmask.name    = 'Brainmask';
brainmask.filter  = 'image';
brainmask.ufilter = 'brainmask';
brainmask.def     = @(val)cat_get_defaults('extopts.brainmask', val{:});
brainmask.num     = [1 1];
brainmask.help    = {
  'Initial brainmask.'
};

%------------------------------------------------------------------------

T1         = cfg_files;
T1.tag     = 'T1';
T1.name    = 'T1';
T1.filter  = 'image';
T1.ufilter = 'T1';
T1.def     = @(val)cat_get_defaults('extopts.T1', val{:});
T1.num     = [1 1];
T1.help    = {
  'Affine registration template.'
};

%---------------------------------------------------------------------

% removed 20161121 because it did not work and their was no reason to use it in the last 5 years
% however it is maybe interesting to create own templates
%{
bb         = cfg_entry;
bb.tag     = 'bb';
bb.name    = 'Bounding box';
bb.strtype = 'r';
bb.num     = [2 3];
bb.def     = @(val)cat_get_defaults('extopts.bb', val{:});
bb.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).'
''
};
%}

%---------------------------------------------------------------------

if expert==0
  regstr        = cfg_menu;
  regstr.labels = {
    'Dartel'
    'Default Shooting'
    'Optimized Shooting'
  };
  regstr.values = {0 4 0.5};
  regstr.help   = [regstr.help; { ...
    'For spatial registration CAT offers the use of the Dartel (Ashburner, 2008) and Shooting (Ashburner, 2011) registrations to an existing template. Furthermore, an optimized shooting approach is available that uses an adaptive threshold and lower initial resolutions to obtain a good tradeoff between accuracy and calculation time.  The CAT default templates were obtained by standard Dartel/Shooting registration of 555 IXI subjects between 20 and 80 years. '
    'The registration time is typically about 3, 10, and 5 minutes for Dartel, Shooting, and optimized Shooting for the default registration resolution. '
    ''
  }];
elseif expert==1
  regstr        = cfg_menu;
  regstr.labels = {
    'Dartel (0)'
    'Default Shooting (4)'
    'Optimized Shooting - vox (5)'
    'Optimized Shooting - fast (eps)'
    'Optimized Shooting - standard (0.5)'
    'Optimized Shooting - fine (1.0)'
    'Optimized Shooting - strong (11)'
    'Optimized Shooting - medium (12)'
    'Optimized Shooting - soft (13)'
  };
  regstr.values = {0 4 5 eps 0.5 1.0 11 12 13};
  regstr.help = [regstr.help; { ...
    'The strength of the optimized Shooting registration depends on the stopping criteria (controlled by the "extopts.regstr" parameter) and by the final registration resolution that can be given by the template (fast,standard,fine), as fixed value (hard,medium,soft), or (iii) by the output resolution (vox).   In general the template resolution is the best choice to allow an adaptive normalization depending on the individual anatomy with some control of the calculation time. Fixed resolution allows to roughly define the degree of normalization for all images with 2.0 mm for smoother and 1.0 mm for stronger deformations.  For special cases the registration resolution can also be set by the output resolution controlled by the "extopts.vox" parameter. '
    ''
    '  0   .. "Dartel"'
    '  4   .. "Default Shooting"'
    '  5   .. "Optimized Shooting - vox"        .. vox/2:vox/4:vox'
    ''
    '  eps .. "Optimized Shooting - fast"       .. TR/2:TR/4:TR (avg. change rate)'
    '  0.5 .. "Optimized Shooting - standard"   .. TR/2:TR/4:TR (avg. change rate)'
    '  1.0 .. "Optimized Shooting - fine"       .. TR/2:TR/4:TR (small change rate)'
    ''
    '  11  .. "Optimized Shooting - stong"      .. max( 1.0 , [3.0:0.5:1.0] )'
    '  22  .. "Optimized Shooting - medium"     .. max( 1.5 , [3.0:0.5:1.0] )'
    '  23  .. "Optimized Shooting - soft"       .. max( 2.0 , [3.0:0.5:1.0] )'
   }];
else
  % allow different registrations settings by using a matrix
  regstr         = cfg_entry;
  regstr.strtype = 'r';
  regstr.num     = [1 inf];
  regstr.help = [regstr.help; { ...
    '"Default Shooting" runs the original Shooting approach for existing templates and takes about 10 minutes per subject for 1.5 mm templates and about 1 hour for 1.0 mm. '
    'The "Optimized Shooting" approach uses lower spatial resolutions in the first iterations and an adaptive stopping criteria that allows faster processing of about 6 minutes for 1.5 mm and 15 minutes for 1.0 mm. '
    ''
    'In the development modus the deformation levels are set by the following values (TR=template resolution) ...'
    '  0         .. "Use Dartel" '                                     
    '  eps - 1   .. "Optimized Shooting" with lower (eps; fast) to higher quality (1; slow; default 0.5)'
    '  2         .. "Optimized Shooting"      .. 3:(3-TR)/4:TR'
    '  3         .. "Optimized Shooting"      .. TR/2:TR/4:TR'
    '  4         .. "Default   Shooting"      .. only TR'
    '  5         .. "Optimized vox Shooting " .. vox/2:vox/4:vox'
    ''
    '  10        .. "Stronger Shooting"       .. max( 0.5 , [2.5:0.5:0.5] )'
    '  11        .. "Strong Shooting"         .. max( 1.0 , [3.0:0.5:1.0] )'
    '  12        .. "Medium Shooting"         .. max( 1.5 , [3.0:0.5:1.0] )'
    '  13        .. "Soft   Shooting"         .. max( 2.0 , [3.0:0.5:1.0] )'
    '  14        .. "Softer Shooting"         .. max( 2.5 , [3.0:0.5:1.0] )'
    '  15        .. "Supersoft Shooting"      .. max( 3.0 , [3.0:0.5:1.0] )'
    ''
    '  10        .. "Stronger Shooting TR"    .. max( max( 0.5 , TR ) , [2.5:0.5:0.5] )'
    '  21        .. "Strong Shooting TR"      .. max( max( 1.0 , TR ) , [3.0:0.5:1.0] )'
    '  22        .. "Medium Shooting TR"      .. max( max( 1.5 , TR ) , [3.0:0.5:1.0] )'
    '  23        .. "Soft   Shooting TR"      .. max( max( 2.0 , TR ) , [3.0:0.5:1.0] )'
    '  24        .. "Softer Shooting TR"      .. max( max( 2.5 , TR ) , [3.0:0.5:1.0] )'
    '  25        .. "Softer Shooting TR"      .. max( max( 3.0 , TR ) , [3.0:0.5:1.0] )'
    ''
    'Double digit variants runs only for a limited resolutions and produce softer maps.  The cases with TR are further limited by the template resolution and to avoid additional interpolation. '
    ''
    'For each given value a separate deformation process is started in inverse order and saved in subdirectories.  The first given value that runs last will be used in the following CAT processing. ' 
    ''
    }]; 
end
regstr.tag    = 'regstr';
regstr.name   = 'Spatial registration';
regstr.def    = @(val)cat_get_defaults('extopts.regstr', val{:});

%---------------------------------------------------------------------

registration        = cfg_branch;
registration.tag    = 'registration';
registration.name   = 'Spatial Registration';
if expert<2
  registration.val  = {darteltpm shootingtpm regstr};
else
  registration.val  = {T1 brainmask cat12atlas darteltpm shootingtpm regstr}; 
end
registration.help   = {
  'For spatial registration CAT offers to use the classical Dartel (Ashburner, 2008) and Shooting (Ashburner, 2011) registrations to a existing template. Furthermore, an optimized shooting approach is available that use adaptive threshold and lower initial resolution to improve accuracy and calculation time at once.  The CAT default templates were obtained by standard Dartel/Shooting registration of 555 IXI subjects between 20 and 80 years. '
  'The registration time is typically about 3, 10, and 5 minutes for Dartel, Shooting, and optimized Shooting for the default registration resolution. '
  ''
}; 

%---------------------------------------------------------------------

pbtres         = cfg_entry;
pbtres.tag     = 'pbtres';
pbtres.name    = 'Voxel size for thickness estimation';
pbtres.strtype = 'r';
pbtres.num     = [1 1];
pbtres.def     = @(val)cat_get_defaults('extopts.pbtres', val{:});
pbtres.help    = {
  'Internal isotropic resolution for thickness estimation in mm.'
  ''
};


%------------------------------------------------------------------------
% special expert and developer options 
%------------------------------------------------------------------------

lazy         = cfg_menu;
lazy.tag     = 'lazy';
lazy.name    = 'Lazy processing';
lazy.labels  = {'yes','No'};
lazy.values  = {1,0};
lazy.val     = {0};
lazy.help    = {
  'Do not process data if result already exist. '
};

experimental        = cfg_menu;
experimental.tag    = 'experimental';
experimental.name   = 'Use experimental code';
experimental.labels = {'No','Yes'};
experimental.values = {0 1};
experimental.def    = @(val)cat_get_defaults('extopts.experimental', val{:});
experimental.help   = {
  'Use experimental code and functions.'
  ''
  'WARNING: This parameter is only for developer and will call functions that are not safe and may change in future versions!'
  ''
};

ignoreErrors        = cfg_menu;
ignoreErrors.tag    = 'ignoreErrors';
ignoreErrors.name   = 'Ignore errors';
ignoreErrors.labels = {'No','Yes'};
ignoreErrors.values = {0 1};
ignoreErrors.def    = @(val)cat_get_defaults('extopts.ignoreErrors', val{:});
ignoreErrors.help   = {
  'Catch preprocessing errors and move on with the next subject'
};

verb         = cfg_menu;
verb.tag     = 'verb';
verb.name    = 'Verbose processing level';
verb.labels  = {'none','default','details'};
verb.values  = {0 1 2};
verb.def     = @(val)cat_get_defaults('extopts.verb', val{:});
verb.help    = {
  'Verbose processing.'
};


print         = cfg_menu;
print.tag     = 'print';
print.name    = 'Create CAT report';
print.labels  = {'No','Yes (volume only)','Yes (volume and surfaces)'};
print.values  = {0 1 2};
print.def     = @(val)cat_get_defaults('extopts.print', val{:});
print.help    = {
  'Create final CAT report that requires Java.'
};


%---------------------------------------------------------------------
% Resolution
%---------------------------------------------------------------------

resnative        = cfg_branch;
resnative.tag    = 'native';
resnative.name   = 'Native resolution ';
resnative.help   = {
    'Preprocessing with native resolution.'
    ''
    'Examples:'
    '  native resolution       internal resolution '
    '   0.95 0.95 1.05     >     0.95 0.95 1.05'
    '   0.45 0.45 1.70     >     0.45 0.45 1.70'
    '   2.00 2.00 2.00     >     2.00 2.00 2.00'
    '' 
  }; 

resbest        = cfg_entry;
resbest.tag    = 'best';
resbest.name   = 'Best native resolution';
resbest.def    = @(val)cat_get_defaults('extopts.resval', val{:});
resbest.num    = [1 2];
resbest.help   = {
    'Preprocessing with the best (minimal) voxel dimension of the native image. The first parameters defines the lowest spatial resolution for every dimension, while the second defines a tolerance range to avoid tiny interpolations for almost correct resolutions. '
    ''
    'Examples:'
    '  Parameters    native resolution       internal resolution'
    '  [1.00 0.10]    0.95 1.05 1.25     >     0.95 1.05 1.00'
    '  [1.00 0.10]    0.95 1.05 1.05     >     0.95 1.05 1.05'
    '  [1.00 0.20]    0.45 0.45 1.50     >     0.45 0.45 1.00'
    '  [0.75 0.20]    0.45 0.45 1.50     >     0.45 0.45 0.75'  
    '  [0.75 0.00]    0.45 0.45 0.80     >     0.45 0.45 0.80'  
    ''
  }; 

resfixed        = cfg_entry;
resfixed.tag    = 'fixed';
resfixed.name   = 'Fixed resolution';
resfixed.val    = {[1.0 0.1]};
resfixed.num    = [1 2];
resfixed.help   = {
    'This option sets an isotropic voxel size that is controlled by the first parameter, whereas the second parameter defines a tolerance range to avoid tiny interpolations for almost correct resolutions. The fixed resolution option can also be used to improve preprocessing stability and speed of high resolution data, for instance protocols with high in-plane resolution and large slice thickness (e.g. 0.5x0.5x1.5 mm) and atypical spatial noise pattern. ' 
    ''
    'Examples: '
    '  Parameters     native resolution       internal resolution'
    '  [1.00 0.10]     0.45 0.45 1.70     >     1.00 1.00 1.00'
    '  [1.00 0.10]     0.95 1.05 1.25     >     0.95 1.05 1.00'
    '  [1.00 0.02]     0.95 1.05 1.25     >     1.00 1.00 1.00'
    '  [0.75 0.10]     0.75 0.95 1.25     >     0.75 0.75 0.75'
    ''
  }; 


restype        = cfg_choice;
restype.tag    = 'restypes';
restype.name   = 'Internal resampling for preprocessing';
switch cat_get_defaults('extopts.restype')
  case 'native', restype.val = {resnative};
  case 'best',   restype.val = {resbest};
  case 'fixed',  restype.val = {resfixed};
end
if ~expert
  restype        = cfg_menu;
  restype.tag    = 'restypes';
  restype.name   = 'Internal resampling for preprocessing';
  restype.labels = {
    'Fixed 1.0 mm'
    'Fixed 0.8 mm'
    'Best native'
  };
  restype.values = {struct('fixed',[1.0 0.1]) ...
                    struct('fixed',[0.8 0.1]) ...
                    struct('best', [0.5 0.1])};
  restype.val    = {struct('fixed',[1.0 0.1])};
  restype.help   = [regstr.help; { ...
    'The default fixed image resolution offers a good trade-off between optimal quality and preprocessing time and memory demands. Standard structural data with a voxel resolution around 1 mm or even data with high in-plane resolution and large slice thickness (e.g. 0.5x0.5x1.5 mm) will benefit from this setting. If you have higher native resolutions the highres option "Fixed 0.8 mm" will sometimes offer slightly better preprocessing quality with an increase of preprocessing time and memory demands. In case of even higher resolutions and high signal-to-noise ratio (e.g. for 7 T data) the "Best native" option will process the data on the highest native resolution. I.e. a resolution of 0.4x0.7x1.0 mm will be interpolated to 0.4x0.4x0.4 mm. A tolerance range of 0.1 mm is used to avoid interpolation artifacts, i.e. a resolution of 0.95x1.01x1.08 mm will not be interpolated in case of the "Fixed 1.0 mm"!  '
    ''
  }];
else
  restype.values = {resnative resbest resfixed};
  restype.help   = {
    'The default fixed image resolution offers a good trade-off between optimal quality and preprocessing time and memory demands. Standard structural data with a voxel resolution around 1mm or even data with high in-plane resolution and large slice thickness (e.g. 0.5x0.5x1.5 mm) will benefit from this setting. If you have higher native resolutions a change of the fixed resolution to smaller values will sometimes offer slightly better preprocessing quality with a increase of preprocessing time and memory demands. In case of even higher resolutions and high signal-to-noise ratio (e.g. for 7T data) the "Best native" option will process the data on the highest native resolution. I.e. a resolution of 0.4x0.7x1.0 mm will be interpolated to 0.4x0.4x0.4 mm. A tolerance range of 0.1 mm is used to avoid interpolation artifacts, i.e. a resolution of 0.95x1.01x1.08 mm will not be interpolated in case of the "Fixed 1.0 mm"!  '
    ''
  }; 
end



%------------------------------------------------------------------------
% AMAP MRF Filter (expert)
%------------------------------------------------------------------------
mrf         = cfg_menu; %
mrf.tag     = 'mrf';
mrf.name    = 'Strength of MRF noise correction';
mrf.labels  = {'none','light','medium','strong','auto'};
mrf.values  = {0 0.1 0.2 0.3 1};
mrf.def     = @(val)cat_get_defaults('extopts.mrf', val{:});
mrf.help    = {
  'Strength of the MRF noise correction of the AMAP segmentation. '
  ''
};


%------------------------------------------------------------------------
% Cleanup
%------------------------------------------------------------------------
cleanupstr         = cfg_menu;
cleanupstr.tag     = 'cleanupstr';
cleanupstr.name    = 'Strength of Final Clean Up';
cleanupstr.def     = @(val)cat_get_defaults('extopts.cleanupstr', val{:});
if ~expert
  cleanupstr.labels  = {'none','light','medium','strong'};
  cleanupstr.values  = {0 0.25 0.50 0.75};
  cleanupstr.help    = {
    'Strength of tissue cleanup after AMAP segmentation. The cleanup removes remaining meninges and corrects for partial volume effects in some regions. If parts of brain tissue were missing then decrease the strength.  If too many meninges are visible then increase the strength. '
    ''
  };
else
  cleanupstr.labels  = {'none (0)','light (0.25)','medium (0.50)','strong (0.75)','heavy (1.00)'};
  cleanupstr.values  = {0 0.25 0.50 0.75 1.00};
  cleanupstr.help    = {
    'Strength of tissue cleanup after AMAP segmentation. The cleanup removes remaining meninges and corrects for partial volume effects in some regions. If parts of brain tissue were missing then decrease the strength.  If too many meninges are visible then increase the strength. '
    ''
    'The strength changes multiple internal parameters: '
    ' 1) Size of the correction area'
    ' 2) Smoothing parameters to control the opening processes to remove thin structures '
    ''
  };
end
if expert==2
  cleanupstr.labels = [cleanupstr.labels 'SPM (2.00)'];
  cleanupstr.values = [cleanupstr.values 2.00]; 
end


%------------------------------------------------------------------------
% Skull-stripping
%------------------------------------------------------------------------
gcutstr           = cfg_menu;
gcutstr.tag       = 'gcutstr';
gcutstr.name      = 'Skull-Stripping';
gcutstr.def       = @(val)cat_get_defaults('extopts.gcutstr', val{:});
gcutstr.help      = {
  'Method of skull-stripping before AMAP segmentation. The SPM approach works quite stable for the majority of data. However, in some rare cases parts of GM (i.e. in frontal lobe) might be cut. If this happens the GCUT approach is a good alternative. GCUT is a graph-cut/region-growing approach starting from the WM area.'
  ''
};
if ~expert
  gcutstr.labels  = {'SPM approach' 'GCUT approach'};
  gcutstr.values  = {0 0.50};
else
  gcutstr.labels  = {'SPM approach (0)','GCUT medium (0.50)','SPM+ approach (2)','SPM+ and GCUT medium (3)'};
  gcutstr.values  = {0 0.50 2 3};
end


%------------------------------------------------------------------------
% Noise correction (expert)
%------------------------------------------------------------------------

% expert only
NCstr        = cfg_menu;
NCstr.tag    = 'NCstr';
NCstr.name   = 'Strength of Noise Corrections';
if expert
  NCstr.help    = {
    'Strength of the spatial adaptive (sub-resolution) non local means (SANLM) noise correction. Please note that the filter strength is automatically estimated. Change this parameter only for specific conditions. Typical values are: none (0), classic (1), light (2), medium (3|-inf), and strong (4). The "classic" option use the ordinal SANLM filter without further adaptions. The "light" option applies half of the filter strength of the adaptive "medium" cases, whereas the "strong" option uses the full filter strength, force sub-resolution filtering and applies an additional iteration. Sub-resolution filtering is only used in case of high image resolution below 0.8 mm or in case of the "strong" option. '
    ''
  };
  NCstr.labels = {'none (0)','classic (1)','light (2)','medium (3|-inf)','strong (4)'};
  NCstr.values = {0 1 2 -inf 4};
else
  NCstr.labels = {'none','light','medium','strong'};
  NCstr.values = {0 2 -inf 4};
  NCstr.help   = {
    'Strength of the (sub-resolution) spatial adaptive  non local means (SANLM) noise correction. Please note that the filter strength is automatically estimated. Change this parameter only for specific conditions. The "light" option applies only half of the filter strength of the adaptive "medium" cases and no sub-resolution filtering. The "medium" case use the full adaptive filter strength and sub-resolution filtering in case of high image resolution below 0.8 mm. The "strong" option uses the full filter strength without adaption, forces the sub-resolution filtering and applies an additional iteration. All cases used an anatomical depending filter strength adaption, i.e. full (adaptive) filter strength for 1 mm data and no filtering for 2.5 mm data. '
    ''
  };
end
NCstr.def    = @(val)cat_get_defaults('extopts.NCstr', val{:});


%------------------------------------------------------------------------
% Blood Vessel Correction (expert)
%------------------------------------------------------------------------

BVCstr         = cfg_menu;
BVCstr.tag     = 'BVCstr';
BVCstr.name    = 'Strength of Blood Vessel Corrections';
BVCstr.labels  = {'none (0)','light (eps)','medium (0.50)','strong (1.00)'};
BVCstr.values  = {0 eps 0.50 1.00};
BVCstr.def     = @(val)cat_get_defaults('extopts.BVCstr', val{:});
BVCstr.help    = {
  'Strength of the Blood Vessel Correction (BVC).'
  ''
};


%------------------------------------------------------------------------
% Local Adaptive Segmentation
%------------------------------------------------------------------------
LASstr         = cfg_menu;
LASstr.tag     = 'LASstr';
LASstr.name    = 'Strength of Local Adaptive Segmentation';
if ~expert 
  LASstr.labels  = {'none','light','medium','strong'};
  LASstr.values  = {0 0.25 0.50 0.75};
else
  LASstr.labels  = {'none (0)','ultralight (eps)','light (0.25)','medium (0.50)','strong (0.75)','heavy (1.00)'};
  LASstr.values  = {0 eps 0.25 0.50 0.75 1.00};
end
LASstr.def     = @(val)cat_get_defaults('extopts.LASstr', val{:});
LASstr.help    = {
  'Additionally to WM-inhomogeneities, GM intensity can vary across different regions such as the motor cortex, the basal ganglia, or the occipital lobe. These changes have an anatomical background (e.g. iron content, myelinization), but are dependent on the MR-protocol and often lead to underestimation of GM at higher intensities and overestimation of CSF at lower intensities. Therefore, a local intensity transformation of all tissue classes is used to reduce these effects in the image. This local adaptive segmentation (LAS) is applied before the final AMAP segmentation.'
  ''
};

%------------------------------------------------------------------------
% XASL quality
%------------------------------------------------------------------------

xasl_quality    = cfg_menu;
xasl_quality.tag = 'xasl_quality';
xasl_quality.name = 'Run xASL on high quality';
xasl_quality.labels = {'no','yes'};
xasl_quality.values = {0 1};
xasl_quality.def = @(val) 1; 
xasl_quality.help    = {'Runs xASL on high quality, if 0 then run at faster settings with lower resolution and fewer iterations'};

%------------------------------------------------------------------------
% XASL use lesions for masking
%------------------------------------------------------------------------

xasl_lesion         = cfg_files;
xasl_lesion.tag     = 'xasl_lesion';
xasl_lesion.name    = 'xASL Lesion';
xasl_lesion.def     = @(val)cat_get_defaults('extopts.xasl_lesion', val{:});
xasl_lesion.num     = [1 1];
xasl_lesion.filter  = 'image';
xasl_lesion.ufilter = 'LesionCAT'; 
xasl_lesion.help    = {
  'Select the lesion file that will be used to mask the registration. '
  ''
  'This should be a single file - a summary of all lesions.'
  ''
};


%------------------------------------------------------------------------
% XASL save steps
%------------------------------------------------------------------------

xasl_savesteps    = cfg_menu;
xasl_savesteps.tag = 'xasl_savesteps';
xasl_savesteps.name = 'Save intermediate registration steps for xASL';
xasl_savesteps.labels = {'no','yes'};
xasl_savesteps.values = {0 1};
xasl_savesteps.def = @(val) 0; 
xasl_savesteps.help    = {'Do not save only the final DARTEL registration but also all the intermediate steps'};

%------------------------------------------------------------------------
% XASL disable DARTEL
%------------------------------------------------------------------------

xasl_disabledartel    = cfg_menu;
xasl_disabledartel.tag = 'xasl_disabledartel';
xasl_disabledartel.name = 'Disables the final full DARTEL spatial normalization step as a part of the CAT12 segmentation';
xasl_disabledartel.labels = {'no','yes'};
xasl_disabledartel.values = {0 1};
xasl_disabledartel.def = @(val) 0; 
xasl_disabledartel.help    = {'When set to 1, the final DARTEL step is not run, and only the basic non-linear registration is utilized'};

%------------------------------------------------------------------------
% WM Hyperintensities (expert)
%------------------------------------------------------------------------
wmhc        = cfg_menu;
wmhc.tag    = 'WMHC';
wmhc.name   = 'WM Hyperintensity Correction (WMHCs) - in development';
wmhc.def    = @(val)cat_get_defaults('extopts.WMHC', val{:});
wmhc.help   = {
  'WARNING: Please note that the detection of WM hyperintensies is still under development and does not have the same accuracy as approaches that additionally consider FLAIR images (e.g. Lesion Segmentation Toolbox)! '
  'In aging or (neurodegenerative) diseases WM intensity can be reduced locally in T1 or increased in T2/PD images. These so-called WM hyperintensies (WMHs) can lead to preprocessing errors. Large GM areas next to the ventricle can cause normalization problems. Therefore, a temporary correction for normalization is useful if WMHs are expected. CAT allows different ways to handle WMHs: '
  ''
  ' 0) No Correction (handled as GM). '
  ' 1) Temporary (internal) correction as WM for spatial normalization. '
  ' 2) Permanent correction to WM. ' 
  ' 3) Handling as separate class. '
  ''
};
wmhc.values = {0 1 2 3};
if expert>1
  wmhc.labels = { ...
    'no correction (0)' ...
    'set WMH as WM only for normalization (1)' ... 
    'set WMH as WM (2)' ...
    'set WMH as own class (3)' ...
  };
else
  wmhc.labels = { ...
    'no WMH correction' ...
    'set WMH as WM only for normalization' ... 
    'set WMH as WM' ...
    'set WMH as own class' ...
  };
end

% deactivated 20180714 because the WMHC in cat_vol_partvol did not support 
% user modification yet
%{
WMHCstr         = cfg_menu;
WMHCstr.tag     = 'WMHCstr';
WMHCstr.name    = 'Strength of WMH Correction';
WMHCstr.labels  = {'none (0)','light (eps)','medium (0.50)','strong (1.00)'};
WMHCstr.values  = {0 eps 0.50 1.00};
WMHCstr.def     = @(val)cat_get_defaults('extopts.WMHCstr', val{:});
WMHCstr.help    = {
  'Strength of the modification of the WM Hyperintensity Correction (WMHC).'
  ''
};
%}

%------------------------------------------------------------------------
% stroke lesion handling (expert)
%------------------------------------------------------------------------
slc        = cfg_menu;
slc.tag    = 'SLC';
slc.name   = 'Stroke Lesion Correction (SLC) - in development';
slc.def    = @(val)cat_get_defaults('extopts.SLC', val{:});
slc.help   = {
  'WARNING: Please note that the handling of stroke lesion is still under development. '
  'Without further correction, stroke lesions will be handled by their most probable tissue class, i.e. typically as CSF or GM. Because the spatial registration tries to normalize these regions, the normalization of large regions lead to storng inproper deformations. '
  'To avoid poor deformations, we created a work-around by manually defined lesion maps. The ... tool can be used to set the tissue intensity to zeros to avoid normalization of stroke lesions. '
  ''
  ' 0) No Correction. '
  ' 1) Correction of manually defined regions that were set to zeros. '
};
if expert>1
  slc.values = {0 1 2};
  slc.labels = { ...
    'no SL handling (0)' ...
    'manual SL handling (1)' ... 
    'manual & automatic handling (2)' ...
  };
  slc.help   = [slc.help;{
    ' 2) Correction automatic detected regions. ' 
    ''}];
else
  slc.values = {0 1};
  slc.labels = { ...
    'no SL handling' ...
    'manual SL handling' ... 
  };
  slc.help   = [slc.help;{
    ''}];
end


%------------------------------------------------------------------------


app        = cfg_menu;
app.tag    = 'APP';
app.name   = 'Affine Preprocessing (APP)';
app.help   = { ...
    'Affine registration and SPM preprocessing can fail in some subjects with deviating anatomy (e.g. other species/neonates) or in images with strong signal inhomogeneities, or untypical intensities (e.g. synthetic images). An initial bias correction can help to reduce such problems (see details below). Recommended are the "rough" and "full" option.' 
    ''
    ' none   - no additional bias correction.' 
    ' rough  - rough APP bias correction (r1070)' 
    ' light  - iterative SPM bias correction on different resolutions' 
    ' full   - iterative SPM bias correction on different resolutions and final high resolution bias correction' 
    ''
  };
  
app.labels = {'none','rough','light','full'};
app.values = {0 1070 1 2};

if expert==2
  app.labels = {'none','light','full','rough','rough (new)','fine (new)'};
  app.values = {0 1 2 1070 3 4};
  app.help   = [app.help;{ 
    ' none       - no additional bias correction.' 
    ' light      - iterative SPM bias correction on different resolutions' 
    ' full       - iterative SPM bias correction on different resolutions and high resolution bias correction' 
    ' rough       - rough APP bias correction (R1070)' 
    ' rough (new) - rough APP bias correction - in development' 
    ' fine  (new) - rough (new) and fine APP bias correction - in development'    
    ''
  }];
end  
app.def    = @(val)cat_get_defaults('extopts.APP', val{:});
app.help   = [app.help; { ...
    'rough: Fast correction (~60s) that identifies large homogeneous areas to estimate the intensity inhomogeneity. A maximum-filter is used to reduce the partial volume effect in T1 data. Moreover, gradient and divergence maps were used to avoid side effects by high intensity tissues (e.g. blood vessels or head tissue). '
    'light: This approach focuses on an iterative application of the standard SPM preprocessing with different bias-correction options from low (samp=6 mm, biasfwhm=120 mm) to high frequency corrections  (samp=4.5 mm, biasfwhm=45 mm). However, the iterative calls require a lot of additional processing time (~500s) and is normally not required in data with low intensity inhomogeneity. '
    'full: In addition to the "light" approach a final maximum-based filter (similar to the ''rough'' method that needs about additional 60s) is used to remove remaining local inhomogeneities. '
    ''
}];
if expert==2
    app.help   = [app.help; { ...
    'rough (new): New version of the ''rough'' approach with improved handling of T2/PD data that is still in development. '
    'fine  (new): Additional fine processing after the ''rough'' processing that incorporates the different brain and head tissues but is also still in development.'
    ''
    }];
end
app.help   = [app.help; { ...
    'In conclusion, we recommend to use the "rough" option or alternatively the "full" option or in case of further problems the "light", "rough (new)" or "fine (new)" options.'
    ''
    }];


%------------------------------------------------------------------------

scale_cortex         = cfg_entry;
scale_cortex.tag     = 'scale_cortex';
scale_cortex.name    = 'Modify cortical surface creation';
scale_cortex.strtype = 'r';
scale_cortex.num     = [1 1];
scale_cortex.def     = @(val)cat_get_defaults('extopts.scale_cortex', val{:});
scale_cortex.help    = {
  'Scale intensity values for cortex to start with initial surface that is closer to GM/WM border to prevent that gyri/sulci are glued if you still have glued gyri/sulci (mainly in the occ. lobe).  You can try to decrease this value (start with 0.6).  Please note that decreasing this parameter also increases the risk of an interrupted parahippocampal gyrus.'
  ''
};

add_parahipp         = cfg_entry;
add_parahipp.tag     = 'add_parahipp';
add_parahipp.name    = 'Modify parahippocampal surface creation';
add_parahipp.strtype = 'r';
scale_cortex.num     = [1 1];
add_parahipp.def     = @(val)cat_get_defaults('extopts.add_parahipp', val{:});
add_parahipp.help    = {
  'Increase values in the parahippocampal area to prevent large cuts in the parahippocampal gyrus (initial surface in this area will be closer to GM/CSF border if the parahippocampal gyrus is still cut.  You can try to increase this value (start with 0.15).'
  ''
};

close_parahipp         = cfg_menu;
close_parahipp.tag     = 'close_parahipp';
close_parahipp.name    = 'Initial morphological closing of parahippocampus';
close_parahipp.labels  = {'No','Yes'};
close_parahipp.values  = {0 1};
close_parahipp.def     = @(val)cat_get_defaults('extopts.close_parahipp', val{:});
close_parahipp.help    = {
  'Apply initial morphological closing inside mask for parahippocampal gyrus to minimize the risk of large cuts of parahippocampal gyrus after topology correction. However, this may also lead to poorer quality of topology correction for other data and should be only used if large cuts in the parahippocampal areas occur.'
  ''
};

%------------------------------------------------------------------------
% special subbranches for experts and developer to cleanup the GUI 
%------------------------------------------------------------------------

segmentation      = cfg_branch;
segmentation.tag  = 'segmentation';
segmentation.name = 'Segmentation Options';
if expert==1
  segmentation.val  = {app,NCstr,xasl_quality,xasl_disabledartel,xasl_lesion,xasl_savesteps,LASstr,gcutstr,cleanupstr,wmhc,slc,restype};
elseif expert==2
  segmentation.val  = {app,NCstr,xasl_quality,xasl_disabledartel,xasl_lesion,xasl_savesteps,LASstr,gcutstr,cleanupstr,BVCstr,wmhc,slc,mrf,restype}; % WMHCstr,
end
segmentation.help = {'CAT12 parameter to control the tissue classification.';''};


admin      = cfg_branch;
admin.tag  = 'admin';
admin.name = 'Administration Options';
if expert==1
  admin.val  = {ignoreErrors verb print};
else
  admin.val  = {experimental lazy ignoreErrors verb print};
end
admin.help = {'CAT12 parameter to control the behaviour of the preprocessing pipeline.';''};

%------------------------------------------------------------------------

surface       = cfg_branch;
surface.tag   = 'surface';
surface.name  = 'Surface Options';
surface.val   = {pbtres scale_cortex add_parahipp close_parahipp};
surface.help  = {'CAT12 parameter to control the surface processing.';''};


%------------------------------------------------------------------------
% main extopts branch .. in order of their call in cat_main
%------------------------------------------------------------------------

extopts       = cfg_branch;
extopts.tag   = 'extopts';
extopts.name  = 'Extended options for CAT12 preprocessing';
if ~spm
  if expert>0 % experimental expert options
    extopts.val   = {segmentation,registration,vox,surface,admin}; 
  else
    extopts.val   = {app,xasl_quality,xasl_disabledartel,xasl_lesion,xasl_savesteps,LASstr,gcutstr,registration,vox,restype}; % NCstr?
  end
else
  % SPM based surface processing and thickness estimation
  if expert>0 % experimental expert options
    extopts.val   = {registration,vox,surface,admin}; 
  else
    extopts.val   = {registration,vox}; 
  end 
end
extopts.help  = {'Using the extended options you can adjust special parameters or the strength of different corrections ("0" means no correction and "0.5" is the default value that works best for a large variety of data).'};
