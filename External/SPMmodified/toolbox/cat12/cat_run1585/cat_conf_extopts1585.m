function extopts = cat_conf_extopts1585(expert,spm)
% Configuration file for extended CAT options
%
% Christian Gaser
% $Id: cat_conf_extopts.m 1552 2020-01-17 10:19:24Z dahnke $
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
vox.def     = @(val)cat_get_defaults1585('extopts.vox', val{:});
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
darteltpm.def     = @(val)cat_get_defaults1585('extopts.darteltpm', val{:});
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
shootingtpm.def     = @(val)cat_get_defaults1585('extopts.shootingtpm', val{:});
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
cat12atlas.def     = @(val)cat_get_defaults1585('extopts.cat12atlas', val{:});
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
brainmask.def     = @(val)cat_get_defaults1585('extopts.brainmask', val{:});
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
T1.def     = @(val)cat_get_defaults1585('extopts.T1', val{:});
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
bb.def     = @(val)cat_get_defaults1585('extopts.bb', val{:});
bb.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).'
''
};
%}

%---------------------------------------------------------------------

if expert==0
  regstr        = cfg_menu;
  regstr.labels = {
    'Optimized Shooting'
    'Default Shooting'
  };
  regstr.values = {0.5 4}; % special case 0 = 0.5 due to Dartel default seting
  regstr.help   = [regstr.help; { ...
    'For spatial registration CAT offers the use of the Dartel (Ashburner, 2008) and Shooting (Ashburner, 2011) registrations to an existing template. Furthermore, an optimized shooting approach is available that uses an adaptive threshold and lower initial resolutions to obtain a good tradeoff between accuracy and calculation time.  The CAT default templates were obtained by standard Dartel/Shooting registration of 555 IXI subjects between 20 and 80 years. '
    'The registration time is typically about 3, 10, and 5 minutes for Dartel, Shooting, and optimized Shooting for the default registration resolution. '
    ''
  }];
elseif expert==1
  regstr        = cfg_menu;
  regstr.labels = {
    'Default Shooting (4)'
    'Optimized Shooting - vox (5)'
    'Optimized Shooting - fast (eps)'
    'Optimized Shooting - standard (0.5)'
    'Optimized Shooting - fine (1.0)'
    'Optimized Shooting - strong (11)'
    'Optimized Shooting - medium (12)'
    'Optimized Shooting - soft (13)'
  };
  regstr.values = {4 5 eps 0.5 1.0 11 12 13}; % special case 0 = 0.5 due to Dartel default seting
  regstr.name   = 'Method';
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
  regstr.name   = 'Spatial registration';
  regstr.help    = [regstr.help; { ...
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
    '  6         .. "Optimized Shooting - hydrocephalus (6)"'
    '                Use many iterations! Very slow! Use k-means AMAP as initial Segmentation!'
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
if cat_get_defaults1585('extopts.regstr')>0
  regstr.def    = @(val)cat_get_defaults1585('extopts.regstr',val{:});
else
  regstr.val    = {0.5};
end

%---------------------------------------------------------------------

dartel        = cfg_branch;
dartel.tag    = 'dartel';
dartel.name   = 'Dartel Registration';
dartel.val    = {darteltpm};
dartel.help   = {
  'Classical Dartel (Ashburner, 2008) registrations to a existing template. The CAT default templates were obtained by standard Dartel registration of 555 IXI subjects between 20 and 80 years. '
  ''
};
 
shooting        = cfg_branch;
shooting.tag    = 'shooting';
shooting.name   = 'Shooting Registration';
shooting.val    = {shootingtpm regstr};
shooting.help   = {
  'Shooting (Ashburner, 2011) registrations to a existing template. Furthermore, an optimized shooting approach is available that use adaptive threshold and lower initial resolution to improve accuracy and calculation time at once. The CAT default templates were obtained by standard Shooting registration of 555 IXI subjects between 20 and 80 years. '
  ''
};

if expert<2
  registration        = cfg_choice;
  registration.tag    = 'registration';
  registration.name   = 'Spatial Registration';
  registration.values = {dartel shooting};
  if cat_get_defaults1585('extopts.regstr')==0
    registration.val  = {dartel};
  else
    registration.val  = {shooting};
  end
else
  registration      = cfg_branch;
  registration.tag  = 'registration';
  registration.name = 'Spatial Registration';
  registration.val  = {T1 brainmask cat12atlas darteltpm shootingtpm regstr}; 
end
registration.help   = {
  'For spatial registration CAT offers to use the classical Dartel (Ashburner, 2008) and Shooting (Ashburner, 2011) registrations to a existing template. Furthermore, an optimized shooting approach is available that use adaptive threshold and lower initial resolution to improve accuracy and calculation time at once.  The CAT default templates were obtained by standard Dartel/Shooting registration of 555 IXI subjects between 20 and 80 years. '
  'The registration time is typically about 3, 10, and 5 minutes for Dartel, Shooting, and optimized Shooting for the default registration resolution. '
  ''
}; 

%---------------------------------------------------------------------

% This version is not ready right now and I packed all working improvments
% (such as the Laplace-based blood-vessel-correction) into cat_surf_createCS2. 
pbtver         = cfg_menu;
pbtver.tag     = 'pbtmethod';
pbtver.name    = 'Projection-based thickness';
pbtver.labels  = {'PBT','PBTx','PBT2'};
pbtver.values  = {'pbt2','pbt2x','pbt3'};
pbtver.def     = @(val) 'pbt2x';
pbtver.help    = {
 ['Version of the projection-based thickness (PBT) thickness and surface reconstruction approach (Dahnke et al., 2013).  ' ...
  'Version 2 first estimates a temporary central surface by utilizing higher CSF and lower WM boundaries to utilize the partial volume effect.  ' ...
  'This surface divides the GM into a lower and upper GM area where PBT is used with sulcal reconstruction in the lower and gyrus reconstruction in the upper part. ' ...
  'The estimated thickness values of each part were projected over the whole GM and summed up to obtain the full thickness.  ' ...
  'Similar to PBT, the project-based thickness and the direct thickness (of the distance maps without sulcus/gyrus reconstruction) are combined by using the minimum. '] 
  ''
  'Experimental development parameter - do not change! '
  ''
};

pbtres         = cfg_entry;
pbtres.tag     = 'pbtres';
pbtres.name    = 'Voxel size for thickness estimation';
pbtres.strtype = 'r';
pbtres.num     = [1 1];
pbtres.def     = @(val)cat_get_defaults1585('extopts.pbtres', val{:});
pbtres.help    = {
  'Internal isotropic resolution for thickness estimation in mm.'
  ''
};

pbtlas         = cfg_menu;
pbtlas.tag     = 'pbtlas';
pbtlas.name    = 'Use correction for cortical myelination';
pbtlas.labels  = {'No','Yes'};
pbtlas.values  = {0 1};
pbtlas.def     = @(val)cat_get_defaults1585('extopts.pbtlas', val{:});
pbtlas.help    = {
  'Apply correction for cortical myelination by local intensity adaption to improve the description of the GM/WM boundary (added in CAT12.7).'
  'Experimental parameter, not yet working properly!'
  ''
};

% currently only for developer
collcorr         = cfg_menu;
collcorr.tag     = 'collcorr';
collcorr.name    = 'Correction for surface collisions';
if expert
collcorr.labels  = {...
  'No (createCS1; 0)',...
  'No (createCS2; 20)',...
  'PBT Self-Intersect (createCS2; 23)',... 
  'PBT + CAT Self-Intersect on Surface Normals (createCS2; 25)',... 
  };
  collcorr.values  = {0 20 23 25};
end
collcorr.help    = {
  ['The creation of the white and pial surface by adding/removing half thickness from the central surface requires further optimization to avoid self-intersections of the surfaces.  ' ...
   'These self-intersections occur for thin structures in strongly folded regions and were caused by different distance metrics.  ' ...
   'Most self-intersections can be fast corrected (~2 minutes) by a new correction added in CAT12.7 (201911) using the percentage position map of the PBT approach (default). ' ... 
   'Although this works quite good in most cases tiny self-intersections still remain and a more efficient but much slower correction would be necessary (plus ~30 minutes). ']
   ''
};
if expert>1
  collcorr.labels  = {...
    'No (createCS; 0)',...
    'CAT Selfintersect surface deformation approach (createCS; 1)',... not working
    'No (createCS2; 20)',...
    'CAT self-intersect surface deformation approach (createCS2; 21)',... 
    'CAT self-intersect surface normal approach (createCS2; 22)',... 
    'PBT self-intersect without optimization (createCS2; 23)',... 
    'PBT self-intersect with optimization (createCS2; 24)',... 
    'PBT + CAT Selfintersect (createCS2; 25)',... 
    ...'Delaunay Approach with Intensity Optimization (createCS2; 26)', ... not working at all
    };
  collcorr.values  = {0 1 20 21 22 23 24 25};
  collcorr.help    = [collcorr.help 
    { ...
   ['There is also an older surface deformation approach (collcorr = 1 | 21) that currently stops too early and results in underestimation of the cortical thickness. ' ...
    'Moreover, there is a faster version of the "CAT self-intersect surface normal approach" (collcorr = 22) and an optimized version of the PBT (collcorr = 24). '] 
    ''} ]; 
end
collcorr.def     = @(val)cat_get_defaults1585('extopts.collcorr', val{:});


reduce_mesh         = cfg_menu;
reduce_mesh.tag     = 'reduce_mesh';
reduce_mesh.name    = 'Reduce Mesh';
reduce_mesh.labels  = { ...
  'No reduction, PBT resolution (0)',...  
  'No reduction, optimal resolution (1)',...  
  'No reduction, interal resolution (2)',...
  'SPM approach init (3)',...
  'SPM approach full (5)',...
  'MATLAB approach init (4)',...
  'MATLAB approach full (6)'
};
reduce_mesh.values  = {0 1 2 3 5 4 6};
reduce_mesh.def     = @(val)cat_get_defaults1585('extopts.reduce_mesh', val{:});
reduce_mesh.help    = {
  ['Limitation of the surface resolution is essential for fast processing and acurate and equaly distributed meshes. ' ...
   'Mesh resolution depends in general on the voxel resolution used for surface creation and can be modified afterwards by refinment and reduction. ' ...
   'However, surface mesh reduction is not trivial and we observered fatal MATLAB errors (full uncatchable crash) and freezing of the following spherical registration on some computers. ' ...
   'This variable therefor controls multiple ways to handle mesh resolution in the surface creation process. '] 
   ''
  ['The first setting (0) uses no reduction at all, creating the intial surface at the PBT resolution and also use no mesh reduction and is very slow. ' ...
   'In general, PBT is processed at 0.5 mm and surface creation result in about 1.200k faces with a quadratic increase of processing time. ' ...
   'However, this resolution is not necessary for nearly all anylsis that often takes place at meshes with 160k (FreeSurfer) or 32k (CIVIT). ']
   ...
  ['Option (1) and (2) use volume reduction to created intial meshes on an optimal (1, depending on the final mesh resolution) or ' ...
   'the internal voxel-resolution (2, depending on your image resolution). ' ...
   'In both cases the maps are refined and further adapted to the PBT position map with a final mesh resolution of about 300k. '];    
   ''
  ['Surface-based reduction by SPM (3,5) or MATLAB (4,6) are used to optimize the initial surface, supporing a fast but still accurate topology correction. ' ...
   'After topology correction the resolution of the mesh is increased again and adapted for PBT position map.  ' ...
   'In option 3 and 4, a supersampling with following reduction is used to obtain an optimal equally distributed sampling. ' ...
   'However, some systems showed problems in the spherical registration (freezing) that seamed to depend on these severe modifications of the mesh. ' ...
   'Hence, option (1) and (2) only use a normal refinement without supersampling.']  
   ''
   'These settings are still in developent!'
   ''
};
if expert
  reduce_mesh.labels  = [reduce_mesh.labels(1:2) {'MATLAB approach external (7)'} reduce_mesh.labels(3:end)];
  reduce_mesh.values  = [reduce_mesh.values(1:2) 7                                reduce_mesh.values(3:end)];
  reduce_mesh.help    = [reduce_mesh.help;{'Option (7) opens an external MATLAB to apply the mesh reduction and only worked on ower systems.';''}]; 
end



% This is just an developer parameter (maybe for experts later) 
% I expect that 300k surfaces support best quality in relation to processing
% time, because surface reconstruction times are dominated by the registration
% that depends on the mesh resolution of the individual and template brain. 
% However a fast and a refined accurate version could be interesting.
% The fast version could be interesting for fast clinical processing but needs also low mesh average surfaces.  
% The more accurate version could be useful for high quality surface rendering.

vdist         = cfg_menu;
vdist.tag     = 'vdist';
vdist.name    = 'Mesh resolution';
vdist.labels  = {'low','optimal','fine'};
vdist.values  = {27/3 4/3 1/6};   % this is the square of the refinement distance 
vdist.def     = @(val)cat_get_defaults1585('extopts.vdist', val{:});
vdist.help    = {
 ['Higher mesh resolution may support more accurate surface reconstruction.  However, this is only useful for high resolution data (<0.8 mm). ' ...
  'For each level, the resolution is doubled and accuracy is increased slightly (square root), resulting in a linear increase of processing time for surface creation. ' ...
  'Because the processing time of the surface registration depends on the resolution of the individual and the template surface (about 300k faces), ' ...
  'processing speed does not benefit from lower individual resolution. ']
  ''
  '  minimal:    maximal vertex distance 4.24 mm > ~100k faces'
  '  optimal:    maximal vertex distance 1.61 mm > ~300k faces'
  '  super-fine: maximal vertex distance 0.57 mm > ~900k faces'
  ''
  'Experimental development parameter that only works for the "createCS2" options of "Correct for surface collisions" (added in CAT12.7, 201909)!'
  ''
};

%------------------------------------------------------------------------
% special expert and developer options 
%------------------------------------------------------------------------

lazy         = cfg_menu;
lazy.tag     = 'lazy';
lazy.name    = 'Lazy processing';
lazy.labels  = {'Yes','No'};
lazy.values  = {1,0};
lazy.val     = {0};
lazy.help    = {
    'Do not process data if the result already exists. '
};

experimental        = cfg_menu;
experimental.tag    = 'experimental';
experimental.name   = 'Use experimental code';
experimental.labels = {'No','Yes'};
experimental.values = {0 1};
experimental.def    = @(val)cat_get_defaults1585('extopts.experimental', val{:});
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
ignoreErrors.def    = @(val)cat_get_defaults1585('extopts.ignoreErrors', val{:});
ignoreErrors.help   = {
  'Catch preprocessing errors and move on with the next subject'
};

verb         = cfg_menu;
verb.tag     = 'verb';
verb.name    = 'Verbose processing level';
verb.labels  = {'none','default','details'};
verb.values  = {0 1 2};
verb.def     = @(val)cat_get_defaults1585('extopts.verb', val{:});
verb.help    = {
  'Verbose processing.'
};


print         = cfg_menu;
print.tag     = 'print';
print.name    = 'Create CAT report';
print.labels  = {'No','Yes (volume only)','Yes (volume and surfaces)'};
print.values  = {0 1 2};
print.def     = @(val)cat_get_defaults1585('extopts.print', val{:});
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
resbest.def    = @(val)cat_get_defaults1585('extopts.resval', val{:});
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

resopt        = cfg_entry;
resopt.tag    = 'optimal';
resopt.name   = 'Optimal resolution';
resopt.def    = @(val)cat_get_defaults1585('extopts.resval', val{:});
resopt.num    = [1 2];
resopt.help   = {
    'Preprocessing with an "optimal" voxel dimension that utilize the median and the volume of the voxel size for special handling of anisotropic images.  In many cases, untypically high slice-resolution (e.g. 0.5 mm for 1.5 Tesla) comes along with higher slice-thickness and increased image interferences.  Our tests showed that a simple interpolation to the best voxel resolution not only resulted in much longer calculation times but also in a worste segmenation (and surface reconstruction) compared to the fixed option with e.g. 1 mm.  Hence, this option tries to incooperate the voxel volume and its isotropy to balance the internal resolution.  E.g., an image with 0.5x0.5x1.5 mm will resampled at a resolution of 0.7x0.7x0.7 mm. ' 
    'The first parameters defines the lowest spatial resolution, while the second defines a tolerance range to avoid tiny interpolations for almost correct resolutions. '
    ''
    'Examples:'
    '  Parameters    native resolution       internal resolution'
    '  [1.00 0.10]    0.95 1.05 1.25     >     0.95 1.05 1.00'
    '  [1.00 0.10]    0.80 0.80 1.00     >     0.80 0.80 0.80'
    '  [1.00 0.10]    0.50 0.50 2.00     >     1.00 1.00 1.00'
    '  [1.00 0.10]    0.50 0.50 1.50     >     0.70 0.70 0.70'
    '  [1.00 0.10]    0.80 1.00 1.00     >     1.00 1.00 1.00'
    ''
  };

restype        = cfg_choice;
restype.tag    = 'restypes';
restype.name   = 'Internal resampling for preprocessing';
switch cat_get_defaults1585('extopts.restype')
  case 'native',  restype.val = {resnative};
  case 'best',    restype.val = {resbest};
  case 'fixed',   restype.val = {resfixed};
  case 'optimal', restype.val = {resopt};
end

if ~expert
  restype        = cfg_menu;
  restype.tag    = 'restypes';
  restype.name   = 'Internal resampling for preprocessing';
  restype.labels = {
    'Optimal'
    'Fixed 1.0 mm'
    'Fixed 0.8 mm'
    'Best native'
  };
  restype.values = {struct('optimal', [1.0 0.1]) ...
                    struct('fixed',   [1.0 0.1]) ...
                    struct('fixed',   [0.8 0.1]) ...
                    struct('best',    [0.5 0.1])};
  restype.val    = {struct('optimal', [1.0 0.1])};
  restype.help   = {
    'The default fixed image resolution offers a good trade-off between optimal quality and preprocessing time and memory demands. Standard structural data with a voxel resolution around 1 mm or even data with high in-plane resolution and large slice thickness (e.g. 0.5x0.5x1.5 mm) will benefit from this setting. If you have higher native resolutions the highres option "Fixed 0.8 mm" will sometimes offer slightly better preprocessing quality with an increase of preprocessing time and memory demands. In case of even higher resolutions and high signal-to-noise ratio (e.g. for 7 T data) the "Best native" option will process the data on the highest native resolution. I.e. a resolution of 0.4x0.7x1.0 mm will be interpolated to 0.4x0.4x0.4 mm. A tolerance range of 0.1 mm is used to avoid interpolation artifacts, i.e. a resolution of 0.95x1.01x1.08 mm will not be interpolated in case of the "Fixed 1.0 mm"!  '
    'This "optimal" option prefers an isotropic voxel size with at least 1.1 mm that is controlled by the median voxel size and a volume term that penalizes highly anisotropic voxels.'
    ''
  };
else
  restype.values = {resopt resnative resbest resfixed};
  restype.help   = {
    'The default fixed image resolution offers a good trade-off between optimal quality and preprocessing time and memory demands. Standard structural data with a voxel resolution around 1mm or even data with high in-plane resolution and large slice thickness (e.g. 0.5x0.5x1.5 mm) will benefit from this setting. If you have higher native resolutions a change of the fixed resolution to smaller values will sometimes offer slightly better preprocessing quality with a increase of preprocessing time and memory demands. In case of even higher resolutions and high signal-to-noise ratio (e.g. for 7T data) the "Best native" option will process the data on the highest native resolution. I.e. a resolution of 0.4x0.7x1.0 mm will be interpolated to 0.4x0.4x0.4 mm. A tolerance range of 0.1 mm is used to avoid interpolation artifacts, i.e. a resolution of 0.95x1.01x1.08 mm will not be interpolated in case of the "Fixed 1.0 mm"!  '
    'This "optimal" option prefers an isotropic voxel size with at least 1.1 mm that is controlled by the median voxel size and a volume term that penalizes highly anisotropic voxels.'
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
mrf.def     = @(val)cat_get_defaults1585('extopts.mrf', val{:});
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
cleanupstr.def     = @(val)cat_get_defaults1585('extopts.cleanupstr', val{:});
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
gcutstr.def       = @(val)cat_get_defaults1585('extopts.gcutstr', val{:});
gcutstr.help      = {
  'Method of initial skull-stripping before AMAP segmentation. The SPM approach works quite stable for the majority of data. However, in some rare cases parts of GM (i.e. in frontal lobe) might be cut. If this happens the GCUT approach is a good alternative. GCUT is a graph-cut/region-growing approach starting from the WM area. '
  'APRG (adaptive probability region-growing) is a new method that refines the probability maps of the SPM approach by region-growing techniques of the gcut approach with a final surface-based optimization strategy. This is currently the method with the most accurate and reliable results. '
  'If you use already skull-stripped data you can turn off skull-stripping although this is automaticaly detected in most cases. '
  'Please note that the choice of the skull-stripping method will also influence the estimation of TIV, because the methods mainly differ in the handling of the outer CSF around the cortical surface. '
  ''
};
if ~expert
  gcutstr.labels  = {'none (already skull-stripped)' 'SPM approach' 'GCUT approach' 'APRG approach'};
  gcutstr.values  = {-1 0 0.50 2};
else
  gcutstr.labels  = {'none (already skull-stripped) (-1)','SPM approach (0)','GCUT medium (0.50)','APRG approach (2)'};
  gcutstr.values  = {-1 0 0.50 2};
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
NCstr.def    = @(val)cat_get_defaults1585('extopts.NCstr', val{:});


%------------------------------------------------------------------------
% Blood Vessel Correction (expert)
%------------------------------------------------------------------------

BVCstr         = cfg_menu;
BVCstr.tag     = 'BVCstr';
BVCstr.name    = 'Strength of Blood Vessel Corrections';
BVCstr.labels  = {'none (0)','light (eps)','medium (0.50)','strong (1.00)'};
BVCstr.values  = {0 eps 0.50 1.00};
BVCstr.def     = @(val)cat_get_defaults1585('extopts.BVCstr', val{:});
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
LASstr.def     = @(val)cat_get_defaults1585('extopts.LASstr', val{:});
LASstr.help    = {
  'Additionally to WM-inhomogeneities, GM intensity can vary across different regions such as the motor cortex, the basal ganglia, or the occipital lobe. These changes have an anatomical background (e.g. iron content, myelinization), but are dependent on the MR-protocol and often lead to underestimation of GM at higher intensities and overestimation of CSF at lower intensities. Therefore, a local intensity transformation of all tissue classes is used to reduce these effects in the image. This local adaptive segmentation (LAS) is applied before the final AMAP segmentation.'
  ''
};


%------------------------------------------------------------------------
% WM Hyperintensities (expert)
%------------------------------------------------------------------------
wmhc        = cfg_menu;
wmhc.tag    = 'WMHC';
wmhc.name   = 'WM Hyperintensity Correction (WMHCs) - in development';
wmhc.def    = @(val)cat_get_defaults1585('extopts.WMHC', val{:});
wmhc.help   = {
  'WARNING: Please note that the detection of WM hyperintensies is still under development and does not have the same accuracy as approaches that additionally consider FLAIR images (e.g. Lesion Segmentation Toolbox)! '
  'In aging or (neurodegenerative) diseases WM intensity can be reduced locally in T1 or increased in T2/PD images. These so-called WM hyperintensies (WMHs) can lead to preprocessing errors. Large GM areas next to the ventricle can cause normalization problems. Therefore, a temporary correction for normalization is useful if WMHs are expected. CAT allows different ways to handle WMHs: '
  ''
  ' 0) No Correction (handled as GM). '
  ' 1) Temporary (internal) correction as WM for spatial normalization and estimation of cortical thickness. '
  ' 2) Permanent correction to WM. ' 
  ' 3) Handling as separate class. '
  ''
};
wmhc.values = {0 1 2 3};
if expert>1
  wmhc.labels = { ...
    'no correction (0)' ...
    'set WMH as WM only for normalization (1) and thickness estimation' ... 
    'set WMH as WM (2)' ...
    'set WMH as own class (3)' ...
  };
else
  wmhc.labels = { ...
    'no WMH correction' ...
    'set WMH as WM only for normalization and thickness estimation' ... 
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
WMHCstr.def     = @(val)cat_get_defaults1585('extopts.WMHCstr', val{:});
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
slc.def    = @(val)cat_get_defaults1585('extopts.SLC', val{:});
slc.help   = {
  'WARNING: Please note that the handling of stroke lesion is still under development. '
  'Without further correction, stroke lesions will be handled by their most probable tissue class, i.e. typically as CSF or GM. Because the spatial registration tries to normalize these regions, the normalization of large regions will lead to strong inproper deformations. '
  'To avoid poor deformations, we created a work-around by manually defined lesion maps. The "Manual image (lesion) masking" tool can be used to set the image intensity to zeros to avoid normalization of stroke lesions. '
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
% Currently there are to much different strategies and this parameter needs 
% revision. There a three basic APP functions that each include an initial 
% rough and a following fine method. The first is the SPM appraoch that 
% is a simple iterative call of the Unified segmentation with following 
% maximum-based bias correction. It is relatively stable but slow and can be 
% combined with the other APP routines. The second one is the classical 
% APP approach with default and fine processing (1070), followed by further 
% developed version that should be more correct with monitor variables and
% T2/PD compatibility but finally worse results. 
%
% So we need more test to find out which strategies will survive to support 
% an alternative if the standard failed with a fast standard and slow but 
% more powerfull other routines. Hence APP1070 (init) or it successor
% should be the standard. The SPM routines are a good alternative due to 
% their differnt concept. 
%------------------------------------------------------------------------


app        = cfg_menu;
app.tag    = 'APP';
app.name   = 'Affine Preprocessing (APP)';
% short help text
app.help   = { ...
    'Affine registration and SPM preprocessing can fail in some subjects with deviating anatomy (e.g. other species/neonates) or in images with strong signal inhomogeneities, or untypical intensities (e.g. synthetic images). An initial bias correction can help to reduce such problems (see details below). Recommended are the "default" and "full" option.' 
    ''
    ' none    - no additional bias correction' 
    ' light   - iterative SPM bias correction on different resolutions' 
    ' full    - iterative SPM bias correction on different resolutions and final high resolution bias correction' 
    ' default - default APP bias correction (r1070)' 
  };
app.def    = @(val)cat_get_defaults1585('extopts.APP', val{:});
app.labels = {'none','light','full','default'};
app.values = {0 1 2 1070};
if expert
  app.labels = [app.labels, {'rough (new)'}];
  app.values = [app.values {1144}];
  app.help   = [app.help;{ 
    ' rough (new) - rough APP bias correction (r1144) - in development' 
  }];
end  
% long help text
app.help   = [app.help; { ...
    ''
    'light: This approach focuses on an iterative application of the standard SPM preprocessing with different bias-correction options from low (samp=6 mm, biasfwhm=120 mm) to high frequency corrections  (samp=4.5 mm, biasfwhm=45 mm). However, the iterative calls require a lot of additional processing time (~500s) and is normally not required in data with low intensity inhomogeneity. '
    'full:  In addition to the "light" approach a final maximum-based filter (similar to the ''default'' method that needs about additional 60s) is used to remove remaining local inhomogeneities. '
    'default: Fast correction (~60s) that identifies large homogeneous areas to estimate the intensity inhomogeneity. A maximum-filter is used to reduce the partial volume effects in T1-weighted data. Moreover, gradient and divergence maps were used to avoid side effects by high intensity tissues (e.g. blood vessels or head tissue). '
}];
if expert
    app.help   = [app.help; { ...
    'rough (new): New version of the ''rough'' approach with improved handling of T2/PD data that is still in development. '
     }];
end
app.help   = [app.help;{''}];


%------------------------------------------------------------------------

if expert>1
  new_release        = cfg_menu;
  new_release.tag    = 'new_release';
  new_release.name   = 'New release functions';
  new_release.help   = { ...
      'Use new rather then standard functions. '
    };
  new_release.val    = {0};  
  new_release.labels = {'No','Yes'};
  new_release.values = {0 1};
end

%------------------------------------------------------------------------
  
if expert>1
  % different Affine registations ... not implemented yet
  %{
  spm_affreg        = cfg_menu;
  spm_affreg.tag    = 'spm_affreg'; 
  spm_affreg.name   = 'Affine registration approach';
  spm_affreg.help   = { ...
      'The affine registion is highly important for the whole pipeline. Failures result in low overlap to the TPM that troubles the Unified Segmenation and all following steps. Therefore, CAT uses different routines to obtain the best solution. However, this can fail especial in atypical subjects (very young/old) and we deside that it is maybe usefull to test the steps separately. Brain or head masks can imrove the results in some but also lead to problems in other cases.'
      ''
      '  The affreg routine process a affine registration based the orignal input (T1) image and a similar weighted scan . '
      '  The maffreg routine use'
    };
  %spm_affreg.def    = @(val)cat_get_defaults1585('extopts.spm_affreg', val{:}); 
  spm_affreg.labels = {
    'no affine registration' ... 
    'only affreg' ...
    'only maffreg' ...
    'affreg + maffreg' ...
    'affreg + maffreg supervised' ...
    };
  spm_affreg.values = {0 1 2 3 4};
  spm_affreg.val    = 3;
  %}
  
  % AMAP rather than SPM segmentation 
  spm_kamap        = cfg_menu;
  spm_kamap.tag    = 'spm_kamap';
  spm_kamap.name   = 'Initial segmentation';
  spm_kamap.help   = { ...
      'In rare cases the Unified Segmentation can fail in highly abnormal brains, where e.g. the cerebrospinal fluid of superlarge ventricles (hydrocephalus) were classified as white matter. However, if the affine registration is correct, the AMAP segmentation with an prior-independent k-means initialization can be used to replace the SPM brain tissue classification. ' 
      'Moreover, if the default Dartel and Shooting registrations will fail then the "Optimized Shooting - superlarge ventricles" option for "Spatial registration" is required! '
      ''
      ' SPM Unified Segmentation - use SPM Unified Segmentation segmentation (default) ' 
      ' k-means AMAP - k-means AMAP approach ' 
      ''
    };
  spm_kamap.def    = @(val)cat_get_defaults1585('extopts.spm_kamap', val{:});  
  spm_kamap.labels = {'SPM Unified Segmentation','k-means AMAP'};
  spm_kamap.values = {0 2};
end

%------------------------------------------------------------------------

scale_cortex         = cfg_entry;
scale_cortex.tag     = 'scale_cortex';
scale_cortex.name    = 'Modify cortical surface creation';
scale_cortex.strtype = 'r';
scale_cortex.num     = [1 1];
scale_cortex.def     = @(val)cat_get_defaults1585('extopts.scale_cortex', val{:});
scale_cortex.help    = {
  'Scale intensity values for cortex to start with initial surface that is closer to GM/WM border to prevent that gyri/sulci are glued if you still have glued gyri/sulci (mainly in the occ. lobe).  You can try to decrease this value (start with 0.6).  Please note that decreasing this parameter also increases the risk of an interrupted parahippocampal gyrus.'
  ''
};

add_parahipp         = cfg_entry;
add_parahipp.tag     = 'add_parahipp';
add_parahipp.name    = 'Modify parahippocampal surface creation';
add_parahipp.strtype = 'r';
scale_cortex.num     = [1 1];
add_parahipp.def     = @(val)cat_get_defaults1585('extopts.add_parahipp', val{:});
add_parahipp.help    = {
  'Increase values in the parahippocampal area to prevent large cuts in the parahippocampal gyrus (initial surface in this area will be closer to GM/CSF border if the parahippocampal gyrus is still cut.  You can try to increase this value (start with 0.15).'
  ''
};

close_parahipp         = cfg_menu;
close_parahipp.tag     = 'close_parahipp';
close_parahipp.name    = 'Initial morphological closing of parahippocampus';
close_parahipp.labels  = {'No','Yes'};
close_parahipp.values  = {0 1};
close_parahipp.def     = @(val)cat_get_defaults1585('extopts.close_parahipp', val{:});
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
  segmentation.val  = {app,NCstr,LASstr,gcutstr,cleanupstr,wmhc,slc,restype};
elseif expert==2
  segmentation.val  = {app,NCstr,spm_kamap,LASstr,gcutstr,cleanupstr,BVCstr,wmhc,slc,mrf,restype}; % WMHCstr,
end
segmentation.help = {'CAT12 parameter to control the tissue classification.';''};


admin      = cfg_branch;
admin.tag  = 'admin';
admin.name = 'Administration Options';
if expert==1
  admin.val  = {lazy ignoreErrors verb print};
elseif expert==2
  admin.val  = {experimental new_release lazy ignoreErrors verb print};
end
admin.help = {'CAT12 parameter to control the behaviour of the preprocessing pipeline.';''};

%------------------------------------------------------------------------

surface       = cfg_branch;
surface.tag   = 'surface';
surface.name  = 'Surface Options';
if expert>1
  surface.val   = {pbtres pbtver pbtlas collcorr reduce_mesh vdist scale_cortex add_parahipp close_parahipp}; % pbtver
else
  surface.val   = {pbtres pbtlas scale_cortex add_parahipp close_parahipp};
end
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
    extopts.val   = {app,LASstr,gcutstr,registration,vox,restype}; % NCstr?
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
