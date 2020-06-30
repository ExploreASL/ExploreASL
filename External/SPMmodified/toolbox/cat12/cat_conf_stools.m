function stools = cat_conf_stools(expert)
%_______________________________________________________________________
% wrapper for calling CAT surface utilities
%_______________________________________________________________________
% Robert Dahnke and Christian Gaser
% $Id: cat_conf_stools.m 1615 2020-05-09 10:19:40Z gaser $
%_______________________________________________________________________

  % try to estimate number of processor cores
  try
    numcores = max(cat_get_defaults('extopts.nproc'),1);
  catch
    numcores = 0;
  end
  
  if isdeployed, numcores = 0; end
  
  if ~exist('expert','var')
    expert = 0; % switch to de/activate further GUI options
  end

  
  % parallelize
  % ____________________________________________________________________
  nproc         = cfg_entry;
  nproc.tag     = 'nproc';
  nproc.name    = 'Split job into separate processes';
  nproc.strtype = 'w';
  nproc.val     = {numcores};
  nproc.num     = [1 1];
  nproc.help    = {
      'In order to use multi-threading the CAT12 segmentation job with multiple subjects can be split into separate processes that run in the background. You can even close Matlab, which will not affect the processes that will run in the background without GUI. If you do not want to run processes in the background then set this value to 0.'
      ''
      'Please note that no additional modules in the batch can be run except CAT12 segmentation. Any dependencies will be broken for subsequent modules.'
    };
 
  
  % do not process, if result already exists
  % ____________________________________________________________________
  lazy         = cfg_menu;
  lazy.tag     = 'lazy';
  lazy.name    = 'Lazy processing';
  lazy.labels  = {'Yes','No'};
  lazy.values  = {1,0};
  lazy.val     = {0};
  lazy.help    = {
    'Do not process data if the result already exists. '
  };


  % merge hemispheres
  % ____________________________________________________________________
  merge_hemi         = cfg_menu;
  merge_hemi.tag     = 'merge_hemi';
  merge_hemi.name    = 'Merge hemispheres';
  merge_hemi.labels  = {
    'No -   save resampled data for each hemisphere',...
    'Yes -  merge hemispheres'
  };
  merge_hemi.values  = {0,1};
  merge_hemi.val     = {1};
  merge_hemi.help    = {
    'Meshes for left and right hemisphere can be merged to one single mesh. This simplifies the analysis because only one analysis has to be made for both hemispheres and this is the recommended approach.'
    'However, this also means that data size is doubled for one single analysis which might be too memory demanding for studies with several hundreds or even more files. If your model cannot be estimated due to memory issues you should not merge the resampled data.'
  };


  % default mesh
  % ____________________________________________________________________
  mesh32k         = cfg_menu;
  mesh32k.tag     = 'mesh32k';
  mesh32k.name    = 'Resample Size';
  mesh32k.labels  = {
    '32k  mesh (HCP)',...
    '164k mesh (Freesurfer)'
  };
  mesh32k.values  = {1,0};
  mesh32k.val     = {1};
  mesh32k.help    = {
    'Resampling can be done either to a higher resolution 164k mesh that is compatible to Freesurfer data or to a lower resolution 32k mesh (average vertex spacing of ~2 mm) that is compatible to the Human Connectome Project (HCP).'
    'The HCP mesh has the advantage of being processed and handled much faster and with less memory demands. Another advantage is that left and right hemispheres are aligned to optionally allow a direct comparison between hemispheres.'
  };



%% import GUIs
%  ------------------------------------------------------------------------
  [check_mesh_cov,check_mesh_cov2] = check_mesh_cov_GUI;
  [surfresamp,surfresamp_fs]       = cat_surf_resamp_GUI(expert,nproc,merge_hemi,mesh32k,lazy);
  [vol2surf,vol2tempsurf]          = cat_surf_vol2surf_GUI(expert,merge_hemi,mesh32k);
  [surfcalc,surfcalcsub]           = cat_surf_calc_GUI(expert);
  renderresults                    = cat_surf_results_GUI(expert);
  surfextract                      = cat_surf_parameters_GUI(expert,nproc,lazy);
  flipsides                        = cat_surf_flipsides_GUI;
  surf2roi                         = cat_surf_surf2roi_GUI(expert,nproc);
  roi2surf                         = cat_roi_roi2surf_GUI; 
  surfstat                         = cat_surf_stat_GUI; 
  surfcon                          = cat_surf_spm_cfg_con;
  
  
%% Toolset
%  ---------------------------------------------------------------------
  stools = cfg_choice;
  stools.name   = 'Surface Tools';
  stools.tag    = 'stools';
  if expert==2
    stools.values = { ...
      check_mesh_cov, ...     .cat.stat
      check_mesh_cov2, ...    .cat.stat
      surfextract, ...        .cat.stools
      surfresamp, ...         .cat.stools
      surfresamp_fs,...       .cat.stools
      vol2surf, ...           .cat.stools
      vol2tempsurf, ...       .cat.stools
      surfcalc, ...           .cat.stools
      surfcalcsub, ...        .cat.stools
      surf2roi, ...           .cat.rtools?
      roi2surf, ...           .cat.rtools?
      flipsides, ...          .cat.stools
      surfstat ...
      surfcon ...
      renderresults, ...      .cat.disp/plot/write?
      };    
  elseif expert==1
    stools.values = { ...
      check_mesh_cov, ...
      check_mesh_cov2, ...
      surfextract, ...
      surfresamp, ...
      surfresamp_fs,...
      vol2surf, ...
      vol2tempsurf, ...
      surfcalc, ...
      surfcalcsub, ...
      surf2roi, ...
      surfstat ...
      surfcon ...
      renderresults, ...      .cat.disp/plot/write?
      };
  else
    stools.values = { ...
      check_mesh_cov, ...
      check_mesh_cov2, ...
      surfextract, ...
      surfresamp, ...
      surfresamp_fs,...
      surf2roi, ...
      vol2surf, ...
      vol2tempsurf, ...
      surfcalc, ...
      surfcalcsub, ...
      surfstat ...
      surfcon ...
      renderresults ...
      };
  end

%==========================================================================
% subfunctions of batch GUIs
%==========================================================================
function con = cat_surf_spm_cfg_con
  con      = spm_cfg_con;
  con.name = ['Surface ' con.name];
  con.vout = @vout_stats; % gifti rather than nifti output

function dep = vout_stats(varargin)
  dep(1)            = cfg_dep;
  dep(1).sname      = 'SPM.mat File';
  dep(1).src_output = substruct('.','spmmat');
  dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
  dep(2)            = cfg_dep;
  dep(2).sname      = 'All Con Surfaces';
  dep(2).src_output = substruct('.','con');
  dep(2).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
  dep(3)            = cfg_dep;
  dep(3).sname      = 'All Stats Surfaces';
  dep(3).src_output = substruct('.','spm');
  dep(3).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});

function surfstat = cat_surf_stat_GUI

  spmmat          = cfg_files;
  spmmat.tag      = 'spmmat';
  spmmat.name     = 'Select SPM.mat files';
  spmmat.filter   = {'mat'};
  spmmat.ufilter  = '^SPM\.mat$';
  spmmat.num      = [1 inf];
  spmmat.help     = {'Select the SPM.mat file that contains the design specification.'};

  surfstat        = cfg_exbranch;
  surfstat.tag    = 'SPM';
  surfstat.name   = 'Estimate Surface Model';
  surfstat.val    = {spmmat};
  surfstat.prog   = @cat_stat_spm;
  surfstat.vout   = @vout_cat_stat_spm;
  surfstat.help   = {
    ''};

function surf2roi = cat_surf_surf2roi_GUI(expert,nproc)
%% surface to ROI (in template space)
%  ---------------------------------------------------------------------
%  * CxN  cdata files [ thickness , curvature , ... any other file ] in template space [ resampled ]
%  * M    atlas files [ choose files ]
%  * average measures [ mean , std , median , max , min , mean95p ]
% 
%  - csv export (multiple files) > catROIs_ATLAS_SUBJECT
%  - xml export (one file)       > catROIs_ATLAS_SUBJECT
%  ---------------------------------------------------------------------
%  surface ROIs have sidewise index?!
%  ---------------------------------------------------------------------

% set of cdata files
  if expert && 0 % this is not ready now
    cdata         = cfg_files;
    cdata.tag     = 'cdata';
    cdata.name    = '(Left) Surface Data Files';
    cdata.filter  = 'any';
    cdata.ufilter = 'lh.(?!cent|pial|white|sphe|defe|hull|pbt).*';
    cdata.num     = [1 Inf];
    cdata.help    = {'Surface data sample. Both sides will be processed'};
  else % only smoothed/resampled
    cdata         = cfg_files;
    cdata.tag     = 'cdata';
    cdata.name    = '(Left) Surface Data Files';
    cdata.filter  = 'any';
    cdata.ufilter = '^lh.(?!cent|pial|white|sphe|defe|hull|pbt).*';
    cdata.num     = [1 Inf];
    cdata.help    = {'Surface data sample. Both sides will be processed'};
  end
  
  cdata_sample         = cfg_repeat;
  cdata_sample.tag     = 'cdata_sub.';
  cdata_sample.name    = 'Surface Data Sample';
  cdata_sample.values  = {cdata};
  cdata_sample.num     = [1 Inf];
  cdata_sample.help = {[...
    'Specify data for each sample (i.e. thickness, gyrification, ...). ' ...
    'All samples must have the same size and same order. ' ...
    ''
  ]};

% ROI files
  ROIs         = cfg_files;
  ROIs.tag     = 'rdata';
  ROIs.name    = '(Left) ROI atlas files';
  ROIs.filter  = 'any';
  if expert 
    ROIs.ufilter = 'lh.aparc.*';
  else
    ROIs.ufilter = 'lh.aparc_(a2009s|DK40|HCP_MMP1).*'; % not yet working for all atlases
  end
  ROIs.dir     = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces'); 
  ROIs.num     = [1 Inf];
  ROIs.help    = {'These are the ROI atlas files. Both sides will be processed.'};

  
% not used
%{ 
% ROI area
  area         = cfg_menu;
  area.tag     = 'area';
  area.name    = 'Estimate ROI Area';
  area.labels  = {'No','Yes'};
  area.values  = {0,1};
  area.val     = {1}; 
  area.help    = {'Estimate area of each ROI.'};
 
% ROI area
  vernum         = cfg_menu;
  vernum.tag     = 'vernum';
  vernum.name    = 'Count ROI Vertices';
  vernum.labels  = {'No','Yes'};
  vernum.values  = {0,1};
  vernum.val     = {1}; 
  vernum.help    = {'Count vertices of each ROI.'};
%}
  
% average mode within a ROI  
  % mean
  avg.mean         = cfg_menu;
  avg.mean.tag     = 'mean';
  avg.mean.name    = 'Mean Estimation';
  avg.mean.labels  = {'No','Yes'};
  avg.mean.values  = {0,1};
  avg.mean.val     = {1}; 
  avg.mean.help    = {'Set mean value estimation per ROI.'}; 
  % std
  avg.std         = cfg_menu;
  avg.std.tag     = 'std';
  avg.std.name    = 'STD Estimation';
  avg.std.labels  = {'No','Yes'};
  avg.std.values  = {0,1};
  avg.std.val     = {1}; 
  avg.std.help    = {'Set standard deviation estimation per ROI.'}; 
  % min
  avg.min         = cfg_menu;
  avg.min.tag     = 'min';
  avg.min.name    = 'Minimum Estimation';
  avg.min.labels  = {'No','Yes'};
  avg.min.values  = {0,1};
  avg.min.val     = {0}; 
  avg.min.help    = {'Set minimum estimation per ROI.'};   
  % max
  avg.max         = cfg_menu;
  avg.max.tag     = 'max';
  avg.max.name    = 'Maximum Estimation';
  avg.max.labels  = {'No','Yes'};
  avg.max.values  = {0,1};
  avg.max.val     = {0}; 
  avg.max.help    = {'Set maximum estimation per ROI.'};   
  % median
  avg.median         = cfg_menu;
  avg.median.tag     = 'median';
  avg.median.name    = 'Median Estimation';
  avg.median.labels  = {'No','Yes'};
  avg.median.values  = {0,1};
  avg.median.val     = {0}; 
  avg.median.help    = {'Set median estimation per ROI.'};   
  % all functions
  avg.main         = cfg_branch;
  avg.main.tag     = 'avg';
  avg.main.name    = 'ROI Average Functions';
  avg.main.val     = {
    avg.mean ...
    avg.std ...
    avg.min ...
    avg.max ...
    avg.median ...
  };

%% main function
  surf2roi      = cfg_exbranch;
  surf2roi.tag  = 'surf2roi';
  surf2roi.name = 'Extract ROI-based surface values';
  switch expert
  case 2
    surf2roi.val  = {
      cdata_sample ...
      ROIs ...
      nproc ... 
      avg.main};
  case {0, 1}
    surf2roi.val  = {cdata_sample};
  end
  surf2roi.prog = @cat_surf_surf2roi;
  surf2roi.vout = @vout_surf_surf2roi;
  surf2roi.help = {
    'While ROI-based values for VBM (volume) data are automatically saved in the label folder as XML file it is necessary to additionally extract these values for surface data. This has to be done after preprocessing the data and creating cortical surfaces. '
    ''
    'You can extract ROI-based values for cortical thickness but also for any other surface parameter that was extracted using the "Extract Additional Surface Parameters" function.'
    ''
    'Please note that these values are extracted from data in native space without any smoothing. As default the mean inside a ROI is calculated and saved as XML file in the label folder.'
  };

%==========================================================================
function [check_mesh_cov,check_mesh_cov2] = check_mesh_cov_GUI

%% Surface correlation and quality check
%-----------------------------------------------------------------------
  data_surf_cov         = cfg_files;
  data_surf_cov.tag     = 'data_surf';
  data_surf_cov.name    = 'Sample';
  data_surf_cov.filter  = 'gifti';
  data_surf_cov.ufilter = 'resampled';
  data_surf_cov.num     = [3 Inf];
  data_surf_cov.help    = {'Select resampled surfaces parameter files.'};
  
  data_xml              = cfg_files;
  data_xml.name         = 'Quality measures (optional)';
  data_xml.tag          = 'data_xml';
  data_xml.filter       = 'xml';
  data_xml.ufilter      = '^cat_.*';
  data_xml.val          = {{''}};
  data_xml.num          = [0 Inf];
  data_xml.help         = {...
'Select optional the quality measures that are saved during segmentation as xml-files in the report folder. This additionally allows to analyze image quality parameters such as noise, and bias. Please note, that the order of the xml-files should be the same as the other data files.'};

  sample_cov            = cfg_repeat;
  sample_cov.tag        = 'sample';
  sample_cov.name       = 'Data';
  sample_cov.values     = {data_surf_cov};
  sample_cov.num        = [1 Inf];
  sample_cov.help       = {...
  'Specify data for each sample. If you specify different samples the mean correlation is displayed in separate boxplots for each sample.'};

  c         = cfg_entry;
  c.tag     = 'c';
  c.name    = 'Vector';
  c.help    = {'Vector of nuisance values'};
  c.strtype = 'r';
  c.num     = [Inf Inf];

  nuisance         = cfg_repeat;
  nuisance.tag     = 'nuisance';
  nuisance.name    = 'Nuisance variable';
  nuisance.values  = {c};
  nuisance.num     = [0 Inf];
  nuisance.help    = {...
  'This option allows for the specification of nuisance effects to be removed from the data. A potential nuisance parameter can be age. In this case the variance explained by age will be removed prior to the calculation of the correlation.'};

  check_mesh_cov      = cfg_exbranch;
  check_mesh_cov.tag  = 'check_mesh_cov';
  check_mesh_cov.name = 'Check sample homogeneity of surfaces';
  check_mesh_cov.val  = {sample_cov,data_xml,nuisance};
  check_mesh_cov.prog = @cat_stat_check_cov;
  check_mesh_cov.help = {
  'In order to identify surfaces with poor image quality or even artefacts you can use this function. Surfaces measures have to be first resampled to the template space (e.g. resampled and smoothed data) and can be then check for each hemisphere separately. Please do not mix data from both hemisphere.'
  ''
  'The idea of this tool is to check the correlation of all files across the sample. The correlation is calculated between all surfaces measures and the mean for each surface measure is plotted using a boxplot and optionally the indicated filenames. The smaller the mean correlation the more deviant is this surface measures from the sample mean. In the plot outliers from the sample are usually isolated from the majority of surface measures which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the order of your measures.'};

%-----------------------------------------------------------------------  
 
  check_mesh_cov2      = check_mesh_cov; 
  check_mesh_cov2.tag  = 'check_mesh_cov2';
  check_mesh_cov2.val  = {sample_cov,nuisance};
  check_mesh_cov2.name = 'Check sample homogeneity of surfaces (new exp. version)';
  check_mesh_cov2.prog = @cat_stat_check_cov2;

%==========================================================================
function [vol2surf,vol2tempsurf] = cat_surf_vol2surf_GUI(expert,merge_hemi,mesh32k)
%% map volumetric data
%-----------------------------------------------------------------------  
  datafieldname         = cfg_entry;
  datafieldname.tag     = 'datafieldname';
  datafieldname.name    = 'Output Name';
  datafieldname.strtype = 's';
  datafieldname.num     = [1 Inf];
  datafieldname.val     = {'intensity'};
  datafieldname.help    = {
    'Name that is prepended to the filename of the mapped volume.'
    ''
    };
 
  interp         = cfg_menu;
  interp.tag     = 'interp';
  interp.name    = 'Interpolation Type';
  interp.labels  = {'Nearest neighbour','Linear','Cubic'};
  interp.values  = {{'nearest_neighbour'},{'linear'},{'cubic'}};
  interp.val     = {{'linear'}};
  interp.help    = {
    'Volume extraction interpolation type. '
    ' -linear:            Use linear interpolation (default).'
    ' -nearest_neighbour: Use nearest neighbour interpolation.'
    ' -cubic:             Use cubic interpolation.'
    ''
  };
  
  % sample function 
  sample         = cfg_menu;
  sample.tag     = 'sample';
  sample.name    = 'Sample Function';
  sample.labels  = {'Mean','Weighted mean','Maximum','Minimum','Absolute maximum','Multi-values'};
  sample.values  = {{'avg'},{'weighted_avg'},{'max'},{'min'},{'maxabs'},{'multi'}};
  sample.val     = {{'maxabs'}};
  sample.help    = {
    'Sample function to combine the values of the grid along the surface normals.'
    ' Mean:          Use average for mapping along normals.'
    ' Weighted mean: Use weighted average with gaussian kernel for mapping along normals. The kernel is so defined that values at the boundary are weighted with 50% while center is weighted with 100% (useful for (r)fMRI data.'
    ' Maximum:       Use maximum value for mapping along normals.'
    ' Minimum:       Use minimum value for mapping along normals.'
    ' Absolute maximum: Use absolute maximum value for mapping along normals (useful for mapping contrast images from 1st-level fMRI analysis).'
    ' Multi-values:  Map data for each grid step separately and save files with indicated grid value. Please note that this option is intended for high-resolution (f)MRI data only (e.g. 0.5mm voxel size).'
    ''
  };

  %% -- sampling points and average function
 
  % startpoint
  abs_startpoint         = cfg_entry;
  abs_startpoint.tag     = 'startpoint';
  abs_startpoint.name    = 'Startpoint';
  abs_startpoint.strtype = 'r';
  abs_startpoint.val     = {-0.5};
  abs_startpoint.num     = [1 1];
  abs_startpoint.help    = {
    'Absolute position of the start point of the grid along the surface normals in mm according to the surface. Give negative value for a start point outside of the surface (CSF direction, outwards). '
  };
  rel_startpoint = abs_startpoint;
  rel_startpoint.val     = {-0.5};
  rel_startpoint.help    = {
    'Relative position of the start point of the grid along the surface normals according to a tissue class. A value of "-0.5" begins at the GM/CSF border and even lower values define a start point outside of the tissue class (CSF direction, outwards). A value of "0" means that the central surface is used as starting point and "0.5" is related to the GM/WM border.'
  };
  
  % steps
  abs_steps         = cfg_entry;
  abs_steps.tag     = 'steps';
  abs_steps.name    = 'Steps';
  abs_steps.strtype = 'w';
  abs_steps.val     = {7};
  abs_steps.num     = [1 1];
  abs_steps.help    = {
    'Number of grid steps. '
  };
  rel_steps = abs_steps; 

  % endpoint
  abs_endpoint         = cfg_entry;
  abs_endpoint.tag     = 'endpoint';
  abs_endpoint.name    = 'Endpoint';
  abs_endpoint.strtype = 'r';
  abs_endpoint.val     = {0.5};
  abs_endpoint.num     = [1 1];
  abs_endpoint.help    = {
    'Absolute position of the end point of the grid along the surface normals (pointing inwards) in mm according to the surface. '
  };
  rel_endpoint = abs_endpoint;
  rel_endpoint.val     = {0.5};
  rel_endpoint.help    = {
    'Relative position of the end point of the grid along the surface normals (pointing inwards) according to a tissue class. A value of "0.5" ends at the GM/WM border and values > 0.5 define an end point outside of the tissue class (WM direction, inwards). A value of "0" ends at the central surface.'
  };

  % tissue class
  rel_class         = cfg_menu;
  rel_class.tag     = 'class';
  rel_class.name    = 'Tissue Class';
  rel_class.labels  = {'GM'};
  rel_class.values  = {'GM'};
  rel_class.val     = {'GM'};
  rel_class.help    = {
    'Tissue class for which the relative positions are estimated.'
  };

  % tissue class
  abs_class         = cfg_menu;
  abs_class.tag     = 'surface';
  abs_class.name    = 'Surface';
  abs_class.labels  = {'WM Surface','Central Surface','Pial Surface'};
  abs_class.values  = {'WM','Central','Pial'};
  abs_class.val     = {'Central'};
  abs_class.help    = {
    'Surface (or tissue boundary) for which the absolute positions are estimated.'
  };

  % absolute position
  abs_mapping         = cfg_branch;
  abs_mapping.tag     = 'abs_mapping';
  abs_mapping.name    = 'Absolute Grid Position From a Surface';
  abs_mapping.val   = {
    abs_class ...
    abs_startpoint ...
    abs_steps ...
    abs_endpoint ...
  }; 
  abs_mapping.help    = {
    'Map volumetric data from abolute grid positions from a surface (or tissue boundary).'
  };
  
  %% relative mapping with equi-distance approach
  rel_mapping         = cfg_branch;
  rel_mapping.tag     = 'rel_mapping';
  rel_mapping.name    = 'Relative Grid Position Within a Tissue Class (Equi-distance Model)';
  rel_mapping.val   = {
    rel_class ...
    rel_startpoint ...
    rel_steps ...
    rel_endpoint ...
  };
  rel_mapping.help    = {
    'Map volumetric data from relative grid positions within a tissue class using equi-distance approach. Here, the grid lines have equal distances in between the tissue.'
  };

  %% relative mapping with equi-volume approach
  rel_equivol_mapping         = cfg_branch;
  rel_equivol_mapping.tag     = 'rel_equivol_mapping';
  rel_equivol_mapping.name    = 'Relative Grid Position Within a Tissue Class (Equi-volume Model)';
  rel_equivol_mapping.val   = {
    rel_class ...
    rel_startpoint ...
    rel_steps ...
    rel_endpoint ...
  };
  rel_equivol_mapping.help    = {
    'Map volumetric data from relative positions within a tissue class using equi-volume approach. '
    'This option is using the approach by Bok (Z. Gesamte Neurol. Psychiatr. 12, 682-750, 1929). '
    'Here, the volume between the grids is constant. The correction is based on Waehnert et al. (NeuroImage, 93: 210-220, 2014).'
    'Please note that this option is intended for high-resolution (f)MRI data only.'
    '' 
  };

  %% -- Mapping function

  mapping         = cfg_choice;
  mapping.tag     = 'mapping';
  mapping.name    = 'Mapping Function';
  mapping.values  = {
    abs_mapping ...
    rel_mapping ...
  }; 
  mapping.val = {rel_mapping};
  mapping.help    = {
    'Volume extraction type. '
    '  Absolute Grid Position From a Surface (or Tissue Boundary):'
    '    Extract values around a surface or tissue boundary with a specified absolute sample '
    '    distance and either combine these values or save values separately.'
    '  Relative Grid Position Within a Tissue Class (Equi-distance approach):' 
    '    Extract values within a tissue class with a specified relative sample distance'
    '    with equally distributed distances and either combine these values or save values separately.'
    '' 
  };

  mapping_native = mapping;
  mapping_native.values{3} = rel_equivol_mapping;
  mapping_native.help    = {
    'Volume extraction type. '
    '  Absolute Grid Position From a Surface (or Tissue Boundary):'
    '    Extract values around a surface or tissue boundary with a specified absolute sample '
    '    distance and either combine these values or save values separately.'
    '  Relative Grid Position Within a Tissue Class (Equi-distance approach):' 
    '    Extract values within a tissue class with a specified relative sample distance'
    '    with equally distributed distances and either combine these values or save values separately.'
    '  Relative Grid Position Within a Tissue Class (Equi-volume approach):' 
    '    Extract values within a tissue class with a specified relative sample distance'
    '    that is corrected for constant volume between the grids and either combine these values or save values separetely.'
    '' 
  };



% extract volumetric data in individual space 
%-----------------------------------------------------------------------  

  data_surf_sub_lh         = cfg_files;
  data_surf_sub_lh.tag     = 'data_mesh_lh';
  data_surf_sub_lh.name    = '(Left) Individual Surfaces';
  data_surf_sub_lh.filter  = 'gifti';
  data_surf_sub_lh.ufilter = '^lh.central.(?!nofix).*';
  data_surf_sub_lh.num     = [1 Inf];
  data_surf_sub_lh.help    = {
    'Select left subject surface files (do not select the *.nofix.* surface).'
    'Right side will be automatically processed.'
    };
   
  data_sub         = cfg_files; 
  data_sub.tag     = 'data_vol';
  data_sub.name    = '(Co-registered) Volumes in Native Space';
  data_sub.filter  = 'image';
  data_sub.ufilter = '^(?!wm|wp|m0wp|mwp|wc).*'; % no normalized images
  data_sub.num     = [1 Inf];
  data_sub.help    = {
    'Select volumes in native (subject) space.'
    'Please note that these images have to be in the same space as the T1-image that was used to extract the cortical surface. An optional co-registration might be necessary if you have functional or structural data that are not yet aligned to the T1-image.'
  };

  vol2surf      = cfg_exbranch;
  vol2surf.tag  = 'vol2surf';
  vol2surf.name = 'Map Volume (Native Space) to Individual Surface';
  if expert
    vol2surf.val = {
      data_sub ...
      data_surf_sub_lh ...
      sample ...
      interp ...
      datafieldname ...
      mapping_native ...
      };
  else
    vol2surf.val = {
      data_sub ...
      data_surf_sub_lh ...
      sample ...
      datafieldname ...
      mapping_native ...
    };
  end
  vol2surf.prog = @cat_surf_vol2surf;
  vol2surf.vout = @vout_vol2surf;
  vol2surf.help = {
    'Map volume (native space) to individual surface. These mapped volumes have to be finally resampled and smoothed before any statistical analysis.'
    ''
    'The output will be named:' 
    '  [rh|lh].OutputName_VolumeName' 
    ''
  };
 
  data_surf_avg_lh         = cfg_files; 
  data_surf_avg_lh.tag     = 'data_mesh_lh';
  data_surf_avg_lh.name    = '(Left) Template Hemisphere';
  data_surf_avg_lh.filter  = 'gifti';
  data_surf_avg_lh.ufilter = '^lh.*';
  data_surf_avg_lh.num     = [1 1];
  data_surf_avg_lh.val{1}  = {fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.central.Template_T1_IXI555_MNI152_GS.gii')};
  data_surf_avg_lh.dir     = fullfile(spm('dir'),'toolbox','cat12');
  data_surf_avg_lh.help    = {
    'Select left template surface file. '
    'Right hemisphere will be automatically processed.'
    };
    
  data_norm         = cfg_files; 
  data_norm.tag     = 'data_vol';
  data_norm.name    = 'Spatially Normalized Volumes';
  data_norm.filter  = 'image';
  data_norm.ufilter = '.*';
  data_norm.num     = [1 Inf];
  data_norm.help    = {
    'Select spatially normalized volumes (in template space).'
    ''
  };

  vol2tempsurf      = cfg_exbranch;
  vol2tempsurf.tag  = 'vol2surftemp';
  vol2tempsurf.name = 'Map Normalized Volume to Template Surface';
  vol2tempsurf.val  = {
    data_norm ...
    merge_hemi ...
    mesh32k ...
    sample ...
    interp ...
    datafieldname ...
    mapping ...
  };
  vol2tempsurf.prog = @cat_surf_vol2surf;
  vol2tempsurf.vout = @vout_vol2surf;
  vol2tempsurf.help = {
    'Map spatially normalized data (in template space) to template surface.'
    'The template surface was generated by CAT12 surface processing of the average of 555 Dartel-normalized images of the IXI database that were also used to create the IXI Dartel template.   '
    ''
    'The ouput will be named:' 
    '  [rh|lh|mesh].OutputName_VolumeName.Template_T1_IXI555_MNI152_GS.gii' 
    ''
    '  WARNING: This function is primarily thought in order to project statistical results '
    '           to the template surface. Do not use the output of this function'
    '           for any statistical analysis. Statistical analysis should only use'
    '           resampled data of individual volume data (in native space) that were mapped'
    '           to individual surfaces.'
    ''
  };

%==========================================================================
function [surfcalc,surfcalcsub] = cat_surf_calc_GUI(expert)
%% surface calculations 
%  ---------------------------------------------------------------------
%  estimation per subject (individual and group sampling):
%  g groups with i datafiles and i result datafile
 
  cdata               = cfg_files;
  cdata.tag           = 'cdata';
  cdata.name          = 'Surface Data Files';
  cdata.filter        = 'any';
  cdata.ufilter       = '(lh|rh|mesh).*';
  cdata.num           = [1 Inf];
  cdata.help          = {'These are the surface data files that are used by the calculator.  They are referred to as s1, s2, s3, etc in the order they are specified.'};
  
  cdata_sub           = cfg_files;
  cdata_sub.tag       = 'cdata';
  cdata_sub.name      = 'Surface Data Files';
  cdata_sub.filter    = 'gifti';
  cdata_sub.ufilter   = '(lh|rh|mesh).(?!cent|pial|white|sphe|defe|hull).*gii';
  cdata_sub.num       = [1 Inf];
  cdata_sub.help      = {'These are the surface data files that are used by the calculator.  They are referred to as s1, s2, s3, etc in the order they are specified.'};
   
  cdata_sample        = cfg_repeat;
  cdata_sample.tag    = 'cdata_sub.';
  cdata_sample.name   = 'Surface Data Sample';
  cdata_sample.values = {cdata_sub};
  cdata_sample.num    = [1 Inf];
  cdata_sample.help   = {...
    'Specify data for each sample. All samples must have the same size and same order.'};
 
  outdir              = cfg_files;
  outdir.tag          = 'outdir';
  outdir.name         = 'Output Directory';
  outdir.filter       = 'dir';
  outdir.ufilter      = '.*';
  outdir.num          = [0 1];
  outdir.val{1}       = {''};
  outdir.help         = {
    'Files produced by this function will be written into this output directory.  If no directory is given, images will be written  to current working directory.  If both output filename and output directory contain a directory, then output filename takes precedence.'
  };
  
  dataname            = cfg_entry;
  dataname.tag        = 'dataname';
  dataname.name       = 'Output Filename';
  dataname.strtype    = 's';
  dataname.num        = [1 Inf];
  dataname.val        = {'output'};
  dataname.help       = {'The output surface data file is written to current working directory unless a valid full pathname is given.  If a path name is given here, the output directory setting will be ignored.'};

  datahandling         = cfg_menu;
  datahandling.name   = 'Sample data handling';
  datahandling.tag    = 'datahandling';
  datahandling.labels = {'subjectwise','datawise'};
  datahandling.values = {1,2};
  datahandling.val    = {1};
  datahandling.help   = {...
    'There are two ways to hand surface calculation for many subjects: (i) subjectwise or (ii) datawise.'
    ''
   ['Subject-wise means that you have multiple surface measures that should ' ...
    'be processed for multiple subjects in the same way within the native subject space. ' ...
    'E.g., To multiply thickness with area as simple approximation of cortical ' ...
    'volume, resulting in a new surface measure or texture for each subject. ' ...
    'So you have to create one sample for each surface measure and select for each of them all subject. ' ...
    'In our example, we have two samples: (i) one with thickness that include the thickness data of all subjects and ' ...
    '(ii) another one that includes the surface area data of the same subjects. ' ...
    'It is therefore important that you select exactly the same subjects in the same order! ' ... 
    'Hence, the name entry specify a the name of the the new surface measure (MEASUREENAME) in the resulting file: ']
    '  [rh|lh].TEXTURNAME[.resampled].subjectname[.gii]' 
    ''
   ['Data-wise means can be used to process each sample by the imcalc operation. ' ...
    'Hence, it require that each surface measure has the same data size - being from one subject or being resampled.' ...
    'This can be useful if you want to average the data of multiple subjects (e.g., mean(S) or std(S)) or ' ...
    'if you want to process multiple new measure of the same subject but with varying number of inputs (e.g., for multiple time points. ' ...
    'If each sample contain different subject than the SUBJECTNAME is replaced otherwise the TEXTURENAME: ']
    '  [rh|lh].TEXTURNAME[.resampled].SUBJECTNAME[.gii]'
    ''
    };
  
  expression         = cfg_entry;
  expression.tag     = 'expression';
  expression.name    = 'Expression';
  expression.strtype = 's';
  expression.num     = [1 Inf];
  expression.val     = {'s1'};
  expression.help    = {
    'Example expressions (f):'
    '  * Mean of six surface textures (select six texture files)'
    '    f = ''(s1+s2+s3+s4+s5+s6)/6'''
    '  * Make a binary mask texture at threshold of 100'
    '    f = ''(s1>100)'''
    '  * Make a mask from one texture and apply to another'
    '    f = ''s2.*(s1>100)'''
    '        - here the first texture is used to make the mask, which is applied to the second texture'
    '  * Sum of n textures'
    '    f = ''s1 + s2 + s3 + s4 + s5 + ...'''
    '  * Sum of n textures (when reading data into a data-matrix - use dmtx arg)'
    '    f = mean(S)'
    ''
  };

  dmtx         = cfg_menu;
  dmtx.tag     = 'dmtx';
  dmtx.name    = 'Data Matrix';
  dmtx.labels  = {
    'No - don''t read images into data matrix',...
    'Yes -  read images into data matrix'
  };
  dmtx.values  = {0,1};
  dmtx.val     = {0};
  dmtx.help    = {
    'If the dmtx flag is set, then textures are read into a data matrix S (rather than into separate variables s1, s2, s3,...). The data matrix should be referred to as S, and contains textures in rows. Computation is vertex by vertex, S is a NxK matrix, where N is the number of input textures, and K is the number of vertices per plane.'
  };


  % main functions
  % ----------------------------------------------------------------------

  % any image with similar properties
  surfcalc      = cfg_exbranch;
  surfcalc.tag  = 'surfcalc';
  surfcalc.name = 'Surface Calculator';
  surfcalc.val  = {
    cdata ...
    dataname ...
    outdir ...
    expression ...
    dmtx ...
  };
  surfcalc.prog = @cat_surf_calc;
  surfcalc.help = {
    'Mathematical operations for surface data (textures).'
    'It works similar to "spm_imcalc".  The input surface data must have the same number of entries (e.g. data of the same hemisphere of a subject or resampled data).'
  };
  surfcalc.vout = @(job) vout_cat_surf_calc(job);

  
  % subject- / group-w
  surfcalcsub      = cfg_exbranch;
  surfcalcsub.tag  = 'surfcalcsub';
  surfcalcsub.name = 'Surface Calculator (subject-wise)';
  if expert>1
    surfcalcsub.val  = {
      cdata_sample ...
      datahandling ...
      dataname ...
      outdir ...
      expression ...
      dmtx ...
    };
    surfcalcsub.help = {
      'Mathematical operations for surface data sets (textures).'
      'In contrast to the "Surface Calculator" it allows to apply the same expression to multiple subjects. Please note that a fixed name structure is expected: [rh|lh].TEXTURENAME[.resampled].subjectname[.gii]. Here TEXTURENAME will be replaced by the output name.'
    };
  else
    surfcalcsub.val  = {
      cdata_sample ...
      dataname ...
      outdir ...
      expression ...
      dmtx ...
    };
    surfcalcsub.help = {
      'Mathematical operations for surface data sets (textures).'
      'In contrast to the "Surface Calculator" it allows to apply the same expression to multiple subjects. Please note that a fixed name structure is expected: [rh|lh].TEXTURENAME[.resampled].subjectname[.gii]. Here TEXTURENAME will be replaced by the output name.'
    };
  end
  surfcalcsub.prog = @cat_surf_calc;
  surfcalcsub.vout = @(job) vout_cat_surf_calcsub(job);
 
%==========================================================================
function flipsides = cat_surf_flipsides_GUI
% Flipsides
%-----------------------------------------------------------------------
  cdata               = cfg_files;
  cdata.tag           = 'cdata';
  cdata.name          = 'Surface Data Files';
  cdata.filter        = 'gifti';
  cdata.ufilter       = '^s.*mm\.lh.*';
  cdata.num           = [1 Inf];
  cdata.help          = {'Texture maps that should be flipped/mirrored from right to left.' ''};
  
  cdata_sample        = cfg_repeat;
  cdata_sample.tag    = 'cdata_sub.';
  cdata_sample.name   = 'Surface Data Sample';
  cdata_sample.values = {cdata};
  cdata_sample.num    = [1 Inf];
  cdata_sample.help   = {'Specify data for each sample.'};
 
  flipsides           = cfg_exbranch;
  flipsides.tag       = 'flipsides';
  flipsides.name      = 'Flip right to left hemisphere';
  flipsides.val       = {cdata}; % replace later by cdata_sample ... update vout ...
  flipsides.prog      = @cat_surf_flipsides;
  flipsides.help      = {'This function flip the right hemisphere to the left side, to allow side' ''};
  flipsides.vout      = @(job) vout_cat_surf_flipsides(job); 

%==========================================================================
function [surfresamp,surfresamp_fs] = cat_surf_resamp_GUI(expert,nproc,merge_hemi,mesh32k,lazy)
%% Resample and smooth surfaces 
%  ------------------------------------------------------------------------

if 0 % old 
  data_surf            = cfg_files;
  data_surf.tag        = 'data_surf';
  data_surf.name       = '(Left) Surfaces Data';
  data_surf.num        = [1 Inf];
  data_surf.filter     = 'any';
  if expert > 1
    data_surf.ufilter  = '^lh.';
  else
    data_surf.ufilter  = '^lh.(?!cent|pial|white|sphe|defe|hull|pbt).*';
  end
  data_surf.help       = {'Select surfaces data files for left hemisphere for resampling to template space.'};

else
  data_surf         = cfg_files;
  data_surf.tag     = 'data_surf';
  data_surf.name    = '(Left) Surfaces Data';
  data_surf.filter  = 'any';
  if expert > 1
    data_surf.ufilter = '^lh.';
  else
    data_surf.ufilter = '^lh.(?!cent|pial|white|sphe|defe|hull|pbt).*';
  end
  data_surf.num     = [1 Inf];
  data_surf.help    = {'Select surfaces data files for left hemisphere for resampling to template space.'};

  sample_surf            = cfg_repeat;
  sample_surf.tag        = 'sample';
  sample_surf.name       = 'Data';
  sample_surf.values     = {data_surf};
  sample_surf.num        = [1 Inf];
  sample_surf.help       = {...
  'Specify data for each sample. If you specify different samples the mean correlation is displayed in separate boxplots for each sample.'};
end

  fwhm_surf         = cfg_entry;
  fwhm_surf.tag     = 'fwhm_surf';
  fwhm_surf.name    = 'Smoothing Filter Size in FWHM';
  fwhm_surf.strtype = 'r';
  fwhm_surf.num     = [1 1];
  fwhm_surf.val     = {12};
  fwhm_surf.help    = {
    'Select filter size for smoothing. For cortical thickness a good starting value is 12mm, while other surface parameters based on cortex folding (e.g. gyrification, cortical complexity) need a larger filter size of about 20-25mm. For no filtering use a value of 0.'};

  surfresamp      = cfg_exbranch;
  surfresamp.tag  = 'surfresamp';
  surfresamp.name = 'Resample and Smooth Surface Data';
  if expert
    surfresamp.val  = {sample_surf,merge_hemi,mesh32k,fwhm_surf,lazy,nproc};
  else
    surfresamp.val  = {data_surf,merge_hemi,mesh32k,fwhm_surf,nproc};
  end
  surfresamp.prog = @cat_surf_resamp;
  surfresamp.vout = @(job) vout_surf_resamp(job);
  surfresamp.help = {
  'In order to analyze surface parameters all data have to be resampled into template space and the resampled data have to be finally smoothed. Resampling is done using the warped coordinates of the resp. sphere.'};




%% Resample and smooth FreeSurfer surfaces
%  ------------------------------------------------------------------------

  data_fs         = cfg_files;
  data_fs.tag     = 'data_fs';
  data_fs.name    = 'Freesurfer Subject Directories';
  data_fs.filter  = 'dir';
  data_fs.ufilter = '.*';
  data_fs.num     = [1 Inf];
  data_fs.help    = {'Select subject folders of Freesurfer data to resample data (e.g. thickness).'};

  measure_fs         = cfg_entry;
  measure_fs.tag     = 'measure_fs';
  measure_fs.name    = 'Freesurfer Measure';
  measure_fs.strtype = 's';
  measure_fs.num     = [1 Inf];
  measure_fs.val     = {'thickness'};
  measure_fs.help    = {
    'Name of surface measure that should be resampled and smoothed.'
    ''
    };

  outdir         = cfg_files;
  outdir.tag     = 'outdir';
  outdir.name    = 'Output Directory';
  outdir.filter  = 'dir';
  outdir.ufilter = '.*';
  outdir.num     = [0 1];
  outdir.help    = {'Select a directory where files are written.'};

  surfresamp_fs      = cfg_exbranch;
  surfresamp_fs.tag  = 'surfresamp_fs';
  surfresamp_fs.name = 'Resample and Smooth Existing FreeSurfer Surface Data';
  surfresamp_fs.val  = {data_fs,measure_fs,merge_hemi,mesh32k,fwhm_surf,outdir};
  surfresamp_fs.prog = @cat_surf_resamp_freesurfer;
  surfresamp_fs.help = {
  'If you have existing Freesurfer data (e.g. thickness) this function can be used to resample these data, smooth the resampled data, and convert Freesurfer data to gifti format.'};

%==========================================================================
function roi2surf = cat_roi_roi2surf_GUI
%% roi to surface  
%  ------------------------------------------------------------------------

% ROI files
  r2s.ROIs         = cfg_files;
  r2s.ROIs.tag     = 'rdata';
  r2s.ROIs.name    = 'ROI atlas files';
  r2s.ROIs.filter  = 'xml';
  r2s.ROIs.ufilter = '.*';
  r2s.ROIs.dir     = fullfile(spm('dir'),'toolbox','cat12','templates_volumes'); 
  r2s.ROIs.num     = [1 Inf];
  r2s.ROIs.help    = {'These are the indiviudal ROI atlas files from the label directory. Choose XML files.'};

 % atlas used for extraction .. if xml 
  r2s.atlas         = cfg_entry;
  r2s.atlas.tag     = 'atlas';
  r2s.atlas.name    = 'Atlas maps';
  r2s.atlas.help    = {'Enter the name of the atlas maps, e.g. "LPBA40", or "all".'};
  r2s.atlas.strtype = 's';
  r2s.atlas.val     = {'all'};
  r2s.atlas.num     = [1 Inf];
  
% datafield in the CSV or XML file eg the GM thickness
  r2s.data         = cfg_entry;
  r2s.data.tag     = 'fields';
  r2s.data.name    = 'Datafields';
  r2s.data.help    = {'Enter the name of the data fields, e.g. "Vgm", "Igm", "Tgm", "mySurfData" or "all"'};
  r2s.data.strtype = 's';
  r2s.data.val     = {'all'};
  r2s.data.num     = [1 Inf];
  
% surface 
  r2s.surf        = cfg_menu;
  r2s.surf.name   = 'Surface type for mapping';
  r2s.surf.tag    = 'surf';
  r2s.surf.labels = {'FS average','Dartel average','subject'};
  r2s.surf.values = {'freesurfer','dartel','subject'};
  r2s.surf.val    = {'freesurfer'};
  r2s.surf.help   = {'Surface type for value projection. '};
 
% average function? - not required 
 
% main function
  roi2surf      = cfg_exbranch;
  roi2surf.tag  = 'roi2surf';
  roi2surf.name = 'Map ROI data to the surface';
  roi2surf.val  = {
    r2s.ROIs ...
    r2s.atlas ...
    r2s.data ...
    r2s.surf ...
    };
  roi2surf.prog = @cat_roi_roi2surf;
  roi2surf.help = {
    ''
  };

%==========================================================================
function surfextract = cat_surf_parameters_GUI(expert,nproc,lazy)
%% surface measures
%  ------------------------------------------------------------------------  
  
  data_surf_extract         = cfg_files;
  data_surf_extract.tag     = 'data_surf';
  data_surf_extract.name    = 'Central Surfaces';
  data_surf_extract.filter  = 'gifti';
  data_surf_extract.ufilter = '^lh.central';
  data_surf_extract.num     = [1 Inf];
  data_surf_extract.help    = {'Select left surfaces to extract values.'};
  
  % absolute mean curvature
  GI        = cfg_menu;
  if expert 
    GI.name   = 'Gyrification index (absolute mean curvature)';
    GI.labels = {'No','Unscaled','Scaled','Both'};
    GI.values = {0,1,2,3};
  else
    GI.name = 'Gyrification index';
    GI.labels = {'No','Yes'};
    GI.values = {0,1};
  end
  GI.tag    = 'GI';
  GI.val    = {1};
  GI.help   = {
    'Extract gyrification index (GI) based on absolute mean curvature. The method is described in Luders et al. NeuroImage, 29: 1224-1230, 2006.'
  };

  if expert>1
    % ---------------------------------------------------------------------
    % DFG project measures
    % ---------------------------------------------------------------------
    GIL        = cfg_menu;
    GIL.name   = 'Laplacian gyrification indices';
    GIL.tag    = 'GIL';
    GIL.labels = {'No','inward','outward','generalized','all'};
    GIL.values = {0,1,2,3,4};
    GIL.val    = {0}; 
    GIL.help   = {[ ...
      'WARNING: This GI measures is still in development and not verified yet!\n\n ' ...
      'Extraction of Laplacian-based gyrification indices as local area relation between the individual folded and the unfolded surface. ' ...
      'The Laplacian approach supports an optimal mapping between the surfaces but the classical definition of the outer hull result in a measure that focuses on inward folding. ' ...
      '\n\n' ...
      'The "inward" option use only the hull surface resulting in high sulcal values (Dahnke et al. 2010, Li et al. 2014) and is a local 3D representation of Zilles gyrification index (Zilles et al. 1989). ' ...
      'The "outward" option use only the core surface and result in high gyral values. However, the core definition is expected to be less robust for very strong tissue atrophy. ' ...
      'The "generalized" option combine both models resulting in an independent folding model that is expected to require less smoothing for surface-based analysis (15 mm). ' ...
      '\n\n' ...
      'This research part of the DFG project DA 2167/1-1 (2019/04 - 2012/04).\n\n' ...
    ]};
  
    if expert>1
      % Developer/expert parameters? 
      GILtype            = GIL; 
    
      % to small values may lead to problems in high resolution data?
      % > compensation term in cat_surf_gyrification 
      laplaceerr         = cfg_menu;
      laplaceerr.name    = 'Laplacian filter accuracy';
      laplaceerr.tag     = 'laplaceerr';
      laplaceerr.labels  = {'0.001','0.0001','0.00001'};
      laplaceerr.values  = {0.001,0.0001,0.00001};
      laplaceerr.val     = {0.0001};
      laplaceerr.help    = {'Filter accuracy of the voxel-base Laplace filter applied before streamline estimation. Smaller values are more accurate but increase processing time. \n\n'};

      % check if this compensate for voxel size!
      GIstreamopt        = cfg_menu;
      GIstreamopt.name   = 'Streamline accuracy';
      GIstreamopt.tag    = 'GIstreamopt';
      GIstreamopt.labels = {'0.01','0.001','0.0001'};
      GIstreamopt.values = {0.01,0.001,0.0001};
      GIstreamopt.val    = {0.01};
      GIstreamopt.help   = {'Stepsize of the Laplacian streamline estimation. Small values are more accurate but increase processing and memory demands. \n\n'};

      % normalization function
      GInorm             = cfg_menu;
      GInorm.name        = 'Normalization function';
      GInorm.tag         = 'GInorm';
      GInorm.labels      = {'none','2nd root','3rd root','6th root','10th root','log2','log','log10'};
      GInorm.values      = {1,2,3,6,10,'log2','log','log10'};
      GInorm.val         = {'log10'};
      GInorm.help        = {[ ...
        'Normalization function (after area filtering) to avoid the exponential increase of the GI values and to have more gaussian-like value distribution.' ...
        'To avoid negative values in case of logarithmic functions the values were corrected by the basis of the logarithmic function, e.g. +10 for log10. ' ...
        'The logarithmic function supports stronger compensation of high values. This is especially important for the inward and outward but not generalized GI. ' ...
        ]};
      
      % Relative (brain size depending) or absolute (in mm) filtering?
      % For animals relative should be more correct! 
      %
      % Comment RD20190409
      % ---------------------------------------------------------------------
      % Smoothing has to be done on the original surfaces to guarantee that 
      % the local GI is equal to the global on. 
      % We got similar problems for "resample and smooth" as for the "area"
      % measure. 
      % Maybe an "resample and smooth" (automatic) normalization option can 
      % help that guarantee the right kind of handling, e.g. global mean or
      % sum?
      % ---------------------------------------------------------------------
    
      GILfs              = cfg_entry;
      GILfs.tag          = 'GIpresmooth';
      GILfs.name         = 'Final filter size';
      GILfs.strtype      = 'r';
      GILfs.val          = {0.1}; 
      GILfs.num          = [1 1];
      GILfs.help         = {[ ...
        'Filter of the Laplacian based GI has to be done for the area variables of the folded and unfolded surfaces. ' ...
        'Values between 0 and 1 describe relative smoothing depending on the brain size, whereas values larger/equal than one describes the filter size in mm.']}; 

      % save temporary surfaces 
      GIwritehull        = cfg_menu;
      GIwritehull.name   = 'Laplacian gyrification indices hull/core surface output';
      GIwritehull.tag    = 'GIwritehull';
      GIwritehull.labels = {'No','Yes'};
      GIwritehull.values = {0,1};
      GIwritehull.val    = {1};
      GIwritehull.help   = {'Write hull/core surface use in the Laplacian gyrification index estimation.'};
      
      % hull model - This is not yet implemented! 
      % The idea behind is that the hemisphere-based hull is artificial too
      % and that are more natural definition is given by the full intra-cranial
      % volume. However, this have to include the cerebellum and 
      % brainstem as well and needs to utilize the Yp0 map!
      %
      %{
      GIhullmodel        = cfg_menu;
      GIhullmodel.name   = 'Hull model';
      GIhullmodel.tag    = 'GIhullmodel';
      GIhullmodel.labels = {'Hemisphere-based','Skull-based'};
      GIhullmodel.values = {0,1};
      GIhullmodel.val    = {0};
      GIhullmodel.help   = {'There are two different hull definitions: (i) the classical using a semi-convex hull for each hemisphere and (ii) the full intracranial volume with both hemispheres and the cerebellum. ' ''};
      %}
      
      % core model 
      GIcoremodel        = cfg_menu;
      GIcoremodel.name   = 'Core model';
      GIcoremodel.tag    = 'GIcoremodel';
      GIcoremodel.labels = {'1','2','3'};
      GIcoremodel.values = {1,2,3};
      GIcoremodel.val    = {1};
      GIcoremodel.help   = {[ ...
        'There are multiple ways to define the central core of the brain structure and it is unclear which one fits best. ' ...
        'It should have some anatomical meaning (so the ventricles would be nice as source of cortical development) ' ...
        'but being at once independent of aging (so ventricles are not real good). ' ...
        'It should be defined on an individual level (to support the scaling/normalization of head size within and between species) ' ...
        'to remove rather than increasing individual effects. ' ...
        ] '' };
      
      % threshold for core estimation 
      GIcoreth           = cfg_entry;
      GIcoreth.name      = 'Core threshold';
      GIcoreth.tag       = 'GIcoreth';
      GIcoreth.strtype   = 'r';
      GIcoreth.num       = [1 1];
      GIcoreth.val       = {0.2};
      GIcoreth.help      = {[ ...
        'Lower thresholds remove less gyri whereas high thresholds remove more gyri depending on the used core model. ' ...
        'Initial test range between 0.05 (remove only very high frequency structures) and 0.25 (remove everything that is not close to the insula). ' ...
        ] ''};
      
      % add own suffix to support multiple GI estimations
      GIsuffix         = cfg_entry;
      GIsuffix.tag     = 'GIsuffix';
      GIsuffix.name    = 'Suffix';
      GIsuffix.strtype = 's';
      GIsuffix.num     = [0 Inf];
      GIsuffix.val     = {''};
      GIsuffix.help    = {[ ...
        'Specify the string to be appended to the filenames of the filtered image file(s). ' ...
        'Default suffix is "".' ...
        ... 'Use "PARA" to add input parameters, e.g. "_lerr%1d_sacc%1d_fs%03d_smodel%d[_coreth%02d]" with ... . '
        ] ''};
    
      % main note
      GIL      = cfg_branch;
      GIL.name = GILtype.name;
      GIL.tag  = GILtype.tag;
      GIL.val  = {GILtype,GIcoremodel,GIcoreth,laplaceerr,GIstreamopt,GILfs,GInorm,GIsuffix,GIwritehull}; %GIhullmodel,
    end
  end
  if expert
    %% Thickness measures
    Tfs        = cfg_menu;
    Tfs.name   = 'Freesurfer thickness (Tfs)';
    Tfs.tag    = 'Tfs';
    Tfs.labels = {'No','Yes'};
    Tfs.values = {0,1};
    Tfs.val    = {0};
    Tfs.help   = {
      'Freesurfer thickness metric: Tfs = mean( [ Tnear(IS) Tnear(OS) ] )'
    };
  
    Tmin        = cfg_menu;
    Tmin.name   = 'Minimum thickness (Tmin)';
    Tmin.tag    = 'Tmin';
    Tmin.labels = {'No','Yes'};
    Tmin.values = {0,1};
    Tmin.val    = {0};
    Tmin.help   = {
      'Minimum thickness metric: Tmin = min( [ Tnear(IS) Tnear(OS) ] )'
    };
      
    Tmax        = cfg_menu;
    Tmax.name   = 'Maximum thickness (Tmax)';
    Tmax.tag    = 'Tmax';
    Tmax.labels = {'No','Yes'};
    Tmax.values = {0,1};
    Tmax.val    = {0};
    Tmax.help   = {
      'Maximum thickness metric: Tmax = max( [ Tnear(IS) Tnear(OS) ] )'
      'This metric is only for error diagnostic!'
    };
      
    % main note
    thickness      = cfg_branch;
    thickness.name = 'Additional thickness metrics';
    thickness.tag  = 'thickness';
    thickness.val  = {Tfs,Tmin,Tmax}; 
    thickness.help = {
      'For comparison of different thickness metrics in general see (MacDonalds et al. 1999, Lerch et al. 2005).'
      '  Tnear      .. closest point from a surface to another one'
      '  Tnormal    .. distance measured by following the surface normal (not really used)'
      '  Tfs        .. Freesurfer distance metric that is the mean of the Tnear metric of '
      '                (i) the white to the pial surface and (ii) the pial to the white surface (Fischl et al. 2000)'      
      '  Tpbt       .. voxel-based thickness metrics (Dahnke et al., 2013)'
      '  Tlaplacian .. Laplacian based thickness metric (Jones et al., 2000, Lerch et al. 2005)'
      '  Tlink      .. distance between surface with the same surface structure that is in general the result of a deformation'
      '  Tmin       .. minimum Tnear distance between two surfaces'
      '  Tmax       .. maximum Tnear distance between two surfaces'
    };
    
    %% Inner and outer surfaces
    %  -----------------------------------------------------------------
    OS        = cfg_menu;
    OS.name   = 'Pial surface';
    OS.tag    = 'OS';
    OS.labels = {'No','Yes'};
    OS.values = {0,1};
    OS.val    = {0};
    OS.help   = {
      'Creates the pial surface (outer cortical surface) by moving each vertices of the central surface by the half local thickness along the surface normal.'
    };
  
    IS        = cfg_menu;
    IS.name   = 'White matter surface';
    IS.tag    = 'IS';
    IS.labels = {'No','Yes'};
    IS.values = {0,1};
    IS.val    = {0};
    IS.help   = {
      'Creates the white matter surface (inner cortical surface) by moving each vertices of the central surface by the half local thickness along the surface normal.'
    };
  
    % main note
    surfaces      = cfg_branch;
    surfaces.name = 'Additional surfaces';
    surfaces.tag  = 'surfaces';
    surfaces.val  = {IS,OS}; 
    
  end

  % Rachels fractal dimension by spherical harmonics
  FD        = cfg_menu;
  FD.name   = 'Cortical complexity (fractal dimension)';
  FD.tag    = 'FD';
  FD.labels = {'No','Yes'};
  FD.values = {0,1};
  FD.val    = {0};
  FD.help   = {
    'Extract cortical complexity (fractal dimension) which is described in Yotter et al. Neuroimage, 56(3): 961-973, 2011.'
    ''
    'Warning: Estimation of cortical complexity is very slow!'
    ''
  };


  % sulcal depth with a hull surface that based on surface-inflation.
  % ADD cite! - VanEssen?
  SD        = cfg_menu;
  SD.name   = 'Sulcus depth';
  SD.tag    = 'SD';
  SD.val    = {1};
  SD.help   = {
    'Extract sqrt-transformed sulcus depth based on the euclidean distance between the central surface and its convex hull.'
    ''
    'Transformation with sqrt-function is used to render the data more normally distributed.'
    ''
  };
  if expert>1
    SD.labels = {'No','Unscaled','Scaled','Both'};
    SD.values = {0,1,2,3};
    SD.help   = [SD.help;{''; ['The "unscaled" version may require normalization by brain size, ' ...
      'whereas the scaled version does not.']}]; 
  else
    SD.labels = {'No','Yes'};
    SD.values = {0,1};
  end
  

  % affine normalized measures
  if expert>1
    tGI        = cfg_menu;
    tGI.name   = 'Toro''s gyrification index';
    tGI.tag    = 'tGI';
    tGI.labels = {'No','Unscaled','Scaled','Scaled extended','All'};
    tGI.values = {0,1,2,3,4};
    tGI.val    = {0};
    tGI.help   = {
      'Extract affine normalized version of Toro''s gyrification index (ntGI). The method is described in Toro et al., 2008.'
    };
  end
  
  
  % surface area ... 
  if expert>1
    area        = cfg_menu;
    area.name   = 'Surface area';
    area.tag    = 'area';
    area.labels = {'No','Unscaled','Scaled','Both'}; 
    area.values = {0,1,2,3}; 
    area.val    = {0}; 
    area.help   = {
      'WARNING: IN DEVELOPMENT!'
      'This method requires a sum-based mapping rather than the mean-based interpolation. The mapping utilize the Delaunay graph to transfer the area around a vertex to its nearest neighbor(s). See Winkler et al.,  2017. '}; 
    %{
    % Winklers method is not implemented right now (201904)
    area.help   = {
      'Extract log10-transformed local surface area using re-parameterized tetrahedral surface. The method is described in Winkler et al. NeuroImage, 61: 1428-1443, 2012.'
      ''
      'Log-transformation is used to render the data more normally distributed.'
      ''
    };
    %}
   
    gmv        = cfg_menu;
    gmv.name   = 'Surface GM volume';
    gmv.tag   = 'gmv';
    gmv.labels = {'No','Unscaled','Scaled','Both'}; 
    gmv.values = {0,1,2,3}; 
    gmv.val    = {0}; 
    gmv.help   = {
      'WARNING: NOT WORKING RIGHT NOW!'
      'This method requires a sum-based mapping rather than the mean-based interpolation. The mapping utilize the Delaunay graph to transfer the area around a vertex to its closes neighbor(s) that is also use for the area mapping. '}; 

  end
  


  % main menu
  % 
  % TODO:
  % - add a branch for measures?
  % - add a branch for options?
  surfextract      = cfg_exbranch;
  surfextract.tag  = 'surfextract';
  surfextract.name = 'Extract additional surface parameters';
  if expert == 2
    surfextract.val  = {data_surf_extract, area,gmv, GI, SD, FD, tGI, GIL, thickness, surfaces, nproc, lazy}; % area, 
  elseif expert == 1
    surfextract.val  = {data_surf_extract,GI,FD,SD, surfaces,nproc,lazy};
  else
    surfextract.val  = {data_surf_extract,GI,FD,SD,nproc};
  end
  surfextract.prog = @cat_surf_parameters;
  surfextract.vout = @vout_surfextract;
  surfextract.help = {'Additional surface parameters can be extracted that can be used for statistical analysis.'};

%==========================================================================
function renderresults = cat_surf_results_GUI(expert)
%% Batch mode for cat_surf_results function. 
%  ---------------------------------------------------------------------

  % input data
  cdata                = cfg_files;
  cdata.tag            = 'cdata';
  cdata.name           = 'Result files';
  cdata.filter         = 'any';
  cdata.ufilter        = '.*';
  cdata.num            = [1 Inf];
  cdata.help           = {'Select data for surface overlay. ' ''};
  
  
  % == render options ==
  surface              = cfg_menu;
  surface.name         = 'Render surface';
  surface.tag          = 'surface';
  surface.labels       = {'Central','Inflated','Dartel','Flatmap'};
  surface.values       = {1,2,3,4};
  surface.val          = {1};
  surface.help         = {'Select on which average surface the data should be projected. ' ''};

  texture              = cfg_menu;
  texture.name         = 'Surface underlay texture';
  texture.tag          = 'texture';
  texture.labels       = {'Mean curvature', 'Sulcal depth'};
  texture.values       = {1,2};
  texture.val          = {1};
  texture.help         = {'Select a underlaying surface texture to illustrate the cortical folding pattern by mean curvature or sulcal depth.' ''};
  
  transparency         = cfg_menu;
  transparency.name    = 'Transparency';
  transparency.tag     = 'transparency';
  transparency.labels  = {'No', 'Yes'};
  transparency.values  = {0,1};
  transparency.val     = {1};
  transparency.help    = {'Select a underlaying surface texture to illustrate the cortical folding pattern by mean curvature or sulcal depth.' ''};
  
  view                  = cfg_menu;
  view.name             = 'Render view';
  view.tag              = 'view';
  view.labels           = {'Show top view', 'Show bottom view', 'Show only lateral and medial views'};
  view.values           = {1,-1,2};
  view.val              = {1};
  view.help             = {'Select different types of surface views. The "top view"/"bottom view" shows later and medial view of the left and right hemisphere and the top/bottom view in middle. The view option has no effect in case of flatmaps.' ''};
  
  colormap              = cfg_menu;
  colormap.name         = 'Colormap';
  colormap.tag          = 'colormap';
  colormap.labels       = {'Jet', 'Hot', 'HSV', 'Cold-hot'};
  colormap.values       = {1,2,3,4};
  colormap.val          = {1};
  colormap.help         = {'Select the colormap for data visualization.' ''};
  
  colorbar              = cfg_menu;
  colorbar.name         = 'Colorbar';
  colorbar.tag          = 'colorbar';
  if expert>1
    colorbar.labels     = {'Off', 'On', 'On with histogram'};
    colorbar.values     = {0,1,2};
  else
    colorbar.labels     = {'Off', 'On'};
    colorbar.values     = {0,1};
  end
  colorbar.val          = {1};
  colorbar.help         = {'Print colorbar (with histogram).' ''};
 
  background            = cfg_menu;
  background.name       = 'Background color';
  background.tag        = 'background';
  background.labels     = {'White','Black'};
  background.values     = {1,2};
  background.val        = {1};
  background.help       = {'Select the color of the background.' ''};

  showfilename          = cfg_menu;
  showfilename.name     = 'Show filename';
  showfilename.tag      = 'showfilename';
  showfilename.labels   = {'No','Yes'};
  showfilename.values   = {0,1};
  showfilename.val      = {1};
  showfilename.help     = {'Print the filename at the top of the left and right hemisphere.' ''};
  
  invcolormap           = cfg_menu;
  invcolormap.name      = 'Inverse colormap';
  invcolormap.tag       = 'invcolormap';
  invcolormap.labels    = {'No','Yes'};
  invcolormap.values    = {0,1};
  invcolormap.val       = {0};
  invcolormap.help      = {'Use inverse colormap.' ''};
 
  if expert>1
    clims               = cfg_menu;
    clims.name          = 'Data range';
    clims.tag           = 'clims';
    clims.labels        = {'default','MSD2','MSD4','MSD8','SD2','SD4','SD8','%99.5','%99','min-max','0-max'};
    clims.values        = {'','MSD2','MSD4','MSD8','SD2','SD4','SD8','%99.5','%99','min-max','0-max'};
    clims.val           = {''};
    clims.help          = {'Define normalized datarange' ''};
  end
  
  render                = cfg_exbranch;
  render.tag            = 'render';
  render.name           = 'Render options';
  if expert>1 % ... not working ... colorbar,
    render.val          = {surface,view,texture,transparency,colormap,invcolormap,background,clims,showfilename};  
  else
    render.val          = {surface,view,texture,transparency,colormap,invcolormap,background,showfilename};  
  end
  render.prog           = @cat_surf_results;
  render.help           = {'Rendering options for the surface output.'};  
  
  
  % == stat ==
  threshold             = cfg_menu;
  threshold.name        = 'Threshold';
  threshold.tag         = 'threshold';
  threshold.labels      = {'No threshold', 'P<0.05', 'P<0.01', 'P<0.001'};
  threshold.values      = {0,1.3,2,3};
  threshold.val         = {0};
  threshold.help        = {'Define the threshold of statistical maps.' ''};
  
  hideNegRes            = cfg_menu;
  hideNegRes.name       = 'Hide negative results';
  hideNegRes.tag        = 'hide_neg';
  hideNegRes.labels     = {'No','Yes'};
  hideNegRes.values     = {0,1};
  hideNegRes.val        = {0};
  hideNegRes.help       = {'Do not show negative results.' ''};
  
  stat                  = cfg_exbranch;
  stat.tag              = 'stat';
  stat.name             = 'Statistical options';
  stat.val              = {threshold,hideNegRes};  
  stat.prog             = @cat_surf_results;
  stat.help             = {'Additional options for statistical maps.'};  
  
  
  % == filename ==
  outdir                = cfg_files;
  outdir.tag            = 'outdir';
  outdir.name           = 'Output directory';
  outdir.filter         = 'dir';
  outdir.ufilter        = '.*';
  outdir.num            = [0 1];
  outdir.help           = {'Select a directory where files are written. Empty writes in the directory of the input files.' ''};
  outdir.val{1}         = {''};

  prefix                = cfg_entry;
  prefix.tag            = 'prefix';
  prefix.name           = 'Filename prefix';
  prefix.strtype        = 's';
  prefix.num            = [0 Inf];
  prefix.val            = {'render_'};
  prefix.help           = {'Specify the string to be prepended to the filenames of the filtered image file(s). Default prefix is "".' ''};

  suffix                = cfg_entry;
  suffix.tag            = 'suffix';
  suffix.name           = 'Filename suffix';
  suffix.strtype        = 's';
  suffix.num            = [0 Inf];
  suffix.val            = {''};
  suffix.help           = {'Specify the string to be appended to the filenames of the filtered image file(s). Default suffix is "".' ''};
  
  fparts                = cfg_exbranch;
  fparts.tag            = 'fparts';
  fparts.name           = 'Filename';
  fparts.val            = {outdir,prefix,suffix};  
  fparts.prog           = @cat_surf_results;
  fparts.help           = {'Define output directory and prefix and suffix.'};  
  
  
  % == main ==
  renderresults         = cfg_exbranch;
  renderresults.tag     = 'renderresults';
  renderresults.name    = 'Render result data';
  renderresults.val     = {cdata,render,stat,fparts};  
  renderresults.prog    = @(job) cat_surf_results('batch',job);
  renderresults.vout    = @(job) vout_cat_surf_results(job);
  renderresults.help    = {'CAT result render function for automatic result image export.' ''};

%==========================================================================
% dependency functions
%==========================================================================  
function dep = vout_vol2surf(job)
%%

  if isfield(job,'merge_hemi') && job.merge_hemi
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Mapped values';
    dep(1).src_output = substruct('.','mesh');
    dep(1).tgt_spec   = cfg_findspec({{'filter','mesh','strtype','e'}});
  else
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Left mapped values';
    dep(1).src_output = substruct('.','lh');
    dep(1).tgt_spec   = cfg_findspec({{'filter','mesh','strtype','e'}});
    dep(2)            = cfg_dep;
    dep(2).sname      = 'Right mapped values';
    dep(2).src_output = substruct('.','rh');
    dep(2).tgt_spec   = cfg_findspec({{'filter','mesh','strtype','e'}});
  end
  
%==========================================================================
function dep = vout_cat_surf_results(job)
%% 

  dep(1)            = cfg_dep;
  dep(1).sname      = 'Rendered surface data';
  dep(1).src_output = substruct('()',{1},'.','png','()',{':'});
  dep(1).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});

%==========================================================================  
function dep = vout_surfextract(job)

measures = { % para-field , para-subfield , para-val , output-var, [left|right] dep-var-name 
  'GI'        ''      [1 3]  'GI'             'gyrification'; 
  'GI'        ''      [2 3]  'GIs'            'scaled gyrification'; 
  'FD'        ''       1     'FD'             'fractal dimension';
  'SD'        ''      [1 3]  'SD'             'sulcal depth';
  'SD'        ''      [2 3]  'SDs'            'scaled sulcal depth';
  ...
  'tGI'       ''      [1 4]  'tGI'            'Toro GI'; 
  'tGI'       ''      [2 4]  'tGIs'           'scaled Toro GI'; 
  'tGI'       ''      [3 4]  'tGIsx'          'exscaled Toro GI'; 
  ...
  'GIL'       ''      [1 4]  'iGI'            'inward-folding Laplacian-based GI'; 
  'GIL'       ''      [2 4]  'oGI'            'outward-folding Laplacian-based GI'; 
  'GIL'       ''      [3 4]  'gGI'            'generalized Laplacian-based GI'; 
  ... 
  'area'      ''      [1 3]  'area'           'surface area'; 
  'area'      ''      [2 3]  'areas'          'scaled surface area'; 
  'gmv'       ''      [1 2]  'gmv'            'surface GM volume'; 
  ... scaled GMV? < maybe directly normalized by TIV?
  ...
  'surfaces'  'IS'     1     'white'          'white matter surface';
  'surfaces'  'OS'     1     'pial'           'pial surface';
  'GIL'       'hull'  [1 4]  'hull'           'hull surface';
  'GIL'       'core'  [2 4]  'core'           'core surface';
  };
sides = {
  'l' 'Left';
  ...'r' 'Right';
  }; 

for si=1:size(sides,1)
  for mi=1:size(measures,1)
    if isfield(job,measures{mi,1}) && ...
        (strcmp(measures{mi,2},'hull') || strcmp(measures{mi,2},'core')) && ...
        isfield(job.(measures{mi,1}),'GIwritehull') && job.(measures{mi,1}).GIwritehull
      % special case for the hull and core surface due to the write field
      if any( job.(measures{mi,1}).(measures{mi,1}) == measures{mi,3} )
        if ~exist('dep','var'), dep = cfg_dep; else, dep(end+1) = cfg_dep; end %#ok<AGROW>
        dep(end).sname      = [sides{si,2} ' ' measures{mi,5}];
        dep(end).src_output = substruct('()',{1}, '.',[sides{si,1} 'P' measures{mi,4}],'()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
      end
    else
      if isfield(job,measures{mi,1}) && ...
       ( ( isnumeric( job.(measures{mi,1}) ) && any( job.(measures{mi,1})==measures{mi,3} ) ) || ... % no subfield
         ( isfield(job.(measures{mi,1}),measures{mi,1}) && any( job.(measures{mi,1}).(measures{mi,1})==measures{mi,3} ) ) || ... % with same subfield - GI dev. mode
         ( isfield(job.(measures{mi,1}),measures{mi,2}) && any( job.(measures{mi,1}).(measures{mi,2})==measures{mi,3} ) ) ) % with other subfield
        if ~exist('dep','var'), dep = cfg_dep; else, dep(end+1) = cfg_dep; end %#ok<AGROW>
        dep(end).sname      = [sides{si,2} ' ' measures{mi,5}];
        dep(end).src_output = substruct('()',{1}, '.',[sides{si,1} 'P' measures{mi,4}],'()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
      end
    end
  end
end
if ~exist('dep','var'), dep = cfg_dep; end

%==========================================================================
function dep = vout_surf_surf2roi(job) %#ok<INUSD>

dep(1)            = cfg_dep;
dep(1).sname      = 'Extracted Surface ROIs';
dep(1).src_output = substruct('()',{1}, '.','xmlname','()',{':'});
dep(1).tgt_spec   = cfg_findspec({{'filter','xml','strtype','e'}});

%==========================================================================
function dep = vout_surf_resamp(job)

% First we have to catch possible dependency definition problems.
if isobject(job.data_surf) && numel(job.data_surf)>1 % simple file input
  error('cat_surf_resamp:FileDepError',...
   ['CAT resample and smooth support only one dependeny per file input \n' ...
    'to support further separate processing. ']); 
elseif iscell(job.data_surf)
  nDEP = zeros(numel(job.data_surf));
  for i = 1:numel(job.data_surf)
    if isobject(job.data_surf{i}) && numel(job.data_surf{i})>1
      nDEP(i) = numel(job.data_surf{i});
    end
  end
  msg = ['The CAT batch "Resample and Smooth Surface Data" supports only one \n' ...
    'dependency per input sample to support further separate processing. '];
  for i = 1:numel(job.data_surf)
    if nDEP(i)>1
      msg = [msg sprintf('\n  Sample %d has %d dependencies.',i,nDEP(i))]; %#ok<AGROW>
    end
  end
  if sum(nDEP)>0
    error('cat_surf_resamp:SampleDepError',msg); 
  end
end
%% here we have to extract the texture name to have useful output names
mname = repmat({''},1,numel(job.data_surf));
if isobject(job.data_surf) 
  sep       = strfind(job.data_surf.sname,':');
  mname{1}  = spm_str_manip( cat_io_strrep( job.data_surf.sname(sep+1:end) , {' ','Left'} ,{'' ''}) ); 
elseif iscell(job.data_surf) 
  for i=1:numel(job.data_surf)
    if isobject(job.data_surf{i})
      sep       = strfind(job.data_surf{i}.sname,':');
      mname{i}  = spm_str_manip( cat_io_strrep( job.data_surf{i}.sname(sep+1:end) , {' ','Left'} ,{'' ''}) ); 
    else
      sinfo = cat_surf_info(job.data_surf{i});
      mname{1} = sinfo.texture; 
    end
  end
% else is empty initialization
end
%%
if iscell(job.data_surf) && ( iscell( job.data_surf{1} ) || isobject(job.data_surf{1} ) )
% this differentiation is important otherwise it has problemes with cellstrings
  for di = 1:numel(job.data_surf)
    if ~exist('dep','var'), dep = cfg_dep; else, dep(end+1) = cfg_dep; end %#ok<AGROW>
    if job.merge_hemi
      dep(end).sname      = ['Merged resampled ' mname{di}];
      dep(end).src_output = substruct('.','sample','()',{di},'.','Psdata'); %,'()',{':'});
      dep(end).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
    else
      dep(end).sname      = ['Left resampled ' mname{di}];
      dep(end).src_output = substruct('.','sample','()',{di},'.','lPsdata'); %,'()',{':'});
      dep(end).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
    end
  end
else % file input - just one texture type 
  dep = cfg_dep; 
  if job.merge_hemi
    dep(end).sname      = ['Merged resampled ' mname{1}];
    dep(end).src_output = substruct('.','sample','()',{1},'.','Psdata'); %,'()',{':'});
    dep(end).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
  else
    dep(end).sname      = ['Left resampled ' mname{1}];
    dep(end).src_output = substruct('.','sample','()',{1},'.','lPsdata'); %,'()',{':'});
    dep(end).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
  end
end
%==========================================================================  
function dep = vout_cat_surf_calc(job)
  % improve data description 
  if isobject(job.cdata) 
    % extract dependency information  
    ni = find(job.cdata.sname==':',1,'first');
    if ~isempty(ni)
      ni = ni + find(job.cdata.sname(ni+1:end)~=' ',1,'first');
      depsfnames = ['S=' job.cdata.sname(ni:end) '; ' job.expression]; 
    end
  else
    % direct input
    depsfnames = job.expression;
  end

  % create new dependency object
  dep(1)            = cfg_dep;
  dep(1).sname      = depsfnames;
  dep(1).src_output = substruct('()',{1}, '.','output','()',{':'});
  dep(1).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});

%==========================================================================  
function dep = vout_cat_surf_calcsub(job)
  depsfnames = job.expression; 
  
  for di = 1:numel(job.cdata)
    
    % extract dependency information  
    if isobject(job.cdata{di}) 
      ni = find(job.cdata{di}.sname==':',1,'first');
      if ~isempty(ni)
        ni = ni + find(job.cdata{di}.sname(ni+1:end)~=' ',1,'first');
        depsfnames = ['S=' job.cdata{di}.sname(ni:end) '; ' job.expression]; 
      end
    else
      % direct input
      depsfnames = job.expression;
    end
    
    % create new dependency object
    if ~exist('dep','var'), dep = cfg_dep; else, dep(end+1) = cfg_dep; end %#ok<AGROW>
    dep(end).sname      = depsfnames;
    dep(end).src_output = substruct('()',{1}, '.','output','()',{di});
    dep(end).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
  end
  

%==========================================================================  
function dep = vout_cat_surf_flipsides(job)
 

  dep(1)            = cfg_dep;
  dep(1).sname      = 'Flipside';
  dep(1).src_output = substruct('()',{1}, '.','files','()',{':'});
  dep(1).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});

%==========================================================================  
function dep = vout_cat_stat_spm(job)
  dep(1)            = cfg_dep;
  dep(1).sname      = 'SPM.mat File';
  dep(1).src_output = substruct('.','spmmat');
  dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

  % further files?
  
%==========================================================================  
