function tools = cat_conf_tools(expert)
% wrapper for calling CAT utilities
% 
% tools = cat_conf_tools(expert)
%_______________________________________________________________________
% Christian Gaser & Robert Dahnke
% $Id: cat_conf_tools.m 1594 2020-03-30 14:55:33Z gaser $

% multi use fields
% -----------------------------------------------------------------------

  % just used once 
  data_xml          = cfg_files;
  data_xml.name     = '';
  data_xml.tag      = 'data_xml';
  data_xml.filter   = 'xml';
  data_xml.ufilter  = '^cat_.*';
  data_xml.val      = {{''}};
  data_xml.num      = [0 Inf];
  data_xml.help     = {''};

  outdir            = cfg_files;
  outdir.tag        = 'outdir';
  outdir.name       = 'Output directory';
  outdir.filter     = 'dir';
  outdir.ufilter    = '.*';
  outdir.num        = [0 1];
  outdir.help       = {'Select a directory where files are written.'};
  outdir.val{1}     = {''};

  data              = cfg_files; 
  data.tag          = 'data';
  data.name         = 'Volumes';
  data.filter       = 'image';
  data.ufilter      = '.*';
  data.num          = [1 Inf];
  data.help         = {''};

% CHECK USE
  data_vol          = cfg_files;
  data_vol.name     = 'Sample data';
  data_vol.tag      = 'data_vol';
  data_vol.filter   = 'image';
  data_vol.num      = [1 Inf];
  data_vol.help     = {'These are the (spatially registered) data. They must all have the same image dimensions, orientation, voxel size etc. Furthermore, it is recommended to use unsmoothed files.'};
  
  % also this is a separate function that is used for the results
  spm_type          = cfg_menu; %
  spm_type.tag      = 'spm_type';
  spm_type.name     = 'Data type of the output images';
  if expert>1 
    % developer! there should be no great difference between uint# and int# due to rescaling 
    spm_type.labels = {'same','uint8','int8','uint16','int16','single'};
    spm_type.values = {0 2 256 512 4 16};
  else
    spm_type.labels = {'same','uint8','uint16','single (32 bit)'};
    spm_type.values = {0 2 512 16};
  end
  spm_type.val      = {16};
  spm_type.help     = {
    'SPM data type of the output image. Single precision is recommended, but  uint16 also provides good results. Internal scaling supports a relative high accuracy for the limited number of bits, special values such as NAN and INF (e.g. in the background) will be lost and NAN is converted to 0, -INF to the minimum, and INF to the maximum value. '
    ''
  };

  % also this limit is a separate function that is used for the noise filter
  % and therefore included here
  intlim            = cfg_entry;
  intlim.tag        = 'intlim';
  intlim.name       = 'Global intensity limitation';
  intlim.strtype    = 'r';
  intlim.num        = [1 1];
  intlim.val        = {100};
  intlim.help       = {
    'General intensity limitation to remove strong outliers by using 99.99%% of the original histogram values  before noise correction. '
    ''
  };

  prefix            = cfg_entry;
  prefix.tag        = 'prefix';
  prefix.name       = 'Filename prefix';
  prefix.strtype    = 's';
  prefix.num        = [0 Inf];
  prefix.val        = {''};
  prefix.help       = {''};

  suffix            = cfg_entry;
  suffix.tag        = 'suffix';
  suffix.name       = 'Filename suffix';
  suffix.strtype    = 's';
  suffix.num        = [0 Inf];
  suffix.val        = {''};
  suffix.help       = {''};

  fname             = prefix; 
  fname.name        = 'Filename';
  fname.tag         = 'fname';
  fname.val         = {'CATcheckdesign_'}; 
  fname.help        = {'Basic filename to save figures.'};
   
  save              = cfg_menu;
  save.name         = 'Save & close windows';
  save.tag          = 'save';
  save.labels       = {'Save & close','Save only','No'};
  save.values       = {2,1,0};
  save.val          = {0};
  save.help         = {'Save and close figures for batch processing.'};
  
  colormap              = cfg_menu;
  save.name         = 'Save & close windows';
  save.tag          = 'save';
  save.labels       = {'Save & close','Save only','No'};
  save.values       = {2,1,0};
  save.val          = {0};
  save.help         = {'Save and close figures for batch processing.'};


% get subbatches
% -------------------------------------------------------------------------
  [T2x,T2x_surf,F2x,F2x_surf] = cat_stat_T2x_GUI;
  [check_cov, check_cov2]     = cat_stat_check_cov_GUI(data_xml,outdir,fname,save,expert);
  [defs,defs2]                = cat_vol_defs_GUI;
  nonlin_coreg                = cat_conf_nonlin_coreg;
  headtrimming                = cat_vol_headtrimming_GUI(intlim,spm_type,prefix,suffix,expert);
  check_SPM                   = cat_stat_check_SPM_GUI(outdir,fname,save,expert); 
  showslice                   = cat_stat_showslice_all_GUI(data_vol);
  maskimg                     = cat_vol_maskimage_GUI(data,prefix);
  calcvol                     = cat_stat_TIV_GUI;
  spmtype                     = cat_io_volctype_GUI(data,intlim,spm_type,prefix,suffix,expert);
  calcroi                     = cat_roi_fun_GUI(outdir);
  avg_img                     = cat_vol_average_GUI(data,outdir);
  realign                     = cat_vol_series_align_GUI(data);
  sanlm                       = cat_vol_sanlm_GUI(data,intlim,spm_type,prefix,suffix,expert);
  urqio                       = cat_vol_urqio_GUI;
  long                        = cat_conf_long;
%   long1173                    = cat_conf_long1173; % EXPLOREASL HACK
%   long1445                    = cat_conf_long1445;
  iqr                         = cat_stat_IQR_GUI(data_xml);
  %qa                         = cat_vol_qa_GUI(data);
  
  
% create main batch 
% -------------------------------------------------------------------------
  tools = cfg_choice;
  tools.name   = 'Tools';
  tools.tag    = 'tools';
  tools.values = { ...
    showslice, ...                        cat.stat.pre 
    ... qa, ...                           cat.stat.pre
    check_cov, ...                        cat.stat.pre
    check_cov2, ...                       cat.stat.pre
    check_SPM, ...                        cat.stat.pre
    ...
    calcvol, ...                          cat.stat.pre
    calcroi, ...                          cat.stat.pre
    iqr, ....                             cat.stat.pre
    ...
    T2x, F2x, T2x_surf, F2x_surf, ...     cat.stat.models?
    ...
    ... SPLIT THIS FILE ?!
    ...
    sanlm, ...                            cat.pre.vtools.
    maskimg, ...                          cat.pre.vtools.
    spmtype, ...                          cat.pre.vtools.
    headtrimming, ...                     cat.pre.vtools.
    ...
    realign, ...                          cat.pre.long.?
    long, ...                             cat.pre.long.?
    ...
    nonlin_coreg, ...                     cat.pre.vtools.
    defs, ...                             cat.pre.vtools.
    defs2, ...                            cat.pre.vtools.
    avg_img, ...                          cat.pre.vtoolsexp.
    };
%     long1173, ... % EXPLOREASL HACK
%     long1445, ...
  if expert 
    tools.values = [tools.values,{urqio}]; 
  end
return

function iqr = cat_stat_IQR_GUI(data_xml)
%  ------------------------------------------------------------------------
  iqr_name         = cfg_entry;
  iqr_name.tag     = 'iqr_name';
  iqr_name.name    = 'Output file';
  iqr_name.strtype = 's';
  iqr_name.num     = [1 Inf];
  iqr_name.val     = {'IQR.txt'};
  iqr_name.help    = {'The output file is written to current working directory unless a valid full pathname is given'};

  iqr       = cfg_exbranch;
  iqr.tag   = 'iqr';
  iqr.name  = 'Get Weighted Overall Image Quality';
  iqr.val   = {data_xml,iqr_name};
  iqr.prog  = @cat_stat_IQR;
  iqr.help  = {'This function reads weighted overall image quality from saved xml-files.' ''};
return

function qa = cat_vol_qa_GUI(data) %#ok<DEFNU>
  % update input
  data.help = {'Select images for quality control.'};

  % main
  qa        = cfg_exbranch;
  qa.tag    = 'qa';
  qa.name   = 'CAT quality control';
  qa.val    = {data};
  qa.prog   = @cat_vol_qa;
  qa.vfiles = @vfiles_qa;
  qa.help   = {'CAT Quality Control of T1 images. '};
return
  
function sanlm = cat_vol_sanlm_GUI(data,intlim,spm_type,prefix,suffix,expert)
%% ------------------------------------------------------------------------

  % --- update input variables ---
  data.help         = {'Select images for filtering.'};
  
  prefix.val        = {'sanlm_'};
  prefix.help       = {
    'Specify the string to be prepended to the filenames of the filtered image file(s). Default prefix is "samlm_". Use the keyword "PARA" to add the strength of filtering, e.g. "sanlm_PARA" result in "sanlm_NC#_*.nii".'
    ''
  };

  suffix.val        = {''};
  suffix.help       = {
    'Specify the string to be appended to the filenames of the filtered image file(s). Default suffix is ''''.  Use "PARA" to add input parameters, e.g. "sanlm_*_NC#.##_RN#_RD#_RIA#.##_SR#_FSR#_RNI#_OL#.##_iterm#_iter#.nii" with NC=NCstr, RN=Rician noise, RD=resolution dependency, RIA=relative intensity adaption, SR=sub-resolutions, FSR=force sub-resolution, RNI=replace NAN and INF, and OL=outlier correction.'
    ''
  };
  

  % --- new fields ---
  rician            = cfg_menu;
  rician.tag        = 'rician';
  rician.name       = 'Rician noise';
  rician.labels     = {'Yes' 'No'};
  rician.values     = {1 0};
  rician.val        = {0};
  rician.help       = {
    'MRIs can have Gaussian or Rician distributed noise with uniform or nonuniform variance across the image. If SNR is high enough (>3) noise can be well approximated by Gaussian noise in the foreground. However, for SENSE reconstruction or DTI data a Rician distribution is expected. Please note that the Rician noise estimation is sensitive for large signals in the neighbourhood and can lead to artefacts, e.g. cortex can be affected by very high values in the scalp or in blood vessels.'
    ''
  };

  % remove artifacts
  outlier           = cfg_entry;
  outlier.tag       = 'outlier';
  outlier.name      = 'Strength of outlier correction';
  outlier.strtype   = 'r';
  outlier.num       = [1 1];
  outlier.val       = {1};
  outlier.help      = {
    'Remove strong outliers (salt and pepper noise) with more than n times of the average local correction strength. Larger values will result in stronger corrections, whereas lower values result in less corrections. Changes will be more visible in high quality areas/images.' 
  };

  if expert 
    % developer with matrix values
    NCstr           = cfg_entry;
    NCstr.tag       = 'NCstr';
    NCstr.name      = 'Strength of noise corrections';
    NCstr.strtype   = 'r';
    NCstr.num       = [1 1]; %inf]; % this case did not work with yet
    NCstr.def       = @(val) cat_get_defaults('extopts.NCstr', val{:});
    NCstr.help      = {
     ['Strength of the spatial adaptive (sub-resolution) non-local means (SANLM) noise correction. Please note that the filter strength is automatically estimated. Change this parameter only for specific conditions. ' ...
      'Typical values are: none (0), classic (1), light (2), medium (3|-inf), strong (4), heavy (5). The "classic" option use the ordinal SANLM filter without further adaptions. The "light" option uses the half filter strength of "medium" cases. The "strong" option use 8-times of the "medium" filter strength. Sub-resolution filtering is only used in case of high image resolution below 0.8 mm or in case of the "heavy" option. ' ...
      'For the global modified scheme use smaller values (>0) for less denoising, higher values (<=1) for stronger denoising, and "inf" for an automatic estimated threshold. Negative values control the local adaptive scheme, with the default "-inf"|"-1", that successfully tested on a variety of scans. Use higher values (>-1,<0) for less filtering and lower values "<-1" for stronger filtering. The value 0 will turn off any noise correction.']
      ''
    };
  end

  % noise correction level
  NCstrm            = cfg_menu;
  NCstrm.tag        = 'NCstr';
  NCstrm.name       = 'Strength of Noise Corrections';
  NCstrm.def        = @(val) cat_get_defaults('extopts.NCstr', val{:});
  NCstrm.help       = {
    ['Strength of the (sub-resolution) spatial adaptive  non local means (SANLM) noise correction. Please note that the filter strength is automatically estimated. Change this parameter only for specific conditions. ' ...
     'The "light" option applies half of the filter strength of the adaptive "medium" cases, whereas the "strong" option uses the full filter strength, force sub-resolution filtering and applies an additional iteration. Sub-resolution filtering is only used in case of high image resolution below 0.8 mm or in case of the "strong" option.']
     ''
  };
  NCstrm.values     = {2 -inf 4};
  if expert
    NCstrm.labels   = {'light (2)','medium (3|-inf)','strong (4)'};
  else
    NCstrm.labels   = {'light','medium','strong'};
  end

  addnoise          = cfg_entry;
  addnoise.tag      = 'addnoise';
  addnoise.name     = 'Strength of additional noise in noise-free regions';
  addnoise.strtype  = 'r';
  addnoise.val      = {0.5}; 
  addnoise.num      = [1 1];
  addnoise.help     = {
    'Add minimal amount of noise in regions without any noise to avoid image segmentation problems. This parameter defines the strength of additional noise as percentage of the average signal intensity. '
    ''
  };

  replaceNANandINF         = cfg_menu;
  replaceNANandINF.tag     = 'replaceNANandINF';
  replaceNANandINF.name    = 'Replace NAN and INF';
  replaceNANandINF.labels  = {'Yes' 'No'};
  replaceNANandINF.values  = {1 0};
  replaceNANandINF.val     = {1};
  replaceNANandINF.help    = {
    'Replace NAN by 0, -INF by the minimum and INF by the maximum of the image.'
    ''
    };

  % relative value vs. on/off
  if expert
    relativeFilterStengthLimit          = cfg_entry;
    relativeFilterStengthLimit.tag      = 'relativeFilterStengthLimit';
    relativeFilterStengthLimit.name     = 'Factor of relative filter strength limit';
    relativeFilterStengthLimit.strtype  = 'r';
    relativeFilterStengthLimit.num      = [1 1];
    relativeFilterStengthLimit.val      = {1};
    relativeFilterStengthLimit.help     = {
      'Limit the relative noise correction to avoid over-filtering of low intensity areas. Low values will lead to less filtering in low intensity areas, whereas high values will be closer to the original filter. INF deactivates the filter. '
      ''
    };
  else
    relativeFilterStengthLimit          = cfg_menu;
    relativeFilterStengthLimit.tag      = 'relativeFilterStengthLimit';
    relativeFilterStengthLimit.name     = 'Use relative filter strength';
    relativeFilterStengthLimit.labels   = {'Yes' 'No'};
    relativeFilterStengthLimit.values   = {1 0};
    relativeFilterStengthLimit.val      = {1};
    relativeFilterStengthLimit.help     = {
      'Limit the relative noise correction to avoid over-filtering of low intensities areas.'
      ''
      };
  end

  relativeIntensityAdaption             = cfg_entry;
  relativeIntensityAdaption.tag         = 'relativeIntensityAdaption';
  relativeIntensityAdaption.name        = 'Strength of relative intensity adaption';
  relativeIntensityAdaption.strtype     = 'r';
  relativeIntensityAdaption.num         = [1 1];
  relativeIntensityAdaption.val         = {1};
  relativeIntensityAdaption.help        = {
    'Strength of relative intensity adaption, with 0 for no adaption and 1 for full adaption. The SANLM filter is often very successful in the background and removed nearly all noise. However, routines such as the SPM Unified Segmentation expect Gaussian distribution in all regions and is troubled by regions with too low variance. Hence, a relative limitation of SANLM correction is added here that is based on the bias reduced image intensity. '
    ''
  };


  % very special parameter ...
  if expert
    iter         = cfg_entry;
    iter.tag     = 'iter';
    iter.name    = 'Number of additional sub-resolution iterations';
    iter.strtype = 'r';
    iter.num     = [1 1];
    iter.val     = {0};
    iter.help    = {
      'Choose number of additional iterations that can further reduce sub-resolution noise but also anatomical information, e.g. larger blood vessel or small gyri/sulci.'
      ''
    };

    iterm         = cfg_entry;
    iterm.tag     = 'iterm';
    iterm.name    = 'Number of additional iterations';
    iterm.strtype = 'r';
    iterm.num     = [1 1];
    iterm.val     = {0};
    iterm.help    = {
      'Choose number of additional iterations that can further reduce noise but also anatomical information, e.g. smaller blood-vessels.'
      ''
    };

    relativeIntensityAdaptionTH         = cfg_entry;
    relativeIntensityAdaptionTH.tag     = 'relativeIntensityAdaptionTH';
    relativeIntensityAdaptionTH.name    = 'Strength of smoothing of the relative filter strength limit';
    relativeIntensityAdaptionTH.strtype = 'r';
    relativeIntensityAdaptionTH.num     = [1 1];
    relativeIntensityAdaptionTH.val     = {2};
    relativeIntensityAdaptionTH.help    = {
      'Smoothing of the relative filter strength limitation.'
      ''
    };

    resolutionDependency                = cfg_menu;
    resolutionDependency.tag            = 'resolutionDependency';
    resolutionDependency.name           = 'Resolution depended filtering';
    resolutionDependency.labels         = {'Yes' 'No'};
    resolutionDependency.values         = {1 0};
    resolutionDependency.val            = {0};
    resolutionDependency.help           = {
      'Resolution depending filtering with reduced filter strength in data with low spatial resolution defined by the "Range of resolution dependency".'
      ''
      };

    resolutionDependencyRange           = cfg_entry;
    resolutionDependencyRange.tag       = 'resolutionDependencyRange';
    resolutionDependencyRange.name      = 'Range of resolution dependency';
    resolutionDependencyRange.strtype   = 'r';
    resolutionDependencyRange.num       = [1 2];
    resolutionDependencyRange.val       = {[1 2.5]};
    resolutionDependencyRange.help      = {
      'Definition of the spatial resolution for "full filtering" (first value) and "no filtering" (second value), with [1 2.5] for typical structural data of humans. '
      ''
    };

    resolutionReduction                 = cfg_menu;
    resolutionReduction.tag             = 'red';
    resolutionReduction.name            = 'Low resolution filtering';
    resolutionReduction.labels          = {'Yes (allways)' 'Yes (only highres <0.8 mm)' 'No'};
    resolutionReduction.values          = {11 1 0};
    resolutionReduction.val             = {0};
    resolutionReduction.help            = {
      'Some MR images were interpolated or use a limited frequency spectrum to support higher spatial resolution with acceptable scan-times (e.g., 0.5x0.5x1.5 mm on a 1.5 Tesla scanner). However, this can result in "low-frequency" noise that can not be handled by the standard NLM-filter. Hence, an additional filtering step is used on a reduces resolution. As far as filtering of low resolution data will also remove anatomical information the filter use by default maximal one reduction with a resolution limit of 1.6 mm. I.e. a 0.5x0.5x1.5 mm image is reduced to 1.0x1.0x1.5 mm, whereas a 0.8x0.8x0.4 mm images is reduced to 0.8x0.8x0.8 mm and a 1x1x1 mm dataset is not reduced at all. '
      ''
      };

    verb                                = cfg_menu;
    verb.tag                            = 'verb';
    verb.name                           = 'Verbose output';
    verb.labels                         = {'No' 'Yes'};
    verb.values                         = {0 1};
    verb.val                            = {1};
    verb.help                           = {
      'Be verbose.'
      ''
      };
  end

  
  
  nlm_default           = cfg_branch;
  nlm_default.tag       = 'classic';
  nlm_default.name      = 'Classic SANLM filter';
  nlm_default.val       = {};
  nlm_default.help      = {
      'Classical SANLM filter without further adaptions, i.e. strong filtering on the full resolution.' 
  }; 

  nlm_optimized         = cfg_branch;
  nlm_optimized.tag     = 'optimized';
  nlm_optimized.name    = 'Optimized filter';
  nlm_optimized.val     = {NCstrm};
  nlm_optimized.help    = {
      'Optimized SANLM filter with adaptive predefined parameters for simplified GUI cases.' 
  }; 

  if expert
    nlm_expert          = cfg_branch;
    nlm_expert.tag      = 'expert';
    nlm_expert.name     = 'Optimized filter (expert options)';
    nlm_expert.val      = {NCstr iter iterm outlier relativeIntensityAdaption relativeIntensityAdaptionTH relativeFilterStengthLimit resolutionDependency resolutionDependencyRange resolutionReduction};
    nlm_expert.help     = {
        'Optimized SANLM filter with all parameters.' 
    }; 
  end

  nlmfilter             = cfg_choice;
  nlmfilter.tag         = 'nlmfilter';
  nlmfilter.name        = 'Filter type';
  if expert
    nlmfilter.values    = {nlm_default nlm_optimized nlm_expert};
  else
    nlmfilter.values    = {nlm_default nlm_optimized};
  end
  if cat_get_defaults('extopts.NCstr')>0 && cat_get_defaults('extopts.NCstr')<=1
    nlmfilter.val       = {nlm_default}; 
  elseif expert
    nlmfilter.val       = {nlm_expert}; 
  else
    nlmfilter.val       = {nlm_optimized}; 
  end  
  if expert
    nlmfilter.help      = {
      'Selection between the classical SANLM filter and an optimized SANLM filter with simplyfied and complex setting.' 
      ''
    }; 
  else
    nlmfilter.help      = {
      'Selection between the classical SANLM filter and an optimized SANLM filter.' 
      ''
    };
  end
  sanlm                 = cfg_exbranch;
  sanlm.tag             = 'sanlm';
  sanlm.name            = 'Spatially adaptive non-local means (SANLM) denoising filter';
  if expert>1 % developer
    sanlm.val           = {data spm_type prefix suffix intlim addnoise rician replaceNANandINF nlmfilter};
  elseif expert
    sanlm.val           = {data spm_type prefix suffix addnoise rician replaceNANandINF nlmfilter};
  else
    sanlm.val           = {data spm_type prefix nlmfilter};
  end
  sanlm.prog            = @cat_vol_sanlm;
  sanlm.vout            = @vfiles_sanlm;
  sanlm.help            = {
    'This function applies an spatial adaptive (sub-resolution) non-local means denoising filter to the data. This filter will remove noise while preserving edges. The filter strength is automatically estimated based on the standard deviation of the noise. '
    ''
    'This filter is internally used in the segmentation procedure anyway. Thus, it is not necessary (and not recommended) to apply the filter before segmentation.'
    ''
  };


return

function spmtype = cat_io_volctype_GUI(data,  intlim,  spm_type,prefix,suffix,expert)
  % update variables 
  data.help           = {'Select images for data type conversion';''};
  
  intlim.tag          = 'range';
  
  prefix.val          = {'PARA'};
  prefix.help         = {
    'Specify the string to be prepended to the filenames of the converted image file(s). Default prefix is "PARA" that is replaced by the chosen datatype.'
    ''
  };
  
  suffix.help         = {
    'Specify the string to be prepended to the filenames of the converted image file(s). Default prefix is ''''. Use "PARA" to add the datatype to the filename.'
    ''
  };

  spm_type.labels(1)  = []; % remove native case
  spm_type.values(1)  = []; % remove native case
  spm_type.tag        = 'ctype';
  
  % new
  spmtype             = cfg_exbranch;
  spmtype.tag         = 'spmtype';
  spmtype.name        = 'Image data type converter'; 
  if expert
    spmtype.val       = {data spm_type prefix suffix intlim};
  else
    spmtype.val       = {data spm_type prefix intlim};
  end
  spmtype.prog        = @cat_io_volctype;
  spmtype.vout        = @vfiles_volctype;
  spmtype.help        = {
    'Convert the image data type to reduce disk-space.'
    'Uses 99.99% of the main intensity histogram to avoid problems due to outliers. Although the internal scaling supports a relative high accuracy for the limited number of bits, special values such as NAN and INF will be lost!'
    ''
  };
return

function headtrimming = cat_vol_headtrimming_GUI(intlim,spm_type,prefix,suffix,expert)

  % update input variables
  intlim1               = intlim;
  intlim1.tag           = 'intlim1';
  intlim1.name          = 'Global intensity limitation for masking';
  intlim1.val           = {90};
  intlim1.help          = {'General intensity limitation to remove strong outliers by using 90% of the original histogram values. To high values will include background noise and do not allow trimming, whereas to low values will cut objects with low values (e.g. by image inhomogeneities). ' ''};

  prefix.val            = {'trimmed_'};

  % many subjects
  simages               = cfg_files;
  simages.tag           = 'simages';
  simages.name          = 'Source images';
  simages.help          = {'Select images for trimming (e.g. T1 images).' ''};
  simages.filter        = 'image';
  simages.ufilter       = '.*';
  simages.num           = [1 Inf];

  images1               = cfg_files;
  images1.tag           = 'oimages';
  images1.name          = 'Images';
  images1.help          = {'Select other images that should be trimmed similar to the source images (e.g. coregistrated images).' ''};
  images1.filter        = 'image';
  images1.ufilter       = '.*';
  images1.num           = [1 Inf];

  oimages               = cfg_repeat;
  oimages.tag           = 'oimages';
  oimages.name          = 'Other images';
  oimages.help          = {'Select other images that should be trimmed similar to the source images. For example, the source images are a set of T1 images, whereas the second set may be a set of coregistered images of the same subjects with the same image dimensions.' ''};
  oimages.values        = {images1};
  oimages.val           = {};
  oimages.num           = [0 Inf];

  manysubjects          = cfg_branch;
  manysubjects.tag      = 'manysubjects';
  manysubjects.name     = 'Many subjects';
  manysubjects.val      = {simages oimages};
  manysubjects.help     = {
    'Create stacks of images of one class that include the same number of many subjects:'
    '  { {S1T1, S2T1,...} {S1T2, S2T2, ...} ... }'
    ''
  };

  % manyimages 
  subjectimages         = cfg_files;
  subjectimages.tag     = 'subjectimages';
  subjectimages.name    = 'Subject';
  subjectimages.help    = {
    'Select all images of one subject that are in the same space and should be trimmed together.'
  };
  if expert
    subjectimages.help  = [ subjectimages.help; {
    'The first image is used to estimate the trimming. ' 
    ''
    }];
  else
    subjectimages.help  = [ subjectimages.help; {
      'In general the first image is used to estimate the trimming (see "Average images" option). ' 
      ''
    }];
  end
  subjectimages.filter  = 'image';
  subjectimages.ufilter = '.*';
  subjectimages.num     = [1 Inf];

  manyimages            = cfg_repeat;
  manyimages.tag        = 'manyimages';
  manyimages.name       = 'Many images';
  manyimages.help       = {
    'Collect images of each subject that should be trimmed together and are in the same space.' 
    '  { {S1T1, S1T2,...} {S2T1, S2T2, ...}  {S2T1, S2T2 } ... }'
    ''};
  manyimages.values     = {subjectimages};
  manyimages.val        = {};
  manyimages.num        = [1 Inf];

  % image selection type
  timages               = cfg_choice;
  timages.tag           = 'image_selector';
  timages.name          = 'Select type of image selection'; 
  timages.values        = {manyimages manysubjects}; 
  timages.val           = {manyimages};
  timages.help          = {
    'Select "many images" if you have a small number of subjects with a VARYING number of images.'
    'Select "many subjects" if you have a large number of subject with the SAME number of images.'
  };


  pth                   = cfg_entry;
  pth.tag               = 'pth';
  pth.name              = 'Percentual trimming threshold';
  pth.strtype           = 'r';
  pth.num               = [1 1];
  pth.val               = {0.4};
  pth.help              = {'Percentual treshold for trimming. Lower values will result in a wider mask, ie. more air, whereas higher values will remove more air but maybe also brain regions with very low intensity.' ''};

  avg                   = cfg_entry;
  avg.tag               = 'avg';
  avg.name              = 'Average images';
  avg.strtype           = 'r';
  avg.num               = [1 1];
  avg.val               = {0};
  avg.help              = {'By default, only the source image is used for masking. However sometimes it is helpful to use the average of all additional images (avg=1) or the first n images inclusive th source image (avg>1) of the given set.' ''};

  open                  = cfg_entry;
  open.tag              = 'open';
  open.name             = 'Size of morphological opening of the mask';
  open.strtype          = 'n';
  open.num              = [1 1];
  open.val              = {2};
  open.help             = {'The morphological opening of the mask allows to avoid problems due to noise in the background. However, too large opening will also remove the skull or parts of the brain.' ''};

  addvox                = cfg_entry;
  addvox.tag            = 'addvox';
  addvox.name           = 'Add voxels around mask';
  addvox.strtype        = 'w';
  addvox.num            = [1 1];
  addvox.val            = {2};
  addvox.help           = {'Add # voxels around the original mask to avoid to hard masking.' ''};

  mask                  = cfg_menu;
  mask.name             = 'Final masking with source image';
  mask.tag              = 'mask';
  mask.labels           = {'Yes','No'};
  mask.values           = {1,0};
  mask.val              = {1};
  mask.help             = {'Use source image for trimming and final masking (e.g. for skull-stripping in longitudinal pipeline).'};

  
  % --- main ---
  headtrimming         = cfg_exbranch;
  headtrimming.tag     = 'datatrimming';
  headtrimming.name    = 'Image data trimming'; 
  if expert
    headtrimming.val   = {timages prefix mask suffix intlim1 pth avg open addvox spm_type intlim};
  else
    headtrimming.val   = {timages prefix mask pth spm_type};
  end
  headtrimming.prog    = @cat_vol_headtrimming;
  headtrimming.vout    = @vfiles_headtrimming;
  headtrimming.help    = {
    'Remove air around the head and convert the image data type to save disk-space but also to reduce memory-space and load/save times. Corresponding images have to have the same image dimenions. '
    'Uses 99.99% of the main intensity histogram to avoid problems due to outliers. Although the internal scaling supports a relative high accuracy for the limited number of bits, special values such as NAN and INF will be lost!'
    ''
  };
return

function maskimg = cat_vol_maskimage_GUI(data,prefix)
%------------------------------------------------------------------------
 
  % update input variables
  data.name       = 'Select images';
  data.help       = {'Select images for lesion or brain masking';''};
 
  prefix.val      = {'msk_'};
  prefix.help     = {
    'Specify the string to be prepended to the filenames of the masked image file(s).'
    ''
  };

  % lesion mask
  mask = data; 
  mask.tag        = 'mask';
  mask.name       = 'Select lesion mask images';
  mask.help       = {'Select (additional) lesion mask images that describe the regions that should be set to zero.';''};
  mask.num        = [0 Inf];
  
  % brain mask
  bmask = data; 
  bmask.tag       = 'bmask';
  bmask.name      = 'Optionally select additional brain mask images';
  bmask.help      = {'Select (additional) brain mask images that describe the regions that should remain in the image.';''};
  bmask.num       = [0 Inf];
  bmask.val       = {{''}}; 
  
  % recalc
  recalc          = cfg_menu;
  recalc.tag      = 'recalc';
  recalc.name     = 'Reprocess';
  recalc.help     = {'If an output image already exist then use this image rather than the original input image for additional masking. This allows you to add lesions from other lesion images.'};
  recalc.labels   = {'Yes' 'No'};
  recalc.values   = {1 0};
  recalc.val      = {1};
  
  % main
  maskimg         = cfg_exbranch;
  maskimg.tag     = 'maskimg';
  maskimg.name    = 'Manual image (lesion) masking'; 
  maskimg.val     = {data mask bmask recalc prefix};
  maskimg.prog    = @cat_vol_maskimage;
  maskimg.vout    = @vfiles_maskimg;
  maskimg.help    = {
    'Mask images to avoid segmentation and registration errors in brain lesion. The number of mask images has to be equal to the number of the original images. Voxels inside the lesion mask(s) and outside the brainmask(s) will be set to zero. '
    'If you have multiple lesion masks than add them with the original images, eg. "images = {sub01.nii; sub02.nii; sub01.nii}" and "mask = {sub01_lesion1.nii; sub02_lesion1.nii; sub01_lesion2.nii}". Alternatively, you can choose only one original image and a various number of mask files.'
    ''
  };

return

function [defs,defs2] = cat_vol_defs_GUI()
  %------------------------------------------------------------------------

  field           = cfg_files;
  field.tag       = 'field';
  field.name      = 'Deformation Field';
  field.filter    = 'image';
  field.ufilter   = '^(i)?y_.*\.nii$';
  field.num       = [1 Inf];
  field.help      = {[
    'Deformations can be thought of as vector fields. These can be represented by three-volume images.' ...
    'Use the "y_*.nii" to project data from subject to template space, and the "iy_*.nii" to map data from template to individual space.' ...
    'Both deformation maps can be created in the CAT preprocessing by setting the "Deformation Field" flag (no written by default).' ... 
  ]};

  field1          = cfg_files;
  field1.tag      = 'field1';
  field1.name     = 'Deformation Field';
  field1.filter   = 'image';
  field1.ufilter  = '^(i)?y_.*\.nii$';
  field1.num      = [1 1];
  field1.help     = {[
    'Deformations can be thought of as vector fields. These can be represented by three-volume images.' ...
    'Use the "y_*.nii" to project data from subject to template space, and the "iy_*.nii" to map data from template to individual space.' ...
    'Both deformation maps can be created in the CAT preprocessing by setting the "Deformation Field" flag (no written by default).' ... 
  ]};

  images1         = cfg_files;
  images1.tag     = 'images';
  images1.name    = 'Images';
  images1.help    = {'Select images to be warped. Note that there should be the same number of images as there are deformation fields, such that each flow field warps one image.'};
  images1.filter  = 'image';
  images1.ufilter = '.*';
  images1.num     = [1 Inf];

  images          = cfg_repeat;
  images.tag      = 'images';
  images.name     = 'Images';
  images.help     = {'The flow field deformations can be applied to multiple images. At this point, you are choosing how many images each flow field should be applied to.'};
  images.values   = {images1 };
  images.num      = [1 Inf];

  interp          = cfg_menu;
  interp.name     = 'Interpolation';
  interp.tag      = 'interp';
  interp.labels   = {
    'Nearest neighbour','Trilinear','2nd Degree B-spline',...
    '3rd Degree B-Spline ','4th Degree B-Spline ','5th Degree B-Spline',...
    '6th Degree B-Spline','7th Degree B-Spline','Categorical'};
  interp.values   = {0,1,2,3,4,5,6,7,-1};
  interp.val      = {1};
  interp.help     = {
    'The method by which the images are sampled when being written in a different space.'
    '    Nearest Neighbour:     - Fastest, but not normally recommended.'
    '    Bilinear Interpolation:     - OK for PET, or realigned fMRI.'
    '    B-spline Interpolation:     - Better quality (but slower) interpolation/* \cite{thevenaz00a}*/, especially with higher degree splines.  Can produce values outside the original range (e.g. small negative values from an originally all positive image). Do not use B-splines when there is any region of NaN or Inf in the images. '
    '    Categorical Interpolation:  - Slow (particularly when there are lots of categories). This is intended to warp categorical images such as label maps.'
  }';

  modulate        = cfg_menu;
  modulate.tag    = 'modulate';
  modulate.name   = 'Modulate image (preserve volume)';
  modulate.labels = {'No','Affine + non-linear (SPM12 default)','Non-linear only'};
  modulate.values = {0 1 2};
  modulate.val    = {0};
  modulate.help   = {
    '"Modulation" is to compensate for the effect of spatial normalisation. Spatial normalisation causes volume changes due to affine transformation (global scaling) and non-linear warping (local volume change). The SPM default is to adjust spatially normalised grey matter (or other tissue class) by using both terms and the resulting modulated images are preserved for the total amount of grey matter. Thus, modulated images reflect the grey matter volumes before spatial normalisation. However, the user is often interested in removing the confound of different brain sizes and there are many ways to apply this correction. We can use the total amount of GM, GM+WM, GM+WM+CSF, or manual estimated total intracranial volume (TIV). Theses parameters can be modeled as nuisance parameters (additive effects) in an AnCova model or used to globally scale the data (multiplicative effects): '
    ''
    '% Correction   Interpretation'
    '% ----------   --------------'
    '% nothing      absolute volume'
    '% globals 	    relative volume after correcting for total GM or TIV (multiplicative effects)'
    '% AnCova 	    relative volume that can not be explained by total GM or TIV (additive effects)'
    ''
    'Modulated images can be optionally saved by correcting for non-linear warping only. Volume changes due to affine normalisation will be not considered and this equals the use of default modulation and globally scaling data according to the inverse scaling factor due to affine normalisation. I recommend this option if your hypothesis is about effects of relative volumes which are corrected for different brain sizes. This is a widely used hypothesis and should fit to most data. The idea behind this option is that scaling of affine normalisation is indeed a multiplicative (gain) effect and we rather apply this correction to our data and not to our statistical model. These modulated images are indicated by "m0" instead of "m". '
    ''
  };

  defs            = cfg_exbranch;
  defs.tag        = 'defs';
  defs.name       = 'Apply deformations (many images)';
  defs.val        = {field1,images1,interp,modulate};
  defs.prog       = @cat_vol_defs;
  defs.vfiles     = @vfiles_defs;
  defs.help       = {'This is a utility for applying a deformation field of one subject to many images.'};

  defs2           = cfg_exbranch;
  defs2.tag       = 'defs2';
  defs2.name      = 'Apply deformations (many subjects)';
  defs2.val       = {field,images,interp,modulate};
  defs2.prog      = @cat_vol_defs;
  defs2.vfiles    = @vfiles_defs2;
  defs2.help      = {'This is a utility for applying deformation fields of many subjects to images.'};
return

function realign  = cat_vol_series_align_GUI(data)
  
  data.help       = {
  'Select all images for this subject'};

  tim             = cfg_entry;
  tim.tag         = 'times';
  tim.name        = 'Times';
  tim.strtype     = 'e';
  tim.val         = {NaN};
  tim.num         = [1 Inf];
  tim.help        = {'Specify the times of the scans in years. If you leave the default NaN value the standard warping regularization will be used for all scans.'};

  bparam          = cfg_entry;
  bparam.tag      = 'bparam';
  bparam.name     = 'Bias Regularisation';
  bparam.help     = {
    'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.'
    ''
    'An important issue relates to the distinction between variations in the difference between the images that arise because of the differential bias artifact due to the physics of MR scanning, and those that arise due to shape differences.  The objective is to model the latter by deformations, while modelling the former with a bias field. We know a priori that intensity variations due to MR physics tend to be spatially smooth. A more accurate estimate of a bias field can be obtained by including prior knowledge about the distribution of the fields likely to be encountered by the correction algorithm. For example, if it is known that there is little or no intensity non-uniformity, then it would be wise to penalise large estimates of the intensity non-uniformity.'
    'Knowing what works best should be a matter of empirical exploration, as it depends on the scans themselves.  For example, if your data has very little of the artifact, then the bias regularisation should be increased.  This effectively tells the algorithm that there is very little bias in your data, so it does not try to model it.'
  }';
  bparam.strtype  = 'e';
  bparam.num      = [1 1];
  bparam.val      = {1e7};

  wparam          = cfg_entry;   
  wparam.tag      = 'wparam';
  wparam.name     = 'Warping Regularisation';
  wparam.help     = {
    'Registration involves simultaneously minimising two terms.  One of these is a measure of similarity between the images (mean-squared difference in the current situation), whereas the other is a measure of the roughness of the deformations.  This measure of roughness involves the sum of the following terms:',...
    '* Absolute displacements need to be penalised by a tiny amount.  The first element encodes the amount of penalty on these.  Ideally, absolute displacements should not be penalised, but it is often necessary for technical reasons.',...
    '* The `membrane energy'' of the deformation is penalised (2nd element), usually by a relatively small amount. This penalises the sum of squares of the derivatives of the velocity field (ie the sum of squares of the elements of the Jacobian tensors).',...
    '* The `bending energy'' is penalised (3rd element). This penalises the sum of squares of the 2nd derivatives of the velocity.',...
    '* Linear elasticity regularisation is also included (4th and 5th elements).  The first parameter (mu) is similar to that for linear elasticity, except it penalises the sum of squares of the Jacobian tensors after they have been made symmetric (by averaging with the transpose).  This term essentially penalises length changes, without penalising rotations.',...
    '* The final term also relates to linear elasticity, and is the weight that denotes how much to penalise changes to the divergence of the velocities (lambda).  This divergence is a measure of the rate of volumetric expansion or contraction.',...
    'Note that regularisation is specified based on what is believed to be appropriate for a year of growth.  The specified values are divided by the number of years time difference.' 
  };
  wparam.strtype  = 'e';
  wparam.num      = [1 5];
  wparam.val      = {[0 0 100 25 100]};
  % Change to (eg): wparam.val     = {[0 0 100 25 12]};

  write_rimg          = cfg_menu;
  write_rimg.tag      = 'write_rimg';
  write_rimg.name     = 'Save rigidly registered images';
  write_rimg.help     = {'Do you want to save the rigidly registered images? The resliced images are named the same as the originals, except that they are prefixed by ''r''.'};
  write_rimg.labels   = {'Save','Dont save'};
  write_rimg.values   = { 1 0 };
  write_rimg.val      = {1};

  write_avg           = cfg_menu;
  write_avg.tag       = 'write_avg';
  write_avg.name      = 'Save Mid-point average';
  write_avg.help      = {'Do you want to save the mid-point average template image? This is likely to be useful for groupwise alignment, and is prefixed by ''avg_'' and written out in the same directory of the first time point data. Please note that with rigid registration the median is stored instead of the average image.'};
  write_avg.labels    = {'Save','Dont save'};
  write_avg.values    = { 1 0 };
  write_avg.val       = {1};

  write_jac           = cfg_menu;
  write_jac.tag       = 'write_jac';
  write_jac.name      = 'Save Jacobians';
  write_jac.help      = {'Do you want to save a map of the Jacobian determinants?  Some consider these useful for morphometrics (although the divergences of the initial velocities may be preferable). Each map of Jacobians encodes the relative volume (at each spatial location) between the scan and the median time-point average. Values less than one indicate contraction (over time), whereas values greater than one indicate expansion.  These files are prefixed by ``j_'''' and written out in the same directory of the first time point data.'};
  write_jac.labels    = {'Save','Dont save'};
  write_jac.values    = { 1 0 };
  write_jac.val       = {1};

  write_def           = cfg_menu;
  write_def.tag       = 'write_def';
  write_def.name      = 'Deformation Fields';
  write_def.help      = {'Deformation fields can be saved to disk, and used by the Deformations Utility. Deformations are saved as y_*.nii files, which contain three volumes to encode the x, y and z coordinates.  They are written in the same directory as the corresponding image.'};
  write_def.labels    = {'Save','Dont save'};
  write_def.values    = { 1 0 };
  write_def.val       = {0};

  use_brainmask         = cfg_menu;
  use_brainmask.name    = 'Use Brainmask';
  use_brainmask.tag     = 'use_brainmask';
  use_brainmask.labels  = {'Yes','No'};
  use_brainmask.values  = {1,0};
  use_brainmask.val     = {1};
  use_brainmask.help    = {'Use brainmask at last level of rigid body registration to obtain better registration by considering brain regions only.'};

  reduce                = cfg_menu;
  reduce.name           = 'Reduce Bounding Box';
  reduce.tag            = 'reduce';
  reduce.labels         = {'Yes','No'};
  reduce.values         = {1,0};
  reduce.val            = {1};
  reduce.help           = {
    'Reduce bounding box at final resolution level because usually there is a lot of air around the head after registration of multiple scans. This helps to save memory and time for later use of these registered images.'
    ''
    'Please note that this option can only be used for rigid registration and will be disabled for non-linear registration.'
  };

  nonlin                = cfg_branch;
  nonlin.tag            = 'nonlin';
  nonlin.name           = 'Non-linear registration';
  nonlin.val            = {tim wparam write_jac write_def};
  nonlin.help           = {''};

  rigid                 = cfg_const;
  rigid.tag             = 'rigid';
  rigid.name            = 'Rigid body registration';
  rigid.val             = {1};
  rigid.help            = {'Rigid registration only'};

  reg                   = cfg_choice;
  reg.name              = 'Registration Method';
  reg.tag               = 'reg';
  reg.values            = {rigid nonlin};
  reg.val               = {rigid};
  reg.help              = {'Choose between rigid body and non-linear registration. The non-linear registration is using the methods of the Longitudinal Toolbox and is thought for data over longer periods, where the deformations can then be used to calculate local volume changes, which are then multiplied (modulated) by the segmented mean image. Rigid body registration can be used to detect more subtle effects over shorter periods of time (e.g. brain plasticity or training effects after a few weeks or even shorter times).'};

  noise                 = cfg_entry;
  noise.tag             = 'noise';
  noise.name            = 'Noise Estimate';
  noise.help            = {'.'};
  noise.strtype         = 'e';
  noise.num             = [Inf Inf];
  noise.val             = {NaN};
  noise.help            = {'Specify the standard deviation of the noise in the images.  If a scalar is entered, all images will be assumed to have the same level of noise.  For any non-finite values, the algorithm will try to estimate the noise from fitting a mixture of two Rician distributions to the intensity histogram of each of the images, and assuming that the Rician with the smaller overall intensity models the intensity distribution of air in the background. This works reasonably well for simple MRI scans, but less well for derived images (such as averages) and it fails badly for scans that are skull-stripped.  The assumption used by the registration is that the residuals, after fitting the model, are i.i.d. Gaussian. The assumed standard deviation of the residuals is derived from the estimated Rician distribution of the air.'
  };

  realign               = cfg_exbranch;
  realign.tag           = 'series';
  realign.name          = 'Longitudinal Registration';
  realign.val           = {data noise bparam use_brainmask reduce reg write_rimg write_avg};
  realign.help          = {
    'Longitudinal registration of series of anatomical MRI scans for a single subject.  It is based on inverse-consistent alignment among each of the subject''s scans, and incorporates a bias field correction.  Prior to running the registration, the scans should already be in very rough alignment, although because the model incorporates a rigid-body transform, this need not be extremely precise.  Note that there are a bunch of hyper-parameters to be specified.  If you are unsure what values to take, then the defaults should be a reasonable guess of what works.  Note that changes to these hyper-parameters will impact the results obtained.'
    ''
    'The alignment assumes that all scans have similar resolutions and dimensions, and were collected on the same (or very similar) MR scanner using the same pulse sequence.  If these assumption are not correct, then the approach will not work as well.  There are a number of settings (noise estimate, regularisation etc). Default settings often work well, but it can be very helpful to try some different values, as these can have a large effect on the results.'
    ''
    'The resliced images are named the same as the originals, except that they are prefixed by ''r''.'
  };
  realign.prog          = @cat_vol_series_align;
  realign.vout          = @vout_reslice;

return

function [T2x,T2x_surf,F2x,F2x_surf] = cat_stat_T2x_GUI

  data_T2x          = cfg_files;
  data_T2x.tag      = 'data_T2x';
  data_T2x.name     = 'Data';
  data_T2x.filter   = {'image'};
  data_T2x.ufilter  = '^spmT.*';
  data_T2x.num      = [1 Inf];
  data_T2x.help     = {'Select spmT-data to transform or convert.'};

  sel               = cfg_menu;
  sel.name          = 'Convert t value to';
  sel.tag           = 'sel';
  sel.labels        = {'p','-log(p)','correlation coefficient cc','apply thresholds without conversion'};
  sel.values        = {1,2,3,5};
  sel.val           = {2};
  sel.help          = {'Select conversion of t-value'};

  thresh05          = cfg_entry;
  thresh05.tag      = 'thresh05';
  thresh05.name     = 'Threshold';
  thresh05.help     = {''};
  thresh05.strtype  = 'r';
  thresh05.num      = [1 1];
  thresh05.val      = {0.05};

  thresh001         = cfg_entry;
  thresh001.tag     = 'thresh001';
  thresh001.name    = 'Threshold';
  thresh001.help    = {''};
  thresh001.strtype = 'r';
  thresh001.num     = [1 1];
  thresh001.val     = {0.001};

  kthresh           = cfg_entry;
  kthresh.tag       = 'kthresh';
  kthresh.name      = 'Extent (voxels)';
  kthresh.help      = {'Enter the extent threshold in voxels'};
  kthresh.strtype   = 'r';
  kthresh.val       = {0};
  kthresh.num       = [1 1];

  noniso            = cfg_menu;
  noniso.name       = 'Correct for non-isotropic smoothness';
  noniso.tag        = 'noniso';
  noniso.labels     = {'Yes','No'};
  noniso.values     = {1,0};
  noniso.val        = {1};
  noniso.help       = {'Correct for non-isotropic smoothness for cluster extent thresholds.'};

  none              = cfg_const;
  none.tag          = 'none';
  none.name         = 'None';
  none.val          = {1};
  none.help         = {'No threshold'};

  k                 = cfg_branch;
  k.tag             = 'k';
  k.name            = 'k-value';
  k.val             = {kthresh, noniso };
  k.help            = {''};

  fwe               = cfg_branch;
  fwe.tag           = 'fwe';
  fwe.name          = 'FWE';
  fwe.val           = {thresh05 };
  fwe.help          = {''};

  fdr               = cfg_branch;
  fdr.tag           = 'fdr';
  fdr.name          = 'FDR';
  fdr.val           = {thresh05 };
  fdr.help          = {''};

  fwe2              = cfg_branch;
  fwe2.tag          = 'fwe2';
  fwe2.name         = 'FWE';
  fwe2.val          = {thresh05, noniso };
  fwe2.help         = {''};

  uncorr            = cfg_branch;
  uncorr.tag        = 'uncorr';
  uncorr.name       = 'uncorrected';
  uncorr.val        = {thresh001 };
  uncorr.help       = {''};

  kuncorr           = cfg_branch;
  kuncorr.tag       = 'kuncorr';
  kuncorr.name      = 'uncorrected';
  kuncorr.val       = {thresh05, noniso };
  kuncorr.help      = {''};

  En                = cfg_branch;
  En.tag            = 'En';
  En.name           = 'Expected voxels per cluster';
  En.val            = {noniso };
  En.help           = {''};

  inverse           = cfg_menu;
  inverse.name      = 'Show also inverse effects (e.g. neg. values)';
  inverse.tag       = 'inverse';
  inverse.labels    = {'Yes','No'};
  inverse.values    = {1,0};
  inverse.val       = {0};
  inverse.help      = {'Show also inverse effects (e.g. neg. values). This is not valid if you convert to (log) p-values.'};

  threshdesc        = cfg_choice;
  threshdesc.name   = 'Threshold type peak-level';
  threshdesc.tag    = 'threshdesc';
  threshdesc.values = {none uncorr fdr fwe};
  threshdesc.val    = {uncorr};
  threshdesc.help   = {'Select method for voxel threshold'};

  cluster           = cfg_choice;
  cluster.name      = 'Cluster extent threshold';
  cluster.tag       = 'cluster';
  cluster.values    = {none k En kuncorr fwe2};
  cluster.val       = {none};
  cluster.help      = {'Select method for extent threshold'};

  conversion        = cfg_branch;
  conversion.tag    = 'conversion';
  conversion.name   = 'Conversion';
  conversion.val    = {sel threshdesc inverse cluster};
  conversion.help   = {''};

  atlas             = cfg_menu;
  atlas.name        = 'Atlas Labeling';
  atlas.tag         = 'atlas';
  atlas.labels{1}   = 'None';
  atlas.values{1}   = 'None';
  list              = spm_atlas('List','installed');
  for i=1:numel(list)
    atlas.labels{i+1} = list(i).name;
    atlas.values{i+1} = list(i).name;
  end
  atlas.val         = {'None'};
  atlas.help        = {
    'Select atlas for labeling. The prepending atlas name ''dartel_'' indicates that this atlas was created using Dartel spatial normalization with the Dartel IXI template as default.'
    ''
    'Please note, that you can install additional atlases for CAT12 using the command ''cat_install_atlases''. '
  };



  % T2x volumes
  % -----------------------------------------------------------------------
  T2x      = cfg_exbranch;
  T2x.tag  = 'T2x';
  T2x.name = 'Threshold and transform spmT images';
  T2x.val  = {data_T2x,conversion,atlas};
  T2x.prog = @cat_stat_spm2x;
  T2x.vout = @vout_stat_spm2x;
  T2x.help = {
    'This function transforms t-maps to P, -log(P), r or d-maps.'
    'The following formulas are used:'
    '--------------------------------'
    'correlation coefficient:'
    '          t'
    'r = ------------------'
    '      sqrt(t^2 + df)'
    'p-value'
    'p = 1-spm_Tcdf'
    'log p-value'
    '-log10(1-P) = -log(1-spm_Tcdf)'
    'For the last case of log transformation this means that a p-value of p=0.99 (0.01) is transformed to a value of 2.'
    'Examples:'
    'p-value  -log10(1-P)'
    '0.1      1'
    '0.05     1.3'
    '0.01     2'
    '0.001    3'
    '0.0001   4'
    'All maps can be thresholded using height and extent thresholds and you can also apply corrections for multiple comparisons based on family-wise error (FWE) or false discovery rate (FDR). You can easily threshold and/or transform a large number of spmT-maps using the same thresholds.'
    'Naming convention of the transformed files:'
    '   Type_Contrast_Pheight_Pextent_K_Neg'
    '   Type:      P    - p-value'
    '              logP - log p-value'
    '              R    - correlation coefficient'
    '              T    - t-value'
    '   Contrast:  name used in the contrast manager with replaced none valid'
    '              strings'
    '   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")'
    '              pFWE - p-value with FWE correction in %'
    '              pFDR - p-value with FDR correction in %'
    '   Pextent:   pk    - uncorr. extent p-value in % (p<0.05 coded with "p5")'
    '              pkFWE - extent p-value with FWE correction in %'
    '   K:         extent threshold in voxels'
    '   Neg:       image also shows thresholded inverse effects (e.g. neg. '
    '              values) '
  }';


  % T2x surfaces
  % -----------------------------------------------------------------------
  
  % Do not use 3D atlases for surfaces
  data_T2x.filter = {'gifti'};
  
  % surfaces
  T2x_surf        = T2x;
  T2x_surf.val    = {data_T2x,conversion};
  T2x_surf.tag    = 'T2x_surf';
  T2x_surf.name   = 'Threshold and transform spmT surfaces';
  T2x_surf.vout   = @vout_stat_spm2x_surf;
  
  
  % F2x volumes
  % -----------------------------------------------------------------------

  data_F2x          = cfg_files;
  data_F2x.tag      = 'data_F2x';
  data_F2x.name     = 'Data';
  data_F2x.filter   = {'image'};
  data_F2x.ufilter  = '^spmF.*';
  data_F2x.num      = [1 Inf];
  data_F2x.help     = {'Select spmF-data to select.'};

  sel               = cfg_menu;
  sel.name          = 'Convert F value to';
  sel.tag           = 'sel';
  sel.labels        = {'p','-log(p)','coefficient of determination R^2','apply thresholds without conversion'};
  sel.values        = {1,2,3,4};
  sel.val           = {2};
  sel.help          = {'Select conversion of F-value'};

  none              = cfg_const;
  none.tag          = 'none';
  none.name         = 'None';
  none.val          = {1};
  none.help         = {'No threshold'};

  cluster           = cfg_choice;
  cluster.name      = 'Cluster extent threshold';
  cluster.tag       = 'cluster';
  cluster.values    = {none k En kuncorr fwe2};
  cluster.val       = {none};
  cluster.help      = {'Select method for extent threshold'};

  conversion        = cfg_branch;
  conversion.tag    = 'conversion';
  conversion.name   = 'Conversion';
  conversion.val    = {sel threshdesc cluster};
  conversion.help   = {''};

  F2x               = cfg_exbranch;
  F2x.tag           = 'F2x';
  F2x.name          = 'Threshold and transform spmF images';
  F2x.val           = {data_F2x,conversion,atlas};
  F2x.prog          = @cat_stat_spm2x;
  F2x.vout          = @vout_stat_spm2x;
  F2x.help          = {
    'This function transforms F-maps to P, -log(P), or R2-maps.'
    'The following formulas are used:'
    '--------------------------------'
    'coefficient of determination R2:'
    '             1'
    'R2 = ------------------'
    '      1 + F*(p-1)/n-p)'
    'p-value:'
    'p = 1-spm_Fcdf'
    'log p-value:'
    '-log10(1-P) = -log(1-spm_Fcdf)'
    'For the last case of log transformation this means that a p-value of p=0.99 (0.01) is transformed to a value of 2.'
    'Examples:'
    'p-value  -log10(1-P)'
    '0.1      1'
    '0.05     1.3'
    '0.01     2'
    '0.001    3'
    '0.0001   4'
    'All maps can be thresholded using height and extent thresholds and you can also apply corrections for multiple comparisons based on family-wise error (FWE) or false discovery rate (FDR). You can easily threshold and/or transform a large number of spmT-maps using the same thresholds.'
    'Naming convention of the transformed files:'
    '   Type_Contrast_Pheight_K'
    '   Type:      P    - p-value'
    '              logP - log p-value'
    '              R2   - coefficient of determination'
    '   Contrast:  name used in the contrast manager with replaced none valid'
    '              strings'
    '   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")'
    '              pFWE - p-value with FWE correction in %'
    '              pFDR - p-value with FDR correction in %'
    '   K:         extent threshold in voxels'
  }';


  % F2x surfaces
  % -----------------------------------------------------------------------

  % Do not use 3D atlases for surfaces
  data_F2x.filter  = {'gifti'};
  
  F2x_surf         = F2x;
  F2x_surf.val     = {data_F2x,conversion};
  F2x_surf.tag     = 'F2x_surf';
  F2x_surf.name    = 'Threshold and transform spmF surfaces';
  F2x_surf.vout    = @vout_stat_spm2x_surf;

return

%_______________________________________________________________________
function showslice = cat_stat_showslice_all_GUI(data_vol)
  data_vol.help = {'Select all images. Images have to be in the same orientation with same voxel size and dimension (e.g. normalized images)'};

  scale           = cfg_menu;
  scale.tag       = 'scale';
  scale.name      = 'Proportional scaling?';
  scale.labels    = {'No','Yes'};
  scale.values    = {0 1};
  scale.val       = {0};
  scale.help      = {'This option should be only used if image intensity is not scaled (e.g. T1 images) or if images have to be scaled during statistical analysis (e.g. modulated images).'};

  orient          = cfg_menu;
  orient.tag      = 'orient';
  orient.name     = 'Spatial orientation';
  orient.labels   = {'axial','coronal','sagittal'};
  orient.values   = {3 2 1};
  orient.val      = {3};
  orient.help     = {'Spatial orientation of slice.'};

  slice           = cfg_entry;
  slice.tag       = 'slice';
  slice.name      = 'Selected slice (in mm)?';
  slice.strtype   = 'r';
  slice.num       = [1 1];
  slice.val       = {0};
  slice.help      = {'Choose slice in mm.'};

  showslice       = cfg_exbranch;
  showslice.tag   = 'showslice';
  showslice.name  = 'Display one slice for all images';
  showslice.val   = {data_vol,scale,orient,slice};
  showslice.prog  = @cat_stat_showslice_all;
  showslice.help  = {'This function displays a selected slice for all images and indicates the respective filenames which is useful to check image quality for a large number of files in a circumscribed region (slice).'};

%_______________________________________________________________________
function [check_cov, check_cov2] = cat_stat_check_cov_GUI(data_xml,outdir,fname,save,expert) 
 
  % --- update input data ---
  data_xml.name     = 'Quality measures (optional)';
  data_xml.help     = {'Select optional the quality measures that are saved during segmentation as xml-files in the report folder. This additionally allows to analyze image quality parameters such as noise, and bias. Please note, that the order of the xml-files should be the same as the other data files.'};
  
  % --- further data ---
  c                 = cfg_entry;
  c.tag             = 'c';
  c.name            = 'Vector';
  c.help            = {'Vector of nuisance values'};
  c.strtype         = 'r';
  c.num             = [Inf 1];

  nuisance          = cfg_repeat;
  nuisance.tag      = 'nuisance';
  nuisance.name     = 'Nuisance variable';
  nuisance.values   = {c};
  nuisance.num      = [0 Inf];
  nuisance.help     = {'This option allows for the specification of nuisance effects to be removed from the data. A potential nuisance parameter can be TIV if you check segmented data with the default modulation. In this case the variance explained by TIV will be removed prior to the calculation of the correlation. Another meaningful nuisance effect is age.'};

  gap               = cfg_entry;
  gap.tag           = 'gap';
  gap.name          = 'Separation';
  gap.strtype       = 'n';
  gap.num           = [1 1];
  gap.val           = {3};
  gap.help          = {'To speed up calculations you can define that covariance is estimated only every x voxel. Smaller values give slightly more accurate covariance, but will be much slower.'};

  data_vol          = cfg_files;
  data_vol.name     = 'Sample data';
  data_vol.tag      = 'data_vol';
  data_vol.filter   = 'image';
  data_vol.num      = [1 Inf];
  data_vol.help     = {'These are the (spatially registered) data. They must all have the same image dimensions, orientation, voxel size etc. Furthermore, it is recommended to use unsmoothed files.'};

  sample            = cfg_repeat;
  sample.tag        = 'sample';
  sample.name       = 'Data';
  sample.values     = {data_vol};
  sample.num        = [1 Inf];
  sample.help       = {'Specify data for each sample. If you specify different samples the mean correlation is displayed in separate boxplots (or violin plots) for each sample.'};


  check_cov         = cfg_exbranch;
  check_cov.tag     = 'check_cov';
  check_cov.name    = 'Check sample homogeneity of 3D data';
  if expert>1
    check_cov.val     = {sample,data_xml,gap,nuisance, outdir,fname,save};
  else
    check_cov.val     = {sample,data_xml,gap,nuisance};
  end
  check_cov.prog    = @cat_stat_check_cov;
  check_cov.help    = {
    'In order to identify images with poor image quality or even artefacts you can use this function. Images have to be in the same orientation with same voxel size and dimension (e.g. normalized images without smoothing). The idea of this tool is to check the correlation of all files across the sample.'
    ''
    'The correlation is calculated between all images and the mean for each image is plotted using a spmmat and the indicated filenames. The smaller the mean correlation the more deviant is this image from the sample mean. In the plot outliers from the sample are usually isolated from the majority of images which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the image order.'
  };


  % --- main ---
  check_cov2        = check_cov; 
  check_cov2.tag    = 'check_cov2';
  check_cov2.name   = 'Check sample homogeneity of 3D data (new exp. version)';
  check_cov2.val    = {sample,gap,nuisance};
  check_cov2.prog   = @cat_stat_check_cov2;

%_______________________________________________________________________
function check_SPM = cat_stat_check_SPM_GUI(outdir,fname,save,expert) 

  spmmat                      = cfg_files;
  spmmat.tag                  = 'spmmat';
  spmmat.name                 = 'Select SPM.mat';
  spmmat.filter               = {'mat'};
  spmmat.ufilter              = '^SPM\.mat$';
  spmmat.num                  = [1 1];
  spmmat.help                 = {'Select the SPM.mat file that contains the design specification.'};

  % check_SPM_cov
  use_unsmoothed_data         = cfg_menu;
  use_unsmoothed_data.name    = 'Use unsmoothed data if found';
  use_unsmoothed_data.tag     = 'use_unsmoothed_data';
  use_unsmoothed_data.labels  = {'Yes','No'};
  use_unsmoothed_data.values  = {1,0};
  use_unsmoothed_data.val     = {1};
  use_unsmoothed_data.help    = {'Check for sample homogeneity results in more reliable values if unsmoothed data are used. Unsmoothed data contain more detailed information about differences and similarities between the data.'};

  adjust_data                 = cfg_menu;
  adjust_data.name            = 'Adjust data using design matrix';
  adjust_data.tag             = 'adjust_data';
  adjust_data.labels          = {'Yes','No'};
  adjust_data.values          = {1,0};
  adjust_data.val             = {1};
  adjust_data.help            = {'This option allows to use nuisance and group parameters from the design matrix to obtain adjusted data. In this case the variance explained by these parameters will be removed prior to the calculation of the correlation. Furthermore, global scaling (if defined) is also applied to the data.'};

  do_check_cov                = cfg_branch;
  do_check_cov.tag            = 'do_check_cov';
  do_check_cov.name           = 'Yes';
  if expert>1
    do_check_cov.val            = {use_unsmoothed_data adjust_data ,outdir,fname,save};
  else
    do_check_cov.val            = {use_unsmoothed_data adjust_data};
  end  
  do_check_cov.help           = {''};

  none                        = cfg_const;
  none.tag                    = 'none';
  none.name                   = 'No';
  none.val                    = {1};
  none.help                   = {''};

  check_SPM_cov               = cfg_choice;
  check_SPM_cov.name          = 'Check for sample homogeneity';
  check_SPM_cov.tag           = 'check_SPM_cov';
  check_SPM_cov.values        = {none do_check_cov};
  check_SPM_cov.val           = {do_check_cov};
  check_SPM_cov.help          = {
    'In order to identify images with poor image quality or even artefacts you can use this function. The idea of this tool is to check the correlation of all files across the sample using the files that are already defined in SPM.mat.'
    ''
    'The correlation is calculated between all images and the mean for each image is plotted using a boxplot (or violin plot) and the indicated filenames. The smaller the mean correlation the more deviant is this image from the sample mean. In the plot outliers from the sample are usually isolated from the majority of images which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the image order'
  };

  check_SPM_ortho             = cfg_menu;
  check_SPM_ortho.name        = 'Check for design orthogonality';
  check_SPM_ortho.tag         = 'check_SPM_ortho';
  check_SPM_ortho.labels      = {'Yes','No'};
  check_SPM_ortho.values      = {1,0};
  check_SPM_ortho.val         = {1};
  check_SPM_ortho.help        = {'Review Design Orthogonality.'};

  check_SPM                   = cfg_exbranch;
  check_SPM.tag               = 'check_SPM';
  check_SPM.name              = 'Check design orthogonality and homogeneity';
  check_SPM.val               = {spmmat,check_SPM_cov,check_SPM_ortho};
  check_SPM.prog              = @cat_stat_check_SPM;
  check_SPM.help              = {'Use design matrix saved in SPM.mat to check for sample homogeneity of the used data and for orthogonality of parameters.'};

%_______________________________________________________________________
function calcvol = cat_stat_TIV_GUI
  calcvol_name         = cfg_entry;
  calcvol_name.tag     = 'calcvol_name';
  calcvol_name.name    = 'Output file';
  calcvol_name.strtype = 's';
  calcvol_name.num     = [1 Inf];
  calcvol_name.val     = {'TIV.txt'};
  calcvol_name.help    = {
  'The output file is written to current working directory unless a valid full pathname is given.'};

  calcvol_TIV         = cfg_menu;
  calcvol_TIV.tag     = 'calcvol_TIV';
  calcvol_TIV.name    = 'Save values';
  calcvol_TIV.labels  = {'TIV only' 'TIV & GM/WM/CSF/WMH'};
  calcvol_TIV.values  = {1 0};
  calcvol_TIV.val     = {1};
  calcvol_TIV.help    = {'You can save either the total intracranial volume (TIV) only or additionally also save the global volumes for GM, WM, CSF, and WM hyperintensities.'
  ''
  };

  clear data_xml
  data_xml = cfg_files;
  data_xml.name = 'XML files';
  data_xml.tag  = 'data_xml';
  data_xml.filter = 'xml';
  data_xml.ufilter = '^cat_.*';
  data_xml.num  = [1 Inf];
  data_xml.help   = {...
  'Select xml-files that are saved during segmentation in the report folder.'};

  calcvol       = cfg_exbranch;
  calcvol.tag   = 'calcvol';
  calcvol.name  = 'Estimate TIV and global tissue volumes';
  calcvol.val   = {data_xml,calcvol_TIV,calcvol_name};
  calcvol.prog  = @cat_stat_TIV;
  calcvol.vout  = @vout_stat_TIV;
  calcvol.help  = {
  'This function reads raw volumes for TIV/GM/WM/CSF/WM hyperintensities (WMH) and saves values in a txt-file. These values can be read with the matlab command: vol = spm_load. If you choode to save all values the entries for TIV/GM/WM/CSF/WMH are now saved in vol(:,1) vol(:,2) vol(:,3), vol(:,4), and vol(:,5) respectively.'
  ''
  'You can use TIV either as nuisance in an AnCova model or as user-specified globals with the "global calculation" option depending on your hypothesis. The use of TIV as nuisance or globals is recommended for modulated data where both the affine transformation and the non-linear warping of the registration are corrected for. '
  ''
  };

%_______________________________________________________________________
function calcroi = cat_roi_fun_GUI(outdir)   
  roi_xml               = cfg_files;
  roi_xml.name          = 'XML files';
  roi_xml.tag           = 'roi_xml';
  roi_xml.filter        = 'xml';
  roi_xml.ufilter       = '^catROI.*';
  roi_xml.num           = [1 Inf];
  roi_xml.help          = {'These are the xml-files that are saved in the label folder.'};

  % NOT USED
  %{
  usefolder             = cfg_menu;
  usefolder.tag         = 'folder';
  usefolder.name        = 'Use foldername';
  usefolder.labels      = {'Yes' 'No'};
  usefolder.values      = {1 0};
  usefolder.val         = {0};
  usefolder.help        = {'Use foldername to describe the subject.'};
  %}
  
  point                 = cfg_menu;
  point.tag             = 'point';
  point.name            = 'Decimal point';
  point.labels          = {',','.'};
  point.values          = {',','.'};
  point.val             = {'.'};
  point.help            = {'Decimal point.'}; % that has to be unequal to the column delimiter.'};

  % tab "\t" does not work and so we automatically switch in case of decimal 
  % point "," to delimiter ";".
  %{
  delimiter             = cfg_menu;
  delimiter.tag         = 'delimiter';
  delimiter.name        = 'column delimiter';
  delimiter.labels      = {',',';',' '};
  delimiter.values      = {',',';',' '};
  delimiter.val         = {','};
  delimiter.help        = {'Delimiter between columns.'};
  %}

  calcroi_name          = cfg_entry;
  calcroi_name.tag      = 'calcroi_name';
  calcroi_name.name     = 'Output file';
  calcroi_name.strtype  = 's';
  calcroi_name.num      = [1 Inf];
  calcroi_name.val      = {'ROI'};
  calcroi_name.help     = {'The volume values in mL are written to the current working directory unless a valid full pathname is given. The output file will also include the name of the atlas and the measure (e.g. Vgm). The file is using tabstops to separate values in order to easily import the file into Excel or SPSS or any other software for subsequent analysis.'};

  calcroi               = cfg_exbranch;
  calcroi.tag           = 'calcroi';
  calcroi.name          = 'Estimate volumes inside ROI';
  calcroi.val           = {roi_xml,point,outdir,calcroi_name}; 
  %calcroi.val   = {roi_xml,usefolder,point,outdir,calcroi_name}; % usefolder is never used
  calcroi.prog          = @(job)cat_roi_fun('exportSample',job);
  calcroi.help          = {
    'This function reads volume values in mL inside a ROI from different atlases and saves values for all data in a csv-file. '
    'Missed values were replaced by NaN.'
  };

%_______________________________________________________________________
function urqio = cat_vol_urqio_GUI
%  ------------------------------------------------------------------------
%  Ultra-High Resolution Quantitative Image Optimization
%  ------------------------------------------------------------------------

  % -- Data ---
  r1          = cfg_files; 
  r1.tag      = 'r1';
  r1.name     = 'R1-Volumes';
  r1.filter   = 'image';
  r1.ufilter  = '.*';
  r1.num      = [1 Inf];
  r1.help     = {'Select R1 weighted images.'};

  pd          = cfg_files; 
  pd.tag      = 'pd';
  pd.name     = 'PD-Volumes';
  pd.filter   = 'image';
  pd.ufilter  = '.*';
  pd.num      = [1 Inf];
  pd.help     = {'Select PD weighted images.'};

  r2s         = cfg_files; 
  r2s.tag     = 'r2s';
  r2s.name    = 'R2s-Volumes';
  r2s.filter  = 'image';
  r2s.ufilter = '.*';
  r2s.num     = [1 Inf];
  r2s.help    = {'Select R2s weighted images.'};

  data        = cfg_branch;
  data.tag    = 'data';
  data.name   = 'Input data';
  data.val    = {r1 pd r2s}; 
  data.help   = {
    'Input Images.'
  };

  
  % --- Parameter ---
  spm         = cfg_menu;
  spm.tag     = 'spm';
  spm.name    = 'Use SPM Preprocessing';
  spm.labels  = {'No','Yes'};
  spm.values  = {0 1};
  spm.val     = {1};
  spm.help    = {
    'Use SPM preprocessing if the data is not skull-stripped.'
  };
  
  bc          = cfg_menu;
  bc.tag      = 'bc';
  bc.name     = 'Bias Correction';
  bc.labels   = {'No','light','medium','strong'};
  bc.values   = {0 0.5 1 2};
  bc.val      = {1};
  bc.help     = {
    'Additional bias correction that is important for detection and correction of blood vessels.'
    ''
    'The correction uses a simple tissue classification and local filter approaches to estimate the local signal intensity in the WM and GM segment, e.g. a minimum/maximum filter in the WM for PD/T1 images.  Next, unclassified voxels were approximated and smoothed depending on the defined strength.  '
    ''
  };

  in          = cfg_menu;
  in.tag      = 'in';
  in.name     = 'Intensity Normalization';
  in.labels   = {'No','Yes'};
  in.values   = {0 1};
  in.val      = {1};
  in.help     = {
    'Additional global intensity normalization that is also important for detection and correction of blood vessels.'
    ''
  };

  bvc         = cfg_menu;
  bvc.tag     = 'bvc';
  bvc.name    = 'Blood Vessel Correction';
  bvc.labels  = {'No','Yes'};
  bvc.values  = {0 1};
  bvc.val     = {1};
  bvc.help    = {
    'Correction of blood vessels with high intensity in T1/R1/R2s and low intensity in PD images by CSF-like intensities. '
    ''
  };

  ss          = cfg_menu;
  ss.tag      = 'ss';
  ss.name     = 'Apply Skull-Stripping';
  ss.labels   = {'No','Yes'};
  ss.values   = {0 1};
  ss.val      = {1};
  ss.help     = {
    'Write skull-stripped images. '
    ''
  };

  nc          = cfg_menu;
  nc.tag      = 'nc';
  nc.name     = 'Noise Correction';
  nc.labels   = {'No','Yes'};
  nc.values   = {0 1};
  nc.val      = {1};
  nc.help     = {
    'Noise corrections of the final images.'
    ''
  };

  prefix         = cfg_entry;
  prefix.tag     = 'prefix';
  prefix.name    = 'Filename prefix';
  prefix.strtype = 's';
  prefix.num     = [0 Inf];
  prefix.val     = {'catsyn_'};
  prefix.help    = {
    'Prefix of output files.'};


  opts        = cfg_branch;
  opts.tag    = 'opts';
  opts.name   = 'Parameter';
  opts.val    = {spm bc in bvc ss nc prefix}; 
  opts.help   = {
    'Parameter settings for image correction.'
  };


  % --- Output ---
  pdo         = cfg_menu;
  pdo.tag     = 'pd';
  pdo.name    = 'PD Output';
  pdo.labels  = {'No','Yes'};
  pdo.values  = {0 1};
  pdo.val     = {1}; 
  pdo.help    = {
    'Write PD output images.'
  };

  t1o         = cfg_menu;
  t1o.tag     = 't1';
  t1o.name    = 'T1 Output';
  t1o.labels  = {'No','Yes'};
  t1o.values  = {0 1};
  t1o.val     = {1}; 
  t1o.help    = {
    'Write synthesized T1 output images based on the PD image.'
  };

  r1o         = cfg_menu;
  r1o.tag     = 'r1';
  r1o.name    = 'R1 Output';
  r1o.labels  = {'No','Yes'};
  r1o.values  = {0 1};
  r1o.val     = {1}; 
  r1o.help    = {
    'Write R1 output images.'
  };

  r2so        = cfg_menu;
  r2so.tag    = 'r2s';
  r2so.name   = 'R2s Output';
  r2so.labels = {'No','Yes'};
  r2so.values = {0 1};
  r2so.val    = {1}; 
  r2so.help   = {
    'Write R2s output images.'
  };

  bvco        = cfg_menu;
  bvco.tag    = 'bv';
  bvco.name   = 'Blood Vessel Output';
  bvco.labels = {'No','Yes'};
  bvco.values = {0 1};
  bvco.val    = {0}; 
  bvco.help   = {
    'Write map of blood vessels.'
  };
    
  output      = cfg_branch;
  output.tag  = 'output';
  output.name = 'Output';
  output.val  = {r1o r2so pdo t1o bvco}; 
  output.help = {
    'Output images.'
  };

  
  % --- main ---
  % batch mode - output is undefined!
  urqio       = cfg_exbranch;
  urqio.tag   = 'urqio';
  urqio.name  = 'Ultra-High Resolution Quantitative Image Optimization';
  urqio.val   = {data opts output};
  urqio.prog  = @cat_vol_urqio;
  %urqio.vout  = @vfiles_urqio;
  urqio.help  = {
    'Additional correction of high resolution PD, R1, and R2s weighted images that includes another bias correction, intensity normalization, and blood vessel correction step. '
    ''
    'WARNING: This tool is in development and was just tested on a small set of subjects!'
  };

%_______________________________________________________________________
function avg_img = cat_vol_average_GUI(data,outdir)
% image average
% -------------------------------------------------------------------------

  % update input functions
  data.name       = 'Select images';
  data.help       = {'Select images for calculating average.'};

  outdir.help     = {'Select a directory where files are written otherwise the path of the first image will be used.'};

  % filename
  output          = cfg_entry;
  output.tag      = 'output';
  output.name     = 'Output Filename';
  output.help     = {
    'The output image is written to current working directory unless a valid full pathname is given. If a path name is given here, the output directory setting will be ignored.'
    'If the field is left empty, i.e. set to '''', then the name of the 1st input image, preprended with ''i'', is used (change this letter in the spm_defaults if necessary).'
  };
  output.strtype  = 's';
  output.num      = [0 Inf];
  output.val      = {'avg.nii'};

  % main
  avg_img         = cfg_exbranch;
  avg_img.tag     = 'avg_img';
  avg_img.name    = 'Image Average';
  avg_img.val     = {data output outdir};
  avg_img.help    = {'This function is for calculating the average of a set of images, which should be of same dimension and voxel size (i.e. after spatial registration).'};
  avg_img.prog    = @cat_vol_avg;
  avg_img.vout    = @vout_avg;

%_______________________________________________________________________
function vf = vfiles_defs(job)

PU = job.field1;
PI = job.images;

vf = cell(numel(PI),1);
for i=1:numel(PU),
    for m=1:numel(PI),
        [pth,nam,ext,num] = spm_fileparts(PI{m});
        
        switch job.modulate
        case 2
            filename = fullfile(pth,['m0w' nam ext num]);
        case 1
            filename = fullfile(pth,['mw' nam ext num]);
        case 0
            filename = fullfile(pth,['w' nam ext num]);
        end;
        vf{m} = filename;
    end
end

return;
%_______________________________________________________________________

function vf = vfiles_defs2(job)

  PU = job.field;
  PI = job.images;

  vf = cell(numel(PU),numel(PI));
  for i=1:numel(PU),
      for m=1:numel(PI),
          [pth,nam,ext,num] = spm_fileparts(PI{m}{i});

          switch job.modulate
          case 2
              filename = fullfile(pth,['m0w' nam ext num]);
          case 1
              filename = fullfile(pth,['mw' nam ext num]);
          case 0
              filename = fullfile(pth,['w' nam ext num]);
          end;
          vf{i,m} = filename;
      end
  end

return;
%_______________________________________________________________________
function cdep = vfiles_urqio(job)
  %%
  cdep = cfg_dep;
  if job.output.r1
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'R1 Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'r1_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if job.output.pd
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'PD Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'pd_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if job.output.t1
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'T1 Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 't1_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if job.output.r2s==1 || job.output.r2s==3
    cdep(end+1)           = cfg_dep;
    cdep(end).sname      = 'R2s nobc Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'nobc_r2s_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end 
  if job.output.r2s==2 || job.output.r2s==3
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'R2s bc Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'bc_r2s_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end 
  if job.output.bv
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Blood Vessel Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'bv_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if numel(cdep)>1
    cdep(1)=[];
  end
%%
return;
%_______________________________________________________________________
function dep = vfiles_sanlm(varargin)
  %job.returnOnlyFilename = 1; 
  %vf = cat_vol_sanlm(job); 
  
  dep(1)            = cfg_dep;
  dep(1).sname      = 'SANLM Images';
  dep(1).src_output = substruct('.','files');
  dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return;
%_______________________________________________________________________
function dep = vfiles_maskimg(varargin)
  %job.returnOnlyFilename = 1; 
  %vf = cat_vol_maskimage(job); 
  
  dep            = cfg_dep;
  dep.sname      = 'Masked Images';
  dep.src_output = substruct('.','files');
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return;
%_______________________________________________________________________
function cdep = vfiles_headtrimming(job)
  job.returnOnlyFilename = 1; 
  %vf = cat_vol_headtrimming(job); 
  vf = job; 
  cdep = cfg_dep;

  if isfield(vf.image_selector,'manysubjects')
    cdep(end).sname      = 'source images';
    cdep(end).src_output = substruct('.','image_selector','.','manysubjects','.','simages'); 
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    if isfield(vf.image_selector.manysubjects,'oimages')
      for i=1:numel(vf.image_selector.manysubjects.oimages)
        cdep(end+1)          = cfg_dep;%#ok<AGROW>
        cdep(end).sname      = sprintf('other images %d',i);
        cdep(end).src_output = substruct('.','image_selector','.','manysubjects','.','oimages','{}',{i}); 
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      end
    end
  elseif isfield(vf.image_selector,'subjectimages') 
    % image-wise
    % - first image
    cdep(end).sname      = sprintf('first images of all subjects');
    cdep(end).src_output = substruct('.','image_selector','.','firstimages'); 
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    % - other images
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = sprintf('other images of all subjects');
    cdep(end).src_output = substruct('.','image_selector','.','otherimages'); 
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    % subject-wise ... 
    % the substruct seams to be correct, but it does not work, 
    % probably because each cdep entry has to be unique  
    % RD201810
    %{
    for si=1:numel(vf.image_selector.subjectimages)
      cdep(end+1)          = cfg_dep;%#ok<AGROW>
      cdep(end).sname      = sprintf('all imgages of subject %d',si);
      cdep(end).src_output = substruct('.','image_selector','.','subjectimages','{}',{si}); 
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    % single scans
    for si=1:numel(vf.image_selector.subjectimages)
      for fi=1:numel(vf.image_selector.subjectimages{si})
        cdep(end+1)          = cfg_dep;%#ok<AGROW>
        cdep(end).sname      = sprintf('subject %d image %d',si,fi);
        cdep(end).src_output = substruct('.','image_selector','.','subjectimages','{}',{si},'()',{fi}); 
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      end
    end
    %}
  else
    for i=1:numel(vf.images)
      cdep(i)            = cfg_dep;
      cdep(i).sname      = sprintf('image %d',i);
      cdep(i).src_output = substruct('.','images','{}',{i}); 
      cdep(i).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
  end

return;
%_______________________________________________________________________
function dep = vfiles_volctype(varargin)
 % job.returnOnlyFilename = 1; 
 % vf = cat_io_volctype(job);
    
  dep            = cfg_dep;
  dep.sname      = ['Images'];
  dep.src_output = substruct('.','files');
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return;
%_______________________________________________________________________
function vf = vfiles_qa(job)
  s  = cellstr(char(job.data)); vf = s; 
  for i=1:numel(s),
      [pth,nam,ext,num] = spm_fileparts(s{i});
      vf{i} = fullfile(pth,[job.prefix,nam,ext,num]);
  end;
return;

%------------------------------------------------------------------------
function dep = vout_avg(job)
  dep            = cfg_dep;
  if ~ischar(job.output) || strcmp(job.output, '<UNDEFINED>')
      dep.sname  = 'Average Image';
  else
      dep.sname  = sprintf('Average Image %s', job.output);
  end
  dep.src_output = substruct('.','files');
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return

%------------------------------------------------------------------------
function dep = vout_stat_TIV(varargin)
  dep            = cfg_dep;
  dep.sname      = 'TIV';
  dep.src_output = substruct('.','calcvol','()',{':'});
  dep.tgt_spec   = cfg_findspec({{'strtype','e','strtype','r'}});
return

%------------------------------------------------------------------------
function cdep = vout_reslice(job)
  ind  = 1;
  if job.write_avg
      cdep(ind)            = cfg_dep;
      cdep(ind).sname      = 'Midpoint Average';
      cdep(ind).src_output = substruct('.','avg','()',{':'});
      cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      ind = ind + 1;
  end
  if job.write_rimg
      cdep(ind)            = cfg_dep;
      cdep(ind).sname      = 'Realigned images';
      cdep(ind).src_output = substruct('.','rimg','()',{':'});
      cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      ind = ind + 1;
  end
  if isfield(job.reg,'nonlin') && job.reg.nonlin.write_jac
      cdep(ind)            = cfg_dep;
      cdep(ind).sname      = 'Jacobian Diff';
      cdep(ind).src_output = substruct('.','jac','()',{':'});
      cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      ind = ind + 1;
  end
  if isfield(job.reg,'nonlin') && job.reg.nonlin.write_def
      cdep(ind)            = cfg_dep;
      cdep(ind).sname      = 'Deformation (1)';
      cdep(ind).src_output = substruct('.','def1','()',{':'});
      cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      ind = ind + 1;

      cdep(ind)            = cfg_dep;
      cdep(ind).sname      = 'Deformation (2)';
      cdep(ind).src_output = substruct('.','def2','()',{':'});
      cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end

return
 
%------------------------------------------------------------------------
function dep = vout_stat_spm2x(job)
  dep            = cfg_dep;
  dep.sname      = 'Transform & Threshold spm volumes';
  dep.src_output = substruct('.','Pname');
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return

%------------------------------------------------------------------------
function dep = vout_stat_spm2x_surf(job)
  dep            = cfg_dep;
  dep.sname      = 'Transform & Threshold spm surfaces';
  dep.src_output = substruct('.','Pname');
  dep.tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
return