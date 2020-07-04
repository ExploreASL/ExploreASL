function [output,output_spm,output1173,output1445] = cat_conf_output(expert)
%function [output,output_spm,output1173,output1445,output1585] = cat_conf_output(expert)
% writing options for data
%_______________________________________________________________________
%
% Christian Gaser
% $Id: cat_conf_output.m 1605 2020-04-14 18:44:36Z gaser $
%
%#ok<*AGROW>
 
  try
    expert = cat_get_defaults('extopts.expertgui');
  catch %#ok<CTCH>
    expert = 0; 
  end

  %------------------------------------------------------------------------

  surface        = cfg_menu;
  surface.tag    = 'surface';
  surface.name   = 'Surface and thickness estimation';
  surface.labels = {'No','Yes'};
  surface.values = {0 1};
  surface.def    = @(val)cat_get_defaults('output.surface', val{:});
  surface.help   = {
  'Use projection-based thickness (PBT) (Dahnke et al. 2012) to estimate cortical thickness and to create the central cortical surface for left and right hemisphere. Surface reconstruction includes topology correction (Yotter et al. 2011), spherical inflation (Yotter et al.) and spherical registration. Additionally you can also estimate surface parameters such as gyrification, cortical complexity or sulcal depth that can be subsequently analyzed at each vertex of the surface. '
  ''
  'Please note, that surface reconstruction additionally requires about 20-60 min of computation time.'
  ''
  'You can also estimate thickness for ROI analysis only. This takes much less time, but does not allow to use the advantages of surface-based registration and smoothing and the extraction of additional surface parameters. Here, the analysis is restricted to cortical thickness in atlas-defined ROIs only.'
  ''
  };

  if expert == 2
    surface.labels = {'No','lh + rh','lh + rh + cb',...
      'lh + rh (fast, no registration)','lh + rh + cb (fast, no registration)', ...
      ...'lh + rh (fast registration)', ... 7
      'lh + rh + cb (fast registration)',... 8 
      'Thickness estimation (for ROI analysis only)', 'Full'};
    surface.values = {0 1 2 , 5 6 , 8 , 9 12}; % 7 8 
    surface.help   = [surface.help; {
      'Cerebellar reconstruction is still in development and is strongly limited due to the high frequency of folding and image properties! '
      ''
      'The fast reconstruction allows a _visual_ check of the processing quality by providing a reconstructed cortical surface and thickness values with slightly lower accuracy. It takes only a few minutes for both hemispheres by internally using 0.8 mm (instead of 0.5 mm) spatial resolution and by skipping topology correction and surface registration. Please note that neither the surface nor the measures on the surface (e.g., thickness) can be used for any further statistical analysis! '
      ''
    }];    
  end

  if expert 
    % write specific output surface maps
    % 0 none:       only surfaces (central,white,pial,sphere,sphere.reg)
    % 1 default:    + thickness
    % 2 expert:     + L4myelination, topology defects, 
    % 3 developer:  + WM and CSF thickness, YppRMSEmap?
    % 4 debug:      + substeps in subdirs 
    surf_measures        = cfg_menu;
    surf_measures.tag    = 'surf_measures';
    surf_measures.name   = 'Surface measures';
    surf_measures.labels = {'Default','Expert'};
    surf_measures.values = {1 2};
    surf_measures.val    = {1};
    surf_measures.help   = {
     ['Write additional surface measures that are currently under development. ' ...
      'The defaults setting include only cortical thickness, whereas the expert level also ' ...
      'writes a myelination map (normalized T1 intensity extracted at the layer 4 surface) and ' ...
      'a map of topology defects coding the percentage size of the effect. '];
    };
    if expert == 2
      surf_measures.labels = [ surf_measures.labels {'Developer','Debug'}];
      surf_measures.values = [ surf_measures.values {3,4}];
      surf_measures.val    = {3};
      surf_measures.help   = [ surf_measures.help 
        {'The developer option further write the gyral and sulcal thickness. '}
        ];
    end
  end
  
  %------------------------------------------------------------------------
  native        = cfg_menu;
  native.tag    = 'native';
  native.name   = 'Native space';
  native.labels = {'No','Yes'};
  native.values = {0 1};
  native.help   = {
    'The native space option allows you to save a tissue class image (p*) that is in alignment with the original image.'
  ''
  };

  warped        = cfg_menu;
  warped.tag    = 'warped';
  warped.name   = 'Normalized';
  warped.labels = {'No','Yes'};
  warped.values = {0 1};
  warped.help   = {'Write image in normalized space without any modulation.'
  ''
  };

  dartel        = cfg_menu;
  dartel.tag    = 'dartel';
  dartel.name   = 'DARTEL export';
  if expert
    dartel.labels = {'No','Rigid (SPM12 default)','Affine','Both'};
    dartel.values = {0 1 2 3};
  else
    dartel.labels = {'No','Rigid (SPM12 default)','Affine'};
    dartel.values = {0 1 2};
  end
  dartel.help   = {
  'This option is to export data into a form that can be used with DARTEL. The SPM default is to only apply rigid body transformation. However, a more appropriate option is to apply affine transformation, because the additional scaling of the images requires less deformations to non-linearly register brains to the template.'
  ''
  };

	native.def  = @(val)cat_get_defaults('output.bias.native', val{:});
	warped.def  = @(val)cat_get_defaults('output.bias.warped', val{:});
	dartel.def  = @(val)cat_get_defaults('output.bias.dartel', val{:});
	bias        = cfg_branch;
	bias.tag    = 'bias';
	bias.name   = 'Bias, noise and global intensity corrected T1 image';
	if expert
		bias.val    = {native warped dartel};
	else
		bias.val    = {warped};
	end
	bias.help   = {
		'This is the option to save a bias, noise, and global intensity corrected version of the original T1 image. MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images. The bias corrected version should have more uniform intensities within the different types of tissues and can be saved in native space and/or normalised. Noise is corrected by an adaptive non-local mean (NLM) filter (Manjon 2008, Medical Image Analysis 12).'
	''
	};

  native.def  = @(val)cat_get_defaults('output.las.native', val{:});
  warped.def  = @(val)cat_get_defaults('output.las.warped', val{:});
  dartel.def  = @(val)cat_get_defaults('output.las.dartel', val{:});
  las        = cfg_branch;
  las.tag    = 'las';
  las.name   = 'Bias, noise and local intensity corrected T1 image';
  las.val    = {native warped dartel};
  las.help   = {
    'This is the option to save a bias, noise, and local intensity corrected version of the original T1 image. MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images. The bias corrected version should have more uniform intensities within the different types of tissues and can be saved in native space and/or normalised. Noise is corrected by an adaptive non-local mean (NLM) filter (Manjon 2008, Medical Image Analysis 12).'
  ''
  };


  %------------------------------------------------------------------------

  jacobianwarped      = warped;
  jacobianwarped.tag  = 'jacobianwarped';
  jacobianwarped.name = 'Jacobian determinant';
  jacobianwarped.def  = @(val)cat_get_defaults('output.jacobian.warped', val{:});
  jacobianwarped.help = {
    'This is the option to save the Jacobian determinant, which expresses local volume changes. This image can be used in a pure deformation based morphometry (DBM) design. Please note that the affine part of the deformation field is ignored. Thus, there is no need for any additional correction for different brain sizes using ICV.'
    ''
  };

  %------------------------------------------------------------------------

  if expert
    native.def  = @(val)cat_get_defaults('output.label.native', val{:});
    warped.def  = @(val)cat_get_defaults('output.label.warped', val{:});
    dartel.def  = @(val)cat_get_defaults('output.label.dartel', val{:});

    label       = cfg_branch;
    label.tag   = 'label';
    label.name  = 'PVE label image';
    label.val   = {native warped dartel};
    label.help  = {
    'This is the option to save a labeled version of your segmentations for fast visual comparision. Labels are saved as Partial Volume Estimation (PVE) values with different mix classes for GM-WM (2.5) and GM-CSF (1.5). BG=0, CSF=1, GM=2, WM=3, WMH=4 (if WMHC=3), SL=1.5 (if SLC)'
    ''
    };
  else
    labelnative      = native;
    labelnative.tag  = 'labelnative';
    labelnative.name = 'PVE label image in native space';
    labelnative.def  = @(val)cat_get_defaults('output.label.native', val{:});
    labelnative.help = {
    'This is the option to save a labeled version of your segmentations in native space for fast visual comparision and preprocessing quality control. Labels are saved as Partial Volume Estimation (PVE) values with different mix classes for GM-WM (2.5) and GM-CSF (1.5). BG=0, CSF=1, GM=2, WM=3, WMH=4 (if WMHC=3), SL=1.5 (if SLC)'
    ''
    }; 
  end

  %------------------------------------------------------------------------

  modulated        = cfg_menu;
  modulated.tag    = 'mod';
  modulated.name   = 'Modulated normalized';
  if expert
    modulated.labels = {'No','Affine + non-linear (SPM12 default)','Non-linear only'};
    modulated.values = {0 1 2};
  else
    modulated.labels = {'No','Yes'};
    modulated.values = {0 1};
  end
  modulated.help = {
  '"Modulation" is to compensate for the effect of spatial normalisation. Spatial normalisation causes volume changes due to affine transformation (global scaling) and non-linear warping (local volume change). After modulation the resulting modulated images are preserved for the total amount of grey matter signal in the normalised partitions. Thus, modulated images reflect the tissue volumes before spatial normalisation. However, the user is almost always interested in removing the confound of different brain sizes and there are many ways to apply this correction. In contrast to previous VBM versions I now recommend to use total intracranial volume (TIV) as nuisance parameter in an AnCova model. '
  ''
  'Please note that I do not use the SPM modulation where the original voxels are projected into their new location in the warped images because this method introduces aliasing artifacts. Here, I use the scaling by the Jacobian determinants to generate "modulated" data. '
  ''
  };

  native.def    = @(val)cat_get_defaults('output.GM.native', val{:});
  warped.def    = @(val)cat_get_defaults('output.GM.warped', val{:});
  modulated.def = @(val)cat_get_defaults('output.GM.mod', val{:});
  dartel.def    = @(val)cat_get_defaults('output.GM.dartel', val{:});
  grey          = cfg_branch;
  grey.tag      = 'GM';
  grey.name     = 'Grey matter';
  grey.help     = {'Options to save grey matter images.'
  ''
  };
  grey_spm      = grey;
  if expert
    grey.val      = {native warped modulated dartel};
    grey_spm.val  = {warped modulated dartel};
  else
    grey.val      = {native modulated dartel};
    grey_spm.val  = {modulated dartel};
  end

  native.def    = @(val)cat_get_defaults('output.WM.native', val{:});
  warped.def    = @(val)cat_get_defaults('output.WM.warped', val{:});
  modulated.def = @(val)cat_get_defaults('output.WM.mod',    val{:});
  dartel.def    = @(val)cat_get_defaults('output.WM.dartel', val{:});
  white         = cfg_branch;
  white.tag     = 'WM';
  white.name    = 'White matter';
  white.help    = {'Options to save white matter images.'
  ''
  };
  white_spm     = white;
  if expert
    white.val      = {native warped modulated dartel};
    white_spm.val  = {warped modulated dartel};
  else
    white.val      = {native modulated dartel};
    white_spm.val  = {modulated dartel};
  end


  native.def    = @(val)cat_get_defaults('output.CSF.native', val{:});
  warped.def    = @(val)cat_get_defaults('output.CSF.warped', val{:});
  modulated.def = @(val)cat_get_defaults('output.CSF.mod',    val{:});
  dartel.def    = @(val)cat_get_defaults('output.CSF.dartel', val{:});
  csf           = cfg_branch;
  csf.tag       = 'CSF';
  csf.name      = 'Cerebro-Spinal Fluid (CSF)';
  csf.help      = {'Options to save CSF images.'
  ''
  };
  csf_spm       = csf;
  csf.val       = {native warped modulated dartel};
  csf_spm.val   = {warped modulated dartel};

  % head/background tissue classes
  native.def    = @(val)cat_get_defaults('output.TPMC.native', val{:});
  warped.def    = @(val)cat_get_defaults('output.TPMC.warped', val{:});
  modulated.def = @(val)cat_get_defaults('output.TPMC.mod',    val{:});
  dartel.def    = @(val)cat_get_defaults('output.TPMC.dartel', val{:});
  tpmc          = cfg_branch;
  tpmc.tag      = 'TPMC';
  tpmc.name     = 'Tissue Probability Map Classes';
  tpmc.val      = {native warped modulated dartel};
  tpmc.help     = {'Option to save the SPM tissue class 4 to 6: p#*.img, wp#*.img and m[0]wp#*.img.'
  ''
  };

  % WMH
  native.def    = @(val)cat_get_defaults('output.WMH.native', val{:});
  warped.def    = @(val)cat_get_defaults('output.WMH.warped', val{:});
  modulated.def = @(val)cat_get_defaults('output.WMH.mod',    val{:});
  dartel.def    = @(val)cat_get_defaults('output.WMH.dartel', val{:});
  wmh           = cfg_branch;
  wmh.tag       = 'WMH';
  wmh.name      = 'White matter hyperintensities (WMHs)';
  wmh.val       = {native warped modulated dartel};
  wmh.help      = {'WARNING: Please note that the detection of WM hyperintensies (WMHs) is still under development and does not have the same accuracy as approaches that additionally consider FLAIR images (e.g. Lesion Segmentation Toolbox)!'
  'Options to save WMH images, if WMHC==3: p7*.img, wp7*.img and m[0]wp7*.img.'
  ''
  };

  % stroke lesions
  native.def    = @(val)cat_get_defaults('output.SL.native', val{:});
  warped.def    = @(val)cat_get_defaults('output.SL.warped', val{:});
  modulated.def = @(val)cat_get_defaults('output.SL.mod',    val{:});
  dartel.def    = @(val)cat_get_defaults('output.SL.dartel', val{:});
  sl           = cfg_branch;
  sl.tag       = 'SL';
  sl.name      = 'Stroke lesions (SLs) - in development';
  sl.val       = {native warped modulated dartel};
  sl.help      = {'WARNING: Please note that the handling of stroke lesions (SLs) is still under development! '
  'To save SL images, SLC has to be active and (SLs has to be labeled): p8*.img, wp8*.img and m[0]wp8*.img.'
  ''
  };

  % main structure atlas
  native.def    = @(val)cat_get_defaults('output.atlas.native', val{:});
  warped.def    = @(val)cat_get_defaults('output.atlas.warped', val{:});
  dartel.def    = @(val)cat_get_defaults('output.atlas.dartel', val{:});
  atlas         = cfg_branch;
  atlas.tag     = 'atlas';
  atlas.name    = 'Atlas label maps';
  if expert>1
      atlas.val     = {native warped dartel};
  else
      atlas.val     = {native dartel};
  end    
  atlas.help    = {
    'WARNING: The functions that create these maps are still under development! This is the option to save an atlas map with major structures (a1*). Odd numbers code the left, even numbers the right hemisphere. Furthermore, AAL and Broadman atlas maps were created based on maps from MRIcron that where adapted to the other VBM maps. Other maps are used from the IBASPM toolbox.  http://www.thomaskoenig.ch/Lester/ibaspm.htmAnatomy toolbox:Alexander Hammers brain atlas from the Euripides project:   www.brain-development.org  Hammers A, Allom R, Koepp MJ, Free SL, Myers R, Lemieux L, Mitchell   TN, Brooks DJ, Duncan JS. Three-dimensional maximum probability atlas   of the human brain, with particular reference to the temporal lobe.   Hum Brain Mapp 2003, 19: 224-247.'
  ''
  };

  % cortical thickness maps
  native.def    = @(val)cat_get_defaults('output.ct.native', val{:});
  warped.def    = @(val)cat_get_defaults('output.ct.warped', val{:});
  dartel.def    = @(val)cat_get_defaults('output.ct.dartel', val{:});
  native.val  = {0};
  warped.val  = {0};
  dartel.val  = {0};
  gmt         = cfg_branch;
  gmt.tag     = 'ct';
  gmt.name    = 'Cortical Thickness';
  gmt.val     = {native warped dartel};
  gmt.help    = {
    'Options to save cortical thickess maps (experimental).'
    ''
  };

  % percentual position maps - uses defaults from thickness
  native.def    = @(val)cat_get_defaults('output.pp.native', val{:});
  warped.def    = @(val)cat_get_defaults('output.pp.warped', val{:});
  dartel.def    = @(val)cat_get_defaults('output.pp.dartel', val{:});
  native.val    = {0};
  warped.val    = {0};
  dartel.val    = {0};
  pp            = cfg_branch;
  pp.tag        = 'pp';
  pp.name       = 'Percentage Position';
  pp.val        = {native warped dartel};
  pp.help       = {
    'Options to save percentage position maps (experimental).'
    ''
  };

  warps = cfg_menu;
  warps.tag    = 'warps';
  warps.name   = 'Deformation Fields';
  warps.labels = {
      'No'
      'Image->Template (forward)'
      'Template->Image (inverse)'
      'inverse + forward'};
  warps.values = {[0 0],[1 0],[0 1],[1 1]};
  warps.def    = @(val)cat_get_defaults('output.warps', val{:});
  warps.help   = {
    'Deformation fields can be saved to disk, and used by the Deformations Utility and/or applied to coregistered data from other modalities (e.g. fMRI). For spatially normalising images to MNI space, you will need the forward deformation, whereas for spatially normalising (eg) GIFTI surface files, you''ll need the inverse. It is also possible to transform data in MNI space on to the individual subject, which also requires the inverse transform. Deformations are saved as .nii files, which contain three volumes to encode the x, y and z coordinates.'
  ''
  };

  %% ------------------------------------------------------------------------
  
  [ROI,sROI]  = cat_conf_ROI(expert);       % ROI options
  
  output      = cfg_branch;
  output.tag  = 'output';
  output.name = 'Writing options';
  if expert==2
    output.val  = {surface surf_measures ROI sROI grey white csf gmt pp wmh sl tpmc atlas label bias las jacobianwarped warps}; 
  elseif expert==1
    output.val  = {surface ROI sROI grey white csf wmh sl tpmc atlas label bias las jacobianwarped warps};
  else
    output.val  = {surface ROI grey white labelnative bias jacobianwarped warps};
  end
  output.help = {
  'There are a number of options about what kind of data you like save. The routine can be used for saving images of tissue classes, as well as bias corrected images. The native space option will save a tissue class image (p*) that is in alignment with the original image. You can also save spatially normalised versions - both with (m[0]wp*) and without (wp*) modulation. In the cat toolbox, the voxel size of the spatially normalised versions is 1.5 x 1.5 x 1.5mm as default. The saved images of the tissue classes can directly be used for doing voxel-based morphometry (both un-modulated and modulated). All you need to do is smooth them and do the stats (which means no more questions on the mailing list about how to do "optimized VBM"). Please note that many less-common options are only available in expert mode (e.g. CSF, labels, atlas maps).'
  ''
  'Modulation is to compensate for the effect of spatial normalisation. When warping a series of images to match a template, it is inevitable that volumetric differences will be introduced into the warped images. For example, if one subject''s temporal lobe has half the volume of that of the template, then its volume will be doubled during spatial normalisation. This will also result in a doubling of the voxels labeled grey matter. In order to remove this confound, the spatially normalised grey matter (or other tissue class) is adjusted by multiplying by its relative volume before and after warping. If warping results in a region doubling its volume, then the correction will halve the intensity of the tissue label. This whole procedure has the effect of preserving the total amount of grey matter signal in the normalised partitions.'
  ''
  'A deformation field is a vector field, where three values are associated with each location in the field. The field maps from co-ordinates in the normalised image back to co-ordinates in the original image. The value of the field at co-ordinate [x y z] in the normalised space will be the co-ordinate [x'' y'' z''] in the original volume. The gradient of the deformation field at a co-ordinate is its Jacobian matrix, and it consists of a 3x3 matrix:'
  '%   /                      \%   | dx''/dx  dx''/dy dx''/dz |%   |                       |%   | dy''/dx  dy''/dy dy''/dz |%   |                       |%   | dz''/dx  dz''/dy dz''/dz |%   \                      /'
  ''
  'The value of dx''/dy is a measure of how much x'' changes if y is changed by a tiny amount. The determinant of the Jacobian is the measure of relative volumes of warped and unwarped structures.  The modulation step simply involves multiplying by the relative volumes.'};

  %%
  %------------------------------------------------------------------------
  % R1173
  %------------------------------------------------------------------------
  warped.def    = @(val)cat_get_defaults1173('output.jacobian.warped', val{:});
  jacobian      = cfg_branch;
  jacobian.tag  = 'jacobian';
  jacobian.name = 'Jacobian determinant';
  jacobian.val  = {warped};
  jacobian.help = {
    'This is the option to save the Jacobian determinant, which expresses local volume changes. This image can be used in a pure deformation based morphometry (DBM) design. Please note that the affine part of the deformation field is ignored. Thus, there is no need for any additional correction for different brain sizes using ICV.'
  ''
  };

%   % additional segmentation versions % EXPLOREASLHACK, DISABLE PREVIOUS
%   CAT12 VERSIONS
%   [ROI1173,atlases1173]         = cat_conf_ROI1173(expert);       % ROI options
%   [ROI1445,atlases1445]         = cat_conf_ROI1445(expert);       % ROI options
% %  [ROI1585,atlases1585]         = cat_conf_ROI1585(expert);       % ROI options
%   
%   output1173            = output;
%   if expert==2
%     output1173.val      = {surface ROI1173 atlases1173 grey white csf wmh tpmc atlas label bias las jacobian warps}; 
%   elseif expert==1
%     output1173.val      = {surface ROI1173 atlases1173 grey white csf wmh label bias las jacobian warps};
%   else
%     output1173.val      = {surface ROI1173 grey white jacobian warps};
%   end
%   
%   output1445            = output;
%   if expert==2
%     output1445.val      = {surface ROI1445 atlases1445 grey white csf wmh tpmc atlas label bias las jacobian warps}; 
%   elseif expert==1
%     output1445.val      = {surface ROI1445 atlases1445 grey white csf wmh label bias las jacobian warps};
%   else
%     output1445.val      = {surface ROI1445 grey white jacobian warps};
%   end
% 
% if 0
%   output1585            = output;
%   if expert==2
%     output1585.val      = {surface ROI1585 atlases1585 grey white csf wmh tpmc atlas label bias las jacobian warps}; 
%   elseif expert==1
%     output1585.val      = {surface ROI1585 atlases1585 grey white csf wmh label bias las jacobian warps};
%   else
%     output1585.val      = {surface ROI1585 atlases1585 grey white jacobian warps};
%   end
% end
  
	output_spm  = output; 
	if expert==2
		output_spm.val  = {ROI surface grey_spm white_spm csf_spm label jacobianwarped warps}; 
	elseif expert==1
		output_spm.val  = {ROI surface grey_spm white_spm csf_spm label jacobianwarped warps};
	else % also CSF output because it is requiered as input ...
		output_spm.val  = {ROI surface grey_spm white_spm csf_spm labelnative jacobianwarped warps};
	end

return
%------------------------------------------------------------------------