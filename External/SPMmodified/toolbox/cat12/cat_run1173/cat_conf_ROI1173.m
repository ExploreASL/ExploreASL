function [ROI,atlases] = cat_conf_ROI1173(expert)
%_______________________________________________________________________
% wrapper for calling CAT ROI options
%_______________________________________________________________________
% Robert Dahnke and Christian Gaser
% $Id: cat_conf_ROI1173.m 1571 2020-02-27 14:51:11Z gaser $
%_______________________________________________________________________

%------------------------------------------------------------------------
% Parameter to choose between different ways to extract data by ROIs.
% Inactive due to missing implementation and evaluation of the optimized 
% space ROI analysis. 
% RD 20160112
%------------------------------------------------------------------------
ROI        = cfg_menu;
ROI.tag    = 'ROI';
ROI.name   = 'ROI analyis';
ROI.labels = {'No ROI analyis','ROI analysis'};
ROI.values = {0 1};
ROI.def    = @(val)cat_get_defaults1173('output.ROI', val{:});
ROI.help   = {
  'Export of ROI data of volume, intensity, and thickness to a xml-files. '
  'For further information see atlas specific text files in "templates_volumes" CAT12 subdir. '
  ''
  'For thickness estimation the projection-based thickness (PBT) [Dahnke:2012] is used that average cortical thickness for each GM voxel. '
  ''
  'There are different atlas maps available: '
  '(1) Hammers (68 CSF/GM/[WM] ROIs of 20 subjects, 2003):'
  '    Alexander Hammers brain atlas from the Euripides project (www.brain-development.org).'
  '    Hammers et al. Three-dimensional maximum probability atlas of the human brain, with particular reference to the temporal lobe. Hum Brain Mapp 2003, 19: 224-247.'
  ''
  '(2) Neuromorphometrics (142 GM ROIs of 15 subjects, 2012):'
  '    Maximum probability tissue labels derived from the MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling'
  '    https://masi.vuse.vanderbilt.edu/workshop2012/index.php/Challenge_Details'
  ''
  '(3) LPBA40 (56 GM ROIs of 40 subjects, 2008):'
  '    The LONI Probabilistic Brain Atlas (LPBA40) is a series of maps of brain anatomic regions. These maps were produced from a set of whole-head MRI of 40 human volunteers. Each MRI was manually delineated to identify a set of 56 structures in the brain, most of which are within the cortex. These delineations were then transformed into a common atlas space to produce a set of coregistered anatomical labels. The original MRI data were also transformed into the atlas space. '
  '    Shattuck et al. 2008. Construction of a 3D Probabilistic Atlas of Human Cortical Structures, NeuroImage 39 (3): 1064-1070. DOI:	10.1016/j.neuroimage.2007.09.031'
  ''
  ...'(4) IBSR (32 CSF/GM ROIs of 18 subjects, 2004):'
  ...'    See IBSR terms "http://www.nitrc.org/projects/ibsr"'
  ...''
  ...'(5) AAL (122 GM ROIs of 1 subject, 2002):'
  ...'    Combining probabilistic cytoarchitectonic maps and functional imaging data of a single Brain.'
  ...'    Tzourio-Mazoyer et al., Automated anatomical labelling of activations in spm using a macroscopic anatomical parcellation of the MNI MRI single subject brain. Neuroimage 2002, 15: 273-289.'
  ...'    http://www.fz-juelich.de/inm/inm-1/DE/Forschung/_docs/SPMAnatomyToolbox/SPMAnatomyToolbox_node.html'
  ...''
  ...'(6) MORI (128 GM/WM ROIs of 1 subject, 2009):'
  ...'    Oishi et al. Atlas-based whole brain white matter analysis using large deformation diffeomorphic metric mapping: application to normal elderly and Alzheimer''s disease participants. 2009'
  ... ''
  ... '() Anatomy ():
  ... '() COBRA ():'
  ''
  ''
};

if expert
  %%
  atlas  = cat_get_defaults1173('extopts.atlas'); 
  matlas = {}; mai = 1; 
  for ai = 1:size(atlas,1)
    if atlas{ai,2}<=expert && exist(atlas{ai,1},'file')
      [pp,ff,ee]  = spm_fileparts(atlas{ai,1}); 

      cat_get_defaults1173(['output.atlases.' ff], atlas{ai,4})
      
      matlas{mai}        = cfg_menu;
      matlas{mai}.tag    = ff;
      matlas{mai}.name   = ff; 
      matlas{mai}.labels = {'No','Yes'};
      matlas{mai}.values = {0 1};
      matlas{mai}.def    = eval(sprintf('@(val) cat_get_defaults1173(''output.atlases.%s'', val{:});',ff)); 
        %@(val)cat_get_defaults1173(['output.atlases.' ff], val{:});
      txtfile = fullfile(pp,[ff '.txt']);
      if exist(txtfile,'file')
        fid = fopen(txtfile,'r');
        txt = textscan(fid,'%s','delimiter','\n');
        fclose(fid);
        matlas{mai}.help   = [{ 
          'Processing flag of this atlas map.'
          ''
          }
          txt{1}];
      else
        matlas{mai}.help   = {
          'Processing flag of this atlas map.'
          ''
          ['No atlas readme text file "' textfile '"!']
        };
      end
      mai = mai+1; 
    end
  end
  
  atlases          = cfg_branch;
  atlases.tag      = 'atlases';
  atlases.name     = 'Atlases';
  atlases.val      = matlas;
  atlases.help     = {'Writing options of ROI atlas maps.'
  ''
  };
else
  atlases = struct();
end

end
 

 