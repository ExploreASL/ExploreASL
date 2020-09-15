function [ROI,sROI] = cat_conf_ROI1585(expert)
%_______________________________________________________________________
% wrapper for calling CAT ROI options
%_______________________________________________________________________
% Robert Dahnke and Christian Gaser
% $Id: cat_conf_ROI.m 1571 2020-02-27 14:51:11Z gaser $
%_______________________________________________________________________


if nargin == 0
  try
    expert = cat_get_defaults1585('extopts.expertgui');
  catch %#ok<CTCH>
    expert = 0; 
  end
end

%% ------------------------------------------------------------------------
% Parameter to choose between different ways to extract data by ROIs.
% Inactive due to missing implementation and evaluation of the optimized 
% space ROI analysis. 
% RD 20160112
%------------------------------------------------------------------------

noROI        = cfg_branch;
noROI.tag    = 'noROI';
noROI.name   = 'No ROI processing';
noROI.help   = {'No ROI processing'};

exatlas  = cat_get_defaults1585('extopts.atlas'); 
matlas = {}; mai = 1; atlaslist = {}; 
for ai = 1:size(exatlas,1)
  if exatlas{ai,2}<=expert && exist(exatlas{ai,1},'file')
    [pp,ff]  = spm_fileparts(exatlas{ai,1}); 

    % if output.atlases.ff does not exist then set it by the default file value
    if isempty(cat_get_defaults1585(['output.atlases.' ff]))
      cat_get_defaults1585(['output.atlases.' ff], exatlas{ai,4})
    end
    atlaslist{end+1,1} = ff; 

    license = {'' ' (no commercial use)' ' (free academic use)'}; 
    if size(exatlas,2)>4
      lic = exatlas{ai,5}; 
    else
      switch ff
        case 'hammers', lic = 2; 
        case 'lpba40' , lic = 1; 
        otherwise,      lic = 0; 
      end
    end
    
    matlas{mai}        = cfg_menu;
    matlas{mai}.tag    = ff;
    matlas{mai}.name   = [ff license{lic+1}]; 
    matlas{mai}.labels = {'No','Yes'};
    matlas{mai}.values = {0 1};
    matlas{mai}.def    = eval(sprintf('@(val) cat_get_defaults1585(''output.atlases.%s'', val{:});',ff)); 
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
        ['No atlas readme text file "' txtfile '"!']
      };
    end
    mai = mai+1; 
  else
    [pp,ff]  = spm_fileparts(exatlas{ai,1}); 
    
    if ~isempty(cat_get_defaults1585(['output.atlases.' ff]))
      cat_get_defaults1585(['output.atlases.' ff],'rmfield');
    end
  end
end

ownatlas              = cfg_files;
ownatlas.tag          = 'ownatlas';
ownatlas.name         = 'own atlas maps';
ownatlas.help         = { 
  sprintf([
    'Select images that should be used as atlas maps.  ' ...
    'The maps should only contain positive integers for regions of interest.  ' ...
    'You can use a CSV-file with the same name as the atlas to define region ' ...
    'names similar to the CSV-files of other atlas files in "%s".  ' ...
    'The CSV-file should have an header line containing the number of the ROI "ROIid", ' ...
    'the abbreviation of the ROI "ROIabbr" and the full name of the ROI "ROIname".  ' ...
    'The GM, WM, and CSF values will be extracted for all regions. '], ...
    fullfile( spm('dir'), 'toolbox', 'cat12', 'templates_volumes') ); 
  ''};
ownatlas.filter       = 'image';
ownatlas.ufilter      = '.*';
ownatlas.val{1}       = {''};
ownatlas.num          = [0 Inf];

atlases          = cfg_branch;
atlases.tag      = 'atlases';
atlases.name     = 'Atlases';
atlases.val      = [matlas,{ownatlas}];
atlases.help     = {'Writing options of ROI atlas maps.'
''
};


ROI        = cfg_choice;
ROI.tag    = 'ROImenu';
ROI.name   = 'Process Volume ROIs';
if cat_get_defaults1585('output.ROI')>0
  ROI.val  = {atlases};
else
  ROI.val  = {noROI};
end
ROI.values = {noROI atlases};
ROI.help   = {
'Export of ROI data of volume to a xml-files. '
'For further information see atlas specific text files in "templates_volumes" CAT12 subdir. '
''
'For thickness estimation the projection-based thickness (PBT) [Dahnke:2012] is used that average cortical thickness for each GM voxel. '
''
'There are different atlas maps available: '
}; 

%%
mai = 1; 
for ali=1:numel(atlaslist)
  if any(~cellfun('isempty',strfind(atlaslist(ali),'hammers')))
    ROI.help = [ROI.help; strrep({
        '(MAI) Hammers (68 CSF/GM/[WM] ROIs of 20 subjects, 2003):'
        '    Alexander Hammers brain atlas from the Euripides project (www.brain-development.org).'
        '    Hammers et al. Three-dimensional maximum probability atlas of the human brain, with particular reference to the temporal lobe. Hum Brain Mapp 2003, 19: 224-247.'
        ''},'MAI',num2str(mai,'%d'))]; mai = mai+1; 
  end
  if any(~cellfun('isempty',strfind(atlaslist(ali),'neuromorphometrics')))
    ROI.help = [ROI.help; strrep({
        '(MAI) Neuromorphometrics (142 GM ROIs of 15 subjects, 2012):'
        '    Maximum probability tissue labels derived from the MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling'
        '    https://masi.vuse.vanderbilt.edu/workshop2012/index.php/Challenge_Details'
        ''},'MAI',num2str(mai,'%d'))]; mai = mai+1; 
  end
  if any(~cellfun('isempty',strfind(atlaslist(ali),'lpba40')))
    ROI.help = [ROI.help; strrep({
        '(MAI) LPBA40 (56 GM ROIs of 40 subjects, 2008):'
        '    The LONI Probabilistic Brain Atlas (LPBA40) is a series of maps of brain anatomical regions. These maps were estimated from a set of whole-head MRI of 40 human volunteers. Each MRI was manually delineated to identify a set of 56 structures in the brain, most of which are within the cortex. These delineations were then transformed into a common atlas space to obtian a set of coregistered anatomical labels. The original MRI data were also transformed into the atlas space. '
        '    Shattuck et al. 2008. Construction of a 3D Probabilistic Atlas of Human Cortical Structures, NeuroImage 39 (3): 1064-1070. DOI:	10.1016/j.neuroimage.2007.09.031'
        ''},'MAI',num2str(mai,'%d'))]; mai = mai+1; 
  end
  if any(~cellfun('isempty',strfind(atlaslist(ali),'cobra')))
    ROI.help = [ROI.help; strrep({
        '(MAI) COBRA (1 GM/WM ROI in amgdala, 2 combined GM/WM ROIs in hippocampus and 13 GM/WM ROIs in cerebellum of 5 subjects):'
        '    The Cobra atlas is build from 3 atlases that are provided by the Computational Brain Anatomy Laboratory at the Douglas Institute (CoBra Lab). The 3 atlases are based on high-resolution (0.3mm isotropic voxel size) images of the amygdala, hippocampus and the cerebellum. Some of the hippocampus subfields were merged because of their small size (CA1/CA2/CA3/stratum radiatum/subiculum/stratum lacunosum/stratum moleculare). Please note that the original labels were changed in order to allow a combined atlas. '
        '    Entis JJ, Doerga P, Barrett LF, Dickerson BC. A reliable protocol for the manual segmentation of the human amygdala and its subregions using ultra-high resolution MRI. Neuroimage. 2012;60(2):1226-35.'
        '    Winterburn JL, Pruessner JC, Chavez S, et al. A novel in vivo atlas of human hippocampal subfields using high-resolution 3 T magnetic resonance imaging.  Neuroimage. 2013;74:254-65.'
        '    Park, M.T., Pipitone, J., Baer, L., Winterburn, J.L., Shah, Y., Chavez, S., Schira, M.M., Lobaugh, N.J., Lerch, J.P., Voineskos, A.N., Chakravarty, M.M. Derivation of high-resolution MRI atlases of the human cerebellum at 3T and segmentation using multiple automatically generated templates. Neurimage. 2014; 95: 217-31.'
        ''},'MAI',num2str(mai,'%d'))]; mai = mai+1; 
  end
  if any(~cellfun('isempty',strfind(atlaslist(ali),'ibsr')))
    ROI.help = [ROI.help; strrep({
        '(MAI) IBSR (32 CSF/GM ROIs of 18 subjects, 2004):'
        '    See IBSR terms "http://www.nitrc.org/projects/ibsr"'
        ''},'MAI',num2str(mai,'%d'))]; mai = mai+1; 
  end
  if any(~cellfun('isempty',strfind(atlaslist(ali),'aal3')))
    ROI.help = [ROI.help; strrep({
        '(MAI) AAL3 (170 GM ROIs of 1 subject, 2020):'
        '    Tzourio-Mazoyer et al., Automated anatomical labelling of activations in spm using a macroscopic anatomical parcellation of the MNI MRI single subject brain. Neuroimage 2002, 15: 273-289.'
        '    Rolls et al., Automated anatomical labelling atlas 3. Neuroimage 2020; 206:116189.'
        ''},'MAI',num2str(mai,'%d'))]; mai = mai+1; 
  end
  if any(~cellfun('isempty',strfind(atlaslist,'mori')))
    ROI.help = [ROI.help; strrep({
        '(MAI) MORI (128 GM/WM ROIs of 1 subject, 2009):'
        '    Oishi et al. Atlas-based whole brain white matter analysis using large deformation diffeomorphic metric mapping: application to normal elderly and Alzheimer''s disease participants. 2009'
        ''},'MAI',num2str(mai,'%d'))]; mai = mai+1; 
  end
  if any(~cellfun('isempty',strfind(atlaslist,'anatomy')))
    ROI.help = [ROI.help; strrep({
        '(MAI) Anatomy (44 GM/WM ROIs in 10 post-mortem subjects, 2014):'
        '    Eickhoff SB, Stephan KE, Mohlberg H, Grefkes C, Fink GR, Amunts K, Zilles K. A new SPM toolbox for combining probabilistic cytoarchitectonic maps and functional imaging data. NeuroImage 25(4), 1325-1335, 2005'
        ''},'MAI',num2str(mai,'%d'))]; mai = mai+1; 
  end
end
 


%% ------------------------------------------------------------------------
% Surface Atalases
% RD 20190416
%------------------------------------------------------------------------

nosROI        = cfg_branch;
nosROI.tag    = 'noROI';
nosROI.name   = 'No surface ROI processing';
nosROI.help   = {'No surface ROI processing'};

exatlas  = cat_get_defaults1585('extopts.satlas'); 
matlas = {}; mai = 1; atlaslist = {}; 
for ai = 1:size(exatlas,1)
  if exatlas{ai,3}<=expert && ~isempty(exatlas{ai,2})
    [pp,ff]  = spm_fileparts( exatlas{ai,2} ); 
    name = exatlas{ai,1}; 

    % if output.atlases.ff does not exist then set it by the default file value
    if isempty(cat_get_defaults1585(['output.satlases.' name]))
      cat_get_defaults1585(['output.atlases.' ff], exatlas{ai,4})
    end
    atlaslist{end+1,1} = name; 

    if cat_get_defaults1585('extopts.expertgui') 
      if strcmp(spm_str_manip(pp,'t'),'atlases_surfaces_32k')
        addname = ' (32k)';
      elseif strcmp(spm_str_manip(pp,'t'),'atlases_surfaces')
        addname = ' (160k)';
      else
        addname = '';
      end
    else
      addname = '';
    end
    
    matlas{mai}        = cfg_menu;
    matlas{mai}.tag    = name;
    matlas{mai}.name   = [name addname]; 
    matlas{mai}.labels = {'No','Yes'};
    matlas{mai}.values = {0 1};
    matlas{mai}.def    = eval(sprintf('@(val) cat_get_defaults1585(''output.atlases.%s'', val{:});',ff)); 
    
    txtfile = fullfile(pp,[name '.txt']);
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
        ['No atlas readme text file "' txtfile '"!']
      };
    end
    mai = mai+1; 
  else
    name = exatlas{ai,1}; 
    
    if ~isempty(cat_get_defaults1585(['output.atlases.' name]))
      cat_get_defaults1585(['output.atlases.' name],'rmfield');
    end
  end
end

ownsatlas          = ownatlas;
ownsatlas.filter   = '';
ownsatlas.ufilter  = '.*'; 
ownsatlas.help     = { 
  'Select FreeSurfer surface annotation files (*.annot), FreeSurfer CURV-files, or GIFTI surfaces with positve integer with 32k or 160k faces. ';
  ''};

satlases          = cfg_branch;
satlases.tag      = 'satlases';
satlases.name     = 'Surface atlases';
satlases.val      = [matlas,{ownatlas}];
satlases.help     = {'Writing options for surface ROI atlas maps.'
''
};


sROI        = cfg_choice;
sROI.tag    = 'sROImenu';
sROI.name   = 'Process Surface ROIs';
if cat_get_defaults1585('output.surface')>0 && cat_get_defaults1585('output.ROI')>0
  sROI.val  = {satlases};
else
  sROI.val  = {nosROI};
end
sROI.values = {nosROI satlases};
sROI.help   = {
'Export of ROI data of volume to a xml-files. '
'For further information see atlas specific text files in "templates_volumes" CAT12 subdir. '
''
'For thickness estimation the projection-based thickness (PBT) [Dahnke:2012] is used that average cortical thickness for each GM voxel. '
''
'There are different atlas maps available: '
}; 
