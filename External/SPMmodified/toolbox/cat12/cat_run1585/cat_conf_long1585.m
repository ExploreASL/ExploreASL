function varargout = cat_conf_long1585(varargin)
% Configuration file for longitudinal data
%
% Christian Gaser
% $Id: cat_conf_long.m 1519 2019-11-19 10:48:29Z gaser $

newapproach = 0; 

if newapproach && nargin>0 
  [dep,varargout{1},varargout{2}] = vout_long(varargin{1});
  return
end

try
  expert = cat_get_defaults1585('extopts.expertgui');
catch %#ok<CTCH>
  expert = 0; 
end

% try to estimate number of processor cores
try
  numcores = cat_get_defaults1585('extopts.nproc');
  % because of poor memory management use only half of the cores for windows
  if ispc
    numcores = round(numcores/2);
  end
  numcores = max(numcores,1);
catch
  numcores = 0;
end

% force running in the foreground if only one processor was found or for compiled version
if numcores == 1 | isdeployed, numcores = 0; end

%------------------------------------------------------------------------
nproc         = cfg_entry;
nproc.tag     = 'nproc';
nproc.name    = 'Split job into separate processes';
nproc.strtype = 'w';
nproc.val     = {numcores};
nproc.num     = [1 1];
nproc.help    = {
    'In order to use multi-threading the CAT12 segmentation job with multiple subjects can be split into separate processes that run in the background. You can even close Matlab, which will not affect the processes that will run in the background without GUI. If you do not want to run processes in the background then set this value to 0.'
    ''
    'Keep in mind that each process needs about 1.5..2GB of RAM, which should be considered to choose the appropriate  number of processes.'
    ''
    'Please further note that no additional modules in the batch can be run except CAT12 segmentation. Any dependencies will be broken for subsequent modules.'
  };
%------------------------------------------------------------------------

mov = cfg_files;
mov.name = 'Longitudinal data for this subject';
mov.tag  = 'mov';
mov.filter = 'image';
mov.num  = [2 Inf];
mov.help   = {...
'These are the data of the same subject.'};
%------------------------------------------------------------------------

subj = cfg_branch;
subj.name = 'Subject';
subj.tag = 'subj';
subj.val = {mov};
subj.help = {...
'Images of the same subject.'};

%------------------------------------------------------------------------

esubjs         = cfg_repeat;
esubjs.tag     = 'esubjs';
esubjs.name    = 'Data';
esubjs.values  = {subj};
esubjs.num     = [1 Inf];
esubjs.help = {...
'Specify data for each subject.'};


%------------------------------------------------------------------------
delete_temp        = cfg_menu;
delete_temp.tag    = 'delete_temp';
delete_temp.name   = 'Delete temporary files';
delete_temp.labels = {'No','Yes'};
delete_temp.values = {0 1};
delete_temp.val    = {1};
delete_temp.help = {
'Temporary files such as the native segmentations or deformation fields are usually removed after preprocessing. However, if you like to keep these files you can use this option.'
''
};

%------------------------------------------------------------------------
extopts = cat_conf_extopts1585(expert);
opts    = cat_conf_opts1585(expert);
output  = cat_conf_output(expert); 
%------------------------------------------------------------------------

long = cfg_exbranch;
long.name = 'CAT12: Segment longitudinal data';
long.tag  = 'long';
if newapproach % new way - not working
  
  % remove major output fields
  clear FN; for vi = 1:numel(output.val), FN{vi} = output.val{vi}.tag; end
  removefields = {'warps','jacobianwarped'};
  for vi = 1:numel(removefields)
    output.val(find(cellfun('isempty',strfind(FN,removefields{vi}))==0)) = []; 
  end
  
  % remove subfields
  removefields = {'native'};
  for vim = 1:numel(output.val)
    clear FN; for vi = 1:numel(output.val{vim}.val), FN{vi} = output.val{vim}.val{vi}.tag; end
    if numel(output.val{vim}.val)
      for vi = 1:numel(removefields)
        output.val{vim}.val(find(cellfun('isempty',strfind(FN,removefields{vi}))==0)) = []; 
      end
    end
  end

  if expert
    output.val = [output.val, delete_temp]; 
  end
  long.val  = {esubjs,nproc,opts,extopts,output}; 
  long.vout = @vout_long2;
else
  % old appraoch
  %------------------------------------------------------------------------
  modulate        = cfg_menu;
  modulate.tag    = 'modulate';
  modulate.name   = 'Modulated GM/WM segmentations';
  modulate.labels = {'No','Yes'};
  modulate.values = {0 1};
  modulate.val    = {1};
  modulate.help = {
  '"Modulation" is to compensate for the effect of spatial normalisation. Spatial normalisation causes volume changes due to affine transformation (global scaling) and non-linear warping (local volume change). After modulation the resulting modulated images are preserved for the total amount of grey matter signal in the normalised partitions. Thus, modulated images reflect the tissue volumes before spatial normalisation. However, the user is almost always interested in removing the confound of different brain sizes and there are many ways to apply this correction. In contrast to previous VBM versions I now recommend to use total intracranial volume (TIV) as nuisance parameter in an AnCova model. '
  ''
  'Please note that I do not use the SPM modulation where the original voxels are projected into their new location in the warped images because this method introduces aliasing artifacts. Here, I use the scaling by the Jacobian determinants to generate "modulated" data. '
  ''
  'For longitudinal data the modulation is actually not necessary because normalization estimates for one subject are the same for all time points and thus modulation will be also the same for all time points. However, modulation might be useful if you want to compare the baseline images in a cross-sectional design in order to test whether there are any differences between the groups at the beginning of the longitudinal study. '
  ''
  };


  %------------------------------------------------------------------------
  dartel        = cfg_menu;
  dartel.tag    = 'dartel';
  dartel.name   = 'DARTEL export of average image';
  if expert
    dartel.labels = {'No','Rigid (SPM12 default)','Affine','Both'};
    dartel.values = {0 1 2 3};
  else
    dartel.labels = {'No','Rigid (SPM12 default)','Affine'};
    dartel.values = {0 1 2};
  end
  dartel.val    = {0};
  dartel.help   = {
  'This option is to export data into a form that can be used with DARTEL. The SPM default is to only apply rigid body transformation. However, a more appropriate option is to apply affine transformation, because the additional scaling of the images requires less deformations to non-linearly register brains to the template.'
  ''
  'Please note, that this option is only useful if you intend to create a customized DARTEl template for your longittudinal data. The DARTEL exported segmentations is saved for the the average image of all time points for one subject and can be used in order to create a customized template with the DARTEL toolbox. The resulting flow fields can be finally applied to the respective native segmentations (e.g. p1/p2 images) to obtain normalized segmentations according to the newly created DARTEL template.'
  ''
  };
  
  ROI       = output.val{2}; 
  
  output.val(2:end) = [];
  
  % use multithreading only if availabe
  if numcores > 1 & ~isdeployed
    if expert
      long.val  = {esubjs,nproc,opts,extopts,output,ROI,modulate,dartel,delete_temp};
    else
      long.val  = {esubjs,nproc,opts,extopts,output,ROI,modulate,dartel};
    end
  else
    long.val  = {esubjs,opts,extopts,output,ROI,modulate,dartel};
  end
  
  long.vout = @vout_long;
end
long.prog = @cat_long_multi_run;

long.help = {
'This option provides customized processing of longitudinal data. Please note that this processing pipeline was optimized for processing and detecting small changes over time as response to short-time plasticity effects (e.g. due to learning and training). This pipelines will not work properly for large longitudinal changes where large parts of the brain will change over time (e.g. atropy due to Alzheimers disease or ageing). This is due to the effect that the spatial normalization parameters are estimated using a mean image of all time points and subsequently applied to all time points. If large atrophy occurs between the time points this can lead to a shift of tissue borders and might result in areas of decreased volumes over time that are surrounded by areas of increased volumes due to this shifting issues. For data with large volume changes over time I would recommend to use the cross-sectional pipeline or the longitudinal toolbox in SPM12.'
''
};

%------------------------------------------------------------------------
varargout{1} = long; 
return;
%------------------------------------------------------------------------
 

%------------------------------------------------------------------------
function dep = vout_long(job)
for k=1:numel(job.subj)
    cdep            = cfg_dep;
    cdep.sname      = sprintf('Segmented longitudinal data (Subj %d)',k);
    cdep.src_output = substruct('.','sess','()',{k},'.','files');
    cdep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    if k == 1
        dep = cdep;
    else
        dep = [dep cdep];
    end
end;
% add all surface/thickness files! of all timepoints and all subjects
if job.output.surface
    for k=1:numel(job.subj)
        dep(end+1)          = cfg_dep;
        dep(end).sname      = 'mwp1 Images';
        dep(end).src_output = substruct('.','mwp1','()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    for k=1:numel(job.subj)
        dep(end+1)          = cfg_dep;
        dep(end).sname      = 'Left Central Surfaces';
        dep(end).src_output = substruct('.','surf','()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
    end
    for k=1:numel(job.subj)
        dep(end+1)          = cfg_dep;
        dep(end).sname      = 'Left Thickness';
        dep(end).src_output = substruct('.','thick','()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
    end
    for k=1:numel(job.subj)
        dep(end+1)          = cfg_dep;
        dep(end).sname      = 'CAT Report';
        dep(end).src_output = substruct('.','catreport','()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','xml','strtype','e'}});
    end
    for k=1:numel(job.subj)
        dep(end+1)          = cfg_dep;
        dep(end).sname      = 'ROI XML File';
        dep(end).src_output = substruct('.','catroi','()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','xml','strtype','e'}});
    end
end


%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [dep,out,inputs] = vout_long2(job)
    inputs = cell(1, numel(job.subj));

    if cat_get_defaults1585('extopts.subfolders')
      mrifolder  = 'mri';
      surffolder = 'surf';
    else
      mrifolder  = '';
      surffolder = ''; 
    end
    
    for i=1:numel(job.subj),
      %%
        out.subj(i).warps = cell(1,1);
        if iscell(job.subj(i).mov)
            [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{1});
        else
            [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov);
        end
        out.subj(i).warps{1} = fullfile(pth,mrifolder,['avg_y_', nam, ext, num]);

        out.subj(i).files = cell(numel(cellstr(job.subj(i).mov)),1);
        m = numel(cellstr(job.subj(i).mov)); % number of scans of this subject
%%
        data = cell(m,1);
        for j=1:m
            if iscell(job.subj(i).mov)
              [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{j});
            else
              [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov);
            end
            % for output (DEP)
            volumes = {
              'GM'  1; 
              'WM'  2; 
              'CSF' 3; 
              'WMH' 7;
              };

            for txi = 1:size(volumes,1)
              %%
              tissue = { 
                volumes{txi,1}  'warped' 1      sprintf('wp%d',volumes{txi,2})   sprintf('wp%dr',volumes{txi,2})    ''; 
                volumes{txi,1}  'mod'    [1 3]  sprintf('mwp%d',volumes{txi,2})  sprintf('mwp%dr',volumes{txi,2})   ''; 
                volumes{txi,1}  'mod'    [2 3]  sprintf('m0wp%d',volumes{txi,2}) sprintf('m0wp%dr',volumes{txi,2})  ''; 
                volumes{txi,1}  'dartel' [1 3]  sprintf('rp%da',volumes{txi,2})  sprintf('rp%dr',volumes{txi,2})    '_affine'; 
                volumes{txi,1}  'dartel' [2 3]  sprintf('rp%dr',volumes{txi,2})  sprintf('rp%dr',volumes{txi,2})    '_affine'; 
              };
              for ti = 1:size(tissue,1)
                if isfield(job.output,tissue{ti,1}) && isfield(job.output.(tissue{ti,1}),tissue{ti,2}) && ...
                  any( job.output.(tissue{ti,1}).(tissue{ti,2}) == tissue{ti,3} )
                  out.subj(i).(tissue{ti,4}){j,1} = fullfile(pth,mrifolder,[(tissue{ti,5}), nam, (tissue{ti,6}), ext, num]);
                end
              end
            end
            %%
            if isfield(job.output,'labelnative') && job.output.labelnative
              out.subj(i).p0{j,1}   = fullfile(pth,mrifolder,['p0r', nam, ext, num]);
            end
            if isfield(job.output.bias,'warped') && job.output.bias.warped
              out.subj(i).wm{j,1}   = fullfile(pth,mrifolder,['wmr', nam, ext, num]);
            end
            if job.output.surface
              out.subj(i).surface{j,1}   = fullfile(pth,surffolder,['lh.central.'  , nam, ext, num]);
              out.subj(i).thickness{j,1} = fullfile(pth,surffolder,['lh.thickness.', nam, ext, num]);
            end

          % for input
          if iscell(job.subj(i).mov)
            data{j} = job.subj(i).mov{j};
          else
            data{j} = cellstr(job.subj(i).mov);
          end
        end
        
        inputs{1,i} = data;
    end
    
    %%
    maps = {...
      ... 'files','warps', ... internal 
      'wp1','wp2','wp3','wp7',...             % unmodulated (warped==1)
      'mwp1','mwp2','mwp3','mwp7',...         % modulated (mod==1 | mod==3)
      'm0wp1','m0wp2','m0wp3','m0wp7', ...    % modulated (mod==2 | mod==3)
      'rp1a','rp2a','rp3a','rp7a',...         % dartel affine (dartel==1 |?dartel==3)
      'rp1r','rp2r','rp3r','rp7r', ...        % dartel rigid  (dartel==2 |?dartel==3)
      'p0','wm', ...                          % 
      'surface','thickness' ...               % surface
      };
    
    for mi=1:numel(maps)
        for k=1:numel(job.subj)
          if isfield(out.subj(k),maps{mi})
            cdep            = cfg_dep;
            cdep.sname      = sprintf('%s files of Subject %d',maps{mi},k);
            cdep.src_output = substruct('.','subj','()',{k},'.',maps{mi});
            cdep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
            if ~exist('dep','var'); 
                dep = cdep;
            else
                dep = [dep cdep];
            end
          end
        end
    end
    %%
return
%------------------------------------------------------------------------
