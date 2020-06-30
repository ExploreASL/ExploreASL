function out = cat_long_multi_run(job)
% Call cat_long_main for multiple subjects
%
% Christian Gaser
% $Id: cat_long_multi_run.m 1515 2019-11-09 23:01:28Z gaser $

global opts extopts output modulate dartel delete_temp ROImenu surfaces

warning off;

if isdeployed, job.nproc = 0; end

% use some options from GUI or default file
opts        = job.opts;
extopts     = job.extopts;
output      = job.output;
modulate    = job.modulate;
dartel      = job.dartel;
ROImenu     = job.ROImenu;
surfaces    = job.output.surface;

if isfield(job,'delete_temp')  
  delete_temp = job.delete_temp;
else
  delete_temp = 1;
end

jobs = repmat({'cat_long_main.m'}, 1, numel(job.subj));
inputs = cell(1, numel(job.subj));

if cat_get_defaults('extopts.subfolders')
  mrifolder    = 'mri';
  surffolder   = 'surf';
  labelfolder  = 'label';
  reportfolder = 'report';
else
  mrifolder    = '';
  surffolder   = ''; 
  labelfolder  = '';
  reportfolder = ''; 
end

out.surf = cell(''); out.thick = cell(''); out.mwp1 = cell('');
out.catreport = cell(''); out.catroi = cell('');
for i=1:numel(job.subj)
		out.sess(i).warps = cell(1,1);
		[pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{1});
		out.sess(i).warps{1} = fullfile(pth,mrifolder,['avg_y_', nam, ext, num]);

    out.sess(i).files = cell(numel(job.subj(i).mov),1);
    m = numel(job.subj(i).mov);
    data = cell(m,1);
    for j=1:m
        [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{j});
        switch modulate
        case 0
          out.sess(i).files{j} = fullfile(pth,mrifolder,['wp1r', nam, ext, num]);
        case 1
          out.sess(i).files{j} = fullfile(pth,mrifolder,['mwp1r', nam, ext, num]);
        case 2
          out.sess(i).files{j} = fullfile(pth,mrifolder,['m0wp1r', nam, ext, num]);
        end
        data{j} = job.subj(i).mov{j};
    
        out.mwp1      = [out.mwp1       fullfile(pth,mrifolder   ,['mwp1r'          nam ext num])]; 
        out.surf      = [out.surf       fullfile(pth,surffolder  ,['lh.central.r'   nam '.gii'])]; 
        out.thick     = [out.thick      fullfile(pth,surffolder  ,['lh.thickness.r' nam])]; 
        out.catreport = [out.catreport  fullfile(pth,reportfolder,['cat_r'          nam '.xml'])]; 
        out.catroi    = [out.catroi     fullfile(pth,labelfolder ,['catROI_r'       nam '.xml'])]; 
        
    end
    inputs{1,i} = data;
    
end


% split job and data into separate processes to save computation time
if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
  if nargout==1
    varargout{1} = cat_parallelize(job,mfilename,'subj');
  else
    cat_parallelize(job,mfilename,'subj');
  end
  return
else
  % inputs is used twice because the first one is only used for command line display
  spm_jobman('run',jobs,inputs{:}); 
end
warning on;
