function out = cat_long_multi_run(job)
% Call cat_long_main for multiple subjects
%
% Christian Gaser
% $Id: cat_long_multi_run1445.m 1470 2019-05-29 11:49:16Z gaser $

global opts extopts output modulate dartel delete_temp ROImenu surfaces

warning off;

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
  mrifolder = 'mri';
else
  mrifolder = '';
end

for i=1:numel(job.subj),
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
    end
    inputs{1,i} = data;
end;


% split job and data into separate processes to save computation time
if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
  if nargout==1
    varargout{1} = cat_parallelize(job,mfilename,'subj');
  else
    cat_parallelize(job,mfilename,'subj');
  end
  return
else
  spm_jobman('run',jobs,inputs{:});
end
warning on;
