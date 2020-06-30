function out = cat_long_multi_run(job)
% Call cat_long_main for multiple subjects
%
% Christian Gaser
% $Id: cat_long_multi_run1173.m 1389 2018-11-11 10:39:41Z dahnke $

global opts extopts output modulate dartel warps

warning off;

% use some options from GUI or default file
opts     = job.opts;
extopts  = job.extopts;
output   = job.output;
modulate = job.modulate;
dartel   = job.dartel;
warps    = job.warps;

jobs = repmat({'cat_long_main.m'}, 1, numel(job.subj));
inputs = cell(1, numel(job.subj));

if cat_get_defaults('extopts.subfolders')
  mrifolder = 'mri';
else
  mrifolder = '';
end

for i=1:numel(job.subj),
    if warps
      out(i).warps = cell(1,1);
      [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{1});
      out(i).warps{1} = fullfile(pth,mrifolder,['y_avg_', nam, ext, num]);
    end

    out(i).files = cell(numel(job.subj(i).mov),1);
    m = numel(job.subj(i).mov);
    data = cell(m,1);
    for j=1:m
        [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{j});
        switch modulate
        case 0
          out(i).files{j} = fullfile(pth,mrifolder,['wp1r', nam, ext, num]);
        case 1
          out(i).files{j} = fullfile(pth,mrifolder,['mwp1r', nam, ext, num]);
        case 2
          out(i).files{j} = fullfile(pth,mrifolder,['m0wp1r', nam, ext, num]);
        end
        data{j} = job.subj(i).mov{j};
    end
    inputs{1,i} = data;
end;

spm_jobman('run',jobs,inputs{:});

warning on;
