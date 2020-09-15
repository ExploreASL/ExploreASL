function out = cat_vol_avg(job)
% Calculate average of images of same dimension and voxel size
% This function also works for 5D images (e.g. deformations).
%
% FORMAT out = cat_vol_avg(job)
% job.data      - data array of input files
% job.output    - output name
% job.outdir    - output directory
%
% out           - output name
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_vol_avg.m 1426 2019-02-04 11:51:40Z gaser $
%

[p,nam,ext] = spm_fileparts(job.output);
if isempty(p)
    if isempty(job.outdir{1})
        p = spm_fileparts(job.data{1});
    else
        p = job.outdir{1};
    end
end
if isempty(nam)
    nam = ['avg_' spm_file(job.data{1},'basename')];
    ext = ['.' spm_file(job.data{1},'ext')];
end
if isempty(ext)
    ext = spm_file_ext;
end
out.files = { fullfile(p,[nam ext]) };

N = nifti(char(job.data));

if length(N)>1 && any(any(diff(cat(1,N.dat.dim),1,1),1))
	error('images don''t all have same dimensions')
end
if max(max(max(abs(diff(cat(3,N.mat),1,3))))) > 1e-8
	error('images don''t all have same orientation & voxel size')
end

d = N(1).dat.dim;

% extend dimensions if necessary
d = [d ones(1, 5-length(d))];

avg = zeros(d);

for i = 1:length(N)
  for j = 1:d(4)
    for k = 1:d(5)
      avg(:,:,:,j,k) = avg(:,:,:,j,k) + N(i).dat(:,:,:,j,k);
    end
  end
end

avg = avg/length(N);

Nout = N(1);
Nout.dat.fname = out.files{1};
create(Nout);
Nout.dat(:,:,:,:,:) = avg;

