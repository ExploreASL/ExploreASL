%-----------------------------------------------------------------------
% Job for longitudinal batch
% Christian Gaser
% $Id: cat_long_main1445.m 1470 2019-05-29 11:49:16Z gaser $
%-----------------------------------------------------------------------

global opts extopts output modulate dartel delete_temp ROImenu surfaces cat

write_CSF = cat_get_defaults('output.CSF.mod') > 0;

if write_CSF
  cat12('expert')
end

warning('off','MATLAB:DELETE:FileNotFound');
matlabbatch{1}.spm.tools.cat.tools.series.data = '<UNDEFINED>';

if exist('opts','var')
	matlabbatch{2}.spm.tools.cat.estwrite.opts = opts;
end
if exist('extopts','var')
	matlabbatch{2}.spm.tools.cat.estwrite.extopts = extopts;
end
if exist('output','var')
	matlabbatch{2}.spm.tools.cat.estwrite.output = output;
end
if exist('ROImenu','var')
  matlabbatch{2}.spm.tools.cat.estwrite.output.ROImenu = ROImenu;
end

% modulation option for applying deformations
if modulate
  matlabbatch{4}.spm.tools.cat.tools.defs.modulate = modulate;
end

% surface estimation
if surfaces
  matlabbatch{2}.spm.tools.cat.estwrite.output.surface = 1;
else
  matlabbatch{2}.spm.tools.cat.estwrite.output.surface = 0;
end

matlabbatch{2}.spm.tools.cat.estwrite.nproc = 0;

% also write CSF?
if write_CSF
  matlabbatch{2}.spm.tools.cat.estwrite.output.CSF.native = 1;
end

% dartel export option
matlabbatch{2}.spm.tools.cat.estwrite.output.GM.dartel = dartel;
matlabbatch{2}.spm.tools.cat.estwrite.output.WM.dartel = dartel;

% longitudinal rigid registration with final masking
matlabbatch{1}.spm.tools.cat.tools.series.bparam = 1e6;
matlabbatch{1}.spm.tools.cat.tools.series.use_brainmask = 1;
matlabbatch{1}.spm.tools.cat.tools.series.reduce = 1;

% cat12 segmentation of realigned images 
matlabbatch{2}.spm.tools.cat.estwrite.data(1) = cfg_dep('Longitudinal Rigid Registration: Realigned images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rimg', '()',{':'}));
matlabbatch{2}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{2}.spm.tools.cat.estwrite.output.GM.mod = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{2}.spm.tools.cat.estwrite.output.WM.mod = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.bias.warped = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.warps = [1 0];

% averaging deformations
matlabbatch{3}.spm.tools.cat.tools.avg_img.data(1) = cfg_dep('CAT12: Segmentation: Deformation Field', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{3}.spm.tools.cat.tools.avg_img.output = '';
matlabbatch{3}.spm.tools.cat.tools.avg_img.outdir = {''};

% applying deformations to native segmentations
matlabbatch{4}.spm.tools.cat.tools.defs.field1(1) = cfg_dep('Image Average: Average Image: ', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{4}.spm.tools.cat.tools.defs.images(1) = cfg_dep('CAT12: Segmentation: p1 Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','p', '()',{':'}));
matlabbatch{4}.spm.tools.cat.tools.defs.images(2) = cfg_dep('CAT12: Segmentation: p2 Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','p', '()',{':'}));
if write_CSF
  matlabbatch{4}.spm.tools.cat.tools.defs.images(3) = cfg_dep('CAT12: Segmentation: p3 Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','p', '()',{':'}));
end
matlabbatch{4}.spm.tools.cat.tools.defs.interp = 1;

% applying deformations to average T1 image
matlabbatch{5}.spm.tools.cat.tools.defs.field1(1) = cfg_dep('Image Average: Average Image: ', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{5}.spm.tools.cat.tools.defs.images(1) = cfg_dep('Longitudinal Registration: Midpoint Average', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avg', '()',{':'}));
matlabbatch{5}.spm.tools.cat.tools.defs.interp = 1;
matlabbatch{5}.spm.tools.cat.tools.defs.modulate = 0;

% delete temporary files
if delete_temp
  matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('CAT12: Segmentation: p1 Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','p', '()',{':'}));
  matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('CAT12: Segmentation: p2 Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','p', '()',{':'}));
  matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('CAT12: Segmentation: Deformation Field', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
  if write_CSF
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(4) = cfg_dep('CAT12: Segmentation: p3 Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','p', '()',{':'}));
  end

  matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
end


