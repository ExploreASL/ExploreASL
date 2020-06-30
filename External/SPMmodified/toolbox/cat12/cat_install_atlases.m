function cat_install_atlases
% Convert CAT12 atlas files (csv) and add Dartel atlas labels to spm12 
% atlas folder (xml)
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_install_atlases.m 1571 2020-02-27 14:51:11Z gaser $

spm_dir = spm('dir');
atlas_dir = fullfile(spm_dir,'atlas');
[ST, RS] = mkdir(atlas_dir);

if ST
  [csv_files, n] = cat_vol_findfiles(fullfile(spm_dir,'toolbox','cat12','templates_volumes'), '*.csv');
  for i = 1:n
    csv_file = deblank(csv_files{i});
    csv = cat_io_csv(csv_file,'','',struct('delimiter',';'));
    [pth,nam] = spm_fileparts(csv_file);
    xml_file = fullfile(atlas_dir, ['labels_cat12_' nam '.xml']);
    old_xml_file = fullfile(atlas_dir, ['labels_dartel_' nam '.xml']);
    create_spm_atlas_xml(xml_file, csv);
    atlas_file = fullfile(pth,[nam '.nii']);
    new_atlas_name = ['cat12_' nam '.nii'];
    try
      copyfile(atlas_file,fullfile(atlas_dir,new_atlas_name));
      delete(old_xml_file);
      fprintf('Install %s\n',xml_file);
    catch
      disp('Writing error: Please check file permissions.');
    end
  end
else
  error(RS);
end

% this is maybe not enough, to update the file in SPM functions
% you may need to remove the old files and finish SPM, update and restart SPM 
spm_atlas('list','installed','-refresh');

fprintf('\nUse atlas function in SPM Results or context menu in orthogonal view (via right mouse button): Display|Labels\n');

function create_spm_atlas_xml(fname,csvx,opt)
% create an spm12 compatible xml version of the csv data
if ~exist('opt','var'), opt = struct(); end

[pp,ff,ee] = spm_fileparts(fname); 

% remove prepending name part
if ~isempty(strfind(ff,'labels_cat12_'))
  ff = ff(length('labels_cat12_')+1:end);
end

def.name   = ff;
def.desc   = '';
def.url    = '';
def.ver    = cat_version;
def.lic    = 'CC BY-NC';
def.cor    = 'MNI DARTEL'; 
def.type   = 'Label';
def.images = ['cat12_' ff '.nii'];

opt = cat_io_checkinopt(opt,def);

xml.header = [...
  '<?xml version="1.0" encoding="ISO-8859-1"?>\n' ...
  '<!-- This is a temporary file format -->\n' ...
  '  <atlas version="' opt.ver '">\n' ...
  '    <header>\n' ...
  '      <name>' opt.name '</name>\n' ...
  '      <version>' opt.ver '</version>\n' ...
  '      <description>' opt.desc '</description>\n' ...
  '      <url>' opt. url '</url>\n' ...
  '      <licence>' opt.lic '</licence>\n' ...
  '      <coordinate_system>' opt.cor '</coordinate_system>\n' ...
  '      <type>' opt.type '</type>\n' ...
  '      <images>\n' ...
  '        <imagefile>' opt.images '</imagefile>\n' ...
  '      </images>\n' ...
  '    </header>\n' ...
  '  <data>\n' ...
  ];
xml.data = '';

for di = 2:size(csvx,1);
  % index      = label id
  % name       = long name SIDE STRUCTURE TISSUE 
  % short_name = short name 
  % RGBA       = RGB color
  % XYZmm      = XYZ coordinate
  if size(csvx,2) < 8
    xml.data = sprintf('%s%s\n',xml.data,sprintf(['    <label><index>%d</index>'...
      '<short_name>%s</short_name><name>%s</name>' ...
      '<RGBA></RGBA><XYZmm></XYZmm></label>'],...
      csvx{di,1}, csvx{di,2},csvx{di,3}));
  else
    xml.data = sprintf('%s%s\n',xml.data,sprintf(['    <label><index>%d</index>'...
      '<short_name>%s</short_name><name>%s</name>' ...
      '<RGBA></RGBA><XYZmm>%s</XYZmm></label>'],...
      csvx{di,1}, csvx{di,2},csvx{di,3},csvx{di,8}));
  end
end

xml.footer = [ ...
  '  </data>\n' ...
  '</atlas>\n' ...
  ];

fid = fopen(fname,'w');
if fid >= 0
  fprintf(fid,[xml.header,xml.data,xml.footer]);
  fclose(fid);
else
  fprintf('Error while writing %s. Check file permissions.\n',fname);
end

