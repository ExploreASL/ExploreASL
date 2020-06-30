function cat_io_rmdir(dirs)
% _________________________________________________________________________
% Recusively remove empty subdirectories, whereas the original MATLAB rmdir  
% removes only the deepest empty subdirectory or all subdirectories incl.
% all files ('s' option). 
%
% cat_io_rmdir(dirs)
%  
% dirs = char or cellstr of directories
% _________________________________________________________________________
% Robert Dahnke
% $Id: cat_io_rmdir.m 1371 2018-10-04 09:21:11Z gaser $

  dirs    = cellstr(dirs); 
  
  subdirs = cell(0);
  for s=1:numel(dirs)
    subdirs = [subdirs;cat_vol_findfiles(dirs{s},'*',struct('dirs',1))]; %#ok<AGROW>
  end
  subdirs = unique(subdirs);
  clear dirs;
  
  % sort subdirs by length
  dsubdirs = cellfun('length',subdirs);
  [tmp,si] = sort(dsubdirs,'descend'); %#ok<ASGLU>
  subdirs = subdirs(si);
  clear dsubdirs tmp si;
  
  % remove empty dirs step by step
  for si=1:numel(subdirs)
    if exist(subdirs{si},'dir')
      try  %#ok<TRYNC>
        rmdir(subdirs{si}); 
      end
    end
  end
  
end