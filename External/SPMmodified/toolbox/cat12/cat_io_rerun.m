function run = cat_io_rerun(files,filedates)
%cat_io_rerun. Test if a file is newer than another file.  
% This function is used to estimated if a file is newer than another given 
% file or date. For instance file is the result of anther file that was 
% changed in the meantime, it has to be reprocessed. 
%
%  run = cat_io_rerun(files,filedates)
% 
%  run      .. logical vector with the number of given files
%              cell if directories or wildcards are used
%  files    .. filenames (cellstr or char)
%  filedat  .. filenames (cellstr or char) or datetimes or datenum
%
% Examples: 
%  1) Is the working directory younger than the SPM dir?
%     cat_io_lazy(pwd,spm('dir'); 
%
%  2) Is the working directory younger than one month?
%     cat_io_lazy(pwd,clock - [0 1 0 0 0 0]) 
%   
%  3) Is this function younger than one year?
%     cat_io_lazy(which('cat_io_lazy'),clock - [1 0 0 0 0 0]) 
%
% _________________________________________________________________________
% Robert Dahnke
% $Id: cat_io_rerun.m 1578 2020-03-10 09:19:19Z gaser $

  files = cellstr(files);
  if iscellstr(filedates) || ischar(filedates)
    filedates = cellstr(filedates);
    if numel(filedates) == 1
      filedates = repmat(filedates,numel(files),1);
    else
      if ~isempty(filedates) && numel(files) ~= numel(filedates)
        error('ERROR:cat_io_lazy:inputsize','Number of files and filedates has to be equal.\n')
      end
    end  
  else 
    if size(filedates,1)
      filedates = repmat(filedates,numel(files),1);
    end
  end
  
  run = ones(size(files)); 
  for fi = 1:numel(files)
    if ~exist(files{fi},'file')
      run(fi) = 1; 
    else 
      fdata = dir(files{fi});
      if numel(fdata)>1
        run = num2cell(run); 
      end
      if exist('filedates','var') && iscellstr(filedates) && exist(filedates{fi},'file')
        fdata2 = dir(filedates{fi});
        if numel(fdata)>1
          run{fi} = [fdata(:).datenum] < fdata2.datenum;
        else
          run(fi) = fdata.datenum < fdata2.datenum;
        end
      elseif ~isempty(filedates) 
        if numel(fdata)>1
          run{fi} = [fdata(:).datenum] < datenum( filedates(fi,:) );
        else
          run(fi) = fdata.datenum < datenum( filedates(fi,:) );
        end
      end
    end
  end
end
