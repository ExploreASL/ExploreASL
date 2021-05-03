function varargout = file_utils(str, varargin)
  %
  % Character array (or cell array of strings) handling facility
  %
  %
  % To list files or directories (with fullpath if necessary)
  %
  % FORMAT [files, dirs] = bids.internal.file_utils('List',   directory,        regexp)
  % FORMAT [files, dirs] = bids.internal.file_utils('FPList', directory,        regexp)
  % FORMAT [dirs]        = bids.internal.file_utils('List',   directory, 'dir', regexp)
  % FORMAT [dirs]        = bids.internal.file_utils('FPList', directory, 'dir', regexp)
  %
  %
  %
  % To get a certain piece of information from a file.
  %
  % FORMAT str = bids.internal.file_utils(str, option)
  %
  % str        - character array, or cell array of strings
  %
  % option     - string of requested item - one among:
  %              {'path', 'basename', 'ext', 'filename', 'cpath', 'fpath'}
  %
  %
  %
  % To set a certain piece of information from a file.
  %
  % FORMAT str = bids.internal.file_utils(str, opt_key, opt_val,...)
  %
  % str        - character array, or cell array of strings
  %
  % opt_key    - string of targeted item - one among:
  %              {'path', 'basename', 'ext', 'filename', 'prefix', 'suffix'}
  %
  % opt_val    - string of new value for feature
  %
  % __________________________________________________________________________
  %
  % Based on spm_file.m and spm_select.m from SPM12.
  % __________________________________________________________________________

  % Copyright (C) 2011-2018 Guillaume Flandin, Wellcome Centre for Human Neuroimaging

  %#ok<*AGROW>

  if ismember(lower(str), {'list', 'fplist'})
    [varargout{1:nargout}] = listfiles(str, varargin{:});
    return
  end

  needchar = ischar(str);
  options = varargin;

  str = cellstr(str);

  % -Get item
  % ==========================================================================
  if numel(options) == 1
    for n = 1:numel(str)
      [pth, nam, ext] = fileparts(deblank(str{n}));
      switch lower(options{1})
        case 'path'
          str{n} = pth;
        case 'basename'
          str{n} = nam;
        case 'ext'
          str{n} = ext(2:end);
        case 'filename'
          str{n} = [nam ext];
        case 'cpath'
          str(n) = canonicalise_path(str(n));
        case 'fpath'
          str{n} = fileparts(char(canonicalise_path(str(n))));
        otherwise
          error('Unknown option: ''%s''', options{1});
      end
    end
    options = {};
  end

  % -Set item
  % ==========================================================================
  while ~isempty(options)
    for n = 1:numel(str)
      [pth, nam, ext] = fileparts(deblank(str{n}));
      switch lower(options{1})
        case 'path'
          pth = char(options{2});
        case 'basename'
          nam = char(options{2});
        case 'ext'
          ext = char(options{2});
          if ~isempty(ext) && ext(1) ~= '.'
            ext = ['.' ext];
          end
        case 'filename'
          nam = char(options{2});
          ext = '';
        case 'prefix'
          nam = [char(options{2}) nam];
        case 'suffix'
          nam = [nam char(options{2})];
        otherwise
          warning('Unknown item ''%s'': ignored.', lower(options{1}));
      end
      str{n} = fullfile(pth, [nam ext]);
    end
    options([1 2]) = [];
  end

  if needchar
    str = char(str);
  end
  varargout = {str};

end

function t = canonicalise_path(t, d)
  % ==========================================================================
  % -Canonicalise paths to full path names
  % ==========================================================================
  %
  % canonicalise paths to full path names, removing xxx/./yyy and xxx/../yyy
  % constructs
  %
  % - t must be a cell array of (relative or absolute) paths
  %
  % - d must be a single cell containing the base path of relative paths in t
  %

  mch = '^/';
  if ispc % valid absolute paths
    % Allow drive letter or UNC path
    mch = '^([a-zA-Z]:)|(\\\\[^\\]*)';
  end

  if (nargin < 2) || isempty(d)
    d = {pwd};
  end

  % Find partial paths, prepend them with d
  ppsel    = cellfun(@isempty, regexp(t, mch, 'once'));
  t(ppsel) = cellfun(@(t1)fullfile(d{1}, t1), t(ppsel), 'UniformOutput', false);

  % Break paths into cell lists of folder names
  pt = pathparts(t);

  % Remove single '.' folder names
  sd = cellfun(@(pt1)strcmp(pt1, '.'), pt, 'UniformOutput', false);
  for cp = 1:numel(pt)
    pt{cp} = pt{cp}(~sd{cp});
  end

  % Go up one level for '..' folders, don't remove drive letter/server name from PC path
  ptstart = 1;

  for cp = 1:numel(pt)

    tmppt = {};

    for cdir = ptstart:numel(pt{cp})
      if strcmp(pt{cp}{cdir}, '..')
        tmppt = tmppt(1:end - 1);
      else
        tmppt{end + 1} = pt{cp}{cdir};
      end
    end

    pt{cp} = tmppt;

  end

  % Assemble paths
  t = cellfun(@(pt1)fullfile(filesep, pt1{:}), pt, 'UniformOutput', false);
  if ispc
    t = cellfun(@(pt1)fullfile(pt1{:}), pt, 'UniformOutput', false);
  end

end

function pp = pathparts(p)
  % ==========================================================================
  % -Parse paths
  % ==========================================================================
  %
  % parse paths in cellstr p
  %
  % returns cell array of path component cellstr arrays
  %
  % For PC (WIN) targets, both '\' and '/' are accepted as filesep, similar
  % to MATLAB fileparts
  %

  file_separator = filesep;
  if ispc
    file_separator = '\\/';
  end

  pp = cellfun(@(p1)textscan(p1, '%s', ...
                             'delimiter', file_separator, ...
                             'MultipleDelimsAsOne', 1), p);

  if ispc
    for k = 1:numel(pp)
      if ~isempty(regexp(pp{k}{1}, '^[a-zA-Z]:$', 'once'))
        pp{k}{1} = strcat(pp{k}{1}, filesep);
      elseif ~isempty(regexp(p{k}, '^\\\\', 'once'))
        pp{k}{1} = strcat(filesep, filesep, pp{k}{1});
      end
    end
  end

end

function [files, dirs] = listfiles(action, directory, varargin)
  % ==========================================================================
  % -List files and directories
  % ==========================================================================
  %
  % FORMAT [files, dirs] = listfiles('List',   directory,        regexp)
  % FORMAT [files, dirs] = listfiles('FPList', directory,        regexp)
  % FORMAT [dirs]        = listfiles('List',   directory, 'dir', regexp)
  % FORMAT [dirs]        = listfiles('FPList', directory, 'dir', regexp)
  %

  files = '';
  dirs = '';

  switch lower(action)
    case 'list'
      fp = false;
    case 'fplist'
      fp = true;
    otherwise
      error('Invalid action: ''%s''.', action);
  end

  directory = bids.internal.file_utils(directory, 'cpath');
  if nargin < 2
    directory = pwd;
  end

  dd = dir(directory);
  if isempty(dd)
    return
  end

  % set if we work on directory or files
  % set regular expression to use
  dirmode = false;
  expr = '.*';

  if nargin >= 3
    expr = varargin{1};
  end
  if nargin == 3 && strcmpi(varargin{1}, 'dir')
    dirmode = true;
  end
  if nargin >= 3 && strcmpi(varargin{1}, 'dir')
    dirmode = true;
    expr = varargin{2};
  end

  files = sort({dd(~[dd.isdir]).name})';
  dirs = sort({dd([dd.isdir]).name})';
  dirs = setdiff(dirs, {'.', '..'});

  if dirmode

    t = regexp(dirs, expr);

    if numel(dirs) == 1 && ~iscell(t)
      t = {t};
    end
    dirs = dirs(~cellfun(@isempty, t));
    files = dirs;

  else

    t = regexp(files, expr);

    if numel(files) == 1 && ~iscell(t)
      t = {t};
    end
    files = files(~cellfun(@isempty, t));

  end

  if fp
    files = cellfun(@(x) fullfile(directory, x), files, 'UniformOutput', false);
    dirs = cellfun(@(x) fullfile(directory, x), dirs, 'UniformOutput', false);
  end

  files = char(files);
  dirs = char(dirs);

end
