function meta = get_metadata(filename, pattern)
  % Read a BIDS's file metadata according to the inheritance principle
  % FORMAT meta = bids.internal.get_metadata(filename, pattern)
  % filename    - name of file following BIDS standard
  % pattern     - regular expression matching metadata file
  % meta        - metadata structure
  % __________________________________________________________________________

  % Copyright (C) 2016-2018, Guillaume Flandin, Wellcome Centre for Human Neuroimaging
  % Copyright (C) 2018--, BIDS-MATLAB developers

  if nargin == 1
    pattern = '^.*%s\\.json$';
  end

  pth = fileparts(filename);
  p = bids.internal.parse_filename(filename);

  meta = struct();

  N = 3;

  % -There is a session level in the hierarchy
  if isfield(p, 'ses') && ~isempty(p.ses)
    N = N + 1;
  end

  % -Loop from the directory where the file of interest is back to the
  % top level of the BIDS hierarchy
  for n = 1:N

    % -List the potential metadata files associated with this file suffix type
    % Default is to assume it is a JSON file
    metafile = bids.internal.file_utils('FPList', pth, sprintf(pattern, p.type));

    if isempty(metafile)
      metafile = {};
    else
      metafile = cellstr(metafile);
    end

    % -For all those files we find which one is potentially associated with
    % the file of interest
    for i = 1:numel(metafile)

      p2 = bids.internal.parse_filename(metafile{i});
      fn = setdiff(fieldnames(p2), {'filename', 'ext', 'type'});

      % -Check if this metadata file contains the same entity-label pairs as its
      % data file counterpart
      ismeta = true;
      for j = 1:numel(fn)
        if ~isfield(p, fn{j}) || ~strcmp(p.(fn{j}), p2.(fn{j}))
          ismeta = false;
          break
        end
      end

      % -Read the content of the metadata file if it is a JSON file and update
      % the metadata concerning the file of interest otherwise store the filename
      if ismeta
        if strcmp(p2.ext, '.json')
          meta = update_metadata(meta, bids.util.jsondecode(metafile{i}), metafile{i});
        else
          meta.filename = metafile{i};
        end
      end

    end

    % -Go up to the parent folder
    pth = fullfile(pth, '..');

  end

  % ==========================================================================
  % -Inheritance principle
  % ==========================================================================
function s1 = update_metadata(s1, s2, file)
  if isempty(s2)
    return
  elseif ~isstruct(s2)
    error('Metadata file contents were neither struct nor empty. File: %s', file);
  end
  fn = fieldnames(s2);
  for i = 1:numel(fn)
    if ~isfield(s1, fn{i})
      s1.(fn{i}) = s2.(fn{i});
    end
  end
