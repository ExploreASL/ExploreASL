function metalist = get_meta_list(filename, pattern)
  %
  % Read a BIDS's file metadata according to the inheritance principle
  %
  % USAGE::
  %
  %    meta = bids.internal.get_metadata(filename, pattern = '^.*%s\\.json$')
  %
  % :param filename: fullpath name of file following BIDS standard
  % :type  filename: string
  % :param pattern:  Regular expression matching the metadata file (default is ``'^.*%s\\.json$'``)
  %                  If provided, it must at least be ``'%s'``.
  % :type  pattern:  string
  %
  %
  % metalist    - list of paths to metafiles
  % __________________________________________________________________________

  % Copyright (C) 2016-2018, Guillaume Flandin, Wellcome Centre for Human Neuroimaging
  % Copyright (C) 2018--, BIDS-MATLAB developers

  if nargin == 1
    pattern = '^.*%s\\.json$';
  end

  pth = fileparts(filename);
  p = bids.internal.parse_filename(filename);
  metalist = {};

  % Default assumes we are dealing with a file in the root directory
  % like "participants.tsv"
  % If the file has underscore separated entities ("sub-01_T1w.nii")
  % then we look through the hierarchy for potential metadata file associated
  % with queried file.
  N = 1;
  if isfield(p, 'entities')
    N = 3;
    % -There is a session level in the hierarchy
    if isfield(p.entities, 'ses') && ~isempty(p.entities.ses)
      N = N + 1;
    end
  end

  for n = 1:N

    % List the potential metadata files associated with this file suffix type
    % Default is to assume it is a JSON file
    metafile = bids.internal.file_utils('FPList', pth, sprintf(pattern, p.suffix));

    if isempty(metafile)
      metafile = {};
    else
      metafile = cellstr(metafile);
    end

    % For all those files we find which one is potentially associated with
    % the file of interest
    for i = 1:numel(metafile)

      p2 = bids.internal.parse_filename(metafile{i});
      entities = {};
      if isfield(p2, 'entities')
        entities = fieldnames(p2.entities);
      end

      % Check if this metadata file contains
      %   - the same entity-label pairs
      %   - same suffix
      % as its data file counterpart
      ismeta = true;
      if ~strcmp(p.suffix, p2.suffix)
        ismeta = false;
      end
      for j = 1:numel(entities)
        if ~isfield(p.entities, entities{j}) || ...
                ~strcmp(p.entities.(entities{j}), p2.entities.(entities{j}))
          ismeta = false;
          break
        end
      end

      % append path to list
      if ismeta
        metalist{end + 1, 1} = metafile{i}; %#ok<AGROW>
      end

    end

    % Go up to the parent folder
    pth = fullfile(pth, '..');
  end
end
