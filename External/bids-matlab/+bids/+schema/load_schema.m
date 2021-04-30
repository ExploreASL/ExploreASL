function schema = load_schema(use_schema)
  % Loads a json schema by recursively looking through a folder structure.
  %
  % The nesting of the output structure reflects a combination of the folder structure and
  % any eventual nesting within each json.
  %
  %
  % Copyright (C) 2021--, BIDS-MATLAB developers

  % TODO:
  %  - folders that do not contain json files themselves but contain
  %  subfolders that do, are not reflected in the output structure (they are
  %  skipped). This can lead to "name conflicts". See "silenced" unit tests
  %  for more info.

  if nargin < 1
    use_schema = true();
  end

  if ~use_schema
    schema = [];
    return
  end

  if ischar(use_schema)
    schema_dir = use_schema;
  else
    schema_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'schema');
  end

  if ~exist(schema_dir, 'dir')
    error('The schema directory %s does not exist.', schema_dir);
  end

  schema = struct();

  [json_file_list, dirs] = bids.internal.file_utils('FPList', schema_dir, '^.*.json$');

  schema = append_json_content_to_structure(schema, json_file_list);

  schema = inspect_subdir(schema, dirs);

end

function structure = append_json_content_to_structure(structure, json_file_list)

  for iFile = 1:size(json_file_list, 1)

    file = deblank(json_file_list(iFile, :));

    field_name = bids.internal.file_utils(file, 'basename');

    structure.(field_name) = bids.util.jsondecode(file);
  end

end

function structure = inspect_subdir(structure, subdir_list)
  % recursively inspects subdirectory for json files and reflects folder
  % hierarchy in the output structure.

  for iDir = 1:size(subdir_list, 1)

    directory = deblank(subdir_list(iDir, :));

    [json_file_list, dirs] = bids.internal.file_utils('FPList', directory, '^.*.json$');

    if ~isempty(json_file_list)
      field_name = bids.internal.file_utils(directory, 'basename');
      structure.(field_name) = struct();
      structure.(field_name) = append_json_content_to_structure(structure.(field_name), ...
                                                                json_file_list);
    end

    structure = inspect_subdir(structure, dirs);

  end

end
