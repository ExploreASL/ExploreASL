function tsvwrite(f, var)
  % Save text and numeric data to .tsv file
  % FORMAT tsvwrite(f, var)
  % f     - filename
  % var   - data array or structure
  %
  % Adapted from spm_save.m
  % __________________________________________________________________________
  % Copyright (C) 2018, Guillaume Flandin, Wellcome Centre for Human Neuroimaging
  % Copyright (C) 2018--, BIDS-MATLAB developers

  delim = sprintf('\t');

  % If the input is a MATLAB table, the built-in functionality is used
  % otherwise export is performed manually
  if isstruct(var) || iscell(var) || isnumeric(var) || islogical(var)

    % Convert input to a common format we will use for writing
    % var will be a cell array where the first row is the header and all
    % following rows contains the values to write
    if isstruct(var)

      fn = fieldnames(var);
      var = struct2cell(var)';

      for i = 1:numel(var)

        if ~ischar(var{i})
          var{i} = var{i}(:);
        end

        if ~iscell(var{i})
          var{i} = cellstr(num2str(var{i}, 16));
          var{i} = strtrim(var{i});
          var{i}(cellfun(@(x) strcmp(x, 'NaN'), var{i})) = {'n/a'};
        end

      end

      var = [fn'; var];

    elseif iscell(var) || isnumeric(var) || islogical(var)

      if isnumeric(var) || islogical(var)
        var = num2cell(var);
      end

      var = cellfun(@(x) num2str(x, 16), var, 'UniformOutput', false);
      var = strtrim(var);
      var(cellfun(@(x) strcmp(x, 'NaN'), var)) = {'n/a'};

    end

    % Actually write to file
    fid = fopen(f, 'Wt');

    if fid == -1
      error('Unble to write file %s.', f);
    end

    for i = 1:size(var, 1)

      for j = 1:size(var, 2)

        to_print = var{i, j};

        if iscell(to_print)
          to_print = to_print{1};
        end

        if isempty(to_print)
          to_print = 'n/a';

        elseif any(to_print == delim)
          to_print = ['"' to_print '"'];

        end

        fprintf(fid, '%s', to_print);

        if j < size(var, 2)
          fprintf(fid, delim);
        end

      end

      fprintf(fid, '\n');
    end

    fclose(fid);

  elseif isa(var, 'table')
    writetable(var, f, ...
               'FileType', 'text', ...
               'Delimiter', delim);

  else
    error('Unknown data type.');

  end
