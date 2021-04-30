function res = ends_with(str, pattern)
  %
  % Checks id character array 'str' ends with 'pat'
  %
  % USAGE  res = bids.internal.endsWith(str, pat)
  %
  % str        - character array
  % pat        - character array
  %
  % __________________________________________________________________________
  %
  % Based on spm_file.m and spm_select.m from SPM12.
  % __________________________________________________________________________

  % Copyright (C) 2011-2018 Guillaume Flandin, Wellcome Centre for Human Neuroimaging
  res = false;
  l_pat = length(pattern);
  if l_pat > length(str)
    return
  end
  res = strcmp(str(end - l_pat + 1:end), pattern);

end
