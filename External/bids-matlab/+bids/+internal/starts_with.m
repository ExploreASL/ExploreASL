function res = starts_with(str, pattern)
  %
  % Checks id character array 'str' starts with 'pat'
  %
  % USAGE  res = bids.internal.startsWith(str, pat)
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
  res = strcmp(str(1:l_pat), pattern);

end
