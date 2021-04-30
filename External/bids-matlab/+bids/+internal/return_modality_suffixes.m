function suffixes = return_modality_suffixes(modality)

  suffixes = '_(';

  % For CI
  if iscell(modality)
    modality = modality{1};
  end

  for iExt = 1:numel(modality(:).suffixes)
    suffixes = [suffixes,  modality.suffixes{iExt}, '|']; %#ok<AGROW>
  end

  % Replace final "|" by a "){1}"
  suffixes(end:end + 3) = '){1}';

end
