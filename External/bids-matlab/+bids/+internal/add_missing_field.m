function structure = add_missing_field(structure, field)
  if ~isfield(structure, field)
    structure(1).(field) = '';
  end
end
