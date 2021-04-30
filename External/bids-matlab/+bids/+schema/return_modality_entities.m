function entities = return_modality_entities(suffix_group, schema)

  % for CI
  if iscell(suffix_group)
    suffix_group = suffix_group{1};
  end

  entity_names = fieldnames(suffix_group.entities);

  for i = 1:size(entity_names, 1)
    entities{1, i} = schema.entities.(entity_names{i}).entity; %#ok<*AGROW>
  end

end
