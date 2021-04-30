function modality_groups = return_modality_groups(schema)

  % dummy variable if we go schema less
  modality_groups = {nan()};
  if ~isempty(schema)
    modality_groups = fieldnames(schema.modalities);
  end

end
