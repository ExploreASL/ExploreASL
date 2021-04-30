function file_info = return_file_info(BIDS, fullpath_filename)

  file_info.path = bids.internal.file_utils(fullpath_filename, 'path');
  file_info.modality = bids.internal.file_utils(file_info.path, 'filename');
  file_info.filename = bids.internal.file_utils(fullpath_filename, 'filename');

  file_info.sub_idx = bids.internal.return_subject_index(BIDS, file_info.filename);
  file_info.file_idx = bids.internal.return_file_index(BIDS, ...
                                                       file_info.modality, ...
                                                       file_info.filename);

end
