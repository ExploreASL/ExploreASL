function xASL_bids_MergeNifti_RenameParms(Fpath,Fname)
%xASL_bids_MergeNifti_RenameParms Find *_parms.m files in directory and shorten to provided name

FileList = xASL_adm_GetFileList(Fpath, '^.*_parms\.mat$', 'List', [], false);

if ~isempty(FileList)
	xASL_Move(fullfile(Fpath, FileList{1}), fullfile(Fpath,[Fname '_parms.mat']));
end


end
