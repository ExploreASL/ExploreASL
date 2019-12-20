function EPAD_CopyDirsTrackProgress(InputPath, OutputPath)
%EPAD_CopyDirsTrackProgress Wrapper of xASL_Copy, that allows tracking progress
% for copying large datasets

fprintf('%s\n', ['Copying ' InputPath ' to ' OutputPath]);
fprintf('%s', 'Progress:  0%');

% Get list of subdirectories
DirList = xASL_adm_GetFileList(InputPath, '.*', 'FPListRec', [0 Inf], 1);
DirList = DirList(end:-1:1); % start with subdirectories

xASL_adm_CreateDir(OutputPath);

for iD=1:length(DirList)
    xASL_TrackProgress(iD, length(DirList)+1); % track progress
    DirIn = DirList{iD};
    DirOut = fullfile(OutputPath, DirList{iD}(length(InputPath)+1:end));
    xASL_adm_CreateDir(DirOut);
    xASL_Copy(DirIn, fileparts(DirOut), true, false);
end
xASL_Copy(InputPath, fileparts(OutputPath), true, false);
xASL_TrackProgress(1, 1); % track progress
fprintf('%s\n', 'Done');


end

