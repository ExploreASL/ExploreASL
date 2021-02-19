function xASL_adm_RemoveLogFilesForRerun(rootDir)
%xASL_adm_RemoveLogFilesForRerun Removes all log files from xASL directory
%
% FORMAT: xASL_adm_RemoveLogFilesForRerun(rootDir);
%
% INPUT:
%        rootDir            - Case root directory
%
% OUTPUT:
%        n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Removes all log files from xASL directory.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          rootDir = '.\Test_Runs\TestDataSet';
%                   xASL_adm_RemoveLogFilesForRerun(rootDir);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2021 ExploreASL

    %% Input Check
    if nargin < 1 || isempty(rootDir)
        error('Missing xASL directory...');
    end
    
    
    %% Get all log files within root directory
    fileList = xASL_adm_GetFileList(rootDir, '^.+\.log$', 'FPListRec');
    
    %% Iterate over log files
    for iFile=1:numel(fileList)
        % Get current file path
        curFile = fileList{iFile};
        % Delete file
        fprintf('Delete: %s\n',curFile);
        xASL_delete(curFile)
    end

end







