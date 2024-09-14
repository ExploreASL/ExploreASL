function xASL_adm_RemoveLogFilesForRerun(rootDir)
%xASL_adm_RemoveLogFilesForRerun Removes all log files from a directory containing .log files.
%
% FORMAT: xASL_adm_RemoveLogFilesForRerun(rootDir);
%
% INPUT:
%        rootDir            - Case root directory (REQUIRED)
%
% OUTPUT:
%        n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:      Removes all log files from any directory containing .log files.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% EXAMPLE:          rootDir = '.\Test_Runs\TestDataSet';
%                   xASL_adm_RemoveLogFilesForRerun(rootDir);
%
% REFERENCES:       ...
% __________________________________
% Copyright @ 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    %% Input Check
    if nargin < 1 || isempty(rootDir)
        error('Missing directory...');
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
