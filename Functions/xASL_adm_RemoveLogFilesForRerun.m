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
    if nargin < 1
        error('Missing xASL directory...');
    end
    
    
    %% Get all log files within root directory
    fileList = dir(fullfile(rootDir, '**\*.log'));  % Get list of log files and folders in any subfolder
    fileList = fileList(~[fileList.isdir]);         % Remove folders from list
    
    
    %% Iterate over log files
    for iFile=1:numel(fileList)
        % Get current file path
        curFile = fullfile(fileList(iFile).folder,fileList(iFile).name);
        % Delete file
        fprintf('Delete: %s\n',fileList(iFile).name);
        delete(curFile)
    end

end







