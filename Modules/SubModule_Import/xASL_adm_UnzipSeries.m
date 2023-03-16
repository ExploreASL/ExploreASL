function [unpackedFiles] = xASL_adm_UnzipSeries(sourcedataDir, tempDerivativeDir, bRecurse)
% Function to unzip the the series of a folder for further processing
%
% FORMAT: x = xASL_adm_UnzipSeries(x[,false, false])
%                     
%
% INPUT:
%   x              x structure containing all input parameters (REQUIRED)
%   bUseGzip       When TRUE then uses GUNZIP 
%                  When FALSE uses UNZIP command (DEFAULT)
%   bRecurse       When TRUE it will recursively unzip all files in the archive (DEFAULT) 
%                  When FALSE it will only unzip the highest layer of the archive 
% OUTPUT:
%   x              x structure containing all input parameters (REQUIRED)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Zip the series of a 
% EXAMPLE: x = xASL_adm_UnzipSeries(x, false, false); Non-recursively zips all matching files
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2023 ExploreASL

if nargin < 2 
    error('Not enough input parameters');
end

if nargin < 3 || isempty(bRecurse)
	bRecurse = true;
end

if nargin > 3
	error('Maximum number of parameters is 3.');
end

if strcmpi(sourcedataDir, tempDerivativeDir)
    error('Input and output directory cannot be the same!');
end

xASL_delete(tempDerivativeDir, true);

%% Unpacking
unpackedFiles = xASL_adm_UnzipRecursive(sourcedataDir, tempDerivativeDir, bRecurse, false);

if isempty(unpackedFiles)
    warning('No files were unpacked')
end

end


% Main function that does the unzipping. 
function unpackedFiles = xASL_adm_UnzipRecursive(dirIn, dirOut, bRecurse, bRemoveSource)

    xASL_adm_CreateDir(dirOut);
    unpackedFiles = xASL_sub_Unzip(dirIn, dirOut, 1, bRemoveSource);

    if ~bRecurse
        return
    end

    directoryList = xASL_adm_GetFileList(dirOut, [], [], [], true);

    for iDir=1:length(directoryList)
        filePath = directoryList{iDir};
        xASL_adm_UnzipRecursive(filePath, filePath, 1, 1);
    end

end


function unpackedFiles = xASL_sub_Unzip(srcDir, destDir, bOverwrite, bRemoveSource)
    if nargin<2 
        warning('not enough input variables')
    end
    
    if nargin<3 || isempty(bOverwrite)
        bOverwrite = false;
    end

    if nargin<4 || isempty(bRemoveSource)
        bRemoveSource = false;
    end

    % Search for all files
    unpackedFiles = [];

    % Looks for case insensitive .zip or .gz files excluding .nii.gz (and in theory .nii.zip). 
    fileList = xASL_adm_GetFileList(srcDir, '(?i).+(?<!\.nii)(\.zip|\.gz)$', [], [], false);

    for iFile=1:length(fileList)
        [ ~, name, ext ] = fileparts(fileList{iFile});
        srcpath = fileList{iFile};
        X = {};

        if strcmpi(ext,'.gz')
            destpath = fullfile(destDir, name );           
            if bOverwrite || ~exist(destpath,'file') 
                X = gunzip(srcpath, destDir);
            end
            if bRemoveSource
                xASL_delete(srcpath);
            end
        elseif strcmpi(ext,'.zip')
            destpath = fullfile(destDir, name );
            if bOverwrite || ~exist(destpath,'file')
                X = unzip(srcpath, destDir);
            end
            if bRemoveSource
                xASL_delete(srcpath);
            end
        end
        if ~isempty(X)
            unpackedFiles{end+1} = X;
        end
    end
end
