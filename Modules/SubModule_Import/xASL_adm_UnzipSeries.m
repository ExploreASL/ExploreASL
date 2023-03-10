function [x] = xASL_adm_UnzipSeries(x, bRecurse)
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

if nargin < 2 || isempty(bRecurse)
	bRecurse = true;
end

if nargin > 2
	error('Maximum number of parameters is 3.');
end

%% Check x fields
if ~isfield(x.D, 'sourcedata')
    warning('x.D.sourcedata missing, skipping enpacking');
    return;
elseif ~isfield(x.D, 'derivatives')
    warning('x.D.derivatives missing, skipping unpacking');
    return;
end

x.D.tempOutputDir = fullfile(x.D.derivatives, 'tempSourceUnzip');

xASL_delete(x.D.tempOutputDir, true);
%% Unpacking
unpackedFiles = xASL_adm_UnzipRecursive(x.D.sourcedata, x.D.tempOutputDir, bRecurse, 0);

if isempty(unpackedFiles)
    warning('No files were unpacked')
end

end


% Main function that does the unzipping. 
function unpackedFiles = xASL_adm_UnzipRecursive(dirIn, dirOut, bRecurse, bRemove)

    xASL_adm_CreateDir(dirOut);
    unpackedFiles = xASL_sub_Unzip(dirIn, dirOut, 1, bRemove);

    if ~bRecurse
        return
    end

    directoryList = xASL_adm_GetFileList(dirOut,[],[],[],true);

    for iDir=1:length(directoryList)
        filePath = directoryList{iDir};
        xASL_adm_UnzipRecursive(filePath, filePath, 1, 1);
    end

end


function unpackedFiles = xASL_sub_Unzip(srcDir, destDir, bOverwrite, bRemove)
    if nargin<2 
        warning('not enough input variables')
    end
    
    if nargin<3 || isempty(bOverwrite)
        bOverwrite = false;
    end

    if nargin<4 || isempty(bRemove)
        bRemove = false;
    end

    % Search for all files
    unpackedFiles = [];

    % Looks for case insensitive .zip or .gz files excluding .nii.gz (and in theory .nii.zip). 
    fileList = xASL_adm_GetFileList(srcDir,'(?i).+(?<!\.nii)(\.zip|\.gz)$',[],[],false);

    for iFile=1:length(fileList)
        [ ~, name, ext ] = fileparts(fileList{iFile});
        srcpath = fileList{iFile};
        X = {};

        if strcmpi(ext,'.gz')
            destpath = fullfile(destDir, name );           
            if bOverwrite || ~exist(destpath,'file') 
                X = gunzip(srcpath, destDir);
            end
            if bRemove
                xASL_delete(srcpath);
            end
        elseif strcmpi(ext,'.zip')
            destpath = fullfile(destDir, name );
            if bOverwrite || ~exist(destpath,'file')
                X = unzip(srcpath, destDir);
            end
            if bRemove
                xASL_delete(srcpath);
            end
        end
        if ~isempty(X)
            unpackedFiles{end+1} = X;
        end
    end
end
