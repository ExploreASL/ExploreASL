function [FSLdir, x, RootWSLdir] = xASL_fsl_SetFSLdir(x, bAutomaticallyDetectFSL)
%xASL_fsl_SetFSLdir Find the FSLdir from Matlab (ExploreASL)
%
% FORMAT: [FSLdir[, x, RootWSLdir]] = xASL_adm_SetFSLdir(x, bUseLatestVersion)
%
% INPUT:
%   x                       - structure containing fields with all information required to run this submodule (OPTIONAL)
%   bAutomaticallyDetectFSL - Boolean to automatically detect the FSL version
%                             if disabled, this function will try to use the system-initialized FSL 
%                             and throw an error if FSL is not initialized
%                             (OPTIONAL, DEFAULT = disabled)
% OUTPUT:
%   FSLdir    - path to FSL (REQUIRED)
%   x         - as input, outputting FSL dir (OPTIONAL)
%   RootWSLdir - if emulated linux (e.g. WSL) this differs from FSLdir
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function finds the FSLdir & puts it out, also in
%              x.FSLdir to allow repeating this function without having to repeat
%              searching.
%              If the FSLdir & RootFSLdir are already defined in x.FSLdir & x.RootFSLdir, this function
%              is skipped.
%              Supports Linux, MacOS & Windows (WSL), & several different
%              default installation folders for different Linux
%              distributions
% 
% EXAMPLE: FSLdir = xASL_fsl_SetFSLdir(x);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL



%% Admin
if nargin<1
    x = struct;
end
if nargin<2 || isempty(bAutomaticallyDetectFSL)
    if isfield(x,'external') && isfield(x.external,'bAutomaticallyDetectFSL')
        bAutomaticallyDetectFSL = x.external.bAutomaticallyDetectFSL;
    else
        bAutomaticallyDetectFSL = false;
    end
end

FSLdir = NaN;
RootWSLdir = NaN;
RootFSLdir = '';

if isfield(x,'FSLdir') && isfield(x,'RootFSLdir') && ~isempty(x.FSLdir) && ~isempty(x.RootFSLdir)
    % if we already have an FSL dir, skip this function
    FSLdir = x.FSLdir;
    RootWSLdir = x.RootFSLdir;
    return;
end

%% Detect OS
if isunix % check for linux (also used for macOS)
    fprintf('Running FSL from Matlab on Linux\n');
elseif ismac
    fprintf('Running FSL from Matlab on macOS\n');
elseif ispc
    [status, ~] = system('wsl ls;'); % leave status2 here, otherwise system will produce output
    if status~=0
        fprintf('Detected windows PC without WSL, needs to be installed first if you want to use FSL');
        fprintf('This warning can be ignored if you dont need to run FSL-specific processing, e.g. TopUp\n');
        return;
    end
end

%% 1) Check if FSL is initialized by system
if ispc 
    [bSuccess, result2] = system('wsl which fsl');
else
    [bSuccess, result2] = system('which fsl');
end 
bSuccess = bSuccess==0;

if bSuccess
    RootFSLdir{end+1} = fileparts(result2);
end

if bAutomaticallyDetectFSL
    %% 2) Try searching at different ROOT paths, first layer subfolder
    PathApps = {'/data/usr/local' '/usr/local' '/opt/amc' '/usr/local/bin' '/usr/local/apps' '/usr/lib'};

    if ispc % for Windows Subsystem of Linux
        [~, result] = system('echo %LOCALAPPDATA%');
        SearchDir = fullfile(strtrim(result), 'Packages');
        % Check distros to save time
        Distros = {'CanonicalGroupLimited.Ubuntu' 'WhitewaterFoundryLtd' 'Ubuntu' 'SUSE' 'Kali' 'Debian'};
        DistroDir = {};
        for iD=1:length(Distros)
            if isempty(DistroDir)
                DistroDir = xASL_adm_GetFileList(SearchDir , ['.*' Distros{iD} '.*'], 'FPList', [0 Inf], true);
                break
            end
        end
        if ~isempty(DistroDir)
            RootWSLdir = xASL_adm_GetFileList(fullfile(DistroDir{1},'LocalState'), '^rootfs$', 'FPList', [0 Inf], true);
            if isempty(RootWSLdir)
                RootWSLdir = xASL_adm_GetFileList(DistroDir{1}, '^rootfs$', 'FPListRec', [0 Inf], true);
            end
        end
        if isempty(RootWSLdir)
            warning('Couldnt find rootfs (filesystem) of WSL, skipping');
            return;
        end
        PathApps{end+1} = fullfile(RootWSLdir,'usr','local');
    end

    %% Search for first & second-layer subfolder
    for iP=1:length(PathApps)
        TempDir = xASL_adm_GetFileList(lower(PathApps{iP}), '(?i)^fsl.*', 'FPList', [0 Inf], true); % we can set this to recursive to be sure, but this will take a long time
        if ~isempty(TempDir)
            for iDir=1:length(TempDir)
                RootFSLdir{end+1} = TempDir{iDir};
                % Also search for second subfolder
                TempDir2 = xASL_adm_GetFileList(RootFSLdir, '(?i)^fsl.*', 'FPList', [0 Inf], true);
                if ~isempty(TempDir2)
                    for iDir2=1:length(TempDir2)
                        RootFSLdir{end+1} = TempDir2{iDir2};
                    end
                end
            end
        end
    end
end

%% Create new list with valid FSL folders
FSLdir = '';
if isempty(RootFSLdir) || isempty(RootFSLdir{end})
    RootFSLdir = NaN;
    FSLdir = NaN;
    RootWSLdir = NaN;
    fprintf('Warning, FSL folder not found!\n');
    return;
end
    
for iDir=1:length(RootFSLdir)
    % First remove 'bin'
    if strcmp(RootFSLdir{iDir}(end-2:end), 'bin')
        RootFSLdir{iDir} = fileparts(RootFSLdir{iDir});
    end
    BinDir = fullfile(RootFSLdir{iDir}, 'bin');
    % Check if root-dir & bin-dir exist
    if exist(RootFSLdir{iDir}, 'dir') && exist(BinDir, 'dir')
        PathBET = fullfile(BinDir, 'bet');
        PathFSL = fullfile(BinDir, 'fsl');
        PathFIT = fullfile(BinDir, 'dtifit');
        % check if few fsl-scripts exist
        if exist(PathBET, 'file') && exist(PathFSL, 'file') && exist(PathFIT, 'file')
            FSLdir{end+1} = RootFSLdir{iDir};
        end
    end
end

% Remove doubles
FSLdir = sort(unique(FSLdir));

%% Pick FSLdir
if length(FSLdir)<1
    warning('Cannot find valid FSL installation dir');
    FSLdir = NaN;
    RootWSLdir = NaN;
    return;
elseif length(FSLdir)>1
    fprintf('%s\n','Found more than 1 FSL version, choosing latest');
else
    % dont say anything
end

FSLdir = FSLdir{end}; % default to most recent version (and convert from cell to char array

if ispc
    if ~isnan(RootWSLDir)
        RootWSLDir = RootWSLdir{end}; % default to most recent version
    end
    FSLdir = FSLdir(length(RootWSLDir)+1:end);
end
FSLdir = strrep(FSLdir,'\','/');

% %% If FSL is installed in a subfolder, find it
% if ~exist(fullfile(FSLdir,'bin'),'dir') || ~exist(fullfile(FSLdir,'bin','fsl'),'file') || ~exist(fullfile(FSLdir,'bin','bet'),'file')
%     fprintf('Found potential FSL folder but not the installation, trying subfolders, this might take a while...\n');
%     FSLsubdir = xASL_adm_GetFileList(FSLdir, '^fsl.*', 'FPListRec', [0 Inf], true); % now we do this recursively
%     % now select the latest folder that has a bin folder in them (assuming
%     % that there can be multiple fsl installations)
%     if ~isempty(FSLsubdir)
%         BinDir = cellfun(@(x) fullfile(x,'bin'), FSLsubdir, 'UniformOutput',false);
%         BetPath = cellfun(@(x) fullfile(x,'bet'), BinDir, 'UniformOutput',false);
%         fslFile = cellfun(@(x) fullfile(x,'fsl'), BinDir, 'UniformOutput',false);
%         ExistBin = cellfun(@(x) logical(exist(x,'dir')), BinDir);
%         ExistBet = cellfun(@(x) logical(exist(x,'file')), BetPath);
%         ExistfslFile = cellfun(@(x) logical(exist(x,'file')), fslFile);
%         
%         DirIndex = max(find(ExistBin & ExistBet & ExistfslFile)); % find the latest dir that has the folders and functions we anticipate it should have
%         FSLdir = FSLsubdir{DirIndex};
%     else
%         warning('Cannot find valid FSL installation dir');
%         FSLdir = NaN;
%         return;        
%     end
% end        


%% Manage RootFSLdir
if ~exist('RootWSLdir','var') || isnan(RootWSLdir)
    RootWSLdir = FSLdir; % default
else
    RootWSLdir = fullfile(RootWSLdir, FSLdir); % for e.g. WSL
end

%% Add to x
x.FSLdir = FSLdir;
x.RootFSLdir = RootWSLdir;


end