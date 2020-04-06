function [FSLdir, x, RootFSLDir] = xASL_fsl_SetFSLdir(x)
%xASL_fsl_SetFSLdir Find the FSLdir from Matlab (ExploreASL)
%
% FORMAT: [FSLdir[, x, RootFSLDir]] = xASL_adm_SetFSLdir(x)
%
% INPUT:
%   x       - structure containing fields with all information required to run this submodule (OPTIONAL)
%
% OUTPUT:
%   FSLdir    - path to FSL (REQUIRED)
%   x         - as input, outputting FSL dir (OPTIONAL)
%   RootFSLDir - if emulated linux (e.g. WSL) this differs from FSLdir
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function finds the FSLdir & puts it out, also in
%              x.FSLdir to allow repeating this function without having to repeat
%              searching.
%              If the FSLdir & RootFSLDir are already defined in x.FSLdir & x.RootFSLDir, this function
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

FSLdir = NaN;
RootFSLDir = {};

if isfield(x,'FSLdir') && isfield(x,'RootFSLDir') && ~isempty(x.FSLdir) && ~isempty(x.RootFSLDir)
    % if we already have an FSL dir, skip this function
    FSLdir = x.FSLdir;
    RootFSLDir = x.RootFSLDir;
    return;
end

%% Detect OS
if isunix % check for linux (also used for macOS)
    fprintf('Running FSL from Matlab on Linux');
elseif ispc
    [status, ~] = system('wsl ls;'); % leave status2 here, otherwise system will produce output
    if status~=0
        fprintf('Detected windows PC without WSL, needs to be installed first if you want to use FSL');
        fprintf('This warning can be ignored if you dont need to run FSL-specific processing, e.g. TopUp\n');
        return;
    end
end

% %% Initialize RootFSLDir
% RootFSLDir = ''; NEW
% 
% %% Check if FSL is initialized by system
% [~, result2] = system('which fsl'); NEW
% RootFSLDir{end+1} = fileparts(result2); NEW

%% Try searching at different ROOT paths
% PathApps = {'/data/usr/local' '/usr/local' '/opt/amc' '/usr/local/bin' '/usr/local/apps' '/usr/lib'}; NEW
PathApps = {'/data/usr/local' '/usr/local' '/opt/amc' '/usr/local/bin' '/usr/local/apps/fsl'};
PathDirect = {'/usr/lib/fsl/5.0'};

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
        RootFSLDir = xASL_adm_GetFileList(fullfile(DistroDir{1},'LocalState'), '^rootfs$', 'FPList', [0 Inf], true);
        if isempty(RootFSLDir)
            RootFSLDir = xASL_adm_GetFileList(DistroDir{1}, '^rootfs$', 'FPListRec', [0 Inf], true);
        end
    end
    if isempty(RootFSLDir)
        warning('Couldnt find rootfs of WSL, skipping');
        return;
    end
    RootFSLDir = RootFSLDir{1};

    PathApps{end+1} = fullfile(RootFSLDir,'usr','local');
end

%% Search for first subfolder
FSLdir = {};
for iP=1:length(PathApps)
    if isempty(FSLdir)
        FSLdir = xASL_adm_GetFileList(lower(PathApps{iP}), '^fsl.*', 'FPList', [0 Inf], true); % we can set this to recursive to be sure, but this will take a long time
    end
end

%% Search for second subfolder


%% Try searching some more - for directories 
for iP=1:length(PathDirect)
    if isempty(FSLdir)
		if exist(PathDirect{iP},'dir')
			FSLdir = {PathDirect{iP}};
		end
    end
end

if length(FSLdir)<1
    warning('Cannot find valid FSL installation dir');
    FSLdir = NaN;
    return;
elseif length(FSLdir)>1
    fprintf('%s\n','Found more than 1 FSL version, choosing latest');
end
FSLdir = FSLdir{end};
if ispc
    FSLdirWin = FSLdir;
    FSLdir = FSLdirWin(length(RootFSLDir)+1:end);
end
FSLdir = strrep(FSLdir,'\','/');

%% If FSL is installed in a subfolder, find it
if ~exist(fullfile(FSLdir,'bin'),'dir') || ~exist(fullfile(FSLdir,'bin','fsl'),'file') || ~exist(fullfile(FSLdir,'bin','bet'),'file')
    fprintf('Found potential FSL folder but not the installation, trying subfolders, this might take a while...\n');
    FSLsubdir = xASL_adm_GetFileList(FSLdir, '^fsl.*', 'FPListRec', [0 Inf], true); % now we do this recursively
    % now select the latest folder that has a bin folder in them (assuming
    % that there can be multiple fsl installations)
    if ~isempty(FSLsubdir)
        BinDir = cellfun(@(x) fullfile(x,'bin'), FSLsubdir, 'UniformOutput',false);
        BetPath = cellfun(@(x) fullfile(x,'bet'), BinDir, 'UniformOutput',false);
        fslFile = cellfun(@(x) fullfile(x,'fsl'), BinDir, 'UniformOutput',false);
        ExistBin = cellfun(@(x) logical(exist(x,'dir')), BinDir);
        ExistBet = cellfun(@(x) logical(exist(x,'file')), BetPath);
        ExistfslFile = cellfun(@(x) logical(exist(x,'file')), fslFile);
        
        DirIndex = max(find(ExistBin & ExistBet & ExistfslFile)); % find the latest dir that has the folders and functions we anticipate it should have
        FSLdir = FSLsubdir{DirIndex};
    else
        warning('Cannot find valid FSL installation dir');
        FSLdir = NaN;
        return;        
    end
end        


%% Manage RootFSLDir
if isnumeric(RootFSLDir) && isnan(RootFSLDir)
    RootFSLDir = FSLdir; % default
else
    RootFSLDir = fullfile(RootFSLDir, FSLdir); % for e.g. WSL
end

%% Add to x
x.FSLdir = FSLdir;
x.RootFSLDir = RootFSLDir;

end