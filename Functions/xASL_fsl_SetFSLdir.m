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
if isunix % check for linux
    fprintf('Running FSL from Matlab on Linux');
elseif ismac
    warning('This Matlab-FSL implementation not yet tested for Mac, skipping');
    return;
elseif ispc
    [status, status2] = system('wsl ls;'); % leave status2 here, otherwise system will produce output
    if status~=0
        warning('Detected windows PC without WSL, needs to be installed first, skipping');
        fprintf('This warning can be ignored if you dont need to run FSL-specific processing, e.g. TopUp\n');
        return;
    end
end

%% Try searching at different ROOT paths
PathApps = {'/data/usr/local' '/usr/local' '/opt/amc' '/usr/local/bin'};
PathDirect = {'/usr/lib/fsl/5.0'};
if ispc
    CurrDir = pwd;
    [~, result] = system('echo %LOCALAPPDATA%');
    SearchDir = fullfile(result(1:end-1), 'Packages');
    % Check distros to save time
    Distros = {'CanonicalGroupLimited.Ubuntu' 'WhitewaterFoundryLtd' 'Ubuntu' 'SUSE' 'Kali' 'Debian'};
    DistroDir = {};
    for iD=1:length(Distros)
        if isempty(DistroDir)
            DistroDir = xASL_adm_GetFileList(SearchDir , ['.*' Distros{iD} '.*'], 'FPList', [0 Inf], true);
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
    cd(CurrDir);
end

%% Try searching some more
FSLdir = {};
for iP=1:length(PathApps)
    if isempty(FSLdir)
        FSLdir = xASL_adm_GetFileList(lower(PathApps{iP}), '^fsl.*', 'FPList', [0 Inf], true); % we can set this to recursive to be sure, but this will take a long time
    end
end

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