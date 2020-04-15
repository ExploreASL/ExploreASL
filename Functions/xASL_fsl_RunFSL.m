function [x, Result1] = xASL_fsl_RunFSL(FSLCommand, x, OutputZipping, NicenessValue, bVerbose)
%xASL_fsl_RunFSL Run FSL from Matlab (ExploreASL)
%
% FORMAT: [x] = xASL_adm_RunFSL(FSLCommand, x[, OutputZipping])
%
% INPUT:
%   FSLCommand      - Command line job for FSL (REQUIRED)
%   x               - structure containing fields with all information required to run this submodule (REQUIRED)
%   OutputZipping   - true for zipping NIfTI output files .nii -> .nii.gz (OPTIONAL, DEFAULT=false)
%   NicenessValue   - the linux nice parameter, a scale with 40 integers 
%                     as index of the priority granted to a process. Lower
%                     = higher priority, higher = lower priority. If no
%                     other processes are running, a lower priority will
%                     still use many resources.
%                     Provide a number between [-20 +19], (OPTIONAL,
%                     DEFAULT=10)
% OUTPUT:
%   x               - as input, outputting FSL dir (OPTIONAL)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function runs an FSL command from ExploreASL:
%              1) Checking the FSL dir
%              2) Manage CUDA/CPU parallelization (currently disabled, WIP)
%              3) Setting up FSL environment
%              4) Running the command
% Supports .nii & .nii.gz, Linux, MacOS & Windows (WSL)
% 
% EXAMPLE: xASL_fsl_RunFSL(FSLCommand, x);
% __________________________________
% Copyright (C) 2015-2019 ExploreASL



%% Admin

if nargin<2 || isempty(x)
    x = struct;
end
if nargin<3 || isempty(OutputZipping)
    OutputZipping = false;
end
if nargin<4 || isempty(NicenessValue)
    NicenessValue = 10;
end
if nargin<5 || isempty(bVerbose)
    bVerbose = true;
end

%% Find FSL directory
[FSLdir, x, RootFSLdir] = xASL_fsl_SetFSLdir(x, x.bAutomaticallyDetectFSL);

if min(isnan(FSLdir))
    warning('No FSL installation found, skipping FSL function');
    return;
end

%% Determine CUDA or openMP (CUDA = multi-core GP, openMP = multi-core CPU)
% if ispc
%     [status, result] = system('wsl ls /usr/local/*cuda*');
% else
%     [status, result] = system('ls /usr/local/*cuda*');
% end
% if status==0
%     ExistCUDA = true;
% else
%     ExistCUDA = false;
% end

%% Define FSL environment script (only declares variables, OK to repeat)
FSLinit0 = ['FSLDIR=' FSLdir ';'];
if exist(fullfile(RootFSLdir,'etc','fslconf','fsl.sh'),'file')
	FSLinit1 = ['. ' FSLdir '/etc/fslconf/fsl.sh;'];
elseif exist(fullfile(RootFSLdir,'etc','fsl','fsl.sh'),'file')
	FSLinit1 = ['. ' FSLdir '/etc/fsl/fsl.sh;'];
else
	warning('Cannot locate fsl.sh, skipping FSL function');
    return;
end
FSLinit2 = 'PATH=${FSLDIR}/bin:${PATH};';
FSLinit3 = 'export FSLDIR PATH;';

FSLinit = [FSLinit0 FSLinit1 FSLinit2 FSLinit3];

if OutputZipping
    FSLoutput = 'FSLOUTPUTTYPE=NIFTI_GZ; export FSLOUTPUTTYPE; ';
    OutputString = '.nii.gz';
else
    FSLoutput = 'FSLOUTPUTTYPE=NIFTI; export FSLOUTPUTTYPE; ';
    OutputString = '.nii';
end

%% Prepend the correct FSL directory if missing
if length(FSLCommand)>5 && strcmp(FSLCommand(1:5),'/bin/')
	if exist(fullfile(RootFSLdir, 'bin'),'dir')
		FSLCommand = [FSLdir FSLCommand];
	elseif exist(RootFSLdir,'dir')
		FSLCommand = [FSLdir FSLCommand(5:end)];
	else
		warning('Cannot locate the command, skipping');
        return;
	end
end

%% Be nice
NiceString = ['nice -' num2str(NicenessValue) ' '];
fprintf('%s\n', ['FSL: NiceNess=' num2str(NicenessValue) ', output=' OutputString]);

%% Run FSL
if ispc
    wslString = 'wsl '; % windows subsystem for linux
else
    wslString = '';
end
if bVerbose
    Result1 = system([wslString FSLinit FSLoutput NiceString FSLCommand], '-echo');
else
    Result1 = system([wslString FSLinit FSLoutput NiceString FSLCommand]);
end
if Result1~=0
    warning('FSL command didnt work nicely:');
end

fprintf('\n');


end