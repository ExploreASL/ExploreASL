function [x, Result1] = xASL_ext_VABYRun(VABYCommand, x, NicenessValue, bVerbose)
%xASL_ext_VABYRun Run VABY from Matlab (ExploreASL)
%
% FORMAT: [x] = xASL_adm_RunVABY(VABYCommand, x[, NicenessValue, bVerbose])
%
% INPUT:
%   VABYCommand     - Command line job for VABY (REQUIRED)
%   x               - structure containing fields with all information required to run this quantification (REQUIRED)
%   NicenessValue   - the linux nice parameter, a scale with 40 integers 
%                     as index of the priority granted to a process. Lower
%                     = higher priority, higher = lower priority. If no
%                     other processes are running, a lower priority will
%                     still use many resources.
%                     Provide a number between [-20 +19], (OPTIONAL,
%                     DEFAULT=10)
%   bVerbose        - verbose output (OPTIONAL, DEFAULT true)
%
% OUTPUT:
%   x               - as input, outputting VABY dir (OPTIONAL)
%   Result1         - Result1 describes if the execution was successful
%                     (0 = successful, NaN = no VABY found, 1 or other = something failed)
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function runs an VABY command from ExploreASL:
%
%              1. Checking the VABY dir
%              2. Running the command
%
% Supports .nii & .nii.gz, Linux, MacOS & Windows
% 
% EXAMPLE: xASL_ext_VABYRun(VABYCommand, x);
% __________________________________
% SPDX-License-Identifier: Apache-2.0
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________



%% Admin
if nargin<2 || isempty(x)
    warning('x input missing');
    x = struct;
end
if nargin<3 || isempty(NicenessValue)
    NicenessValue = 10;
end
if nargin<4 || isempty(bVerbose)
    bVerbose = true;
end

% Defaults
Result1 = NaN;

%% Find VABY directory
[VABYdir, x] = xASL_ext_VABYSetDir(x);

if isempty(VABYdir) || ~xASL_exist(VABYdir, 'dir')
    % Script will return Result1=NaN to show that there is no VABY
    % installation found
    warning('No VABY installation found, skipping VABY function');
    return;
end

% Check if VABYcommand contains the command and parameters, or if it also contains the full path to the vaby command
if strcmp(VABYCommand(1:4),'vaby')
	% Check if the path to the vaby command was provided, if not, add the path
	VABYCommand = fullfile(VABYdir, VABYCommand);
end

%% Be nice
NiceString = ['nice -' num2str(NicenessValue) ' '];
fprintf('%s\n', ['VABY: NiceNess=' num2str(NicenessValue)]);

%% Run VABY
if bVerbose
    Result1 = system([NiceString VABYCommand], '-echo');
else
    Result1 = system([NiceString VABYCommand]);
end
if Result1~=0
    warning('VABY command didnt work nicely:');
end

fprintf('\n');

end
