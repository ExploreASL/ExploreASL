function [VABYdir, x] = xASL_ext_VABYSetDir(x, bAutomaticallyDetectVABY)
%xASL_ext_VABYSetDir Find the VABYdir from Matlab (ExploreASL)
%
% FORMAT: [VABYdir[, x]] = xASL_ext_VABYSetDir(x, bAutomaticallyDetectVABY)
%
% INPUT:
%   x                        - structure containing fields with all information required to run this submodule (OPTIONAL)
%   bAutomaticallyDetectVABY - Boolean to automatically detect the VABY version
%                              if disabled, this function will try to use the system-initialized VABY
%                              and throw an error if VABY is not initialized
%                              (OPTIONAL, DEFAULT = disabled)
% OUTPUT:
%   VABYdir    - path to VABY (REQUIRED)
%   x         - as input, outputting VABY dir (OPTIONAL)
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function finds the VABYdir & puts it out, also in
%              x.VABYdir to allow repeating this function without having to repeat
%              searching. If the VABYdir is already defined in x.VABYdir, this function
%              is skipped. Currently, it only checks the dir, automatic localization does not work yet
% 
% EXAMPLE: VABYdir = xASL_ext_VABYSetDir(x);
% __________________________________
% Copyright (C) 2015-2024 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


%% Admin
if nargin<1
    x = struct;
end
if nargin<2 || isempty(bAutomaticallyDetectVABY)
    if isfield(x,'external') && isfield(x.external, 'bAutomaticallyDetectVABY')
        bAutomaticallyDetectVABY = x.external.bAutomaticallyDetectVABY;
    else
        bAutomaticallyDetectVABY = false;
    end
end

VABYdir = NaN;

if isfield(x,'VABYdir') && ~isempty(x.VABYdir)
    % if we already have an VABY dir, skip this function
    VABYdir = x.VABYdir;
    return;
end

% For VABY quantification, we cannot run automatic detection
if isfield(x, 'external') && isfield(x.external, 'ExternalQuantificationType') && strcmp(x.external.ExternalQuantificationType, 'VABY')
	error('External quantification with VABY requested. You need to provide path in x.VABYdir');
end


%% Detect OS
if ismac
    fprintf('Running VABY from Matlab on macOS\n');
elseif isunix % check for linux (also used for macOS)
    fprintf('Running VABY from Matlab on Linux\n');
elseif ispc
    fprintf('Running VABY from Matlab on Linux\n');
end

%% AUTOMATIC DETECTION NOT YET IMPLEMENTED

end