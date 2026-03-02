function [VABYdir, x] = xASL_ext_VABYSetDir(x)
%xASL_ext_VABYSetDir Find the VABYdir from Matlab (ExploreASL)
%
% FORMAT: [VABYdir[, x]] = xASL_ext_VABYSetDir(x)
%
% INPUT:
%   x                        - structure containing fields with all information required to run this submodule (OPTIONAL)
%
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
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


%% Admin
if nargin<1
    x = struct;
end

if isfield(x,'external') && isfield(x.external, 'bAutomaticallyDetectVABY')
	bAutomaticallyDetectVABY = x.external.bAutomaticallyDetectVABY;
else
	bAutomaticallyDetectVABY = false;
end

VABYdir = ''; % Initialize the output as an empty string

if isfield(x, 'external') && isfield(x.external,'VABYdir') && ~isempty(x.external.VABYdir)
    % if we already have an VABY dir, skip this function
    VABYdir = x.external.VABYdir;
    return;
end

% For VABY quantification, we cannot run automatic detection
if isfield(x, 'modules') && isfield(x.modules, 'asl') && isfield(x.modules.asl, 'ExternalQuantificationType') && strcmp(x.modules.asl.ExternalQuantificationType, 'VABY')
	error('External quantification with VABY requested. You need to provide the path to the VABY command in x.VABYdir');
end

%% AUTOMATIC DETECTION NOT YET IMPLEMENTED

end