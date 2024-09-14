function xASL_adm_RemoveDirectories(root)
%xASL_adm_RemoveDirectories Script to remove all ExploreASL related paths
%
% FORMAT:       xASL_adm_RemoveDirectories(root)
% 
% INPUT:        root - Root directory (REQUIRED, CHAR ARRAY)
%
% OUTPUT:       n/a
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Script to remove all ExploreASL related paths.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      xASL_adm_RemoveDirectories(x.opts.MyPath);
% __________________________________
% Copyright 2015-2021 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________

    %% Check input
    if nargin < 1
        error('Please provide a root directory...');
    end
    
    % ExploreASL
    warning('off', 'MATLAB:rmpath:DirNotFound');
    rmpath(genpath(root));
    warning('on', 'MATLAB:rmpath:DirNotFound');
end
