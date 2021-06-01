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
% EXAMPLE:      xASL_adm_RemoveDirectories(x.MyPath);
% __________________________________
% Copyright 2015-2021 ExploreASL

    %% Check input
    if nargin < 1
        error('Please provide a root directory...');
    end
    
    % ExploreASL
    try
        warning('off', 'MATLAB:rmpath:DirNotFound');
        rmpath(genpath(root));
        warning('on', 'MATLAB:rmpath:DirNotFound');
    catch ME
        fprintf('%s\n', ME.message);
    end

end




