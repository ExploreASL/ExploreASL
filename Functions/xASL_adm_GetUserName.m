function UserName = xASL_adm_GetUserName()
%xASL_adm_GetUserName ...
%
% FORMAT:        UserName = xASL_adm_GetUserName()
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  ...
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL

    if isunix() 
        UserName = getenv('USER'); 
    else 
        UserName = getenv('username'); 
    end

end

