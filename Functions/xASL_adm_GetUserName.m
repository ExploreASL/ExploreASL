function UserName = xASL_adm_GetUserName()
%xASL_adm_GetUserName Get the name of the current user.
%
% FORMAT:       UserName = xASL_adm_GetUserName()
% 
% INPUT:        n/a
%
% OUTPUT:       UserName   - Name of the current user (CHAR ARRAY)
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Get the name of the current user.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      UserName = xASL_adm_GetUserName();
% __________________________________
% Copyright 2015-2020 ExploreASL

    if isunix || ismac 
        UserName = getenv('USER'); 
    else 
        UserName = getenv('username'); 
    end

end

