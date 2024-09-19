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
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________


    if isunix
        UserName = getenv('USER'); 
    else 
        UserName = getenv('username'); 
    end

end

