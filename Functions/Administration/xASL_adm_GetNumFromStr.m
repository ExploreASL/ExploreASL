function num = xASL_adm_GetNumFromStr(str)
%xASL_adm_GetNumFromStr Obtains single number from string. 
% CAVE there should only be one number!
%
% FORMAT:       num = xASL_adm_GetNumFromStr(str)
% 
% INPUT:        ...
%
% OUTPUT:       ...
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Obtains single number from string. 
%               {{CAVE}} there should only be one number!
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      ...
% __________________________________
% Copyright 2015-2020 ExploreASL
% Licensed under Apache 2.0, see permissions and limitations at
% https://github.com/ExploreASL/ExploreASL/blob/main/LICENSE
% you may only use this file in compliance with the License.
% __________________________________



    indices     = regexp(str,'\d');
    num         = str2num(str(indices(1):indices(end)));

end

