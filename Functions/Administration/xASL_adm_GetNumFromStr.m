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
% SPDX-License-Identifier: Apache-2.0
% __________________________________



    indices     = regexp(str,'\d');
    num         = str2num(str(indices(1):indices(end)));

end

