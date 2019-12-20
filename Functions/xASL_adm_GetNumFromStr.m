function num = xASL_adm_GetNumFromStr( str )
%xASL_adm_GetNumFromStr Obtains single number from string. 
% CAVE there should only be one number!


    indices     = regexp(str,'\d');
    num         = str2num(str(indices(1):indices(end)));

end

