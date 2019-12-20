function UserName = xASL_adm_GetUserName()
    if isunix() 
        UserName = getenv('USER'); 
    else 
        UserName = getenv('username'); 
    end

end

