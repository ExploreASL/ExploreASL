function [InPath] = xASL_spm_admin(InPath)
%xASL_spm_admin Force ,1 at end of IMname.
% This is useful for refIM/srcIM in CoregInit, OldNormalizeWrapper etc
% IMname should be a cell

    %% Manage unzipping
    [~, InPath]  = xASL_io_ReadNifti(InPath);

    %% -----------------------------------------------------
    %% Manage iscell
	if ~iscell(InPath)
		tempL = InPath;
		IMnameOut{1} = tempL;
	else
		IMnameOut = InPath;
    end

    InPath = IMnameOut;

    %% -----------------------------------------------------
    %% Manage ',1' suffix
    FoundStr = strfind(InPath{1}, ',1');
    if isempty(FoundStr)
        InPath{1} = [InPath{1} ',1'];
    else
        FoundStr = FoundStr(length(FoundStr));
        if FoundStr+1~=length(InPath{1})
            InPath{1} = [InPath{1} ',1'];
        end
    end

end
