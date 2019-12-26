function [InPath] = xASL_spm_admin(InPath, bPadComma1)
%xASL_spm_admin Force ,1 at end of IMname.
% This is useful for refIM/srcIM in CoregInit, OldNormalizeWrapper etc
% IMname should be a cell
% bPadComma1 is to add the ,1 to the end of the pathstring, which SPM uses
% to assign the first image of a 4D image array (OPTIONAL, DEFAULT = true)

if nargin<2 || isempty(bPadComma1)
    bPadComma1 = true; % default
end


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
    if bPadComma1
        String2Pad = ',1';
    else
        String2Pad = '';
    end
    
    FoundStr = strfind(InPath{1}, ',1');
    if isempty(FoundStr)
        InPath{1} = [InPath{1} String2Pad];
    else
        FoundStr = FoundStr(length(FoundStr));
        if FoundStr+1~=length(InPath{1})
            InPath{1} = [InPath{1} String2Pad];
        end
    end

end
