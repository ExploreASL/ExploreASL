function [ newString ] = xASL_adm_ConvertSlash( StringOriginal,ForceUnix)
%xASL_adm_ConvertSlash Converts Windows forward slashes to backward slashes
% Prevents confusion file separation & regular expression forward slashes
% in Windows

if ~exist('ForceUnix','var')
    ForceUnix   = 0;
end

oldStr  = '\';
oldStr2 = '/';

if  ForceUnix || isunix
    newStr  = '/';
else
	newStr  = '\';
end
	

newString = strrep(StringOriginal,oldStr,newStr);
newString = strrep(newString,oldStr2,newStr);


end

