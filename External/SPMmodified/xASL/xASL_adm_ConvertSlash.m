function [newString] = xASL_adm_ConvertSlash(StringOriginal, ForceUnix)
%xASL_adm_ConvertSlash Converts Windows forward slashes to backward slashes
% Prevents confusion file separation & regular expression forward slashes
% in Windows
%
% FORMAT: [newString] = xASL_adm_ConvertSlash( StringOriginal,ForceUnix)
%
% INPUT:
%   StringOriginal	    - Original string
%   ForceUnix           - Force slash direction
%
% OUTPUT:
%   newString           - Resulting string
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: Converts Windows forward slashes to backward slashes
% Prevents confusion file separation & regular expression forward slashes
% in Windows.
% EXAMPLE:     n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
%
% __________________________________
% Copyright 2015-2021 ExploreASL

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

