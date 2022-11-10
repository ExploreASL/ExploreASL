function px = xASL_adm_FontPtToPx(pt)
% xASL_adm_FontPtToPx converts font point size to pixel size
%
% FORMAT:       px = xASL_adm_FontPtToPx(pt)
% 
% INPUT:        font size in points
%
% OUTPUT:       font size in pixels
% 
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION:  Converts numeric font size in points to pixel height 
% Use           Used in printing of graphics for pixel positioning calculations.
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE:      px = xASL_adm_FontPtToPx(12);
% __________________________________
% Copyright 2015-2022 ExploreASL

if ismac || isunix
    px = floor(pt * 1.2);
elseif ispc
    px = floor(pt * 1.33);
end


end