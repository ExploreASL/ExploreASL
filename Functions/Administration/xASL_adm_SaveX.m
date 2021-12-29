function xASL_adm_SaveX(x, Path_xASL, bOverwrite)
%xASL_adm_SaveX Saves x.mat file that keeps track of QC output
%
% FORMAT: xASL_adm_SaveX(x[, Path_xASL, bOverwrite])
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   Path_xASL   - path to the x.mat that contains the QC output (OPTIONAL, DEFAULT = x.dir.SUBJECTDIR/x.mat)
%   bOverwite   - boolean specifying if we overwrite (OPTIONAL, DEFAULT = true)
%
% OUTPUT:
%   n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function saves the x.mat either to the predefined path or the the subject x.mat
%
% EXAMPLE: outside ExploreASL: xASL_adm_SaveX(x, fullfile(x.D.ROOT,'x.mat'));
%           inside ExploreASL: xASL_adm_SaveX(x);
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% __________________________________
% Copyright (C) 2015-2020 ExploreASL



%% -------------------------------------
%  Admin
if nargin<2 || isempty(Path_xASL)
	if ~isfield(x.dir,'SUBJECTDIR')
		error('Path_xASL not provided and x.dir.SUBJECTDIR not defined.');
	end
	Path_xASL = fullfile(x.dir.SUBJECTDIR,'x.mat');
end
if nargin<3 || isempty(bOverwrite)
    bOverwrite = true;
end

if exist(Path_xASL, 'file') && ~bOverwrite
    error('x.mat already existed, skipping...');
else
    xASL_delete(Path_xASL);
    save(Path_xASL,'x');
end


end
