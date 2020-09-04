function xASL_adm_SaveX(x, Path_xASL)
%xASL_adm_SaveX Saves x.mat file that keeps track of QC output
% FORMAT: xASL_adm_SaveX(x[, Path_xASL])
%
% INPUT:
%   x           - structure containing fields with all information required to run this submodule (REQUIRED)
%   Path_xASL   - path to the x.mat that contains the QC output (OPTIONAL, DEFAULT = x.SUBJECTDIR/x.mat)
%
% OUTPUT:
%   n/a
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This function saves the x.mat either to the predefined path or the the subject x.mat
%
% EXAMPLE: xASL_adm_SaveX(x, fullfile(x.D.ROOT,'x.mat'));
%
% __________________________________
% Copyright (C) 2015-2020 ExploreASL



%% -------------------------------------
%  Admin
if nargin<2 || isempty(Path_xASL)
	if ~isfield(x,'SUBJECTDIR')
		error('Path_xASL not provided and x.SUBJECTDIR not defined.');
	end
	Path_xASL = fullfile(x.SUBJECTDIR,'x.mat');
end

xASL_delete(Path_xASL);
save(Path_xASL,'x'); 

end
